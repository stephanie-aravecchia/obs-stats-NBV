#ifndef INFORMATION_GAIN_H
#define INFORMATION_GAIN_H
#include <vector>
#include <array>
#include <assert.h>
#include <map>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <pcl_ros/point_cloud.h>
#include <pcl_ros/transforms.h>
#include <pcl/point_types.h>
#include <obs_stat/stat_point_type.h>
#include "explo_quality_aware/Candidate.h"
#include "explo_quality_aware/PointRay.h"
#include "explo_quality_aware/Policy.h"
#include "explo_quality_aware/proba_point_type.h"
#include <typeinfo>
#include <list>
#include <functional>
#include <thread>
#include <mutex>

//That class is the core of the NBV on viewpoint stats exploration policy
//Based on a current obs-stat-grid and a costmap, it select the reachable candidates,
//and computes the Information Gain for each candidate
class InformationGain {
    protected:
        
        pcl::PointXYZ robot_position_;
        double resolution_;
        double search_range_;
        double sensor_range_;
        grid_tools::BBox igBBox_;
        double occupancy_threshold_ = OCC_THRESHOLD;
        pcl::PointCloud<PointStat> point_stat_cloud_;
        PolicyChoice policy_;
        bool debug_{false};
        //List of occupied points, with their probability stored
        std::vector<PointXYZP> list_of_occupied_points_;
        //std::vector<Candidate> list_of_candidates_;
        std::list<Candidate> list_of_candidates_;
        std::mutex candidate_mutex_; 
        //VoxelGrid reprensentation of the octomap probabilities
        grid_tools::VoxelGrid occupancy_vox_grid_;

        size_t num_threads_ = 8;	
        
        
        std::multimap<double, pcl::PointXYZ, std::less<double> > ig_per_candidate_;
        std::vector<std::pair<double, pcl::PointXYZ> > ig_per_candidate_vector_;
        std::vector<explo_utils::CandidateResult>  result_per_candidate_;
        std::mutex results_mutex_; 
        
        void setVoxelGridFromPointCloud(const grid_tools::BBox & bbox, const pcl::PointCloud<PointXYZP>& pc);
        
        //gets list of reachable candidates (possible goals to visit) from the costmap
        //within a range r from the robot position
        void setListOfReachableCandidatesFromCostmap(const pcl::PointXYZ&, const explo_utils::CostMap&, double r);
        //get the bbox around the current robot position in which the InformationGain is calculated
        void setBBoxFromRanges();
        //And constrain to a valid bbox        
        void setBBoxFromRangesAndValid(const grid_tools::BBox&);        
        
        //Interate on the point cloud, and stores the occupied Points inside 
        //the BBOx into a vector
        void setListOfOccupiedPointsFromPointCloud(const pcl::PointCloud<PointXYZP>&); 

        static bool point3dLess(const cv::Point3i& v1, const cv::Point3i& v2) 
        {   if (v1.x < v2.x) return true;
            if (v1.x > v2.x) return false;
            // v1.x == v2.x
            if (v1.y < v2.y) return true;
            if (v1.y > v2.y) return false;
            // v1.x == v2.x && v1.y == v2.y
            return v1.z < v2.z;
        };
        
        typedef std::map<cv::Point3i, std::pair<double,obs_stat_tools::SphericalPoint>, decltype(point3dLess)*> GridIntersectMap;
        
        //get the point_ray: 
        //set the origin with the point hit, an occupied point
        //update the vector in point_ray with the proba until the coordinates of the candidate
        void getPointRayFromPointToCandidate(PointRay& point_ray, const pcl::PointXYZ& point, const pcl::PointXYZ& candidate_coord) const;
        
        //Raycast operation
        void getGridSegmentIntersect(const pcl::PointXYZ& orig_point, const pcl::PointXYZ& target_point, GridIntersectMap& grid_intersect) const;
        
        //core
        void updateStatsInListOfCandidates();
        void processListOfCandidates();


    public:
        InformationGain() {};
        InformationGain(const pcl::PointXYZ& p, double resolution, double search_range, double sensor_range, const explo_utils::CostMap& cm, 
                const pcl::PointCloud<PointXYZP>& pc, const pcl::PointCloud<PointStat>& stat_cloud, PolicyChoice policy, const grid_tools::BBox& valid_box, bool debug);

        const grid_tools::BBox& getIgBBox() const { return igBBox_;};
        const pcl::PointXYZ& getRobotPosition() const {return robot_position_;};
        const std::vector<PointXYZP>& getListOfOccupiedPoints() const {return list_of_occupied_points_;};
        const grid_tools::VoxelGrid& getOccupancyVoxGrid() const {return occupancy_vox_grid_;};
        const std::vector<std::pair<double, pcl::PointXYZ> >& getIGPerCandidate() const {return 
            ig_per_candidate_vector_;};
        const std::vector<explo_utils::CandidateResult>& getResultsPerCandidate() const {return 
            result_per_candidate_;};
};

#endif
