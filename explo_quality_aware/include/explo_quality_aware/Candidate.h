#ifndef CANDIDATE_H
#define CANDIDATE_H
#include <vector>
#include <set>
#include <pcl_ros/point_cloud.h>
#include <pcl_ros/transforms.h>
#include <pcl/point_types.h>
#include "obs_stat/stat_point_type.h"
#include "obs_stat/obs_stat_lib.h"
#include "explo_quality_aware/PointRay.h"
#include "explo_quality_aware/Policy.h"


//Class Candidate computes the information gain for a single candidate
class Candidate : public ObsStatGridBase {
    
    protected:
        pcl::PointXYZ coords_;
        pcl::PointXYZ robot_;
        double sensor_range_;
        grid_tools::StatGrid init_stat_grid_;
        double information_gain_;
        std::vector<PointRay> point_ray_list_;
        std::vector<PointXYZP> points_to_raycast_;

        typedef std::map<cv::Point3i, double, decltype(point3dLess)*> VoxPointMap;
        
        int n_times_update_;
        PolicyChoice policy_choice_;
        double cost_;
        
        void update(const cv::Point3i& vox, const obs_stat_tools::SphericalPoint& obs);
        void updateAngularSectors(const cv::Point3i& vox);
        uint32_t getAngularSectorsForVoxFromCandidate(const cv::Point3i& vox) const;
        uint32_t updateAngularSectorsClockwise(uint32_t in_sectors, uint8_t start, uint8_t end) const;
        uint32_t updateAngularSectorsCounterClockwise(uint32_t in_sectors, uint8_t start, uint8_t end) const;

        void updateObsStatWithPointRay(const PointRay&);
        void setBBox(double, const grid_tools::BBox&);
        void voxGridFromPointCloud(const pcl::PointCloud<PointStat>&);
        void initializePointStatsCoords();

    public:
        Candidate() {} ;
        //Candidate(const pcl::PointXYZ& coords): coords_(coords) {};
        //Candidate(const pcl::PointXYZ& coords, double resolution, double sensor_range, pcl::PointXYZ& robot, grid_tools::PolicyChoice policy, double times_per_meter, double cost): 
        //Candidate(const pcl::PointXYZ& coords, double resolution, double sensor_range, const grid_tools::BBox& igbbox, pcl::PointXYZ& robot, PolicyChoice policy, double cost): 
        //    coords_(coords), resolution_(resolution), sensor_range_(sensor_range), policy_choice_(policy), cost_(cost) {
        Candidate(const pcl::PointXYZ& coords, double resolution, double sensor_range, const grid_tools::BBox& igbbox, pcl::PointXYZ& robot, PolicyChoice policy, double cost) {
                coords_ = coords;
                robot_ = robot;
                resolution_ = resolution;
                sensor_range_ = sensor_range;
                policy_choice_ = policy;
                cost_ = cost;
                setBBox(sensor_range, igbbox);
                grid_size_ = grid_tools::getGridSizeFromBBoxAndRes(bbox_, resolution_);
                grid_tools::resizeAndInitializeGrid(stat_grid_, grid_size_, PointStat());
                grid_tools::resizeAndInitializeGrid(init_stat_grid_, grid_size_, PointStat());
                grid_tools::resizeAndInitializeGrid(mask_grid_, grid_size_, false);
                initializePointStatsCoords();

        };
        void updateObsStat(const pcl::PointCloud<PointStat>&);
        const grid_tools::BBox& getBBox() const {return bbox_;};
        const pcl::PointXYZ& getCoords() const {return coords_;};
        const std::vector<PointRay>& getPointRayList() const {return point_ray_list_;};
        //const pcl::PointCloud<PointStat>& getStatCloud() const {return stat_cloud_;};
        double getIG() const ;
        double getCost() const {return cost_;};
        void addPointRay(const PointRay&);
        //set the list of points to raycast, from the cylinder sensor_range + the occupied_points
        void setPointsToRaycast(const std::vector<PointXYZP>&);
        const std::vector<PointXYZP>& getPointsToRaycast() const {return points_to_raycast_;};
        //return the updated_stat_grid, debug function
        const grid_tools::StatGrid& getStatGrid() const {return stat_grid_;    };
        void printDebug() const;
        void clear() {
            point_ray_list_.clear();
            points_to_raycast_.clear();

        } ;
};
#endif