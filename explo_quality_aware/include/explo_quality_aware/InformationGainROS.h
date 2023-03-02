
#ifndef INFORMATION_GAIN_ROS_H
#define INFORMATION_GAIN_ROS_H
#include <ros/ros.h>
#include <nav_msgs/OccupancyGrid.h>
#include <tf/tf.h>
#include <tf/transform_listener.h>
#include <obs_stat/stat_point_type.h>
#include <visualization_msgs/MarkerArray.h>
#include <geometry_msgs/PointStamped.h>
#include "explo_quality_aware/InformationGain.h"
#include <random>

//ROS Wrapper on the InformationGain class,
//Whereas InformationGain computes the IG for all the candidates, the NBV is selected here
class Explorator {
    protected:
        ros::Subscriber cm_sub_;
        ros::Subscriber stat_sub_;
        ros::Subscriber occupancy_sub_;
        ros::Publisher pc_ig_pub_;
        ros::Publisher occupancy_pc_pub_;
        ros::Publisher best_candidate_pub_;
        ros::Publisher debug_position_pub_;
        ros::Publisher debug_ig_pub_;
        tf::TransformListener listener_;
        double resolution_;
        double search_range_;
        double sensor_range_;
        double target_rate_;
        double octomap_res_;
        double alpha_;
        double security_margin_border_;
        PolicyChoice policy_;
        explo_utils::CostMap costmap_;
        pcl::PointCloud<pcl::PointXYZ> testcloud_;
        pcl::PointCloud<PointStat> stat_cloud_;
        pcl::PointCloud<PointXYZP> occupancy_point_cloud_;
        

        grid_tools::BBox exploration_bbox_; 
        bool is_costmap_{false};
        bool is_occ_points_{false};
        bool is_stats_{false};
        
        bool debug_;
        
        ros::NodeHandle nh_;
        std::string map_frame_;
        std::string robot_frame_;

        void cm_callback(const nav_msgs::OccupancyGridConstPtr msg);
        void octomap_pts_callback(const visualization_msgs::MarkerArrayConstPtr msg);
        void stat_pc_callback(const sensor_msgs::PointCloud2ConstPtr msg);
        
        void run();
        size_t getBestCandidateIndex(const std::vector<explo_utils::CandidateResult>& candidates) const; 
        void publishBestCandidate(const std::vector<explo_utils::CandidateResult>& candidates) const; 
        bool isValidCandidate(const explo_utils::CandidateResult) const;
        
        void setPolicy(const std::string&);
        //debug functions, display stuff in RVIZ
        void publishListOfCandidatesAndRobotMarkers(const pcl::PointXYZ&, const std::vector<Candidate>&) const; 
        void publishIgPerCandidateAndRobotMarkers(const pcl::PointXYZ&, const std::vector<std::pair<double, pcl::PointXYZ> >&) const; 
        void publishListOfOccupiedPoints(const std::vector<PointXYZP>&) const; 

    public:
        Explorator();

};

#endif