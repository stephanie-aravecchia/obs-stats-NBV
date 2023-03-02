
#ifndef TRAJECTORY_COSTMAP_H
#define TRAJECTORY_COSTMAP_H
#include <ros/ros.h>
#include <nav_msgs/OccupancyGrid.h>
#include <tf/tf.h>
#include <tf/transform_listener.h>
#include <vector>
#include <array>
#include <assert.h>
#include <map>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>


//Compute the Trajectory Costmap from the robot current position, in the free space
class TrajectoryCostmap {
    protected:
        ros::Subscriber og_sub_;
        ros::Publisher costmap_pub_;
        tf::TransformListener listener_;

        
        ros::NodeHandle nh_;
        std::string map_frame_;
        std::string robot_frame_;
        double cm_res_;
        double robot_radius_;
        double max_update_rate_;
        double cost_scale_factor_;
        bool debug_;
        nav_msgs::MapMetaData og_info_;
        nav_msgs::MapMetaData cm_info_;
        cv::Mat_<uint8_t> og_;
        cv::Mat_<uint8_t> cm_;
        std::vector<int8_t> cm_data_;
        //cv::Point3i og_center_;
        cv::Point og_center_;
        bool first_map_received_{false};

        //debug
        cv::Mat_<cv::Vec3b> og_rgb_, og_rgb_marked_;

        unsigned int neighbourhood_; 
        typedef std::multimap<float, cv::Point> Heap; 
        
        cv::Point P2(const cv::Point3i & P) const {return cv::Point(P.x,P.y);}

        void og_callback(const nav_msgs::OccupancyGridConstPtr msg);
        void publishCostMap();
        void updateCostMap();
        
        void run();
        // Generic test if a point is within the occupancy grid
        bool isInGrid(const cv::Point3i & p, const nav_msgs::MapMetaData& info) const;
        bool isInGrid(const cv::Point & p, const nav_msgs::MapMetaData& info) const;

    public:
        TrajectoryCostmap();

};

#endif