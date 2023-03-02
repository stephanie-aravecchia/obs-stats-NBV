#ifndef OBSSTATGRID_ROS_H
#define OBSSTATGRID_ROS_H
#include <ros/ros.h>
#include <sensor_msgs/PointCloud2.h>
#include <geometry_msgs/Twist.h>
#include <pcl_ros/point_cloud.h>
#include <pcl_ros/transforms.h>
#include <pcl/point_types.h>
#include <tf/tf.h>
#include <tf/transform_listener.h>
#include <obs_stat/obs_stat_lib.h>
#include <fstream>

//This class takes a PointCloud2, compute statistics on the observation point of view in a voxel grid, 
//and publishes a PointCloud2 of those statistics
//This is a ROS wrapper on ObsStat
class ObsStatGridROS {
    protected:
        ros::Subscriber pc_sub_;
        ros::Subscriber cmd_sub_;
        ros::Publisher pc_stat_pub_;
        tf::TransformListener listener_;

        ros::NodeHandle nh_;
        std::string map_frame_;
        double resolution_;
        double range_max_;
        double range_min_;
        double bottom_thres_;
        grid_tools::GridSize grid_size_;
        grid_tools::BBox bbox_;
        bool is_3d_;
        bool filter_negative_z_;
        ObsStatGrid obs_stat_grid_;
        bool on_bag_only_ = false;

        void pc_callback(const sensor_msgs::PointCloud2ConstPtr msg);
        void publishStatPC();
        //utility functions to make sure the bbox starts and end with multiple of res
        double getStartingDistance(double minval) const;
    
    public:
        void processBagFile(const std::string& bagfile, const std::string& topic);
        void saveGridToDisk(const std::string& filename) const;
        void checkFile(std::ofstream & f, const std::string& s) const;


    public:
        ObsStatGridROS();

};

#endif