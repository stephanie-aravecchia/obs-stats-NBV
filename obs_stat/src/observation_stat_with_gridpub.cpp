#include <geometry_msgs/PointStamped.h>
#include <std_msgs/Bool.h>
#include "obs_stat/observation_stat_with_gridpub.h"
#include <random>
#include <rosbag/bag.h>
#include <rosbag/view.h>

using namespace std;
using namespace obs_stat_tools;
using namespace grid_tools;

//ROS Wrapper on ObsStat

ObsStatGridROS::ObsStatGridROS() : nh_("~") {        

    nh_.param<std::string>("map_frame",map_frame_,"/map");
    nh_.param<bool>("filter_negative_z",filter_negative_z_,false);
    nh_.param<double>("range_max",range_max_,20);
    nh_.param<bool>("is_3d",is_3d_,false);
    double xmin; double xmax; double ymin; double ymax; double zmin; double zmax;
    nh_.param<double>("grid_xmin",xmin,0);
    nh_.param<double>("grid_ymin",ymin,0);
    nh_.param<double>("grid_zmin",zmin,0);
    nh_.param<double>("resolution_m_per_pix",resolution_,double(0.1));
    nh_.param<double>("grid_xmax",xmax,0);
    nh_.param<double>("grid_ymax",ymax,0);
    nh_.param<double>("grid_zmax",zmax,0);
    nh_.param<double>("bottom_thres",bottom_thres_,0.2);
    // From the user input, we change the bbox coordinates to start and end with a multiple of res
    // This is mandatory when we compare the statistics of observations to the reconstruction
    // If we are not interested in that step, we can remove that
    xmin = xmin + getStartingDistance(xmin);
    ymin = ymin + getStartingDistance(ymin);
    zmin = zmin + getStartingDistance(zmin);
    size_t grid_nx = floor((xmax - xmin)/resolution_);
    size_t grid_ny = floor((ymax - ymin)/resolution_);
    size_t grid_nz = floor((zmax - zmin)/resolution_);
    xmax = xmin + grid_nx * resolution_;
    ymax = ymin + grid_ny * resolution_;
    zmax = zmin + grid_nz * resolution_;

    grid_size_ = GridSize(grid_nx, grid_ny, grid_nz);
    bbox_ = BBox(xmin, xmax, ymin, ymax, zmin, zmax);
    obs_stat_grid_ = ObsStatGrid(bbox_, resolution_);
    assert(grid_size_ == obs_stat_grid_.getGridSize());

    ROS_INFO_STREAM("grid bbox set to: " << bbox_);
    ROS_INFO_STREAM("grid size set to: " << grid_size_);

    //We set a range_min, we don't want to calculate statistics for observations in the voxel the robot currently is
    //We set that range_min to the diagonal of the voxel
    range_min_ = sqrt(3)*resolution_;
    std::string bag_list; 
    nh_.param<std::string>("bag_list",bag_list,"");
    //If we provide a bag name, process all the messages on the bag and dump to disk.
    //Else run realtime and publish the stats
    if (bag_list == "") {
        ros::Duration(.5).sleep();
        pc_sub_ = nh_.subscribe("input_scan",10,&ObsStatGridROS::pc_callback,this);
        pc_stat_pub_ = nh_.advertise<sensor_msgs::PointCloud2>("stats_pointcloud",1);
    } else {
        on_bag_only_ = true;
        ROS_WARN_STREAM("Running obs_stats on bags " << bag_list);
    }
}


void ObsStatGridROS::processBagFile(const std::string& bagfile, const std::string& topic) {
    rosbag::Bag bag;
    bag.open(bagfile);  // BagMode is Read by default

    for(rosbag::MessageInstance const m: rosbag::View(bag))
    {
      sensor_msgs::PointCloud2::ConstPtr pc = m.instantiate<sensor_msgs::PointCloud2>();
      if ((pc != nullptr) && ((m.getTopic()==topic) || topic.empty())) {
          pc_callback(pc);
          continue;
      }
      tf2_msgs::TFMessage::ConstPtr tf = m.instantiate<tf2_msgs::TFMessage>();
      if (tf != nullptr) {
          for(size_t i=0;i<tf->transforms.size();i++) {
              tf::StampedTransform stf;
              tf::transformStampedMsgToTF(tf->transforms[i],stf);
              listener_.setTransform(stf);
          }
          continue;
      }
    }
    bag.close();
}

double ObsStatGridROS::getStartingDistance(double minval) const {
    return ceil(minval / resolution_) * resolution_ - minval;
}

//What follows supposes the pointcloud we subscribe to is published only when the robot moves
//Else, the code work, but the statistics are less meaningfull: 
//i.e the number of observations will increase rapidly
void ObsStatGridROS::pc_callback(const sensor_msgs::PointCloud2ConstPtr msg) {
    //First, we rest the mask. Each voxel in the map should be updated only once
    obs_stat_grid_.resetMaskGrid();
    pcl::PointCloud<pcl::PointXYZ> cloud;
    pcl::fromROSMsg(*msg, cloud);
    tf::StampedTransform transform_to_map;
    try{
        listener_.lookupTransform(map_frame_,msg->header.frame_id,ros::Time(0),transform_to_map);
    } catch (tf::TransformException &ex){
        ROS_ERROR("%s",ex.what());
        return;
    }
    //we need the position of the laser in map
    tf::Vector3 laser_in_map_tfv(0, 0, 0); 
    laser_in_map_tfv = transform_to_map * laser_in_map_tfv;
    pcl::PointXYZ laser_in_map(laser_in_map_tfv.x(), laser_in_map_tfv.y(), laser_in_map_tfv.z());

    //In simulation, with the gaussian noise on the robot, the laser can be outside the map (specially on z), then we cannot calculate obs_stat
    //If the robot is not in the map, we do not update
    if (!isInBBox(laser_in_map, bbox_)) {
        ROS_INFO("Laser is not in bbox, skipping");
        return;
    }

    //Since we update the stats of a voxel only with the first encountered point inside
    //We shuffle it to avoid taking always the same
    static std::random_device rd;
    static std::default_random_engine rng(rd());
    std::shuffle(cloud.points.begin(), cloud.points.end(), rng);
    
    //Now, we loop on the points
    for (auto it = cloud.points.begin(); it != cloud.points.end(); ++it) {
        //First, we discarded points based on the user choices, on the z value
        if (filter_negative_z_ && (it->z < 0)) {
            continue;
        }
        if ((it->z < (bbox_.zmin + bottom_thres_))) {
            continue;
        }
        //we get the coordinate of the point in the map frame, to compute the observed voxel
        tf::Vector3 point_in_map_tfv(it->x, it->y, it->z); 
        point_in_map_tfv = transform_to_map * point_in_map_tfv;
        pcl::PointXYZ point_in_map(point_in_map_tfv.x(), point_in_map_tfv.y(), point_in_map_tfv.z());
        //if point is not the user input bbox, we discard:
        if (!isInBBox(point_in_map, bbox_)) {
            continue;
        }
        
        //If the point is in the grid, we want to compute the statistics on the voxel it is in
        //First, we get the coordinates of the center of the voxel the point is in
        cv::Point3i vox = metersToVox(point_in_map, bbox_, resolution_);
        pcl::PointXYZ vox_center = voxToMetersCenter(vox, bbox_, resolution_);
        
        //Then, we get the coordinates of the laser in the point frame, to compute the observation viewpoint curr_obs
        pcl::PointXYZ laser_in_point(laser_in_map.x - vox_center.x, laser_in_map.y - vox_center.y, laser_in_map.z - vox_center.z); 
        SphericalPoint curr_obs(cart2spherical(laser_in_point));
        //we discard observations from too far away, or potentially in the voxel
        if (curr_obs.r >= range_max_) {
            continue;
        }
        if (curr_obs.r < range_min_) {
            continue;
        }
        //Now we can compute statisatics on all the voxels between laser_in_map and point_in_map
        std::stringstream warn_stream;
        obs_stat_grid_.addObsStatRayToGrid(point_in_map, laser_in_map, curr_obs, warn_stream);
        //We update NHits for the voxel where the laser returns at least a point
        //(later, we could use that to infer if a voxel is probably empty or occupied)
        obs_stat_grid_.updateNHits(point_in_map);
        if (warn_stream.str().size()>0) {
            ROS_WARN_STREAM(warn_stream.str());
        }
    }
    if (!on_bag_only_){ 
        publishStatPC();
    }
}

void ObsStatGridROS::publishStatPC() {
    pcl::PointCloud<PointStat> cloud;
    cloud.reserve(grid_size_.nx*grid_size_.ny*grid_size_.nz);
    const StatGrid& stat_grid = obs_stat_grid_.getStatGrid();
    for (size_t ix{0}; ix < grid_size_.nx; ix++) {
        for (size_t iy{0}; iy < grid_size_.ny; iy++) {
            for (size_t iz{0}; iz < grid_size_.nz; iz++) {
                cloud.push_back(stat_grid[ix][iy][iz]);
            }
        }
    }
    sensor_msgs::PointCloud2 msg;
    pcl::toROSMsg(cloud, msg);
    msg.header.stamp = ros::Time::now();
    msg.header.frame_id = map_frame_; 
    pc_stat_pub_.publish(msg);
}
void ObsStatGridROS::checkFile(ofstream & f, const string& s) const {
    if (!f) {
        cerr << "Output file cannot be opened: " << s << endl;
        exit(EXIT_FAILURE);
    }
}
void ObsStatGridROS::saveGridToDisk(const std::string& filename) const {
    //then, we save the actual observations
    ROS_INFO("Writing to disk is starting, please wait, it can take a while.");
    ofstream out_file_stat{filename, ios::out};
    checkFile(out_file_stat, filename);
    out_file_stat << "x y z theta_mean phi_mean r_min r_max r_mean sigma_r sigma_angle sigma_angle2d angular_viewpoints_2d n_sectors n_hits n_obs";
    out_file_stat << endl;
    out_file_stat << setprecision(4);
    //we loop on every voxel of our grid and, if there is at least an observation, we store it to the file
    const StatGrid& stat_grid = obs_stat_grid_.getStatGrid();
    for (size_t ix{0}; ix < grid_size_.nx; ix++) {
        for (size_t iy{0}; iy < grid_size_.ny; iy++) {
            for (size_t iz{0}; iz < grid_size_.nz; iz++) {
                if (stat_grid[ix][iy][iz].n_obs > 0) {
                    out_file_stat   << stat_grid[ix][iy][iz].x
                                    << " " << stat_grid[ix][iy][iz].y
                                    << " " << stat_grid[ix][iy][iz].z
                                    << " " << stat_grid[ix][iy][iz].theta_mean
                                    << " " << stat_grid[ix][iy][iz].phi_mean
                                    << " " << stat_grid[ix][iy][iz].r_min
                                    << " " << stat_grid[ix][iy][iz].r_max
                                    << " " << stat_grid[ix][iy][iz].r_mean
                                    << " " << stat_grid[ix][iy][iz].sigma_r
                                    << " " << stat_grid[ix][iy][iz].sigma_angle
                                    << " " << stat_grid[ix][iy][iz].sigma_angle2d
                                    << " " << stat_grid[ix][iy][iz].angular_viewpoints_2d
                                    << " " << grid_tools::getNAnglesFromV(stat_grid[ix][iy][iz].angular_viewpoints_2d)
                                    << " " << stat_grid[ix][iy][iz].n_hits
                                    << " " << stat_grid[ix][iy][iz].n_obs
                                    << endl;
                }
            }
        }
    }
    out_file_stat.close();

}


int main(int argc, char *argv[]) {
    ros::init(argc, argv, "observation_stat");
    ObsStatGridROS obs;
    std::string bag_list(""), topic_name(""), output_file(""); 
    const ros::NodeHandle private_nh("~");
    private_nh.getParam("bag_list", bag_list);
    private_nh.getParam("topic_name", topic_name);
    private_nh.getParam("output_file", output_file);
    //If we provide a bag name, process all the messages on the bag and dump to disk.
    //Else run realtime and publish the stats
    if (bag_list != "") {
        ifstream f(bag_list, ios::in);
        if (!f) {
            ROS_ERROR_STREAM("Output file cannot be opened: " << bag_list);
            return 1;
        }
        string bag_name;
        std::vector<string> nameVec;
        while (f >> bag_name) {
            nameVec.push_back(bag_name);
        }
        f.close();
        for (const std::string& name: nameVec) {
            ROS_WARN_STREAM("bag: " << name);
            obs.processBagFile(name, topic_name);
        }
        obs.saveGridToDisk(output_file);
        ROS_WARN_STREAM("Writing complete, sleeping for 10 seconds");
        ros::Duration(10.).sleep(); 
        ROS_WARN_STREAM("Done, exiting");
        return 0;
    } else {
        ros::spin();
    }
    return 0;
}
