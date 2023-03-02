#include <algorithm>
#include <iostream>
#include "explo_quality_aware/TrajectoryCostmap.h"
#include "explo_quality_aware/Constants.h"

using namespace std;

TrajectoryCostmap::TrajectoryCostmap() : nh_("~") {        

    nh_.param<std::string>("robot_frame",robot_frame_,"/base_link");
    nh_.param<double>("robot_radius",robot_radius_,.5);
    nh_.param<double>("costmap_resolution",cm_res_,.5);
    nh_.param<double>("cost_scale_factor",cost_scale_factor_,1.);
    nh_.param<double>("max_update_rate",max_update_rate_,5.);
    nh_.param<bool>("debug",debug_,false);
    og_sub_ = nh_.subscribe("map",1,&TrajectoryCostmap::og_callback,this);
    costmap_pub_ = nh_.advertise<nav_msgs::OccupancyGrid>("trajectory_costmap",1);
    run();
}

//Convert the msg into a cv::Mat_, erode
void TrajectoryCostmap::og_callback(const nav_msgs::OccupancyGridConstPtr msg) {
    first_map_received_ = true;
    og_info_ = msg->info;
    map_frame_ = msg->header.frame_id;
    og_ = cv::Mat_<uint8_t>(msg->info.height, msg->info.width,FREE);
    og_center_ = cv::Point(-og_info_.origin.position.x/og_info_.resolution,
            -og_info_.origin.position.y/og_info_.resolution);

    for (unsigned int j=0;j<msg->info.height;j++) {
        for (unsigned int i=0;i<msg->info.width;i++) {
            int8_t v = msg->data[j*msg->info.width + i];
            switch (v) {
                case 0: 
                    og_(j,i) = FREE; 
                    break;
                case 100: 
                    og_(j,i) = OCCUPIED; 
                    break;
                case -1: 
                default:
                    og_(j,i) = UNKNOWN; 
                    break;
            }
        }
    }
    // dilatation/erosion part
    int erosion_type = cv::MORPH_RECT ;
    int erosion_size = robot_radius_/og_info_.resolution ;
    cv::Mat element = cv::getStructuringElement(erosion_type,
            cv::Size(2*erosion_size+1,2*erosion_size+1),
            cv::Point( erosion_size, erosion_size));
    cv::erode( og_, og_, element );
    //visualization
    if (debug_) {
        cv::cvtColor(og_, og_rgb_, cv::COLOR_GRAY2RGB);
        cv::imshow( "OccGrid orig", og_rgb_ );
    }
    updateCostMap();

}

bool TrajectoryCostmap::isInGrid(const cv::Point3i & p, const nav_msgs::MapMetaData& info) const {
    if ((p.x < 0) || (p.x >= (signed)info.width)
            || (p.y < 0) || (p.y >= (signed)info.height)) {
        return false;
    }
    return true;
}
bool TrajectoryCostmap::isInGrid(const cv::Point & p, const nav_msgs::MapMetaData& info) const {
    return isInGrid(cv::Point3i(p.x, p.y, 0), info);
    return true;
}

void TrajectoryCostmap::publishCostMap() {
    nav_msgs::OccupancyGrid msg;
    msg.header.stamp = ros::Time::now();
    msg.header.frame_id = map_frame_; 
    msg.info = cm_info_;
    msg.data = cm_data_;
    if (cm_data_.size() != (msg.info.width * msg.info.height)) {
        ROS_ERROR("Size mismatch: %d != (%d * %d)", (int) cm_data_.size(),msg.info.width, msg.info.height);
        return;
    } 
    if (cm_data_.size() ==0 ) {
        ROS_ERROR("Empty costmap");
        return;
    } 
    costmap_pub_.publish(msg);
}

void TrajectoryCostmap::updateCostMap() {

    if (!first_map_received_) {
        ROS_INFO("Waiting for occupancy grid");
        ros::Duration(0.5).sleep();
        return;
    }
    //get the robot pose in the map frame
    tf::StampedTransform transform;
    try { 
        listener_.lookupTransform(map_frame_,robot_frame_, ros::Time(0), transform);
    } catch (tf::TransformException &ex){
        ROS_ERROR("%s",ex.what());
        return;
    }

    //now, put the point in the new costmap
    double scale_factor = og_info_.resolution / cm_res_ ;
    cm_info_.resolution = cm_res_;
    cm_info_.width = (uint32_t) og_info_.width * scale_factor;
    cm_info_.height = (uint32_t) og_info_.height * scale_factor;
    cm_info_.origin = og_info_.origin;
    cv::Point new_center = cv::Point(og_center_.x*scale_factor, og_center_.y* scale_factor);
    cv::Point start = cv::Point(transform.getOrigin().x() / cm_info_.resolution,
                                   transform.getOrigin().y() / cm_info_.resolution)
                    + new_center;
    
    if (!isInGrid(start, cm_info_)) {
        ROS_ERROR("Invalid start point (%.2f %.2f) -> (%d %d)",
            transform.getOrigin().x(), transform.getOrigin().y(), start.x, start.y);
        return;
    }
    //we resize the og to match the new cmap resolution:
    cv::Mat_<uint8_t> og_resized;
    cv::resize(og_,og_resized, cv::Size(cm_info_.width, cm_info_.height), 0,0, cv::INTER_NEAREST);
    
    //And we make sure og_resized still contains only FREE, OCC and UNKNOWN:
    for (size_t j=0;j<static_cast<size_t>(og_resized.size().height);j++) {
        for (size_t i=0;i<static_cast<size_t>(og_resized.size().width);i++) {
            if (og_resized(j,i) < 50) {
                og_resized(j,i) = OCCUPIED;
            } else if (og_resized(j,i) > 200) { 
                og_resized(j,i) = FREE;
            } else {
                og_resized(j,i) = UNKNOWN;
            }
        }
    }

    // The starting point may not be FREE because of the resize, free it
    if (og_resized(start) != FREE) {
        og_resized(start) == FREE;
    }

    //visualization
    if (debug_){
        cv::imshow( "Og_resized", og_resized );
        cv::cvtColor(og_resized, og_rgb_marked_, cv::COLOR_GRAY2RGB);
        cv::circle(og_rgb_marked_,start, 10, cv::Scalar(0,0,255));
        cv::imshow( "OccGrid Marked", og_rgb_marked_ );
    }
    // Breadth-first algorithm 
    int dims[2] = {(int)cm_info_.width, (int)cm_info_.height};
    cv::Mat_<double> cell_value(2,dims, NAN);

    // The neighbour of a given cell in relative coordinates.
    const size_t n_neighbours = 8;
    cv::Point neighbours_0[n_neighbours] = {cv::Point(1,0), cv::Point(0,1), cv::Point(0,-1), cv::Point(-1,0),
                                 cv::Point(1,1), cv::Point(1,-1), cv::Point(-1,1), cv::Point(-1,-1),};

    // Cost of displacement corresponding the neighbours.
    double cost = 1.*cost_scale_factor_;
    double cost_0[n_neighbours] = {cost, cost, cost, cost, cost, cost, cost, cost};
    // The core of the breadth-first, a sorted heap, where the first
    // element the closer to the robot.
    Heap heap;
    heap.insert(Heap::value_type(0, start));
    cell_value(start.x,start.y) = 0;
    while (!heap.empty()) {
        // Select the cell at the top of the heap
        Heap::iterator hit = heap.begin();
        // the cell it contains is this_cell
        cv::Point this_cell = hit->second;
        // and its score is this_cost
        double this_cost = cell_value(this_cell.x,this_cell.y);
        // We can remove it from the heap now.
        heap.erase(hit);
        // Now see where we can go from this_cell
        for (unsigned int i=0;i<n_neighbours;i++) {
            cv::Point dest = this_cell + neighbours_0[i];
            if (!isInGrid(dest, cm_info_)) {
                // outside the grid
                continue;
            }
            uint8_t og = og_resized(dest);
            if (og != FREE) {
                // occupied or unknown
                continue;
            }
            double cv = cell_value(dest.x,dest.y);
            double new_cost = this_cost + cost_0[i];
            if (isnan(cv) || (new_cost < cv)) {
                cell_value(dest.x,dest.y) = new_cost;
                // And insert the selected cells in the map.
                heap.insert(Heap::value_type(new_cost,dest));
            }
        }
    }
    //Now that we have the costs, we build the cost_map cm_ (only for visualisation actually), 
    //and the cm_data_ vector that will be published
    int cm_w = cm_info_.width;
    int cm_h = cm_info_.height;
    cm_ = cv::Mat_<int8_t>(cm_h, cm_w);
    
    cm_data_.clear();
    for (size_t j=0;j<static_cast<size_t>(cm_h);j++) {
        for (size_t i=0;i<static_cast<size_t>(cm_w);i++) {
            assert(i<static_cast<size_t>(cell_value.size().height));
            assert(j<static_cast<size_t>(cell_value.size().width));
            int8_t v;
            if (isnan(cell_value(i,j))) {
                v = CMAP_UNKNOWN;
            } else {
                v = static_cast<uint8_t>(std::max(0., std::min(cell_value(i,j), (double)CMAP_MAXVALID)));
            }
            //if it is an obstacle, mark it lethal:
            if (og_resized(j,i) == OCCUPIED) {
                v = CMAP_LETHAL;
            } 
            cm_(j,i) = v;
            cm_data_.push_back(v);
        }
    }
    if (debug_) {
        cv::imshow( "CostMap", cm_ );
    }
    publishCostMap();
}

void TrajectoryCostmap::run() {
    ros::Rate rate(max_update_rate_);
    while (ros::ok()) {
        ros::spinOnce();
        if (cv::waitKey( 50 )== 'q') {
            ros::shutdown();
        }
        rate.sleep();
    }
}


int main(int argc, char *argv[]) {
    ros::init(argc, argv, "trajectory_costmap");
    TrajectoryCostmap cm;
}
