#include <algorithm>
#include <iostream>
#include "explo_quality_aware/InformationGainROS.h"
#include <random>

using namespace std;

//ROS Wrapper on the InformationGain class,
//Whereas InformationGain computes the IG for all the candidates, the NBV is selected here

//The user needs to make sure the resolution is the same than in pointstat (this is a param in both launchfiles)
Explorator::Explorator() : nh_("~") {        

    double xmin; double ymin; double zmin;
    double xmax; double ymax; double zmax;
    nh_.param<std::string>("robot_frame",robot_frame_,"base_link");
    nh_.param<std::string>("map_frame",map_frame_,"world");
    nh_.param<bool>("debug",debug_,false);
    nh_.param<double>("resolution",resolution_,.5);
    nh_.param<double>("search_range",search_range_,10.);
    nh_.param<double>("sensor_range",sensor_range_,5.);
    nh_.param<double>("target_rate",target_rate_,0.2);
    nh_.param<double>("alpha",alpha_,2.);//policy=max(Ig - alpha * cost)
    nh_.param<double>("xmin",xmin,-29.);
    nh_.param<double>("ymin",ymin,-29.);
    nh_.param<double>("zmin",zmin,-0.);
    nh_.param<double>("xmax",xmax,29.);
    nh_.param<double>("ymax",ymax,29.);
    nh_.param<double>("zmax",zmax,8.);
    nh_.param<double>("octomap_res_in_navigation_stack",octomap_res_,0.5);
    nh_.param<double>("security_margin_border",security_margin_border_,1.);
    std::string policy_choice;
    nh_.param<std::string>("policy",policy_choice,"");
    cm_sub_ = nh_.subscribe("costmap",1,&Explorator::cm_callback,this);
    stat_sub_ = nh_.subscribe("obs_stat",1,&Explorator::stat_pc_callback,this);
    occupancy_sub_ = nh_.subscribe("occupied_points",1,&Explorator::octomap_pts_callback,this);
    best_candidate_pub_ = nh_.advertise<geometry_msgs::PointStamped>("next_goal",1);
    exploration_bbox_ = grid_tools::BBox(xmin, xmax, ymin, ymax, zmin, zmax);

    setPolicy(policy_choice);
    
    if (debug_) {
        std::cout << "Resolution: " << resolution_ << ", search_range: " << search_range_ << ", sensor_range: " << sensor_range_ << std::endl;
        occupancy_pc_pub_ = nh_.advertise<sensor_msgs::PointCloud2>("debug_occupancy_pointcloud",1);
        debug_position_pub_ = nh_.advertise<visualization_msgs::MarkerArray>("debug_position",1);
        debug_ig_pub_ = nh_.advertise<visualization_msgs::MarkerArray>("debug_ig_candidates",1);

    }
    listener_.waitForTransform(map_frame_,robot_frame_, ros::Time(0), ros::Duration(3.0));
    run();
}

//Set the policy, based on what is implemented so far
void Explorator::setPolicy(const std::string& policy_choice) {
    if (policy_choice == std::string("nobs")) {
        policy_ = PolicyChoice::NOBS;
    } else if (policy_choice == std::string("sigma")) {
        policy_ = PolicyChoice::SIGMA;
    } else if (policy_choice == std::string("range")) {
        policy_ = PolicyChoice::RANGE;
    } else if (policy_choice == std::string("angles")) {
        policy_ = PolicyChoice::ANGLES;
    } else if (policy_choice == std::string("all")) {
        policy_ = PolicyChoice::ALL;
    } else {
        ROS_ERROR("Error, the choice of the policy should be in: [nobs, sigma, range, all]");
        exit(EXIT_FAILURE);
    }
}

//set costmap_ 
void Explorator::cm_callback(const nav_msgs::OccupancyGridConstPtr msg){
    int cm_w = msg->info.width;
    int cm_h = msg->info.height;
    cv::Mat_<int8_t> cm = cv::Mat_<int8_t>(cm_h, cm_w);
    if ((cm_h<=0) || (cm_w<=0)) {
        is_costmap_ = false;
        return;
    }
    for (int j=0;j<cm_h;j++) {
        for (int i=0;i<cm_w;i++) {
            assert((j*cm_w+i)<(cm_w*cm_h));
            assert(cm_h-j-1>=0);
            assert(cm_h-j-1<cm_h);
            assert(i>=0);
            assert(i<cm_w);
            cm(cm_h-j-1,i) = msg->data[j*cm_w+i];
        }
    }
    pcl::PointXYZ cm_orig_topleft = pcl::PointXYZ(msg->info.origin.position.x,msg->info.origin.position.y+ cm_h*msg->info.resolution,msg->info.origin.position.z);
    costmap_ = explo_utils::CostMap(cm, msg->info.resolution, cm_orig_topleft);
    if (!is_costmap_) is_costmap_ = true;
}

//We don't use the octomap occupied points because thw points are the center of the leaves and so the size of the leaves is not constant
//We use the marker array, containing the information
//The markers contain in point the center of each voxel, and in id the depth in the octree
void Explorator::octomap_pts_callback(const visualization_msgs::MarkerArrayConstPtr msg) {
    //iterate on the points in the MarkerArray,
    //push_back the points to occupancy_point_cloud_ with probability 0.95
    //We use a small trick: if the points are aligned to our grid, they may be ignored when we calculate the intersection of rays and grid
    //To prevent that, we add a small noise on the points
    static std::random_device rd;
    static std::default_random_engine engine(rd());
    static std::uniform_real_distribution<float> distri(-resolution_/10, resolution_/10);
    occupancy_point_cloud_.clear();
    if (msg->markers.size()>0) {
        for (size_t i = 0; i < msg->markers.size(); i++) {
            //the maxdepth in octomap is 16
            int depth = msg->markers[i].id;
            double vox_size = (17-depth) * octomap_res_;
            size_t npoints = msg->markers[i].points.size();
            if (npoints == 0) continue;
            if (vox_size <= resolution_) {
                for (size_t j = 0 ; j < npoints; j++) {
                    occupancy_point_cloud_.push_back(PointXYZP(msg->markers[i].points[j].x + distri(engine), 
                                                            msg->markers[i].points[j].y + distri(engine),
                                                            msg->markers[i].points[j].z + distri(engine),
                                                            0.95));
                }
            } else {
                //if the size of the voxel is larger than the resolution, we add a point every (resolution/2 around the center of the vox), in the marker volume
                for (size_t j = 0 ; j < npoints; j++) {
                    pcl::PointXYZ center{static_cast<float>(msg->markers[i].points[j].x), static_cast<float>(msg->markers[i].points[j].y), static_cast<float>(msg->markers[i].points[j].z)};
                    grid_tools::BBox marker_bbox{center.x - vox_size/2, center.x + vox_size/2,
                                                center.y - vox_size/2, center.y + vox_size/2,
                                                center.z - vox_size/2, center.z + vox_size/2};
                    grid_tools::GridSize insertion_grid_size = grid_tools::getGridSizeFromBBoxAndRes(marker_bbox, resolution_/2);
                    assert(insertion_grid_size.nx>1);
                    assert(insertion_grid_size.ny>1);
                    assert(insertion_grid_size.nz>1);
                    //we add a resolution/4 offset inside the bbox
                    pcl::PointXYZ min_point{static_cast<float>(marker_bbox.xmin + resolution_/4), static_cast<float>(marker_bbox.ymin + resolution_/4),static_cast<float>(marker_bbox.xmin + resolution_/4)};
                    for (size_t ix{0}; ix < insertion_grid_size.nx-1; ++ix) {
                        for (size_t iy{0}; iy < insertion_grid_size.ny-1; ++iy) {
                            for (size_t iz{0}; iz < insertion_grid_size.nz-1; ++iz) {
                                occupancy_point_cloud_.push_back(PointXYZP(min_point.x + ix*(resolution_/2) + distri(engine),
                                                                            min_point.y + iy*(resolution_/2) + distri(engine),
                                                                            min_point.z + iz*(resolution_/2) + distri(engine),
                                                                            0.95));
                            }
                        }
                    }
                }
            }
        }
    }
    if (!is_occ_points_) is_occ_points_ = true;
    if (msg->markers.size() == 0) is_occ_points_ = false;
}

//set stat_cloud_
void Explorator::stat_pc_callback(const sensor_msgs::PointCloud2ConstPtr msg){
    pcl::fromROSMsg(*msg, stat_cloud_);
    if (!is_stats_) is_stats_ = true;
    if (stat_cloud_.size()==0) is_stats_ = false;
}

//Visualization only, used only in debug mode
//This publishes the list of occupied points 
void Explorator::publishListOfOccupiedPoints(const std::vector<PointXYZP>& points) const {
        pcl::PointCloud<PointXYZP> cloud;
        cloud.resize(points.size());
        for (size_t i=0; i<points.size(); i++) {
            cloud[i] = points[i];
        }
        sensor_msgs::PointCloud2 msg;
        pcl::toROSMsg(cloud, msg);
        msg.header.stamp = ros::Time::now();
        msg.header.frame_id = map_frame_; 
        occupancy_pc_pub_.publish(msg);

}

//check that the candidate is valid, inside the bbox with a security_margin of 1m
bool Explorator::isValidCandidate(const explo_utils::CandidateResult candidate) const {
    if ((candidate.coords.x < exploration_bbox_.xmin + security_margin_border_) ||
        (candidate.coords.x > exploration_bbox_.xmax - security_margin_border_)) {
            return false;
    }
    if ((candidate.coords.y < exploration_bbox_.ymin + security_margin_border_) ||
        (candidate.coords.y > exploration_bbox_.ymax - security_margin_border_)) {
            return false;
    }
    return true;
}

//Here, we publish and choose the best candidate based on
    // the information gain
    // the cost
    // the actual policy max(Ig - alpha * cost)
size_t Explorator::getBestCandidateIndex(const std::vector<explo_utils::CandidateResult>& candidates) const {
    size_t best_index;
    std::vector<double> policy(candidates.size());
    std::vector<double> gains(candidates.size());
    std::vector<double> costs(candidates.size());
    for (size_t i=0; i< candidates.size(); i++) {
        gains[i] = candidates[i].ig; 
        costs[i] = candidates[i].cost; 
    }
    double min_gain = *std::min_element(gains.begin(), gains.end());
    double max_gain = *std::max_element(gains.begin(), gains.end());
    //if the cost is uniform, we probably are the beginning of the exploration, we select one randomly
    //Also, with a bad tuning of the information gain formulation for an element of the grid, the information gain could be somehow uniform
    //We define uniform if (max - min)/max < 5%
    if ((fabs(max_gain/min_gain) - 1) <= 0.05) {
        assert(candidates.size() > 0);
        static std::default_random_engine engine(static_cast<unsigned int>(ros::Time::now().toSec()));
        static std::uniform_int_distribution<unsigned int> distri(0, static_cast<unsigned int>(candidates.size()-1));
        size_t m = distri(engine);
        best_index = m;
        while (!isValidCandidate(candidates[m])) {
            m = distri(engine);
        }
        ROS_INFO("Random selection of the candidate. Beginning of the Exploration ?");
        if (debug_) {
            std::cout << "random selection, m is: " << m << std::endl;
        }
    } else {
        //else we select the best
        for (size_t i=0; i< candidates.size(); i++) {
            policy[i] = gains[i] - alpha_ * costs[i];
        }
        auto max = std::max_element(policy.begin(), policy.end());
        size_t argmax = distance(policy.begin(), max);
        assert(argmax<candidates.size());
        best_index = argmax;
        //if the candidate is too close to the border or outside, we select the second best
        while (!isValidCandidate(candidates[best_index])) {
            policy[best_index] = 0;
            max = std::max_element(policy.begin(), policy.end());
            argmax = distance(policy.begin(), max);
            assert(argmax<candidates.size());
            best_index = argmax;
        }
    }
    return best_index;
}
void Explorator::publishBestCandidate(const std::vector<explo_utils::CandidateResult>& candidates) const {

    size_t best_index = getBestCandidateIndex(candidates);
    explo_utils::CandidateResult best = candidates[best_index];

    geometry_msgs::PointStamped msg;
    msg.header.frame_id = map_frame_;
    msg.header.stamp = ros::Time::now();
    msg.point.x = best.coords.x;
    msg.point.y = best.coords.y;
    msg.point.z = best.coords.z;
    best_candidate_pub_.publish(msg);
    if (debug_) {
        std::cout << "best is: " << best.coords.x << ", " << best.coords.y << std::endl;
    }
}
//visualization in debug only
void Explorator::publishListOfCandidatesAndRobotMarkers(const pcl::PointXYZ& robot, const std::vector<Candidate>& candidates) const {
    visualization_msgs::MarkerArray msg;
    ros::Time now = ros::Time::now();
    size_t n = candidates.size() + 1;
    msg.markers.resize(n);
    std::cout << "number of markers: " <<  n << endl;
    for (size_t i =0; i < n; i++) {
        msg.markers[i].header.frame_id = map_frame_;
        msg.markers[i].header.stamp = now;
        msg.markers[i].ns = "candidates";
        msg.markers[i].id = i;
        msg.markers[i].action = visualization_msgs::Marker::ADD;
        msg.markers[i].type = visualization_msgs::Marker::CUBE;
        if (i == 0) {
            msg.markers[i].pose.position.x = robot.x;
            msg.markers[i].pose.position.y = robot.y;
            msg.markers[i].pose.position.z = robot.z;
        } else {
            msg.markers[i].pose.position.x = candidates[i-1].getCoords().x;
            msg.markers[i].pose.position.y = candidates[i-1].getCoords().y;
            msg.markers[i].pose.position.z = candidates[i-1].getCoords().z;

        }
        msg.markers[i].pose.orientation.x = 0.0;
        msg.markers[i].pose.orientation.y = 0.0;
        msg.markers[i].pose.orientation.z = 0.0;
        msg.markers[i].pose.orientation.w = 1.0;
        msg.markers[i].scale.x = 0.5;
        msg.markers[i].scale.y = 0.5;
        msg.markers[i].scale.z = 0.5;
        msg.markers[i].color.a = 1.0;
        if (i==0) {
            msg.markers[i].color.r = 1.0;
            msg.markers[i].color.g = 0.0;
            msg.markers[i].color.b = 0.0;
        } else {
            msg.markers[i].color.r = 0.0;
            msg.markers[i].color.g = 0.0;
            msg.markers[i].color.b = 1.0;
        }
    }
    debug_position_pub_.publish(msg);

}
//visualization if debug only
void Explorator::publishIgPerCandidateAndRobotMarkers(const pcl::PointXYZ& robot, const std::vector<std::pair<double, pcl::PointXYZ> >& candidates) const {
    std::vector<double> ig_vec;
    for (auto pair : candidates) {
        ig_vec.push_back(pair.first);
    }
    double max_gain; double min_gain;
    max_gain = *std::max_element(ig_vec.begin(), ig_vec.end());
    min_gain = *std::min_element(ig_vec.begin(), ig_vec.end());
    double norm_val; 
    visualization_msgs::MarkerArray msg;
    ros::Time now = ros::Time::now();
    size_t n = candidates.size() + 1;
    msg.markers.resize(n);
    for (size_t i =0; i < n; i++) {
        msg.markers[i].header.frame_id = map_frame_;
        msg.markers[i].header.stamp = now;
        msg.markers[i].ns = "candidates_ig";
        msg.markers[i].id = i;
        msg.markers[i].action = visualization_msgs::Marker::ADD;
        msg.markers[i].type = visualization_msgs::Marker::CUBE;
        if (i == 0) {
            msg.markers[i].pose.position.x = robot.x;
            msg.markers[i].pose.position.y = robot.y;
            msg.markers[i].pose.position.z = robot.z;
        } else {
            msg.markers[i].pose.position.x = candidates[i-1].second.x;
            msg.markers[i].pose.position.y = candidates[i-1].second.y;
            msg.markers[i].pose.position.z = candidates[i-1].second.z;

        }
        msg.markers[i].pose.orientation.x = 0.0;
        msg.markers[i].pose.orientation.y = 0.0;
        msg.markers[i].pose.orientation.z = 0.0;
        msg.markers[i].pose.orientation.w = 1.0;
        msg.markers[i].scale.x = 0.5;
        msg.markers[i].scale.y = 0.5;
        msg.markers[i].scale.z = 0.5;
        if (i==0) {
            msg.markers[i].color.r = 1.0;
            msg.markers[i].color.g = 0.0;
            msg.markers[i].color.b = 0.0;
            msg.markers[i].color.a = 1.0;
        } else {
            msg.markers[i].color.r = 0.0;
            msg.markers[i].color.g = 1.0;
            msg.markers[i].color.b = 0.0;
            if (max_gain - min_gain < 1e-30) {
                msg.markers[i].color.a = 1.;
            } else {
                norm_val = (candidates[i-1].first - min_gain) / (max_gain - min_gain);
                msg.markers[i].color.a = min(max(0.2, norm_val),1.);
            }
            if ((msg.markers[i].color.a <0) || (msg.markers[i].color.a >1)){
                std::cout << msg.markers[i].color.a << " ";
            }
        }
    }
    debug_ig_pub_.publish(msg);

}

void Explorator::run() {
    ros::Rate rate(target_rate_);
   
    while (ros::ok()) {

        //We wait for the costmap, the occupancy point cloud and the statistics on the viewpoints to be available
        //before we calculate the information gain
        if ((is_costmap_ && is_occ_points_) && is_stats_) {

            //get the robot pose in the map frame
            tf::StampedTransform transform;
            try { 
                listener_.lookupTransform(map_frame_,robot_frame_, ros::Time(0), transform);
            } catch (tf::TransformException &ex){
                ROS_ERROR("%s",ex.what());
                continue;
            }
            pcl::PointXYZ robot_position(transform.getOrigin().x(), transform.getOrigin().y(), transform.getOrigin().z());
            
            if (debug_) {
                std::cout << "robot_position: " << robot_position << std::endl;
            }

            InformationGain ig(robot_position, resolution_, search_range_, sensor_range_, 
                                                costmap_, occupancy_point_cloud_, stat_cloud_, policy_, exploration_bbox_, debug_);
            
            const std::vector<explo_utils::CandidateResult>& candidates = ig.getResultsPerCandidate();
            
            if (debug_) {
                //In debug mode: 
                //publish the list of occupied points, should match the octomap markers in a bbox around the robot
                //publish the robot position as a red marker
                //publish the candidates as a markers, the highest the IG, the less transparent it is
                const std::vector<PointXYZP>& list_of_occupied_points = ig.getListOfOccupiedPoints();
                const std::vector<std::pair<double, pcl::PointXYZ> >& ig_per_candidate = ig.getIGPerCandidate();
                std::cout << "number of candidates: " << candidates.size()
                          << ", number of occupied points: " << list_of_occupied_points.size() << std::endl;
                pcl::PointXYZ robot = ig.getRobotPosition();
                std::vector<double> ig_vec;
                for (auto pair : ig_per_candidate) {
                    ig_vec.push_back(pair.first);
                }
                double max_gain; double min_gain;
                max_gain = *std::max_element(ig_vec.begin(), ig_vec.end());
                min_gain = *std::min_element(ig_vec.begin(), ig_vec.end());
                std::cout << "max gain: " << std::scientific << max_gain << ", min gain: " << std::scientific << min_gain << std::fixed << std::endl;
                publishIgPerCandidateAndRobotMarkers(robot, ig_per_candidate);
                publishListOfOccupiedPoints(list_of_occupied_points);
            }

            //Now we can publish the best candidate
            publishBestCandidate(candidates);
        }
        rate.sleep();
        ros::spinOnce();
    }
}


int main(int argc, char *argv[]) {
    ros::init(argc, argv, "explorator");
    Explorator ig;
}
