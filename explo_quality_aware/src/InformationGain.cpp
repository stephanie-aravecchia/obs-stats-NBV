#include <algorithm>
#include <iostream>
#include <random>
#include "explo_quality_aware/InformationGain.h"

using namespace std;

//We calculate the information gain in a BBox around the current robot_position, with a target resolution
//The size of the BBox depends on the range in which we will select candidates and the sensor range (around the candidates)
InformationGain::InformationGain(const pcl::PointXYZ& p, double resolution, double search_range, double sensor_range, const explo_utils::CostMap& costmap, 
    const pcl::PointCloud<PointXYZP>& pc, const pcl::PointCloud<PointStat>& stat_cloud, PolicyChoice policy, const grid_tools::BBox& valid_box, bool debug)
    : robot_position_(p), resolution_(resolution), search_range_(search_range), sensor_range_(sensor_range), policy_(policy), debug_(debug)
    {
        point_stat_cloud_ = stat_cloud;
        setBBoxFromRangesAndValid(valid_box);
        assert(igBBox_.xmax > (igBBox_.xmin + resolution_));
        assert(igBBox_.ymax > (igBBox_.ymin + resolution_));
        assert(igBBox_.zmax > (igBBox_.zmin + resolution_));
        setListOfOccupiedPointsFromPointCloud(pc);
        setVoxelGridFromPointCloud(igBBox_, pc);
        setListOfReachableCandidatesFromCostmap(p, costmap, search_range_);
        updateStatsInListOfCandidates();
}

void InformationGain::setVoxelGridFromPointCloud(const grid_tools::BBox & bbox,const pcl::PointCloud<PointXYZP>& pc) {
    grid_tools::GridSize grid_size = grid_tools::getGridSizeFromBBoxAndRes(igBBox_, resolution_);
    assert(grid_size.nx>0);
    assert(grid_size.ny>0);
    assert(grid_size.nz>0);
    resizeAndInitializeGrid(occupancy_vox_grid_, grid_size, 0.);
    assert(occupancy_vox_grid_.size()>0);
    for (const PointXYZP& point : list_of_occupied_points_) {
        cv::Point3i vox = grid_tools::metersToVox(explo_utils::point3D(point), igBBox_, resolution_);
        assert(grid_tools::isVoxInGrid(vox, occupancy_vox_grid_));
        occupancy_vox_grid_[vox.x][vox.y][vox.z] = point.proba;
    }
}
void InformationGain::setBBoxFromRanges() {
    igBBox_.xmin = robot_position_.x - search_range_ - sensor_range_;
    igBBox_.ymin = robot_position_.y - search_range_ - sensor_range_;
    igBBox_.zmin = robot_position_.z - sensor_range_;
    igBBox_.xmax = robot_position_.x + search_range_ + sensor_range_;
    igBBox_.ymax = robot_position_.y + search_range_ + sensor_range_;
    igBBox_.zmax = robot_position_.z + sensor_range_;
}     
void InformationGain::setBBoxFromRangesAndValid(const grid_tools::BBox& valid) {
    setBBoxFromRanges();
    igBBox_.xmin = max(igBBox_.xmin, valid.xmin);
    igBBox_.ymin = max(igBBox_.ymin, valid.ymin);
    igBBox_.zmin = max(igBBox_.zmin, valid.zmin);
    igBBox_.xmax = min(igBBox_.xmax, valid.xmax);
    igBBox_.ymax = min(igBBox_.ymax, valid.ymax);
    igBBox_.zmax = min(igBBox_.zmax, valid.zmax);
}     

void InformationGain::setListOfOccupiedPointsFromPointCloud(const pcl::PointCloud<PointXYZP>& pc) {
    list_of_occupied_points_.clear();
    size_t counter = 0;
    for (const PointXYZP& point : pc) {
        if (grid_tools::isInBBox(explo_utils::point3D(point), igBBox_)) {
            if (point.proba >= occupancy_threshold_) {
                ++counter;
            }
        }
    }
    list_of_occupied_points_.resize(counter);
    size_t index =0;
    for (const PointXYZP& point : pc) {
        if (grid_tools::isInBBox(explo_utils::point3D(point), igBBox_)) {
            if (point.proba >= occupancy_threshold_) {
                list_of_occupied_points_[index] = point;
                index++;
            }
        }
    }
} 
//set the list of candidates: every free point in the costmap within search_range
void InformationGain::setListOfReachableCandidatesFromCostmap(const pcl::PointXYZ& p, const explo_utils::CostMap& costmap, double search_range){
    double xmin = max(costmap.bbox.xmin, igBBox_.xmin);
    double ymin = max(costmap.bbox.ymin, igBBox_.ymin);
    double xmax = min(costmap.bbox.xmax, igBBox_.xmax);
    double ymax = min(costmap.bbox.ymax, igBBox_.ymax);
    //Min and max pixels in top left origin coords
    int imin = static_cast<int>(max(0., (xmin-costmap.bbox.xmin)/costmap.resolution));
    int jmin = static_cast<int>(max(0., (costmap.bbox.ymax - ymax)/costmap.resolution));
    int imax = static_cast<int>(min(imin + (xmax-xmin)/costmap.resolution, (double)costmap.width));
    int jmax = static_cast<int>(min(jmin + (ymax-ymin)/costmap.resolution, (double)costmap.height));
    assert(imin<costmap.width);
    assert(jmin<costmap.height);
    assert(imin<imax);
    assert(jmin<jmax);
    cv::Point2i pix_cm;
    cv::Point2d point_cm;
    cv::Point2d robot = cv::Point2d(p.x,p.y);
    cv::Point2i robot_pix = explo_utils::metersToPix(robot, costmap);
	{
	    lock_guard<mutex> lock(candidate_mutex_);
	    
	    for (int i=imin; i<imax; i++){
		for (int j=jmin; j<jmax; j++){
		    //We want the coordinates of the pixels in the map to be the coordinate of the bottom-left corner of the pixel in the world
		    //We compare the coordinates in the world of pixel (bottom left) and robot
		    pix_cm = cv::Point2i(i,j);
		    point_cm = explo_utils::pixToMeters(pix_cm, costmap);
		    if (pix_cm == robot_pix) {
			continue;
		    }
		    if (grid_tools::distance(point_cm, robot) <= search_range) {
			assert(i>=0);
			assert(j>=0);
			assert(i<costmap.width);
			assert(j<costmap.height);
			//this still keep candidate in invalid space
			//to get rid of them, we keep candidates with all voxels free around
			if ((costmap.cm(pix_cm) < CMAP_MAXVALID) && (costmap.cm(pix_cm) > CMAP_UNKNOWN)) {
			    bool is_around_free = true;
			    for (int a=-1; a<2; a++) {
				for (int b=-1; b<2; b++) {
				    cv::Point2i q = pix_cm + cv::Point(a,b);
				    if ((q.x<0) || (q.x>=costmap.width)) {
					continue;
				    }
				    if ((q.y<0) || (q.y>=costmap.height)) {
					continue;
				    }
				    if (costmap.cm(q) >= CMAP_MAXVALID) {
					is_around_free = false;
					break;
				    } 
				}
			    }
			    if (is_around_free) {
				double cost = static_cast<double>(costmap.cm(pix_cm))/CMAP_MAXVALID;
				    list_of_candidates_.push_back(Candidate(pcl::PointXYZ(point_cm.x,point_cm.y,0),resolution_,sensor_range_, igBBox_, robot_position_, policy_, cost));
				
			    }
			}
		    }
		}
	    }
	}
    if (debug_) {
        std::cout << "robot " << robot << " robot_pix " << robot_pix << std::endl;
    }
}

void InformationGain::getGridSegmentIntersect(const pcl::PointXYZ& orig_point, const pcl::PointXYZ& target_point, 
        InformationGain::GridIntersectMap& grid_intersect) const {
    //we define the direction vector and the origin of the segment
    //we look for the voxel contaning the start and end of the segment
    cv::Point3i vox_laser = grid_tools::metersToVox(target_point, igBBox_, resolution_);
    cv::Point3i vox_point = grid_tools::metersToVox(orig_point, igBBox_, resolution_);
    //we get the coordinates of the voxel in the map
    pcl::PointXYZ start_vox = grid_tools::voxToMeters(vox_point, igBBox_, resolution_); 
    pcl::PointXYZ end_vox = grid_tools::voxToMeters(vox_laser, igBBox_, resolution_);
    //we get the direction vector from point to laser, that gives the direction of the line
    pcl::PointXYZ laser_in_point(target_point.x - orig_point.x, target_point.y - orig_point.y, target_point.z - orig_point.z); 
    obs_stat_tools::SphericalPoint obs{obs_stat_tools::cart2spherical(laser_in_point.x, laser_in_point.y, laser_in_point.z)};
    pcl::PointXYZ uVec = obs_stat_tools::spherical2unitVector(obs);
    //we define the functions relative to the equation of the line
    double a = uVec.x;
    double b = uVec.y; 
    double c = uVec.z;
    if ((fabs(a)<MIN_POSITIVE_VALUE) || (fabs(b)<MIN_POSITIVE_VALUE) || (fabs(c)<MIN_POSITIVE_VALUE)) {
        cout << "one is too small: " << a << ", " << b << ", " << c << ", target_point: " << target_point << ", orig_point: " << orig_point <<endl;
        return;
    }
    double x1 = orig_point.x;
    double y1 = orig_point.y;
    double z1 = orig_point.z;
    auto fy_x = [a,b,c,x1,y1,z1](double x) {return ((x-x1)/a)*b + y1;};
    auto fz_x = [a,b,c,x1,y1,z1](double x) {return ((x-x1)/a)*c + z1;};
    auto fx_y = [a,b,c,x1,y1,z1](double y) {return ((y-y1)/b)*a + x1;};
    auto fz_y = [a,b,c,x1,y1,z1](double y) {return ((y-y1)/b)*c + z1;};
    auto fx_z = [a,b,c,x1,y1,z1](double z) {return ((z-z1)/c)*a + x1;};
    auto fy_z = [a,b,c,x1,y1,z1](double z) {return ((z-z1)/c)*b + y1;};
    //for each point on the segment that intersect axis x, we calculate coordinates y,z to get the coordinates of the voxel
    //we store the pair (voxel, occupancy probability) in a map to avoid duplicate
    double x_increment = (start_vox.x< end_vox.x) ? resolution_ : -resolution_;
    double y_increment = (start_vox.y< end_vox.y) ? resolution_ : -resolution_;
    double z_increment = (start_vox.z< end_vox.z) ? resolution_ : -resolution_;
    //what is the greater dimension
    double x_diff = fabs(start_vox.x - end_vox.x);
    double y_diff = fabs(start_vox.y - end_vox.y);
    double z_diff = fabs(start_vox.z - end_vox.z);
    double max_diff = std::max({x_diff, y_diff, z_diff});
    if (fabs(max_diff) < resolution_) {
        return;
    }
    if (x_diff == max_diff) {
        for (double x = start_vox.x; ((x - start_vox.x)*x_increment) < ((end_vox.x - start_vox.x)*x_increment);x+=x_increment){
            double y = fy_x(x);
            double z = fz_x(x);
            if (!grid_tools::isInBBox(pcl::PointXYZ(x,y,z), igBBox_)) {
                continue;
            }
            cv::Point3i vox = grid_tools::metersToVox(pcl::PointXYZ(x,y,z), igBBox_, resolution_);
            if (!grid_tools::isVoxInGrid(vox, occupancy_vox_grid_)) {
                continue;
            }
            double proba = occupancy_vox_grid_[vox.x][vox.y][vox.z];
            grid_intersect.insert(std::make_pair(vox, std::make_pair(proba, obs)));
        }
    } else if (y_diff == max_diff) {
        //same for intersection with y axis
        for (double y = start_vox.y; ((y - start_vox.y)*y_increment) < ((end_vox.y - start_vox.y)*y_increment);y+=y_increment){
            double x = fx_y(y);
            double z = fz_y(y);
            if (!grid_tools::isInBBox(pcl::PointXYZ(x,y,z), igBBox_)) {
                continue;
            }
            cv::Point3i vox = grid_tools::metersToVox(pcl::PointXYZ(x,y,z), igBBox_, resolution_);
            if (!grid_tools::isVoxInGrid(vox, occupancy_vox_grid_)) {
                continue;
            }
            double proba = occupancy_vox_grid_[vox.x][vox.y][vox.z];
            grid_intersect.insert(std::make_pair(vox, std::make_pair(proba,obs)));
        }
    } else {
        //same with intersection with z axis
        for (double z = start_vox.z; ((z - start_vox.z)*z_increment) < ((end_vox.z - start_vox.z)*z_increment);z+=z_increment){
            double x = fx_z(z);
            double y = fy_z(z);
            if (!grid_tools::isInBBox(pcl::PointXYZ(x,y,z), igBBox_)) {
                continue;
            }
            cv::Point3i vox = grid_tools::metersToVox(pcl::PointXYZ(x,y,z), igBBox_, resolution_);
            if (!grid_tools::isVoxInGrid(vox, occupancy_vox_grid_)) {
                continue;
            }
            double proba = occupancy_vox_grid_[vox.x][vox.y][vox.z];
            grid_intersect.insert(std::make_pair(vox, std::make_pair(proba,obs)));
        }
    }
}

void InformationGain::getPointRayFromPointToCandidate(PointRay& point_ray, const pcl::PointXYZ& point, const pcl::PointXYZ& candidate_coord) const {

    GridIntersectMap grid_intersect(point3dLess);
    //GridSegmentIntersect does not work if the segment is aligned with one of the axis
    //To avoid that, point is the center of the voxel, candidate also, but with z=0
    //And finally, we add a small random value on the coordinates to avoid that behaviour
    static std::random_device rd;
    static std::default_random_engine engine(rd());
    static std::uniform_real_distribution<float> distri(-resolution_/4, resolution_/4);
    cv::Point3i point_vox = grid_tools::metersToVox(point, igBBox_,resolution_);
    pcl::PointXYZ point_center_vox = grid_tools::voxToMetersCenter(point_vox, igBBox_, resolution_);
    cv::Point3i candidate_vox = grid_tools::metersToVox(candidate_coord, igBBox_,resolution_);
    pcl::PointXYZ candidate_center_vox = grid_tools::voxToMetersCenter(candidate_vox, igBBox_, resolution_);
    
    pcl::PointXYZ npoint{point_center_vox.x+distri(engine), point_center_vox.y+distri(engine), point_center_vox.z+distri(engine)};
    pcl::PointXYZ ncandidate{candidate_center_vox.x+distri(engine), candidate_center_vox.y+distri(engine), candidate_center_vox.z+distri(engine)};
    //If the distance is less than the diagonal of the cube, we discard
    if (grid_tools::distance(npoint, ncandidate) < sqrt(3) * resolution_) {
        return;
    }
    getGridSegmentIntersect(npoint, ncandidate, grid_intersect);
    //now we have the intersection grid / segment = the ray, we store it into the point_ray
    point_ray.setCoords(point);
    PointRay::RayProbaVector ray_v;
    if (grid_intersect.size() == 0) {
        std::cout << "grid intersect is empty. point: " << npoint << "candidate: " << ncandidate << std::endl;
        return;
    }
    for (auto pair : grid_intersect) {
        //The coordinates are in grid_vox, we want them in space
        pcl::PointXYZ p = grid_tools::voxToMeters(pair.first, igBBox_, resolution_);
        ray_v.push_back(std::make_pair(PointXYZP(p.x, p.y, p.z, pair.second.first),pair.second.second));
    }
    point_ray.setRayProbaVector(ray_v);
}

void InformationGain::processListOfCandidates() {
    
    while (true) {

        Candidate candidate;
        //mutex
        {
            lock_guard<mutex> lock(candidate_mutex_);
	    if (list_of_candidates_.empty()) {
		break;
	    }
	    candidate = list_of_candidates_.front();
            list_of_candidates_.pop_front();
        }
        
        candidate.setPointsToRaycast(list_of_occupied_points_);
        for (const PointXYZP& point : candidate.getPointsToRaycast()) {
            if (grid_tools::isInBBox(explo_utils::point3D(point), candidate.getBBox())) {
                PointRay point_ray;
                getPointRayFromPointToCandidate(point_ray, explo_utils::point3D(point), candidate.getCoords());
                if (point_ray.isVisible()) {
                    candidate.addPointRay(point_ray);
                }
            }
        }
        //Initialize and update the stats on the candidate 
        //the update is computed for all visible voxels from candidate (from the previously added pointray)
        candidate.updateObsStat(point_stat_cloud_);
        //finally, we get the IG
        double ig_for_candidate = candidate.getIG();
        double cost_for_candidate = candidate.getCost();
	{
	    lock_guard<mutex> lock(results_mutex_);
	    ig_per_candidate_.insert(std::make_pair(ig_for_candidate, candidate.getCoords()));
        result_per_candidate_.push_back(explo_utils::CandidateResult(ig_for_candidate, cost_for_candidate, candidate.getCoords()));
	}
        candidate.clear();
    }
}

void InformationGain::updateStatsInListOfCandidates() {
    
    {
        lock_guard<mutex> lock(results_mutex_);
        result_per_candidate_.clear();
        ig_per_candidate_.clear();
    }
    vector<shared_ptr<thread> > threads(num_threads_);
    for (size_t i=0; i<num_threads_; i++) {
	threads[i].reset(new thread(&InformationGain::processListOfCandidates, this));
    }
    for (size_t i=0; i<num_threads_; i++) {
	threads[i]->join();
    }
    
    //now it is sorted (from multimap), store to vector
    ig_per_candidate_vector_.clear();
    for (const auto& pair : ig_per_candidate_) {
        ig_per_candidate_vector_.push_back(pair);
    }
}
