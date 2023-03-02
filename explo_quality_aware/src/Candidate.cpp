#include "explo_quality_aware/Candidate.h"
#include <numeric>
#include <random>

//Core operation, we update stat_cloud with a point_ray
void Candidate::updateObsStatWithPointRay(const PointRay& p) {
    for (const std::pair<PointXYZP,obs_stat_tools::SphericalPoint>& point : p.getRayProbaVector()) {
        pcl::PointXYZ current = pcl::PointXYZ(point.first.x, point.first.y, point.first.z);
        if (!grid_tools::isInBBox(current, bbox_)){
            continue;
        }
        cv::Point3i i = grid_tools::metersToVox(current, bbox_, resolution_);
        assert(grid_tools::isVoxInGrid(i, stat_grid_));
        //In the Information Gain, we add pointsRay only if the point is visible
        //We cannot reach here if the voxel is hidden
        update(i, point.second);
    } 
}


//We update the statistics, as in ObsStatGridBase, but we also update which angular sectors this candidate adds
void Candidate::update(const cv::Point3i& vox, const obs_stat_tools::SphericalPoint& obs) {
    ObsStatGridBase::update(vox, obs);
    updateAngularSectors(vox);
}

void Candidate::updateAngularSectors(const cv::Point3i& vox) {
    //We compute the new angular sectors viewpoints the candidate adds to the considered voxel
    //And we add them to the current viewpoints (boolean OR)
    uint32_t new_angular_sectors = getAngularSectorsForVoxFromCandidate(vox);
    stat_grid_[vox.x][vox.y][vox.z].angular_viewpoints_2d = (stat_grid_[vox.x][vox.y][vox.z].angular_viewpoints_2d | new_angular_sectors);
}
        
//If we rot clockwise, we fill between start and end
uint32_t Candidate::updateAngularSectorsClockwise(uint32_t in_sectors, uint8_t start, uint8_t end) const {
    uint32_t angular_sectors = in_sectors;
    for (uint8_t i=start; i< end+1; i++) { 
        angular_sectors = setBitItoV(angular_sectors, i, true);
    }
    return angular_sectors;
}
//If we rot counterclockwise, we fill between end and start
//If max_angle is Pi, end could be 0, we treat that differently
uint32_t Candidate::updateAngularSectorsCounterClockwise(uint32_t in_sectors, uint8_t start, uint8_t end) const {
    uint32_t angular_sectors = in_sectors;
    if (end != 0) {
        for (uint8_t i=end; i < 32; i++) { 
            angular_sectors = setBitItoV(angular_sectors, i, true);
        }
        for (uint8_t i=0; i < start+1; i++) { 
            angular_sectors = setBitItoV(angular_sectors, i, true);
        }
    } else {
        for (uint8_t i=0; i < start+1; i++) { 
            angular_sectors = setBitItoV(angular_sectors, i, true);
        }
    }
    return angular_sectors;

}

uint32_t Candidate::getAngularSectorsForVoxFromCandidate(const cv::Point3i& vox) const {
    pcl::PointXYZ point = grid_tools::voxToMeters(vox, bbox_, resolution_);
    //the obs view point from the robot (theta1), from the candidate (theta3) 
    double theta1 = obs_stat_tools::cart2spherical(robot_.x - point.x, robot_.y - point.y, robot_.z - point.z).theta;
    double theta3 = obs_stat_tools::cart2spherical(coords_.x - point.x, coords_.y - point.y, coords_.z - point.z).theta;

    uint32_t angular_sectors = 0;
    uint8_t start; uint8_t end;
    bool rot_clockwise = false;
    //The observed sector is smaller than pi
    //We always start with the smaller angle. If the difference greater-smaller > Pi, we go clockwise, else counter clockwise
    double min_angle = std::min(theta1, theta3);
    double max_angle = std::max(theta1, theta3);
    if ((max_angle - min_angle)>=M_PI){
        rot_clockwise = false;
    } else {
        rot_clockwise = true; 
    }
    start = thetaToAngularSector(min_angle);
    end = thetaToAngularSector(max_angle);
    static std::random_device rd;
    static std::default_random_engine eng(rd());
    static std::uniform_int_distribution<int> distri(0,1); 
    //Edge case, we have only one sector to fill
    if (start == end) {
        angular_sectors = setBitItoV(angular_sectors, start, true);
    //Edge case, we have exactly half a circle to fill, we select the half circle randomly
    }else if ((end - start) == 16) {
        int j = distri(eng);
        if (j == 0) {
            angular_sectors = updateAngularSectorsClockwise(angular_sectors, start, end);
        } else {
        angular_sectors = updateAngularSectorsCounterClockwise(angular_sectors, start, end);
        }
    //Default case, we fill either clockwise or counterclockwise
    } else if (rot_clockwise) {
        angular_sectors = updateAngularSectorsClockwise(angular_sectors, start, end);
    } else {
        angular_sectors = updateAngularSectorsCounterClockwise(angular_sectors, start, end);
    }
    return angular_sectors;
}

void Candidate::setPointsToRaycast(const std::vector<PointXYZP>& list_occupied_points) {
    //Get the Points on the cylinder around the candidate
    //later on, when we raycast, the points on the border hidden by an occupied point will be ignored
    VoxPointMap point_map(point3dLess);
    
    //First, iterate on list of occupied_points, if point in Candidate BBOX -> add it to points_to_raycast
    for (const PointXYZP& p : list_occupied_points) {
        if (grid_tools::isInBBox(explo_utils::point3D(p), bbox_)) {
            cv::Point3i vox = grid_tools::metersToVox(explo_utils::point3D(p), bbox_, resolution_);
            point_map.insert(std::make_pair(vox, p.proba));
        }
    }
    // Iterate on all the voxels in the bbox, and store them in the same map with occ=0 (if the point is already occupied, it will not be inserted) 
    for (double x = bbox_.xmin+resolution_; x < bbox_.xmax; x += resolution_) {
        for (double y = bbox_.ymin+resolution_; y < bbox_.ymax; y += resolution_) {
            for (double z = bbox_.zmin+resolution_; z < bbox_.zmax; z += resolution_) {
                pcl::PointXYZ pz = pcl::PointXYZ(x, y, z);
                if (grid_tools::isInBBox(pz, bbox_)) {
                    cv::Point3i vox = grid_tools::metersToVox(pz, bbox_, resolution_);
                        point_map.insert(std::make_pair(vox, 0));
                }
            }
        }
    }
    //now we store the point_set in a vector to be used later (this is more convenient to use outside the class as a vector)
    points_to_raycast_.clear();
    points_to_raycast_.resize(point_map.size());
    size_t i=0;
    for (const auto& pair : point_map) {
        pcl::PointXYZ point = grid_tools::voxToMeters(pair.first, bbox_, resolution_);
        points_to_raycast_.at(i) = PointXYZP(point.x, point.y, point.z, pair.second);
        ++i;
    }
    //And now we sort points_to_raycast by occupancy likelihood, so later occupied points are processed first
    std::sort(points_to_raycast_.begin(), points_to_raycast_.end(),
        [](const auto& lhs, const auto& rhs) {return lhs.proba > rhs.proba;});
}


void Candidate::initializePointStatsCoords() {
    double x; double y; double z;
    for (size_t ix{0}; ix < grid_size_.nx; ix++) {
        x = ix * resolution_ + bbox_.xmin;
        for (size_t iy{0}; iy < grid_size_.ny; iy++) {
            y = iy * resolution_ + bbox_.ymin;
            for (size_t iz{0}; iz < grid_size_.nz; iz++) {
                z = iz * resolution_ + bbox_.zmin;
                stat_grid_[ix][iy][iz].x = x;
                stat_grid_[ix][iy][iz].y = y;
                stat_grid_[ix][iy][iz].z = z;
                init_stat_grid_[ix][iy][iz].x = x;
                init_stat_grid_[ix][iy][iz].y = y;
                init_stat_grid_[ix][iy][iz].z = z;
            }
        }
    }
}

//Transform the PointCloud of PointStat into a VoxelGrid ObsStatGrid
void Candidate::voxGridFromPointCloud(const pcl::PointCloud<PointStat>& stat_cloud) {
    for (const auto& point : stat_cloud) {
        //We keep only points inside the candidate bbox
        if (!grid_tools::isInBBox(pcl::PointXYZ(point.x, point.y, point.z), bbox_)) {
            continue;
        }
        cv::Point3i i = grid_tools::metersToVox(pcl::PointXYZ(point.x, point.y, point.z), bbox_, resolution_);
        if (grid_tools::isVoxInGrid(i, stat_grid_)) {
            stat_grid_[i.x][i.y][i.z].n_obs = point.n_obs;
            stat_grid_[i.x][i.y][i.z].sigma_angle = point.sigma_angle;
            stat_grid_[i.x][i.y][i.z].sigma_angle2d = point.sigma_angle2d;
            stat_grid_[i.x][i.y][i.z].r_mean = point.r_mean;
            stat_grid_[i.x][i.y][i.z].r_min = point.r_min;
            stat_grid_[i.x][i.y][i.z].r_max = point.r_max;
            stat_grid_[i.x][i.y][i.z].sigma_r = point.sigma_r;
            stat_grid_[i.x][i.y][i.z].n_hits = point.n_hits;
            stat_grid_[i.x][i.y][i.z].x_mean = point.x_mean;
            stat_grid_[i.x][i.y][i.z].y_mean = point.y_mean;
            stat_grid_[i.x][i.y][i.z].z_mean = point.z_mean;
            stat_grid_[i.x][i.y][i.z].x_sum = point.x_sum;
            stat_grid_[i.x][i.y][i.z].y_sum = point.y_sum;
            stat_grid_[i.x][i.y][i.z].z_sum = point.z_sum;
            stat_grid_[i.x][i.y][i.z].theta_mean = point.theta_mean;
            stat_grid_[i.x][i.y][i.z].phi_mean = point.phi_mean;
            stat_grid_[i.x][i.y][i.z].angular_viewpoints_2d = point.angular_viewpoints_2d;
        }
    }
    
    //stat_grid_ will be updated with the observations, we keep init_stat_grid_ to compute the information gain
    init_stat_grid_ = stat_grid_;
}

void Candidate::updateObsStat(const pcl::PointCloud<PointStat>& stat_cloud) {
    voxGridFromPointCloud(stat_cloud);
    for (const PointRay& p : point_ray_list_) {
        updateObsStatWithPointRay(p);
    }
}

void Candidate::setBBox(double sensor_range, const grid_tools::BBox& valid) {
    bbox_.xmin = std::max(coords_.x - sensor_range, valid.xmin);
    bbox_.ymin = std::max(coords_.y - sensor_range, valid.ymin);
    bbox_.zmin = std::max(coords_.z - sensor_range, valid.zmin); 
    bbox_.xmax = std::min(coords_.x + sensor_range, valid.xmax);
    bbox_.ymax = std::min(coords_.y + sensor_range, valid.ymax);
    bbox_.zmax = std::min(coords_.z + sensor_range, valid.zmax);
}     

void Candidate::addPointRay(const PointRay& point_ray) {
    point_ray_list_.push_back(point_ray);
}

//Compute the actual IG for the current candidate
double Candidate::getIG() const {
    Policy policy(policy_choice_);
    return policy.getPolicyFromGrids(stat_grid_, init_stat_grid_, grid_size_);
}

void Candidate::printDebug() const {
    
    std::vector<int> n_obs_init;
    std::vector<int> n_obs_updated;
    std::vector<double> sigma_init;
    std::vector<double> sigma_updated;
    std::vector<double> range_init;
    std::vector<double> range_updated;
    for (size_t ix{0}; ix < grid_size_.nx; ix++) {
        for (size_t iy{0}; iy < grid_size_.ny; iy++) {
            for (size_t iz{0}; iz < grid_size_.nz; iz++) {
                n_obs_init.push_back(init_stat_grid_[ix][iy][iz].n_obs);
                sigma_init.push_back(init_stat_grid_[ix][iy][iz].sigma_angle);
                range_init.push_back(init_stat_grid_[ix][iy][iz].r_min);
                n_obs_updated.push_back(stat_grid_[ix][iy][iz].n_obs);
                sigma_updated.push_back(stat_grid_[ix][iy][iz].sigma_angle);
                range_updated.push_back(stat_grid_[ix][iy][iz].r_min);
            }
        }
    }
    int max_n; int min_n;
    double max_sigma; double min_sigma;
    double max_range; double min_range;
    max_n = *std::max_element(n_obs_init.begin(), n_obs_init.end());
    min_n = *std::min_element(n_obs_init.begin(), n_obs_init.end());
    max_sigma = *std::max_element(sigma_init.begin(), sigma_init.end());
    min_sigma = *std::min_element(sigma_init.begin(), sigma_init.end());
    max_range = *std::max_element(range_init.begin(), range_init.end());
    min_range = *std::min_element(range_init.begin(), range_init.end());
    
    std::cout << "before the update, min and max vals are :" << std::endl;
    std::cout << "n : " << min_n << ", " << max_n << std::endl;
    std::cout << "sigma : " << std::scientific << min_sigma << ", " << std::scientific << max_sigma << std::endl;
    std::cout << "range : " << std::scientific << min_range << ", " << std::scientific << max_range << std::endl;
    
    max_n = *std::max_element(n_obs_updated.begin(), n_obs_updated.end());
    min_n = *std::min_element(n_obs_updated.begin(), n_obs_updated.end());
    max_sigma = *std::max_element(sigma_updated.begin(), sigma_updated.end());
    min_sigma = *std::min_element(sigma_updated.begin(), sigma_updated.end());
    max_range = *std::max_element(range_updated.begin(), range_updated.end());
    min_range = *std::min_element(range_updated.begin(), range_updated.end());
    
    std::cout << "after the update, min and max vals are :" << std::endl;
    std::cout << "n : " << min_n << ", " << max_n << std::endl;
    std::cout << "sigma : " << std::scientific << min_sigma << ", " << std::scientific << max_sigma << std::endl;
    std::cout << "range : " << std::scientific << min_range << ", " << std::scientific << max_range << std::endl;

    std::cout << "EXITING" << std::endl;
    exit(EXIT_SUCCESS);
}