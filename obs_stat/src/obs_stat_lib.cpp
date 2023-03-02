#include <algorithm>
#include <iostream>
#include <map>
#include "obs_stat/obs_stat_lib.h"

using namespace std;
using namespace obs_stat_tools;
using namespace grid_tools;

//Check that the voxel has not been already updated, if not update it
void ObsStatGridBase::update(const cv::Point3i& vox, const SphericalPoint& obs) {
    assert(isVoxInGrid(vox, stat_grid_));
    if (isMaskTrue(vox, mask_grid_)) return;
    update(stat_grid_[vox.x][vox.y][vox.z], obs);
    setMaskTrue(vox, mask_grid_);
}

//Core of the update function, where the statistics are computed
void ObsStatGridBase::update(PointStat& stat, const SphericalPoint& obs) {
    //we first update the stats on ranges:
    stat.r_mean = (stat.n_obs*stat.r_mean+obs.r)/(stat.n_obs+1);
    stat.sigma_r = (pow(obs.r-stat.r_mean, 2)+stat.n_obs*stat.sigma_r)/(stat.n_obs+1);
    if (stat.n_obs == 0) {
        stat.r_min = obs.r;
        stat.r_max = obs.r;
    }
    if (obs.r < stat.r_min) {
        stat.r_min = obs.r;
    } else if (obs.r > stat.r_max) {
        stat.r_max = obs.r;
    } 
    
    //we then update sums, means, n_obs and finally the stats on the angles
    pcl::PointXYZ uVec = spherical2unitVector(obs);
    stat.x_sum += uVec.x;
    stat.y_sum += uVec.y;
    stat.z_sum += uVec.z;
    stat.x_mean = (stat.n_obs*stat.x_mean+uVec.x)/(stat.n_obs+1);
    stat.y_mean = (stat.n_obs*stat.y_mean+uVec.y)/(stat.n_obs+1);
    stat.z_mean = (stat.n_obs*stat.z_mean+uVec.z)/(stat.n_obs+1);
    
    double R = sqrt(pow(stat.x_sum,2)+pow(stat.y_sum,2)+pow(stat.z_sum,2));
    double Rxy = sqrt(pow(stat.x_sum,2)+pow(stat.y_sum,2));
    stat.n_obs +=1;
    if (stat.n_obs > 1 ) {
        stat.sigma_angle = std::max(0., 1 - (R/stat.n_obs));
        stat.sigma_angle2d = std::max(0., 1 - (Rxy/stat.n_obs));
    } else {
        stat.sigma_angle = 0;
        stat.sigma_angle2d = 0;
    }
    //Finally, we update the viewpoints table
    stat.angular_viewpoints_2d = getUpdatedAngularSector(stat.angular_viewpoints_2d, obs.theta);
}
//This is not really a statistic, n_hits can be used as a proxy to the occupancy likelihood of a voxel
//n_hits is just a "help" variable, use with caution
void ObsStatGridBase::updateNHits(const pcl::PointXYZ& point) {
    cv::Point3i vox = metersToVox(point, bbox_, resolution_); 
    assert(isVoxInGrid(vox, stat_grid_));
    if (isMaskTrue(vox, n_hits_mask_grid_)) return;
    updateNHits(vox);
    setMaskTrue(vox, n_hits_mask_grid_);
}

void ObsStatGridBase::updateNHits(const cv::Point3i& vox) {
    assert(isVoxInGrid(vox, stat_grid_));
    stat_grid_[vox.x][vox.y][vox.z].n_hits +=1;
}
//Initialize the GridStat with PointStats with the correct coordinates (the grid is initialized with default values in the constructor)
void ObsStatGridBase::initializePointStatsCoords() {
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
            }
        }
    }
}
//Specialization of the ObsStatGrid, that call the appropriate function to perform the raycasting
void ObsStatGrid::addObsStatRayToGrid(const pcl::PointXYZ& point_in_map, const pcl::PointXYZ& laser_in_map, 
    const SphericalPoint& obs, std::ostream& warn_stream) {
    GridIntersectMap grid_intersect(point3dLess);
    getGridSegmentIntersect(point_in_map, laser_in_map, obs, grid_intersect, warn_stream);
    for(const auto& pair : grid_intersect) {
        update(pair.first, pair.second);
    }
}

//Raycasting operation: all the voxels in the intersection GRID [point-laser] are updated with the observation viewpoint [vox-laser]
void ObsStatGrid::getGridSegmentIntersect(const pcl::PointXYZ& point_in_map, const pcl::PointXYZ& laser_in_map, 
        const SphericalPoint& obs, GridIntersectMap& grid_intersect, std::ostream& warn_stream) const {
    //we definte the direction vector and the origin of the segment
    //we look for the voxel contaning the start and end of the segment
    static size_t too_small_counter = 0;
    cv::Point3i vox_laser = metersToVox(laser_in_map, bbox_, resolution_);
    cv::Point3i vox_point = metersToVox(point_in_map, bbox_, resolution_);
    //we get the coordinates of the voxel in the map
    pcl::PointXYZ start_vox = voxToMeters(vox_point, bbox_, resolution_); 
    pcl::PointXYZ end_vox = voxToMeters(vox_laser, bbox_, resolution_);
    //we get the direction vector from point to laser, that gives the direction of the line
    pcl::PointXYZ uVec = spherical2unitVector(obs);
    //we define the functions relative to the equation of the line
    double a = uVec.x;
    double b = uVec.y; 
    double c = uVec.z;
    if ((fabs(a)<1e-6) || (fabs(b)<1e-6) || (fabs(c)<1e-6)) {
        too_small_counter++;
        if (too_small_counter > 50) {
            warn_stream << "At least 50 cubes were discarded because parallels to the axes, last coefficients of the unit vector are: " 
                    << a << ", " << b << ", " << c << std::endl;
            warn_stream << "Last observation is: r " << obs.r << ", theta: " << obs.theta << ", phi: " << obs.phi << ". Resetting counter." << std::endl;
            too_small_counter =0;
        }
        return;
    }
    double x1 = point_in_map.x;
    double y1 = point_in_map.y;
    double z1 = point_in_map.z;
    auto fy_x = [a,b,c,x1,y1,z1](double x) {return ((x-x1)/a)*b + y1;};
    auto fz_x = [a,b,c,x1,y1,z1](double x) {return ((x-x1)/a)*c + z1;};
    auto fx_y = [a,b,c,x1,y1,z1](double y) {return ((y-y1)/b)*a + x1;};
    auto fz_y = [a,b,c,x1,y1,z1](double y) {return ((y-y1)/b)*c + z1;};
    auto fx_z = [a,b,c,x1,y1,z1](double z) {return ((z-z1)/c)*a + x1;};
    auto fy_z = [a,b,c,x1,y1,z1](double z) {return ((z-z1)/c)*b + y1;};
    //for each point on the segment that intersect axis x, we calculate coordinates y,z to get the coordinates of the voxel
    //we calculate the Observation (r,theta,phi) from which this voxel was observed
    //we store the pair (voxel, observation) in a map to avoid duplicate
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
            if (!isInBBox(pcl::PointXYZ(x,y,z), bbox_)) {
                continue;
            }
            cv::Point3i vox = metersToVox(pcl::PointXYZ(x,y,z), bbox_, resolution_);
            if (!isVoxInGrid(vox, stat_grid_)) {
                continue;
            }
            pcl::PointXYZ laser_in_point(laser_in_map.x - x, laser_in_map.y - y, laser_in_map.z - z); 
            SphericalPoint curr_obs{cart2spherical(laser_in_point)};
            grid_intersect.insert(std::make_pair(vox, curr_obs));
        }
    } else if (y_diff == max_diff) {
        //same for intersection with y axis
        for (double y = start_vox.y; ((y - start_vox.y)*y_increment) < ((end_vox.y - start_vox.y)*y_increment);y+=y_increment){
            double x = fx_y(y);
            double z = fz_y(y);
            if (!isInBBox(pcl::PointXYZ(x,y,z), bbox_)) {
                continue;
            }
            cv::Point3i vox = metersToVox(pcl::PointXYZ(x,y,z), bbox_, resolution_);
            if (!isVoxInGrid(vox, stat_grid_)) {
                continue;
            }
            pcl::PointXYZ laser_in_point(laser_in_map.x - x, laser_in_map.y - y, laser_in_map.z - z); 
            SphericalPoint curr_obs{cart2spherical(laser_in_point)};
            grid_intersect.insert(std::make_pair(vox, curr_obs));
        }
    } else {
        //same with intersection with z axis
        for (double z = start_vox.z; ((z - start_vox.z)*z_increment) < ((end_vox.z - start_vox.z)*z_increment);z+=z_increment){
            double x = fx_z(z);
            double y = fy_z(z);
            if (!isInBBox(pcl::PointXYZ(x,y,z), bbox_)) {
                continue;
            }
            cv::Point3i vox = metersToVox(pcl::PointXYZ(x,y,z), bbox_, resolution_);
            if (!isVoxInGrid(vox, stat_grid_)) {
                continue;
            }
            pcl::PointXYZ laser_in_point(laser_in_map.x - x, laser_in_map.y - y, laser_in_map.z - z); 
            SphericalPoint curr_obs{cart2spherical(laser_in_point)};
            grid_intersect.insert(std::make_pair(vox, curr_obs));
        }
    }
}
