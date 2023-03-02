#ifndef OBSERVATION_STAT_LIB_H
#define OBSERVATION_STAT_LIB_H
#include <obs_stat/Tools.h>
#include <obs_stat/CoordinateConversions.h>

//This class implements a voxel grid containing PointStats
//And the functions to update the statistics on the observations viewpoints
class ObsStatGridBase {
    protected:
        grid_tools::StatGrid stat_grid_;
        grid_tools::MaskGrid mask_grid_;
        grid_tools::MaskGrid n_hits_mask_grid_;
        grid_tools::BBox bbox_;
        double resolution_;
        grid_tools::GridSize grid_size_;

        //Core of the update function, the statistics are implemented here
        void update(PointStat& stat, const obs_stat_tools::SphericalPoint& obs);
        
        //Utility Functions to manipulate MaskGrid 
        bool isMaskTrue(const cv::Point3i& vox, const grid_tools::MaskGrid& mask) const {
            return mask[vox.x][vox.y][vox.z];
        }
        void setMaskTrue(const cv::Point3i& vox, grid_tools::MaskGrid& mask) {
            mask[vox.x][vox.y][vox.z] = true;
        }
        void initializePointStatsCoords();
        
        static bool point3dLess(const cv::Point3i& v1, const cv::Point3i& v2) 
        {   if (v1.x < v2.x) return true;
            if (v1.x > v2.x) return false;
            // v1.x == v2.x
            if (v1.y < v2.y) return true;
            if (v1.y > v2.y) return false;
            // v1.x == v2.x && v1.y == v2.y
            return v1.z < v2.z;
        };
        
        //update the stats of a voxel, with an observation, if the voxel has not been already updated
        void update(const cv::Point3i& vox, const obs_stat_tools::SphericalPoint& obs);

        //The statistics on the angles are encoded as follows:
        //The theta angle (azimuthal in xy plane) is discretized in 32
        //The 32 angular sectors are encoded with a uint_32t
        //If a voxel has been observed from an angular sector, the corresponding bit is set to 1
        //The following functions implement that:
        uint32_t setBitItoV(uint32_t init_value, uint8_t bit_i, bool v) const {
            assert(bit_i < 32);
            return (init_value & ~(1 << bit_i)) | (v << bit_i);
        }
        bool readBitIFromV(uint32_t value, uint8_t bit_i) const {
            assert(bit_i < 32);
            return (value & (1 << bit_i)) != 0;
        }
        uint8_t thetaToAngularSector(double theta) const {
            //theta is in the range -pi, pi, we discretize it in 32
            return static_cast<uint8_t>(floor(fmod(theta + 3*M_PI, 2*M_PI)*(32/(2*M_PI))));
        }
        uint32_t getUpdatedAngularSector(uint32_t init_sectors, double theta) const {
            uint8_t bit = thetaToAngularSector(theta);
            return setBitItoV(init_sectors, bit, 1);
        }
    
    public:
        ObsStatGridBase() {};
        ObsStatGridBase(const grid_tools::BBox bbox, double resolution) :
            bbox_(bbox), resolution_(resolution) {
                grid_size_ = grid_tools::getGridSizeFromBBoxAndRes(bbox_, resolution_);
                grid_tools::resizeAndInitializeGrid(stat_grid_, grid_size_, PointStat());
                grid_tools::resizeAndInitializeGrid(mask_grid_, grid_size_, false);
                grid_tools::resizeAndInitializeGrid(n_hits_mask_grid_, grid_size_, false);
                initializePointStatsCoords();

            };
        //get the class members
        const grid_tools::BBox& getBBox() const {return bbox_;};
        double getResolution() const {return resolution_;};
        const grid_tools::GridSize& getGridSize() const {return grid_size_;};
        const grid_tools::StatGrid& getStatGrid() const {return stat_grid_;};
        const grid_tools::MaskGrid& getMaskGrid() const {return mask_grid_;};
        

        void resetMaskGrid() {
            grid_tools::resetGrid(mask_grid_, grid_size_, false);
            grid_tools::resetGrid(n_hits_mask_grid_, grid_size_, false);
            };
        
        //This update n_hits, which is not really a stat
        //The value is useful to know if the voxel contains a returned point from the laser, or is just raycasted
        //This does not check the mask
        void updateNHits(const pcl::PointXYZ& point);
        void updateNHits(const cv::Point3i& vox);

};
class ObsStatGrid : public ObsStatGridBase {
    protected:
        //GridIntersectMap and the function below are specialised, grid_intersect depends on what we intend to do
        //Here, we just deal with the stats
        typedef std::map<cv::Point3i, obs_stat_tools::SphericalPoint, decltype(point3dLess)*> GridIntersectMap;
        
        //This function expects an ostream, to display potential warnings
        void getGridSegmentIntersect(const pcl::PointXYZ& point_in_map, const pcl::PointXYZ& laser_in_map, 
            const obs_stat_tools::SphericalPoint& obs, GridIntersectMap& grid_intersect,
            std::ostream& warn_stream = std::cout) const;
    
    public:
        ObsStatGrid() {};
        ObsStatGrid(const grid_tools::BBox bbox, double resolution) : ObsStatGridBase(bbox, resolution) {};

        //This function prepare the raycasting operation to infer the voxels in the segment point-laser
        //And then calls the appropriate update function 
        void addObsStatRayToGrid(const pcl::PointXYZ& point_in_map, const pcl::PointXYZ& laser_in_map, 
            const obs_stat_tools::SphericalPoint& obs, std::ostream& warn_stream = std::cout);
};
//Useful to debug
std::ostream& operator<<(std::ostream& out, const PointStat& stat) {
    out << "x: " << stat.x << " " 
        << "y: " << stat.y << " " 
        << "z: " << stat.z << " " 
        << "n_obs: " << stat.n_obs << " " 
        << ", x_mean: " << stat.x_mean << " "
        << ", y_mean: " << stat.y_mean << " "
        << ", z_mean: " << stat.z_mean << " "
        << ", x_sum: " << stat.x_sum << " "
        << ", y_sum: " << stat.y_sum << " "
        << ", z_sum: " << stat.z_sum << " "
        << ", theta_mean: " << stat.theta_mean << " "
        << ", phi_mean: " << stat.phi_mean << " " 
        << ", r_min: " << stat.r_min << " "
        << ", r_max: " << stat.r_max << " "
        << ", r_mean: " << stat.r_mean << " "
        << ", sigma_r: " << stat.sigma_r << " "
        << ", sigma_angle: " << stat.sigma_angle << " "
        << ", sigma_angle2d: " << stat.sigma_angle2d << " "
        << ", angular_viewpoints_2d: " << stat.angular_viewpoints_2d << " "
        << ", n_hits: " << stat.n_hits << ".";
    return out;
};

#endif