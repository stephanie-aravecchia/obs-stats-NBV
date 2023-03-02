#ifndef POINTRAY_H
#define POINTRAY_H
#include <vector>
#include <pcl_ros/point_cloud.h>
#include <pcl_ros/transforms.h>
#include <pcl/point_types.h>
#include "explo_quality_aware/proba_point_type.h"
#include "explo_quality_aware/Utils.h"
#include "explo_quality_aware/Constants.h"
#include "obs_stat/CoordinateConversions.h"

//This class manipulates PointRay, a point and an associated ray from a raycast operation in a voxelgrid
//Because of raycasting, the ray is expected sorted (from coord, in the ray direction to the robot)
class PointRay {
    protected:
        double occupancy_threshold_ = OCC_THRESHOLD;
    
        pcl::PointXYZ coords_;
        bool is_visible_{true};
        bool is_a_hit_{false};
    public:
        typedef std::vector<std::pair<PointXYZP,obs_stat_tools::SphericalPoint>> RayProbaVector; 
        
        RayProbaVector ray_proba_vector_;

    protected:
        void setIsVisible();
        void setIsAHit();

    public:
        PointRay() {};
        PointRay(double x, double y, double z);
        PointRay(const pcl::PointXYZ& p);
        PointRay(const pcl::PointXYZ& p, const RayProbaVector&);
        void setCoords(const pcl::PointXYZ& p);
        const pcl::PointXYZ& getCoords() const;
        bool isVisible() const;
        bool isAHit() const;
        void setRayProbaVector(const RayProbaVector&);
        const RayProbaVector& getRayProbaVector() const;
        double getOccupancyThreshold() const ;
        void setOccupancyThreshold(double);
};
#endif