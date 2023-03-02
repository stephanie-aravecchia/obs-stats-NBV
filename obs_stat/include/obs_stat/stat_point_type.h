#ifndef STAT_POINT_TYPE_H
#define STAT_POINT_TYPE_H

#define PCL_NO_PRECOMPILE
#include <pcl/pcl_macros.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/io/pcd_io.h>

#include "obs_stat/Constants.h"

struct EIGEN_ALIGN16 _PointStat      // enforce SSE padding for correct memory alignment
{
  PCL_ADD_POINT4D;                  // preferred way of adding a XYZ+padding
  std::uint32_t n_obs ;
  float x_mean;
  float y_mean;
  float z_mean;
  float x_sum;
  float y_sum;
  float z_sum;
  float theta_mean;
  float phi_mean;
  float r_min;
  float r_max;
  float r_mean;
  float sigma_r;
  float sigma_angle;
  float sigma_angle2d;
  std::uint32_t angular_viewpoints_2d;
  std::uint32_t n_hits ;
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW     // make sure our new allocators are aligned
} ;                    
struct EIGEN_ALIGN16 PointStat : public _PointStat
{
    inline PointStat(const _PointStat &p) 
    {
        x = p.x;
        y = p.y;
        z = p.z;
        n_obs = p.n_obs;
        x_mean = p.x_mean;
        y_mean = p.y_mean;
        z_mean = p.z_mean;
        x_sum = p.x_sum;
        y_sum = p.y_sum;
        z_sum = p.z_sum;
        theta_mean = p.theta_mean;
        phi_mean = p.phi_mean;
        r_min = p.r_min;
        r_max = p.r_max;
        r_mean = p.r_mean;
        sigma_r = p.sigma_r;;
        sigma_angle = p.sigma_angle;
        sigma_angle2d = p.sigma_angle2d;
        angular_viewpoints_2d = p.angular_viewpoints_2d;
        n_hits = p.n_hits;
    }
    
    inline PointStat()
    {
        x = 0;
        y = 0;
        z = 0;
        n_obs = 0;
        x_mean = 0;
        y_mean = 0;
        z_mean = 0;
        x_sum = 0;
        y_sum = 0;
        z_sum = 0;
        theta_mean = 0;
        phi_mean = 0;
        r_min = MIN_RANGE_INIT;
        r_max = 0;
        r_mean = 0;
        sigma_r = 0;
        sigma_angle = INVALID_SIGMA;
        sigma_angle2d = INVALID_SIGMA;
        angular_viewpoints_2d = 0;
        n_hits = 0;
    }
    
    inline PointStat(float _x, float _y, float _z, std::uint32_t _n_obs, 
        float _x_mean, float _y_mean, float _z_mean, 
        float _x_sum, float _y_sum, float _z_sum, 
        float _theta_mean, float _phi_mean, float _r_min, float _r_max, float _r_mean, 
        float _sigma_r, float _sigma_angle, float _sigma_angle2d, std::uint32_t _angular_viewpoints_2d,    
        std::uint32_t _n_hits)
    {
        x = _x ;
        y = _y;
        z = _z;
        n_obs = _n_obs;
        x_mean = _x_mean;
        y_mean = _y_mean;
        z_mean = _z_mean;
        x_sum = _x_sum;
        y_sum = _y_sum;
        z_sum = _z_sum;
        theta_mean = _theta_mean;
        phi_mean = _phi_mean;
        r_min = _r_min;
        r_max = _r_max;
        r_mean = _r_mean;
        sigma_r = _sigma_r;;
        sigma_angle = _sigma_angle;
        sigma_angle2d = _sigma_angle2d;
        angular_viewpoints_2d = _angular_viewpoints_2d;
        n_hits = _n_hits;
    }
    inline PointStat(float _x, float _y, float _z)
    {
        x = _x;
        y = _y;
        z = _z;
        n_obs = 0;
        x_mean = 0;
        y_mean = 0;
        z_mean = 0;
        x_sum = 0;
        y_sum = 0;
        z_sum = 0;
        theta_mean = 0;
        phi_mean = 0;
        r_min = MIN_RANGE_INIT;
        r_max = 0;
        r_mean = 0;
        sigma_r = 0;
        sigma_angle = INVALID_SIGMA;
        sigma_angle2d = INVALID_SIGMA;
        angular_viewpoints_2d = 0;
        n_hits = 0;
    }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
} ;

POINT_CLOUD_REGISTER_POINT_STRUCT (PointStat,           
                                   (float, x, x)
                                   (float, y, y)
                                   (float, z, z)
                                   (std::uint32_t, n_obs, n_obs)
                                   (float, x_mean, x_mean) 
                                   (float, y_mean, y_mean)
                                   (float, z_mean, z_mean)
                                   (float, x_sum, x_sum)
                                   (float, y_sum, y_sum)
                                   (float, z_sum, z_sum)
                                   (float, theta_mean, theta_mean)
                                   (float, phi_mean, phi_mean)
                                   (float, r_min, r_min)
                                   (float, r_max, r_max)
                                   (float, r_mean, r_mean)
                                   (float, sigma_r, sigma_r)
                                   (float, sigma_angle, sigma_angle)
                                   (float, sigma_angle2d, sigma_angle2d)
                                   (std::uint32_t, angular_viewpoints_2d, angular_viewpoints_2d)
                                   (std::uint32_t, n_hits, n_hits)
)

POINT_CLOUD_REGISTER_POINT_WRAPPER(PointStat, _PointStat)

#endif // STAT_POINT_TYPE_H
