#ifndef COORDS_TOOLS_H
#define COORDS_TOOLS_H

#include <pcl/point_types.h>
#include <math.h>

namespace obs_stat_tools {
    
    struct SphericalPoint {
        float r;
        float theta;
        float phi;
        SphericalPoint() : r(0), theta(0), phi(0) {}
        SphericalPoint(float r, float theta, float phi):r(r), theta(theta), phi(phi) {}
    };

    //We follow the mathematical convention
    pcl::PointXYZ spherical2unitVector(const SphericalPoint& obs) {
        return pcl::PointXYZ(cos(obs.theta) * sin(obs.phi),
                            sin(obs.theta) * sin(obs.phi),
                            cos(obs.phi));
    };
    //we follow the mathematical convention
    pcl::PointXYZ spherical2cart(const SphericalPoint& obs){
        return pcl::PointXYZ(obs.r * cos(obs.theta) * sin(obs.phi),
                            obs.r * sin(obs.theta) * sin(obs.phi),
                            obs.r * cos(obs.phi));
    };
    //The values are: r, phi (polar angle, measured from z), theta (azimutal angle in the xy plane from the x axis)
    //This is the correct definition of phi and theta in mathematical convention
    SphericalPoint cart2spherical(float x, float y, float z) {
        SphericalPoint spherical_point;
        spherical_point.r = sqrt(x * x + y * y + z * z);
        spherical_point.phi = atan2(hypot(x, y), z);
        spherical_point.theta = atan2(y, x);
        return spherical_point;
    };
    SphericalPoint cart2spherical(const pcl::PointXYZ& p) {
        return cart2spherical(p.x, p.y,p.z);
    };
}
std::ostream& operator<<(std::ostream& out, const obs_stat_tools::SphericalPoint& point) {
    out << "r: " << point.r << " " 
        << "theta: " << point.theta << " " 
        << "phi:" << point.phi; 
    return out;
};
#endif