#include "explo_quality_aware/PointRay.h"

PointRay::PointRay(double x, double y, double z) {
    coords_ = pcl::PointXYZ(x,y,z);
}
PointRay::PointRay(const pcl::PointXYZ& p): coords_(p) {
}

PointRay::PointRay(const pcl::PointXYZ& p, const RayProbaVector& v):
        coords_(p), ray_proba_vector_(v) {
            setIsVisible();
            setIsAHit();
}
void PointRay::setIsVisible() {
    is_visible_ = true;
    for (const auto& point : ray_proba_vector_) {
        if (point.first.proba >= occupancy_threshold_) {
            if (explo_utils::isSamePoint(explo_utils::point3D(point.first),coords_)) {
                continue;
            } else {
                is_visible_ = false;
                return;
            }
        }
    }

}
void PointRay::setIsAHit() {
    is_a_hit_ = false;
    for (const auto& point : ray_proba_vector_) {
        if (explo_utils::isSamePoint(explo_utils::point3D(point.first),coords_)) {
            if (point.first.proba >= occupancy_threshold_) {
                is_a_hit_ = true;
                return;
            }
        }
    }
}
void PointRay::setRayProbaVector(const RayProbaVector& v) {
    ray_proba_vector_ = v;
    setIsVisible();
    setIsAHit();
}
const PointRay::RayProbaVector& PointRay::getRayProbaVector() const {
    return ray_proba_vector_;
}

void PointRay::setCoords(const pcl::PointXYZ& p) {
    coords_ = p;
}
const pcl::PointXYZ& PointRay::getCoords() const {
    return coords_;
}
bool PointRay::isVisible() const {
    return is_visible_;
}
bool PointRay::isAHit() const {
    return is_a_hit_;
}
double PointRay::getOccupancyThreshold() const {
    return occupancy_threshold_;
}
void PointRay::setOccupancyThreshold(double v) {
    occupancy_threshold_ =  v;
}