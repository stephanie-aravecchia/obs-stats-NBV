#ifndef EXPLO_UTILS_H
#define EXPLO_UTILS_H

#include <pcl/point_types.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <opencv2/opencv.hpp>
#include "explo_quality_aware/Constants.h"
#include "obs_stat/Tools.h"

//Librairy of useful functions to manipulate 3D grids and cost maps
namespace explo_utils {
    using namespace grid_tools; 
    
    //Origin of image is coordinates of the (0,0) pixel
    BBox2d get2dBBoxFromImgInfo(double origin_x, double origin_y, int height, int width, double res) {
        double xmin = origin_x;
        double ymin = origin_y - height*res;
        double xmax = origin_x + width*res;
        double ymax = origin_y;
        return BBox2d(xmin, xmax, ymin, ymax);
    };

    struct CostMap {
        cv::Mat_<int8_t> cm;
        double resolution;
        cv::Point2d origin;//coords of top left pixel
        int height;
        int width;
        BBox2d bbox;
        CostMap() {};
        CostMap(const cv::Mat_<int8_t>& costmap, double res, const pcl::PointXYZ& orig):
                resolution(res) {
                    cm = costmap;
                    origin = cv::Point2d(orig.x, orig.y);
                    height = cm.size().height;
                    width = cm.size().width;
                    bbox = get2dBBoxFromImgInfo(origin.x, origin.y, height, width, resolution);
                }
        CostMap(const cv::Mat_<int8_t>& costmap, double res, const cv::Point2d& orig):
                resolution(res), origin(orig) {
                    cm = costmap;
                    height = cm.size().height;
                    width = cm.size().width;
                    bbox = get2dBBoxFromImgInfo(origin.x, origin.y, height, width, resolution);
                }
    };

    struct CandidateResult {
        double ig;
        double cost;
        pcl::PointXYZ coords;
        CandidateResult(): ig(0), cost(1000), coords(pcl::PointXYZ(0,0,0)) {};
        CandidateResult(double ig, double cost, const pcl::PointXYZ& coords) :
            ig(ig), cost(cost), coords(coords) {};

    };

    //convert PointXYZP to PointXYZ
    pcl::PointXYZ point3D(const PointXYZP& p) 
        {return pcl::PointXYZ(p.x,p.y,p.z);};

    //All the world coordinates are expressed with xmin,ymin = bottom left
    //In the world, a pixel coordinate is the bottom left of the pixel
    //Same in the voxelgrid, because it is constructed that way
    //On the contrary, in the images, the pixel coordinate is the top left of the pixel
    //because we work on opencv images

    //Here, we want the coordinates of the pixels, 
    //as we have in the gridvox : origin of frame is bottom left, coords of pixel is bottomleft
    cv::Point2d getCoordsBLPixelBLCoords(const cv::Point2i& pix, const CostMap& costmap) {
        return cv::Point2d{pix.x*costmap.resolution + costmap.bbox.xmin, 
                           pix.y*costmap.resolution + costmap.bbox.ymin}; 
    };
    //In the costmap, the origin is top left, the coords of the bottom left of a pixel is:
    cv::Point2d getCoordsBLPixelTLCoords(const cv::Point2i& pix, const CostMap& costmap) {
        return cv::Point2d{pix.x*costmap.resolution + costmap.bbox.xmin, 
                           (costmap.height-pix.y-1)*costmap.resolution + costmap.bbox.ymin}; 
    };
    //In the costmap, the origin is top left, the coords of the top left of a pixel is:
    cv::Point2d getCoordsTLPixelTLCoords(const cv::Point2i& pix, const CostMap& costmap) {
        return cv::Point2d{pix.x*costmap.resolution + costmap.bbox.xmin, 
                           (costmap.height-pix.y)*costmap.resolution + costmap.bbox.ymin}; 
    };

    cv::Point2i metersToBLPixelBLCoords(const cv::Point2d& p, const CostMap& costmap) {
        cv::Point2i pix;
        pix.x = std::min(std::max(static_cast<int>((p.x-costmap.bbox.xmin)/costmap.resolution), 0), costmap.width-1);
        pix.y = std::min(std::max(static_cast<int>((p.y-costmap.bbox.ymin)/costmap.resolution), 1), costmap.height)-1;
        return pix;
    };
    cv::Point2i metersToTLPixelTLCoords(const cv::Point2d& p, const CostMap& costmap) {
        cv::Point2i pix;
        pix.x = std::min(std::max(static_cast<int>((p.x-costmap.bbox.xmin)/costmap.resolution), 0), costmap.width-1);
        pix.y = costmap.height -(std::min(std::max(static_cast<int>((p.y-costmap.bbox.ymin)/costmap.resolution)+1, 1), costmap.height));
        return pix;
    };
    cv::Point2i metersToBLPixelTLCoords(const cv::Point2d& p, const CostMap& costmap) {
        cv::Point2i pix;
        pix.x = std::min(std::max(static_cast<int>((p.x-costmap.bbox.xmin)/costmap.resolution), 0), costmap.width-1);
        pix.y = std::min(std::max(static_cast<int>((p.y-costmap.bbox.ymin)/costmap.resolution), 1), costmap.height)-1;
        return pix;
    };
    cv::Point2i metersToBLPixelBLCoords(const pcl::PointXYZ& p, const CostMap& costmap) {
        return metersToBLPixelBLCoords(cv::Point2d(p.x,p.y), costmap);
    };
    cv::Point2i metersToTLPixelTLCoords(const pcl::PointXYZ& p, const CostMap& costmap) {
        return metersToTLPixelTLCoords(cv::Point2d(p.x,p.y), costmap);
    };
    cv::Point2i metersToBLPixelTLCoords(const pcl::PointXYZ& p, const CostMap& costmap) {
        return metersToBLPixelTLCoords(cv::Point2d(p.x,p.y), costmap);
    };
    //Within the above functions, only the following are usefull!
    cv::Point2i metersToPix(const cv::Point2d& p, const CostMap& costmap) {
        return metersToTLPixelTLCoords(p, costmap);
    };
    cv::Point2i metersToPix(const pcl::PointXYZ& p, const CostMap& costmap) {
        return metersToTLPixelTLCoords(cv::Point2d(p.x, p.y), costmap);
    };
    cv::Point2d pixToMeters(const cv::Point2i& pix, const CostMap& costmap) {
        return getCoordsBLPixelTLCoords(pix, costmap);
    };

    std::vector<pcl::PointXYZ> getHorizontalPlaneCircleIntersect(const pcl::PointXYZ& center, double R, double res) {
        std::vector<pcl::PointXYZ> circle;
        //Bresenham Algorithm
        double x = 0;
        double y = R;
        double m = 5 -4*R;
        while (x <= y) {
            circle.push_back(pcl::PointXYZ((x  * res) + center.x, ( y * res) + center.y, center.z));
            circle.push_back(pcl::PointXYZ((y  * res) + center.x, ( x * res) + center.y, center.z));
            circle.push_back(pcl::PointXYZ((-x * res) + center.x, ( y * res) + center.y, center.z));
            circle.push_back(pcl::PointXYZ((-y * res) + center.x, ( x * res) + center.y, center.z));
            circle.push_back(pcl::PointXYZ((x  * res) + center.x, (-y * res) + center.y, center.z));
            circle.push_back(pcl::PointXYZ((y  * res) + center.x, (-x * res) + center.y, center.z));
            circle.push_back(pcl::PointXYZ((-x * res) + center.x, (-y * res) + center.y, center.z));
            circle.push_back(pcl::PointXYZ((-y * res) + center.x, (-x * res) + center.y, center.z));
            if ( m > 0 ) {
                y = y - 1;
                m = m - 8 * y;
            }
            x = x + 1;
            m = m + 8*x + 4;
        }
        return circle;
    }



};

#endif