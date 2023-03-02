#ifndef GRID_TOOLS_H
#define GRID_TOOLS_H

#include <opencv2/opencv.hpp>
#include <pcl/point_types.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include "obs_stat/stat_point_type.h"

//Implementation of various useful functions dealing with 3D grids manipulation

namespace grid_tools {
    //Metric BBox
    struct BBox {
        double xmin;
        double xmax;
        double ymin;
        double ymax;
        double zmin;
        double zmax;
        BBox() {}
        BBox(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax):
            xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax), zmin(zmin), zmax(zmax) {}
    };
    
    //Pixel bbox (unsigned int)      
    struct PixBBox {
        unsigned int xmin;
        unsigned int xmax;
        unsigned int ymin;
        unsigned int ymax;
        unsigned int zmin;
        unsigned int zmax;
        PixBBox() {}
        PixBBox(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax):
            xmin(static_cast<unsigned int>(xmin)), xmax(static_cast<unsigned int>(xmax)), 
            ymin(static_cast<unsigned int>(ymin)), ymax(static_cast<unsigned int>(ymax)), 
            zmin(static_cast<unsigned int>(zmin)), zmax(static_cast<unsigned int>(zmax)) {}
        PixBBox(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax):
            xmin(static_cast<unsigned int>(xmin)), xmax(static_cast<unsigned int>(xmax)), 
            ymin(static_cast<unsigned int>(ymin)), ymax(static_cast<unsigned int>(ymax)), 
            zmin(static_cast<unsigned int>(zmin)), zmax(static_cast<unsigned int>(zmax)) {}
    };
    //Metric 2D BBOX
    struct BBox2d {
        double xmin;
        double xmax;
        double ymin;
        double ymax;
        BBox2d() {}
        BBox2d(double xmin, double xmax, double ymin, double ymax):
            xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax) {}
    };
    
    typedef std::vector< std::vector < std::vector <PointStat> > > StatGrid;

    typedef std::vector <std::vector <std::vector <double> > > VoxelGrid;
    
    typedef std::vector< std::vector < std::vector <bool> > > MaskGrid;
    
    struct GridSize {
        size_t nx;
        size_t ny;
        size_t nz;
        GridSize() {}
        GridSize(size_t nx, size_t ny, size_t nz):
            nx(nx), ny(ny), nz(nz) {}
    }; 

    //check if a point is in a BBox, limits included    
    bool isInBBox(const pcl::PointXYZ& p, const BBox& b) {
        if ((p.x < b.xmin) || (p.x) > b.xmax) return false;
        if ((p.y < b.ymin) || (p.y) > b.ymax) return false;
        if ((p.z < b.zmin) || (p.z) > b.zmax) return false;
        return true;
    };
    //check if a voxel is in the voxelgrid
    template<typename T>
    bool isVoxInGrid(const cv::Point3i& p, const T& grid) {
        size_t nx = grid.size();
        assert(nx>0);
        size_t ny = grid[0].size();
        assert(ny>0);
        size_t nz = grid[0][0].size();
        assert(nz>0);
        if ((p.x < 0) || (static_cast<size_t>(p.x) >= nx)) return false;
        if ((p.y < 0) || (static_cast<size_t>(p.y) >= ny)) return false;
        if ((p.z < 0) || (static_cast<size_t>(p.z) >= nz)) return false;
        return true; 
    }; 
    //resize a voxel grid to nx, ny, nz
    template<typename T>
    void resizeGrid(T& grid, size_t nx, size_t ny, size_t nz) {
        grid.resize(nx);
        for (size_t ix{0}; ix < nx; ++ix) {
            grid[ix].resize(ny);
            for (size_t iy{0}; iy < ny; ++iy) {
                grid[ix][iy].resize(nz);
            }
        }
    };
    template<typename T>
    void resizeGrid(T& grid, GridSize sz) {
        resizeGrid(grid, sz.nx, sz.ny, sz.nz);
    };
    template<typename T, typename U>
    void resizeAndInitializeGrid(T& grid, size_t nx, size_t ny, size_t nz, const U& init) {
        grid.resize(nx);
        for (size_t ix{0}; ix < nx; ++ix) {
            grid[ix].resize(ny);
            for (size_t iy{0}; iy < ny; ++iy) {
                grid[ix][iy].resize(nz);
                std::fill(grid[ix][iy].begin(), grid[ix][iy].end(), init);
            }
        }

    }
    template<typename T, typename U>
    void resizeAndInitializeGrid(T& grid, GridSize sz, const U& init) {
        resizeAndInitializeGrid(grid, sz.nx, sz.ny, sz.nz, init);
    }
    template<typename T, typename U>
    void resetGrid(T& grid, size_t nx, size_t ny, size_t nz, const U& init) {
        for (size_t ix{0}; ix < nx; ++ix) {
            for (size_t iy{0}; iy < ny; ++iy) {
                std::fill(grid[ix][iy].begin(), grid[ix][iy].end(), init);
            }
        }

    }
    template<typename T, typename U>
    void resetGrid(T& grid, GridSize sz, const U& init) {
        resetGrid(grid, sz.nx, sz.ny, sz.nz, init);
    }

    //check if the distance between two points is smaller than a threshold
    //consider the two points are the same if true 
    bool isSamePoint(const pcl::PointXYZ& p, const pcl::PointXYZ& q) {
        double norm = sqrt((q.x-p.x)*(q.x-p.x) + (q.y-p.y)*(q.y-p.y) + (q.z-p.z)*(q.z-p.z));
        return (norm <= DIST_THRESHOLD_SAME_POINT) ? true: false;
    };
    double distance(const pcl::PointXYZ& p, const pcl::PointXYZ& q) {
        double norm = sqrt((q.x-p.x)*(q.x-p.x) + (q.y-p.y)*(q.y-p.y) + (q.z-p.z)*(q.z-p.z));
        return norm;
    };
    double distance(const cv::Point3d& p, const cv::Point3d& q) {
        double norm = sqrt((q.x-p.x)*(q.x-p.x) + (q.y-p.y)*(q.y-p.y) + (q.z-p.z)*(q.z-p.z));
        return norm;
    };
    double distance(const cv::Point2d& p, const cv::Point2d& q) {
        double norm = sqrt((q.x-p.x)*(q.x-p.x) + (q.y-p.y)*(q.y-p.y));
        return norm;
    };

    GridSize getGridSizeFromBBoxAndRes(const BBox& bbox, double res) {
        GridSize sz;
        sz.nx = floor((bbox.xmax - bbox.xmin)/res);
        sz.ny = floor((bbox.ymax - bbox.ymin)/res);
        sz.nz = floor((bbox.zmax - bbox.zmin)/res);
        return sz;
    };
    cv::Point3i metersToVox(const pcl::PointXYZ& p, const BBox& bbox, double res) {
        assert(isInBBox(p, bbox));
        cv::Point3i vox;
        GridSize sz  = getGridSizeFromBBoxAndRes(bbox, res);
        vox.x = std::min(std::max(static_cast<int>((p.x-bbox.xmin)/res), 0), static_cast<int>(sz.nx)-1);
        vox.y = std::min(std::max(static_cast<int>((p.y-bbox.ymin)/res), 0), static_cast<int>(sz.ny)-1);
        vox.z = std::min(std::max(static_cast<int>((p.z-bbox.zmin)/res), 0), static_cast<int>(sz.nz)-1);
        return vox;
    };
    pcl::PointXYZ voxToMeters(const cv::Point3i& vox, const BBox& bbox, double res) {
        return pcl::PointXYZ{static_cast<float>(vox.x*res + bbox.xmin), 
                             static_cast<float>(vox.y*res + bbox.ymin), 
                             static_cast<float>(vox.z*res + bbox.zmin)};
    };
    pcl::PointXYZ voxToMetersCenter(const cv::Point3i& vox, const BBox& bbox, double res) {
        return pcl::PointXYZ{static_cast<float>(vox.x*res + bbox.xmin + res/2.), 
                             static_cast<float>(vox.y*res + bbox.ymin + res/2.), 
                             static_cast<float>(vox.z*res + bbox.zmin + res/2.)};
    };
    //Gets the number of angular sectors encoded in the uint32_t
    int getNAnglesFromV(uint32_t v) {
        int n_angles = 0;
        for (uint8_t i=0; i<32; i++) {
            if (v & (1 << i)) {
                ++n_angles;
            }
        }
        return n_angles;
    };

};

std::ostream& operator<<(std::ostream& out, const grid_tools::BBox& bbox) {
    out << "xmin: " << bbox.xmin << " " 
        << "xmax: " << bbox.xmax << " " 
        << "ymin: " << bbox.ymin << " " 
        << "ymax: " << bbox.ymax << " " 
        << "zmin: " << bbox.zmin << " " 
        << "zmax: " << bbox.zmax;
    return out;
};
std::ostream& operator<<(std::ostream& out, const grid_tools::BBox2d& bbox) {
    out << "xmin: " << bbox.xmin << " " 
        << "xmax: " << bbox.xmax << " " 
        << "ymin: " << bbox.ymin << " " 
        << "ymax: " << bbox.ymax; 
    return out;
};
std::ostream& operator<<(std::ostream& out, const grid_tools::GridSize& gz) {
    out << "nx: " << gz.nx << " " 
        << "ny: " << gz.ny << " " 
        << "nz: " << gz.nz; 
    return out;
};

bool operator==(const grid_tools::GridSize& a, const grid_tools::GridSize& b){
    if (a.nx != b.nx) return false;
    if (a.ny != b.ny) return false;
    if (a.nz != b.nz) return false;
    return true;
};

 
#endif