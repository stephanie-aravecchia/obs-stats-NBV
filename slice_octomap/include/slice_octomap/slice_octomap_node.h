#ifndef SLICE_OCTOMAP_H
#define SLICE_OCTOMAP_H
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "octomap/octomap.h"
#include "octomap/OcTree.h"
#include <vector>
#include <array>
#include "comparator/dataset_info.h"

//Takes an octomap as an input, and slice it into images
class SliceOctomap {
    protected:
        octomap::OcTree* tree_;
        std::string fname_;
        std::string outdir_;
        std::vector<cv::Mat_<uint8_t> > mat_vector_;
        std::vector<std::string> img_list_;
        double resolution_;
        double occ_threshold_;
        cv::Size img_size_;

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
            BBox(const octomap::point3d& point, double leafSize): 
                xmin(point.x() - leafSize/2), xmax(point.x() + leafSize/2), 
                ymin(point.y() - leafSize/2), ymax(point.y() + leafSize/2), 
                zmin(point.z() - leafSize/2), zmax(point.z() + leafSize/2) {}
        };
        struct PixBBox {
            unsigned int xmin;
            unsigned int xmax;
            unsigned int ymin;
            unsigned int ymax;
            unsigned int zmin;
            unsigned int zmax;
            PixBBox() {}
            PixBBox(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax):
                xmin(static_cast<unsigned int>(xmin)), xmax(static_cast<unsigned int>(xmax)), 
                ymin(static_cast<unsigned int>(ymin)), ymax(static_cast<unsigned int>(ymax)), 
                zmin(static_cast<unsigned int>(zmin)), zmax(static_cast<unsigned int>(zmax)) {}
        };
        
        BBox bbox_;
        //initialize a vector of cv::Mat, size of the vector depend on the height of the octree and the resolution,
        //size of the cv::Mat depends on the lenght and widht of the octree and the resolution
        void initializeMatVector();
        //iterate over all the leaves, if a voxel is occupied, draw a rectange in the appropriate cv::Mat
        void putOccupiedLeavesInImg();
        //normalize the occupancy to use to color the rectangle
        uint8_t getNormalizedOccupancy(double) const;
        //get the indices of the slice images in which to draw the rectangles, from the height of the voxel
        std::array<size_t, 2> getSliceIndexes(const BBox&) const;
        void bboxToPix(const BBox&, PixBBox&) const;
        //save the created slice images to disk
        void saveSliceImgs();
        void saveBBoxFile() const;
        void saveImgListFile() const;

    public:
        //constructor with the size of the box of the octree read directly from it
        SliceOctomap(double, double, std::string, std::string);
        //constructor with an offset in x and y (space ignored)
        SliceOctomap(double, double, double, std::string, std::string);
        //Constructor with a target BBOX to constrain the slicer to this BBOX
        SliceOctomap(double, double, const ComparatorDatatypes::BBox&, std::string, std::string);
        //slice the octree
        void slice();
};

#endif