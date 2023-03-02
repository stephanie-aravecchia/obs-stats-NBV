#ifndef DATASET_INFO_H
#define DATASET_INFO_H
#include "comparator/datatypes.h"
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

//This class computes and stores general information of the ground truth and recontruction datasets
class DatasetInfo {
    protected:
        double gt_resolution_;
        double ot_resolution_;
        ComparatorDatatypes::PixBBox gtPixBbox_;
        ComparatorDatatypes::PixBBox otPixBbox_;
        ComparatorDatatypes::PixBBox intersectionPixBbox_;
        ComparatorDatatypes::BBox gtBbox_;
        ComparatorDatatypes::BBox otBbox_;
        ComparatorDatatypes::BBox intersectionBbox_;
        ComparatorDatatypes::Offsets gt_offsets_;
        ComparatorDatatypes::Offsets ot_offsets_;
        
        //set the values of the max bbox including both gt and rec dataset
        //both in pixel and meters
        void setIntersectionBBoxes();

        //get the offset in pixel of a dataset with respect to a bounding box
        //first is the outer bbox in meters, second is the bbox from the dataset
        void getOffsets(const ComparatorDatatypes::BBox&, const ComparatorDatatypes::BBox&, ComparatorDatatypes::Offsets&) const;
        
        //get the bounding and the resolution of a dataset from the file containing this information 
        void getBBoxAndResFromDatasetParam(const std::string&, ComparatorDatatypes::BBox&, double&) const;

        //get the pixel bounding box of a dataset from the actual images of the dataset
        void getPixBBoxFromImgList(std::string&, std::string&, ComparatorDatatypes::PixBBox&) const;
        
        //load one image from the dataset to get its size, and load the image list to get the number of images 
        void loadOneImgAndGetNSlices(cv::Mat_<uint8_t>& img, int& nslices, const std::string& imglist, const std::string& dir) const;
    
    public:
        DatasetInfo();
        DatasetInfo(std::string base_dir, std::string xp);

        double getDatasetResolution() const;
        ComparatorDatatypes::PixBBox getDatasetPixBbox() const;
        ComparatorDatatypes::BBox getDatasetMetricsBbox() const;
        ComparatorDatatypes::Offsets getDatasetGTOffsets() const;
        ComparatorDatatypes::Offsets getDatasetOTOffsets() const;
};

#endif