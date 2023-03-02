#include "comparator/dataset_info.h"
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

DatasetInfo::DatasetInfo(std::string base_dir, std::string xp) {
    //First, we define the paths
    std::string gt_base_dir = base_dir + "ground_truth_stl/img/" + xp + "/";
    std::string ot_base_dir = base_dir + "3d_grid/img/" + xp + "/";
    std::string gtlist =  gt_base_dir + "list.txt";
    std::string otlist = ot_base_dir + "list.txt";
    std::string gtbbox =  gt_base_dir + "bbox";
    std::string otbbox = ot_base_dir + "bbox";

    //1. we load the dataset metric parameters in bbox
    getBBoxAndResFromDatasetParam(gtbbox, gtBbox_, gt_resolution_);
    getBBoxAndResFromDatasetParam(otbbox, otBbox_, ot_resolution_);
    //2. we load the dataset pixel info in pixbox
    getPixBBoxFromImgList(gtlist, gt_base_dir, gtPixBbox_);
    getPixBBoxFromImgList(otlist, ot_base_dir, otPixBbox_);
    //we calculate the mapping beween the two datasets,
    //first, we check that the resolutions are the same:
    if (!(fabs(gt_resolution_ - ot_resolution_)  <= 1e-6)){
        cerr << "gt and ot resolutions don't match, aborting" << endl;
        exit(EXIT_FAILURE);
    }
    //We compute the intersection
    setIntersectionBBoxes();
    //and now, we calculate the offset, based on the coordinates of xmin, ymin, zmin in the real box
    //we want to know how many pixels we will have to discard when reading the images to be in the intersctionBBox
    getOffsets(gtBbox_, intersectionBbox_, gt_offsets_);
    getOffsets(otBbox_, intersectionBbox_, ot_offsets_);
    assert((gt_offsets_.x + intersectionPixBbox_.xmax*gt_resolution_) < gtPixBbox_.xmax);
    assert((gt_offsets_.y + intersectionPixBbox_.ymax*gt_resolution_) < gtPixBbox_.ymax);
    assert((gt_offsets_.z + intersectionPixBbox_.zmax*gt_resolution_) < gtPixBbox_.zmax);
    assert((ot_offsets_.x + intersectionPixBbox_.xmax*gt_resolution_) < otPixBbox_.xmax);
    assert((ot_offsets_.y + intersectionPixBbox_.ymax*gt_resolution_) < otPixBbox_.ymax);
    assert((ot_offsets_.z + intersectionPixBbox_.zmax*gt_resolution_) < otPixBbox_.zmax);
    cout << "Dataset Info loaded completely." << endl << endl;
}
void DatasetInfo::getOffsets(const ComparatorDatatypes::BBox& outer, const ComparatorDatatypes::BBox& inner, ComparatorDatatypes::Offsets& offset) const{
    offset.x = round((inner.xmin - outer.xmin)/gt_resolution_); 
    offset.y = round((inner.ymin - outer.ymin)/gt_resolution_); 
    offset.z = round((inner.zmin - outer.zmin)/gt_resolution_);
}
ComparatorDatatypes::BBox DatasetInfo::getDatasetMetricsBbox() const { 
    return intersectionBbox_;
}
ComparatorDatatypes::PixBBox DatasetInfo::getDatasetPixBbox() const { 
    return intersectionPixBbox_;
}
double  DatasetInfo::getDatasetResolution() const {
    return gt_resolution_;
}
ComparatorDatatypes::Offsets DatasetInfo::getDatasetGTOffsets() const {
    return gt_offsets_;
}
ComparatorDatatypes::Offsets DatasetInfo::getDatasetOTOffsets() const {
    return ot_offsets_;
}

void DatasetInfo::setIntersectionBBoxes() {
    //the metric bbox:
    intersectionBbox_.xmin = max(gtBbox_.xmin, otBbox_.xmin);
    intersectionBbox_.ymin = max(gtBbox_.ymin, otBbox_.ymin);
    intersectionBbox_.zmin = max(gtBbox_.zmin, otBbox_.zmin);
    intersectionBbox_.xmax = min(gtBbox_.xmax, otBbox_.xmax);
    intersectionBbox_.ymax = min(gtBbox_.ymax, otBbox_.ymax);
    intersectionBbox_.zmax = min(gtBbox_.zmax, otBbox_.zmax);
    cout << "The metric box is set to: " << intersectionBbox_ << endl;
    //and now, the pixel bbox:
    intersectionPixBbox_.xmin = 0;
    intersectionPixBbox_.ymin = 0;
    intersectionPixBbox_.zmin = 0;
    intersectionPixBbox_.xmax = (intersectionBbox_.xmax - intersectionBbox_.xmin)/gt_resolution_ ;
    intersectionPixBbox_.ymax = (intersectionBbox_.ymax - intersectionBbox_.ymin)/gt_resolution_ ;
    intersectionPixBbox_.zmax = (intersectionBbox_.zmax - intersectionBbox_.zmin)/gt_resolution_ ;
    cout << "The pixel box is set to: " << intersectionPixBbox_ << endl;

}
void DatasetInfo::getBBoxAndResFromDatasetParam(const std::string& boxfile, ComparatorDatatypes::BBox& bbox, double& res) const {
    ifstream file(boxfile, ios::in);
    if (!file) {
        cerr << "cannot open " << boxfile << endl;
        exit(EXIT_FAILURE);
    }
    //discard first line:
    std::string line;
    getline(file, line,'\n');
    vector<double> vals;
    double val;
    while(file >> val) {
        vals.push_back(val);
    }
    bbox.xmin = vals.at(0);
    bbox.ymin = vals.at(1);
    bbox.zmin = vals.at(2);
    bbox.xmax = vals.at(3);
    bbox.ymax = vals.at(4);
    bbox.zmax = vals.at(5);
    res = vals.at(6);
    cout << "dataset loaded: " << bbox << ", res: " << res << endl;

}
void DatasetInfo::getPixBBoxFromImgList(std::string& imglist, std::string& dir, ComparatorDatatypes::PixBBox& pixbox) const{
    cv::Mat_<uint8_t> img;
    int nslices;
    loadOneImgAndGetNSlices(img, nslices, imglist, dir);
    pixbox.xmin = 0;
    pixbox.xmax = img.size().width;
    pixbox.ymin = 0;
    pixbox.ymax = img.size().height;
    pixbox.zmin = 0;
    pixbox.zmax = nslices;
}

void DatasetInfo::loadOneImgAndGetNSlices(cv::Mat_<uint8_t>& img , int& nslices, const std::string& imglist, const std::string& dir) const{
    ifstream file(imglist, ios::in);
    if (!file) {
        cerr << "cannot open " << imglist << endl;
        exit(EXIT_FAILURE);
    }
    string fname;
    //we load the first image only
    file >> fname;
    //we count the number of lines in imglist:
    nslices = 1;
    string tmp;
    while (file >> tmp) {
        nslices++;
    }
    img = cv::imread(dir+fname, cv::IMREAD_UNCHANGED);
    if (tmp.empty()) {
        cerr << dir+fname << " is empty, check the paths." << endl;
        exit(EXIT_FAILURE);
    }
    cout << dir + fname << " loaded." << endl;
    file.close();
}
