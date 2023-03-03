#include "comparator/compare_rec_gt.h"
#include "comparator/compute_metrics.h"
#include "comparator/dataset_info.h"
#include <ros/ros.h>
#include <assert.h>
#include <thread>
#include <algorithm>
#include <signal.h>
#include <unistd.h>

using namespace std;

ComparatorRecGt::ComparatorRecGt(const std::string& dir, const std::string& xp, const std::string& outputdir, double target_res, double img_res,
        ComparatorDatatypes::PixBBox pixbox, ComparatorDatatypes::BBox metricsbox,
        ComparatorDatatypes::Offsets gt_offset, ComparatorDatatypes::Offsets rec_offset, double ground_thres,
        bool is_dist_metrics, bool is_dkl_metrics, bool is_ot_metrics, int nthreads, double ot_reg, int ot_maxiter, double occ_thres, double ot_stopthres,
        double divergence_to_unknown_thres, double non_observed_vi, double non_observed_cov, double non_observed_acc, double non_observed_l1,
        double non_observed_dkl, double non_observed_wd,
        bool unit_test, const std::string& test_filename):
            base_dir_(dir), xp_name_(xp), output_dir_(outputdir), target_res_(target_res), img_res_(img_res), 
            space_pix_bbox_(pixbox), space_bbox_(metricsbox), rec_offset_(rec_offset), gt_offset_(gt_offset), ground_thres_(ground_thres),
            is_dist_metrics_(is_dist_metrics), is_dkl_metrics_(is_dkl_metrics), is_ot_metrics_(is_ot_metrics),n_threads_(nthreads),
            ot_reg_(ot_reg), ot_maxiter_(ot_maxiter), occupancy_thres_(occ_thres), ot_stopthres_(ot_stopthres),
            divergence_to_unknown_thres_(divergence_to_unknown_thres), non_observed_vi_(non_observed_vi), non_observed_cov_(non_observed_cov),
            non_observed_acc_(non_observed_acc), non_observed_l1_(non_observed_l1), non_observed_dkl_(non_observed_dkl),
            non_observed_wd_(non_observed_wd), unit_test_(unit_test) {
    
    //First, we set all we need to construct our 3d grid
    img_width_ = space_pix_bbox_.xmax-space_pix_bbox_.xmin;
    img_height_ = space_pix_bbox_.ymax-space_pix_bbox_.ymin;
    box_size_ = static_cast<int>(target_res_/img_res_);
    //if we cannot calculate results at target_res, abort:
    if (!(fabs(target_res_/img_res_ - box_size_) <= 1e-9)) {
        cout << "target res: " << target_res_ << " cannot be achieved (output resolution " << box_size_*img_res_ << "). Abort. "<< endl;
        exit(EXIT_FAILURE);
    }
    if (target_res_/img_res_ < 1) {
        cerr << "Voxel size set to a value <1, change the target resolution. Aborting." << endl;
        exit(EXIT_FAILURE);
    }
    assert(img_height_>0); 
    assert(img_width_>0); 
    assert(img_height_ -box_size_ >=0);
    assert(img_width_ -box_size_ >=0);
    assert(box_size_>1);
    const int mat_sz[3] = {img_width_,img_height_,box_size_};
    cout << "mat_sz: " << mat_sz[0] << ", "<< mat_sz[1] << ", "<< mat_sz[2] <<endl;
    cv::Mat_<uint8_t> gtMat= cv::Mat_<uint8_t>(3, mat_sz);
    cv::Mat_<uint8_t> recMat= cv::Mat_<uint8_t>(3, mat_sz);
    cout << "3DMat gtMat_ : dim: "<< gtMat.dims << ", width: "<< gtMat.size[0] << ", height: "<< gtMat.size[1] 
            << ", z : "<< gtMat.size[2] << endl;
    cout << "3DMat recMat_ : dim: "<< recMat.dims << ", width: "<< recMat.size[0] << ", height: "<< recMat.size[1] 
            << ", z : "<< recMat.size[2] << endl;
    cout << "target res: " << target_res_ << ", output resolution " << box_size_*img_res_ << endl;
    cout << "box size: " << box_size_<< endl;
   
    //now our space is defined, we prepare all we need to load the pixels in the images of the datasets that correspond to the same coordinates
    std::string gt_base_dir = base_dir_ + "ground_truth_stl/img/" + xp_name_ + "/";
    std::string ot_base_dir = base_dir_ + "3d_grid/img/" + xp_name_ + "/";
    std::string gtlist =  gt_base_dir + "list.txt";
    std::string otlist = ot_base_dir + "list.txt";
    cout << "gt list selected: " << gtlist << endl;
    cout << "ot list selected: " << otlist << endl;
    cout << "output dir: " << output_dir_ << endl;
    //we need to deal with the index at which the two img lists are aligned, so we start loading at the offset of the other one:
    assert(rec_offset_.x>=0);
    assert(rec_offset_.y>=0);
    assert(rec_offset_.z>=0);
    assert(gt_offset_.x>=0);
    assert(gt_offset_.y>=0);
    assert(gt_offset_.z>=0);
    int ot_load_slice_start = 0;
    int gt_load_slice_start = 0;
    //this part is to deal with the z offset of our dataset:
    assert(!((gt_offset_.z > 0) && (rec_offset_.z > 0)));
    getImgList(gt_imglist_, gtlist, gt_base_dir, gt_offset_.z); 
    getImgList(rec_imglist_, otlist, ot_base_dir, rec_offset_.z); 
    cout << "Ot offset: " << rec_offset_.x << ", " << rec_offset_.y << ", " << rec_offset_.z
        << ". Gt offset: " << gt_offset_.x << ", " << gt_offset_.y << ", " << gt_offset_.z << endl;
    n_slices_ = static_cast<int>(min(static_cast<unsigned int>(min(rec_imglist_.size(), gt_imglist_.size())), space_pix_bbox_.zmax-static_cast<int>(ground_thres_/img_res_)))+1;
    cout << "img_width: " << img_width_ << ", " << "img_height: " << img_height_ << ", nslices: " << n_slices_ << endl; 
    cout << "space_bbox_ is set to: " << space_bbox_ << endl;
    assert(n_slices_>0);
    assert(n_slices_ -box_size_ >=0);
    
    //and finally, we set a constrained_bbox  and a starting index to start the comparaison at a coordinate which is a multiple of re
    //i.e if xmin is -5.28 and res .1, we start the comparaison a -5.20
    //this is mandatory for the comparaison to the observation in next stages (with obs_stat)
    istart_ = getStartingIndex(space_bbox_.xmin);
    jstart_ = getStartingIndex(space_bbox_.ymin);
    kstart_ = getStartingIndex(space_bbox_.zmin + ground_thres_);
    cout << "istart_: " << istart_ 
         << ", jstart_: " << jstart_ 
         << ", kstart_: " << kstart_ << endl;
    assert(istart_ >=0);
    assert(jstart_ >=0);
    assert(kstart_ >=0);
    assert(istart_ <=img_width_-box_size_);
    assert(jstart_ <=img_height_-box_size_);
    assert(kstart_ <=n_slices_-box_size_);

    //now, we know how many boxes we will compare:
    nx_ = floor((img_width_-1-istart_)/box_size_);
    ny_ = floor((img_height_-1-jstart_)/box_size_);
    nz_ = floor((n_slices_-1-kstart_)/box_size_);
    tot_box_ = nx_*ny_*nz_;
    cout << "nx " << nx_ << ", ny " << ny_ << ", nz " << nz_ << endl;
    constrained_bbox_.xmin = space_bbox_.xmin + getStartingDistance(space_bbox_.xmin);
    constrained_bbox_.ymin = space_bbox_.ymin + getStartingDistance(space_bbox_.ymin);
    constrained_bbox_.zmin = space_bbox_.zmin +(ot_load_slice_start + gt_load_slice_start)*img_res_ +ground_thres_+ getStartingDistance(space_bbox_.zmin + ground_thres_+ (ot_load_slice_start + gt_load_slice_start)*img_res_);
    //Below may fail with rounding, we need to make sure that min + n*target_res stays in space
    constrained_bbox_.xmax = ((constrained_bbox_.xmin + nx_*target_res_) <= space_bbox_.xmax) ? (constrained_bbox_.xmin + (nx_)*target_res_) : (constrained_bbox_.xmin + (nx_-1)*target_res_ );
    constrained_bbox_.ymax = ((constrained_bbox_.ymin + ny_*target_res_) <= space_bbox_.ymax) ? (constrained_bbox_.ymin + (ny_)*target_res_) : (constrained_bbox_.ymin + (ny_-1)*target_res_ );
    constrained_bbox_.zmax = ((constrained_bbox_.zmin + nz_*target_res_) <= space_bbox_.zmax) ? (constrained_bbox_.zmin + (nz_)*target_res_) : (constrained_bbox_.zmin + (nz_-1)*target_res_ );
    assert(constrained_bbox_.xmax>constrained_bbox_.xmin);
    assert(constrained_bbox_.ymax>constrained_bbox_.ymin);
    assert(constrained_bbox_.zmax>constrained_bbox_.zmin);

    //constrained bbox should be included in space_box
    cout << "constrainded_bbox_ is set to: " << constrained_bbox_ << endl;
    assert(isQueryBoxIncludedInRef(constrained_bbox_, space_bbox_));
    
    cout << "We will start comparing " << tot_box_ << " boxes." << endl;
    
    if (is_ot_metrics_) {
        ot_metrics_ = ComputeOTMetrics(box_size_, box_size_, box_size_, ot_reg_, ot_maxiter_, ot_stopthres_);
        ot_metrics_.setCostMatrixSquareDist();
    }
    if (unit_test_) {
        std::string test_file_path = base_dir_ + "res/obs/" + xp_name_ + "/sample_list_to_test.csv";
    }
}

bool ComparatorRecGt::isQueryBoxIncludedInRef(const ComparatorDatatypes::BBox& q, const ComparatorDatatypes::BBox& r) const {
    if (((q.xmin < r.xmin) || (q.ymin < r.ymin)) || (q.zmin < r.zmin)) {
        return false;
    } else if (((q.xmax > r.xmax) || (q.ymax > r.ymax)) || (q.zmax > r.zmax)) {
        return false;
    }
    return true;
}
bool ComparatorRecGt::isQueryBoxIncludedInRef(const ComparatorDatatypes::PixBBox& q, const ComparatorDatatypes::PixBBox& r) const {
    if (((q.xmin < r.xmin) || (q.ymin < r.ymin)) || (q.zmin < r.zmin)) {
        return false;
    } else if (((q.xmax > r.xmax) || (q.ymax > r.ymax)) || (q.zmax > r.zmax)) {
        return false;
    }
    return true;
}

//Reshape MatVect to contain nz_ 3DMat of shape (img_width_, img_height_, box_size_)
void ComparatorRecGt::reshapeMatVect() {
    gtMatVect_.clear();
    recMatVect_.clear();
    gtMatVect_.resize(nz_);
    recMatVect_.resize(nz_);
    const int mat_sz[3] = {img_width_,img_height_,box_size_};
    for (int i{0}; i < nz_; ++i) {
        gtMatVect_.at(i) = cv::Mat_<uint8_t>(3, mat_sz);
        recMatVect_.at(i) = cv::Mat_<uint8_t>(3, mat_sz);
    }
}
void ComparatorRecGt::loadGTOnlyDataset() {
    reshapeMatVect();
    int thick_slice_number = 0;
    for (int img_index{kstart_}; img_index < n_slices_ - box_size_; img_index+=box_size_) {
        load3DMat(gtMatVect_.at(thick_slice_number), gt_imglist_, gt_offset_, img_index);
        ++thick_slice_number;
    }
}
//For all the elements of MatVect_, fill the 3DMatrix with the values of the corresponding images
void ComparatorRecGt::loadCompleteDataset() {
    reshapeMatVect();
    int thick_slice_number = 0;
    for (int img_index{kstart_}; img_index < n_slices_ - box_size_; img_index+=box_size_) {
        load3DMat(recMatVect_.at(thick_slice_number), rec_imglist_, rec_offset_, img_index);
        load3DMat(gtMatVect_.at(thick_slice_number), gt_imglist_, gt_offset_, img_index);
        ++thick_slice_number;
    }
}
//we loop on the complete dataset to get a list of the boxes to process
void ComparatorRecGt::getListOfBoxes() {
    listOfBoxes_.clear();
    int thick_slice_number = 0;
    assert(thick_slice_number < nz_);
    for (int start_slice{kstart_}; start_slice < n_slices_ - box_size_; start_slice+=box_size_) {
        for (int i{istart_}; i < img_width_-box_size_; i+=box_size_){
            for (int j{jstart_}; j < img_height_-box_size_; j+=box_size_){
                vector<cv::Range> ranges{cv::Range(i, i+box_size_),cv::Range(j, j+box_size_),cv::Range(0,box_size_)};
                assert(ranges[0].start >= 0);
                assert(ranges[0].start < ranges[0].end);
                assert(ranges[0].end <= img_width_);
                assert(ranges[1].start >= 0);
                assert(ranges[1].start < ranges[1].end);
                assert(ranges[1].end <= img_height_);
                assert(ranges[2].start >= 0);
                assert(ranges[2].start < ranges[2].end);
                assert(ranges[2].end <= box_size_);
                ComparatorDatatypes::BoxToProcess box;
                box.index_in_vect = thick_slice_number;
                box.ranges = ranges;
                ComparatorDatatypes::BBox mbox;
                ComparatorDatatypes::PixBBox pixBox{static_cast<double>(ranges[0].start), static_cast<double>(ranges[0].end), static_cast<double>(ranges[1].start), static_cast<double>(ranges[1].end), 
                                                    static_cast<double>(start_slice), static_cast<double>(start_slice+box_size_)};
                bboxPixToMeters(pixBox, mbox);
                box.mbbox = mbox;
                listOfBoxes_.push_back(box);
            }
        }
        ++thick_slice_number;
    }
}
//we loop on the complete dataset to get a list of the boxes to process, the list of boxes is specific to the calculation of the limit of the WD on the dataset
void ComparatorRecGt::getListOfBoxesLimitOnly() {
    listOfBoxesLimitOnly_.clear();
    int thick_slice_number = 0;
    assert(thick_slice_number < nz_);
    assert(n_slices_>0);
    assert(n_slices_ -box_size_ >=0);
    for (int start_slice{kstart_}; start_slice < n_slices_ - box_size_; start_slice+=box_size_) {
        for (int i{istart_}; i < img_width_-box_size_; i+=box_size_){
            for (int j{jstart_}; j < img_height_-box_size_; j+=box_size_){
                vector<cv::Range> ranges{cv::Range(i, i+box_size_),cv::Range(j, j+box_size_),cv::Range(0,box_size_)};
                assert(ranges[0].start >= 0);
                assert(ranges[0].start < ranges[0].end);
                assert(ranges[0].end <= img_width_);
                assert(ranges[1].start >= 0);
                assert(ranges[1].start < ranges[1].end);
                assert(ranges[1].end <= img_height_);
                assert(ranges[2].start >= 0);
                assert(ranges[2].start < ranges[2].end);
                assert(ranges[2].end <= box_size_);
                ComparatorDatatypes::BoxToProcessLimitOnly box;
                box.index_in_vect = thick_slice_number;
                box.ranges = ranges;
                ComparatorDatatypes::BBox mbox;
                ComparatorDatatypes::PixBBox pixBox{static_cast<double>(ranges[0].start), static_cast<double>(ranges[0].end), static_cast<double>(ranges[1].start), static_cast<double>(ranges[1].end), 
                                                    static_cast<double>(start_slice), static_cast<double>(start_slice+box_size_)};
                bboxPixToMeters(pixBox, mbox);
                box.mbbox = mbox;
                listOfBoxesLimitOnly_.push_back(box);
            }
        }
        ++thick_slice_number;
    }
}

int ComparatorRecGt::getStartingIndex(double minval) const {
    return floor((ceil(minval / target_res_) * target_res_ - minval)/img_res_);
}
double ComparatorRecGt::getStartingDistance(double minval) const {
    return ceil(minval / target_res_) * target_res_ - minval;
}

void ComparatorRecGt::getImgList(std::vector<std::string >& imglist_vector, const std::string& imglist, 
        const std::string& dir, int load_slice_start) const {
    ifstream file(imglist, ios::in);
    if (!file) {
        cerr << "cannot open " << imglist << endl;
        exit(EXIT_FAILURE);
    }
    string fname;
    //we discard the first images before we can actually compare:
    for (int i{0}; i < load_slice_start; i++) {
        file >> fname;
    }
    imglist_vector.clear();
    while (file >> fname) {
        imglist_vector.push_back(dir + fname);
    }
    //check that the first path is correct, if not raise an error
    cv::Mat_<uint8_t> tmp = cv::imread(imglist_vector[0], cv::IMREAD_UNCHANGED);
    if (tmp.empty()) {
        cerr << imglist_vector[0] << " is empty, check the paths." << endl;
        exit(EXIT_FAILURE);
    }
    cout << "able to load first img: " << imglist_vector[0] << ", shape: width " << tmp.size().width << ", height " << tmp.size().height 
        << ", dataset contains: " << imglist_vector.size() << " images, and " << load_slice_start << "images were discarded, not comparable because of z offset." << endl;
    file.close();
}

void ComparatorRecGt::bboxPixToMeters(const ComparatorDatatypes::PixBBox& pixbox, ComparatorDatatypes::BBox& mbox) const {
    mbox.xmin = static_cast<double>(pixbox.xmin-istart_)*img_res_ + constrained_bbox_.xmin; 
    mbox.xmax = static_cast<double>(pixbox.xmax-istart_)*img_res_ + constrained_bbox_.xmin;
    mbox.ymin = static_cast<double>(pixbox.ymin-jstart_)*img_res_ + constrained_bbox_.ymin;
    mbox.ymax = static_cast<double>(pixbox.ymax-jstart_)*img_res_ + constrained_bbox_.ymin;
    mbox.zmin = static_cast<double>(pixbox.zmin-kstart_)*img_res_ + constrained_bbox_.zmin;
    mbox.zmax = static_cast<double>(pixbox.zmax-kstart_)*img_res_ + constrained_bbox_.zmin;
}

void ComparatorRecGt::load3DMat(cv::Mat_<uint8_t>& mat, const std::vector<std::string >& imglist, 
        const ComparatorDatatypes::Offsets& offset, int start_slice) const {
    for (int k{start_slice}; k < start_slice + box_size_; k++) {
        assert(k >= 0);
        assert(k < (int) imglist.size());
        assert(k -start_slice >=0);
        cv::Mat_<uint8_t> tmp = cv::imread(imglist[k], cv::IMREAD_GRAYSCALE);
        assert(tmp.size().height*tmp.size().width >0);
        assert(tmp.size().width >= img_width_ + offset.x);
        assert(tmp.size().height >= img_height_ + offset.y);
        for (int i{0}; i < img_width_; i++) {
            for (int j{0}; j < img_height_; j++) {
                assert(i >= 0);
                assert(j >= 0);
                //reading in the dataset:
                assert(i + offset.x< tmp.size().width);
                assert(j + offset.y < tmp.size().height);
                assert(i + offset.x >= 0);
                assert(j + offset.y >= 0);
                //writing in the 3D Mat
                assert(i < img_width_);
                assert(j < img_height_);
                assert(k-start_slice < box_size_);
                //remainder: opencv: tmp() refers to tmp(row, col)
                mat.at<uint8_t>(i, j, k-start_slice) = tmp.at<uint8_t>(j+offset.y,i+offset.x);
            }
        }
    }

}


void ComparatorRecGt::getVectorFromBox(const ComparatorDatatypes::BoxToProcessLimitOnly& box, 
            const cv::Mat_<uint8_t>& refMat, vector<double>& vect) const {
                BoxTools::getVectorFromBox(box.ranges, refMat, vect);
}
void ComparatorRecGt::getVectorFromBox(const ComparatorDatatypes::BoxToProcess& box, 
            const cv::Mat_<uint8_t>& refMat, vector<double>& vect) const {
                BoxTools::getVectorFromBox(box.ranges, refMat, vect);
}

bool ComparatorRecGt::isBoxNotEmptyWithLevel(const ComparatorDatatypes::BoxToProcessLimitOnly& box, const cv::Mat_<uint8_t>& refMat, double level) const {         
    return BoxTools::isBoxNotEmptyWithLevel(box.ranges, refMat, level);
}

bool ComparatorRecGt::isBoxOnlyZeros(const ComparatorDatatypes::BoxToProcessLimitOnly& box, const cv::Mat_<uint8_t>& refMat) const {         
    return BoxTools::isBoxOnlyZeros(box.ranges, refMat);
}

bool ComparatorRecGt::isBoxOnlyZeros(const ComparatorDatatypes::BoxToProcess& box, const cv::Mat_<uint8_t>& refMat) const {         
    return BoxTools::isBoxOnlyZeros(box.ranges, refMat);
}

double ComparatorRecGt::getBoxOccupancyRate(const ComparatorDatatypes::BoxToProcessLimitOnly& box, const cv::Mat_<uint8_t>& refMat) const {         
    return BoxTools::getBoxOccupancyRate(box.ranges, refMat);
}
double ComparatorRecGt::getBoxOccupancyRate(const ComparatorDatatypes::BoxToProcess& box, const cv::Mat_<uint8_t>& refMat) const {         
    return BoxTools::getBoxOccupancyRate(box.ranges, refMat);
}

bool ComparatorRecGt::isSampleDebug(const ComparatorDatatypes::BBox& mbox, double x, double y, double z) const {
    auto is_nearly_equal = [](double a, double b) {return (fabs(a-b)<=1e-9);};
    if ((is_nearly_equal(mbox.xmin,x) && is_nearly_equal(mbox.ymin,y)) && is_nearly_equal(mbox.zmin,z)) {
        return true;
    } else {
        return false;
    }

}
void ComparatorRecGt::unitTestBox(ComparatorDatatypes::BoxToProcess& box){
    //if we are in one of the selected boxes, we display the results step by step
    bool do_test = false;
    if (isSampleDebug(box.mbbox, 28, -1.5, 1.5)) {
        do_test = true;
    } else if (isSampleDebug(box.mbbox, 28, -1.5, 2)){
        do_test = true;
    } else if (isSampleDebug(box.mbbox, 28, -1.5, 2.5)){
        do_test = true;
    } else if (isSampleDebug(box.mbbox, 28, -1.5, 3)){
        do_test = true;
    }

    static int debug_n = 0;    

    if (do_test) {
        calcMetricFromBox(box);
        ComparatorDatatypes::PixBBox pixBox{static_cast<double>(box.ranges[0].start), static_cast<double>(box.ranges[0].end), 
                                        static_cast<double>(box.ranges[1].start), static_cast<double>(box.ranges[1].end), 
                                        static_cast<double>(box.index_in_vect), static_cast<double>(box.index_in_vect+box_size_)};
    
        std::cout << "box: " << box.mbbox.xmin << " " << box.mbbox.ymin << " "  << box.mbbox.zmin 
                << ", wasserstein_occ_remapped: " << box.metrics.wasserstein_occ_remapped
                << ", wasserstein_free_remapped: " << box.metrics.wasserstein_free_remapped
                << ", wasserstein_direct: " << box.metrics.wasserstein_direct
                << ", dkl: " << box.metrics.dkl
                << ", surface_coverage_1: " << box.metrics.surface_coverage_1
                << ", surface_coverage_2: " << box.metrics.surface_coverage_2
                << ", surface_coverage_3: " << box.metrics.surface_coverage_3
                << ", surface_coverage_4: " << box.metrics.surface_coverage_4
                << ", reconstruction_accuracy_1: " << box.metrics.reconstruction_accuracy_1
                << ", reconstruction_accuracy_2: " << box.metrics.reconstruction_accuracy_2
                << ", reconstruction_accuracy_3: " << box.metrics.reconstruction_accuracy_3
                << ", reconstruction_accuracy_4: " << box.metrics.reconstruction_accuracy_4
                << ", volumetric_information: " << box.metrics.volumetric_information
                << ", emptyness_l1: " << box.metrics.emptyness_l1
                << ", n_match_1: " << box.metrics.n_match_1
                << ", n_match_2: " << box.metrics.n_match_2
                << ", n_match_3: " << box.metrics.n_match_3
                << ", n_match_4: " << box.metrics.n_match_4
                << ", n_rec_points: " << box.metrics.n_rec_points
                << ", n_gt_points: " << box.metrics.n_gt_points
                << endl;
        save_sample_debug(box, debug_n);
        ++debug_n;
    }
    
}

void ComparatorRecGt::getNoisyGtVector(const ComparatorDatatypes::BoxToProcessLimitOnly& box, 
        const cv::Mat_<uint8_t>& refMat, vector<double>& vect, double uniform_noise_level, bool use_fixed_sigma, double sigma) const {
    BoxTools::getNoisyGtVector(box.ranges,refMat, vect, uniform_noise_level, use_fixed_sigma, sigma);
} 

//calculate the limit of the Wasserstein distance, with a reconstruction drawned randomly from gt
void ComparatorRecGt::calcLimitMetricFromBoxTestOnly(ComparatorDatatypes::BoxToProcessLimitOnly& box, std::ofstream& outfile) {
    vector<double> gt_vector;
    vector<double> rec_random_vector;
    vector<double> rec_noisygt_vector;
    getVectorFromBox(box, gtMatVect_.at(box.index_in_vect), gt_vector);
    //We will process the SAME vector several times, with different level of noise:
    //And no uniform noise
    std::vector<double> sigmas = {0,.2,.4,.6,.8,1.,1.2,1.4,1.6};
    std::vector<double> levels = {0,.05, .1, .15};
    for (double level : levels) {
        for (double sigma: sigmas) {
            cout << "uniform_noise_level: " << level << ", sigma:" << sigma << endl;
            getNoisyGtVector(box, gtMatVect_.at(box.index_in_vect), rec_noisygt_vector,level,true,sigma);
            ComputeOTMetrics::OccFreeWasserstein wasserstein_norm;
            outfile << setprecision(2);
            wasserstein_norm = ot_metrics_.normalizeAndGetSinkhornDistanceSignedNormalized(rec_noisygt_vector, gt_vector);
            outfile << box.mbbox.xmin << " " << box.mbbox.ymin<< " " << box.mbbox.zmin << " "
                    << sigma << " " << level << " " << wasserstein_norm.occ << " " << wasserstein_norm.free << endl;
            cout << "wd free: " << wasserstein_norm.free << ", wd occ: " << wasserstein_norm.occ << endl;
        }
    } 

}


//calculate the limit of the Wasserstein distance, with a reconstruction drawned randomly from gt
void ComparatorRecGt::calcLimitMetricFromBox(ComparatorDatatypes::BoxToProcessLimitOnly& box) {
    vector<double> gt_vector;
    vector<double> rec_random_vector;
    
    //get the ground-truth vector
    getVectorFromBox(box, gtMatVect_.at(box.index_in_vect), gt_vector);
    //draw a random reconstruction
    drawRandomDistri(gt_vector.size(), rec_random_vector); 

    //ideal_conv1
    //ideal_conv2
    //ideal_conv3 
    //random

    //first, we calculate the common info, and the metrics form the random reconstruction:
    int n_gt_points = ComputeMetrics::getNpoints(gt_vector, 0.1);//in gt, all points >0 are occupied    
    box.random.n_gt_points = n_gt_points;
    bool is_gt_empty = isBoxOnlyZeros(box, gtMatVect_.at(box.index_in_vect));
    if (is_gt_empty) {
        //we calculate L1 norm between random and gt
            box.random.emptyness_l1 = ComputeMetrics::getEmptynessL1(rec_random_vector);
    } else {
    //for each couple sigma / level AND random:
        // get noisyGTVector
        // compute the metrics
        //random
        computeLimitMetricsOnVect(box.random, rec_random_vector, gt_vector);
        vector<double> rec_noisy_gt_vector;
        //conv1:
        getNoisyGtVector(box, gtMatVect_.at(box.index_in_vect), rec_noisy_gt_vector, 
                noise_on_gt_[0].sigma, noise_on_gt_[0].ksize, noise_on_gt_[0].additionnal_uniform_noise);
        computeLimitMetricsOnVect(box.ideal_conv1, rec_noisy_gt_vector, gt_vector);
        rec_noisy_gt_vector.clear();
        //conv2
        getVectorFromBox(box, gtMatVect_.at(box.index_in_vect), gt_vector);
        getNoisyGtVector(box, gtMatVect_.at(box.index_in_vect), rec_noisy_gt_vector, 
                noise_on_gt_[1].sigma, noise_on_gt_[1].ksize, noise_on_gt_[1].additionnal_uniform_noise);
        computeLimitMetricsOnVect(box.ideal_conv2, rec_noisy_gt_vector, gt_vector);
        rec_noisy_gt_vector.clear();
        //conv3:
        getVectorFromBox(box, gtMatVect_.at(box.index_in_vect), gt_vector);
        getNoisyGtVector(box, gtMatVect_.at(box.index_in_vect), rec_noisy_gt_vector, 
                noise_on_gt_[2].sigma, noise_on_gt_[2].ksize, noise_on_gt_[2].additionnal_uniform_noise);
        computeLimitMetricsOnVect(box.ideal_conv3, rec_noisy_gt_vector, gt_vector);
    }

}
void ComparatorRecGt::computeLimitMetricsOnVect(ComparatorDatatypes::Metrics& metrics, std::vector<double>& rec_vector, std::vector<double>& gt_vector) {
    metrics.volumetric_information = ComputeMetrics::getVolumetricInformation(rec_vector);        
    metrics.dkl = ComputeMetrics::getKLDivergence(rec_vector, gt_vector);
    ComputeOTMetrics::OccFreeWasserstein wasserstein;
    wasserstein = ot_metrics_.normalizeAndGetSinkhornDistanceSigned(rec_vector, gt_vector);
    metrics.wasserstein_occ_remapped = wasserstein.occ;
    metrics.wasserstein_free_remapped = wasserstein.free;
}

//Compute the metrics based on the selection: distance and/or dkl and/or wassserstein
void ComparatorRecGt::calcMetricFromBox(ComparatorDatatypes::BoxToProcess& box) {

    vector<double> gt_vector;
    vector<double> rec_vector;
    getVectorFromBox(box, gtMatVect_.at(box.index_in_vect), gt_vector);
    getVectorFromBox(box, recMatVect_.at(box.index_in_vect), rec_vector);
    box.metrics.n_gt_points = ComputeMetrics::getNpoints(gt_vector, occupancy_thres_);    
    box.metrics.n_rec_points = ComputeMetrics::getNpoints(rec_vector, occupancy_thres_);    

    bool is_observed = BoxTools::isBoxObserved(rec_vector, divergence_to_unknown_thres_);
    bool is_gt_empty = isBoxOnlyZeros(box, gtMatVect_.at(box.index_in_vect));
    if (is_dist_metrics_) {
        if (is_observed) {
            box.metrics.volumetric_information = ComputeMetrics::getVolumetricInformation(rec_vector);        
            //We do that for all 4 couples (likelihood, registration distance) defined
            //coverage 1
            int i = 0;
            box.metrics.n_match_1 = ComputeMetrics::getBoxNMatch(box,gtMatVect_.at(box.index_in_vect), recMatVect_.at(box.index_in_vect),
                        cov_thresholds[i].reg_dist, cov_thresholds[i].likelihood, img_width_, img_height_, box_size_, img_res_);        
            box.metrics.surface_coverage_1 = ComputeMetrics::getSurfaceCoverage(box,gtMatVect_.at(box.index_in_vect), recMatVect_.at(box.index_in_vect),
                        cov_thresholds[i].reg_dist, cov_thresholds[i].likelihood, img_width_, img_height_, box_size_, img_res_);        
            box.metrics.reconstruction_accuracy_1 = ComputeMetrics::getReconstructionAccuracy(box,gtMatVect_.at(box.index_in_vect), recMatVect_.at(box.index_in_vect),
                        cov_thresholds[i].reg_dist, cov_thresholds[i].likelihood, img_width_, img_height_, box_size_, img_res_);        
            //Average Hausdorff Distance and Kappa:
            box.metrics.ahd_1 = ComputeMetrics::getAverageHausdorffDistance(box,gtMatVect_.at(box.index_in_vect), recMatVect_.at(box.index_in_vect),
                        cov_thresholds[i].likelihood, img_width_, img_height_, box_size_, img_res_);        
            box.metrics.kappa_1 = ComputeMetrics::getKappa(box,gtMatVect_.at(box.index_in_vect), recMatVect_.at(box.index_in_vect),
                        cov_thresholds[i].likelihood, img_width_, img_height_, box_size_, img_res_);        
            //coverage 2
            i = 1;
            box.metrics.n_match_2 = ComputeMetrics::getBoxNMatch(box,gtMatVect_.at(box.index_in_vect), recMatVect_.at(box.index_in_vect),
                        cov_thresholds[i].reg_dist, cov_thresholds[i].likelihood, img_width_, img_height_, box_size_, img_res_);        
            box.metrics.surface_coverage_2 = ComputeMetrics::getSurfaceCoverage(box,gtMatVect_.at(box.index_in_vect), recMatVect_.at(box.index_in_vect),
                        cov_thresholds[i].reg_dist, cov_thresholds[i].likelihood, img_width_, img_height_, box_size_, img_res_);        
            box.metrics.reconstruction_accuracy_2 = ComputeMetrics::getReconstructionAccuracy(box,gtMatVect_.at(box.index_in_vect), recMatVect_.at(box.index_in_vect),
                        cov_thresholds[i].reg_dist, cov_thresholds[i].likelihood, img_width_, img_height_, box_size_, img_res_);        
            //coverage 3
            i = 2;
            box.metrics.n_match_3 = ComputeMetrics::getBoxNMatch(box,gtMatVect_.at(box.index_in_vect), recMatVect_.at(box.index_in_vect),
                        cov_thresholds[i].reg_dist, cov_thresholds[i].likelihood, img_width_, img_height_, box_size_, img_res_);        
            box.metrics.surface_coverage_3 = ComputeMetrics::getSurfaceCoverage(box,gtMatVect_.at(box.index_in_vect), recMatVect_.at(box.index_in_vect),
                        cov_thresholds[i].reg_dist, cov_thresholds[i].likelihood, img_width_, img_height_, box_size_, img_res_);        
            box.metrics.reconstruction_accuracy_3 = ComputeMetrics::getReconstructionAccuracy(box,gtMatVect_.at(box.index_in_vect), recMatVect_.at(box.index_in_vect),
                        cov_thresholds[i].reg_dist, cov_thresholds[i].likelihood, img_width_, img_height_, box_size_, img_res_);        
            //Average Hausdorff Distance and Kappa:
            box.metrics.ahd_3 = ComputeMetrics::getAverageHausdorffDistance(box,gtMatVect_.at(box.index_in_vect), recMatVect_.at(box.index_in_vect),
                        cov_thresholds[i].likelihood, img_width_, img_height_, box_size_, img_res_);        
            box.metrics.kappa_3 = ComputeMetrics::getKappa(box,gtMatVect_.at(box.index_in_vect), recMatVect_.at(box.index_in_vect),
                        cov_thresholds[i].likelihood, img_width_, img_height_, box_size_, img_res_);        
            //coverage 4
            i = 3;
            box.metrics.n_match_4 = ComputeMetrics::getBoxNMatch(box,gtMatVect_.at(box.index_in_vect), recMatVect_.at(box.index_in_vect),
                        cov_thresholds[i].reg_dist, cov_thresholds[i].likelihood, img_width_, img_height_, box_size_, img_res_);        
            box.metrics.surface_coverage_4 = ComputeMetrics::getSurfaceCoverage(box,gtMatVect_.at(box.index_in_vect), recMatVect_.at(box.index_in_vect),
                        cov_thresholds[i].reg_dist, cov_thresholds[i].likelihood, img_width_, img_height_, box_size_, img_res_);        
            box.metrics.reconstruction_accuracy_4 = ComputeMetrics::getReconstructionAccuracy(box,gtMatVect_.at(box.index_in_vect), recMatVect_.at(box.index_in_vect),
                        cov_thresholds[i].reg_dist, cov_thresholds[i].likelihood, img_width_, img_height_, box_size_, img_res_);


        } else {
            box.metrics.n_match_1 = 0;        
            box.metrics.n_match_2 = 0;        
            box.metrics.n_match_3 = 0;        
            box.metrics.n_match_4 = 0;        
            box.metrics.volumetric_information = non_observed_vi_;        
            box.metrics.surface_coverage_1 = non_observed_cov_;
            box.metrics.surface_coverage_2 = non_observed_cov_;
            box.metrics.surface_coverage_3 = non_observed_cov_;
            box.metrics.surface_coverage_4 = non_observed_cov_;
            box.metrics.reconstruction_accuracy_1 = non_observed_acc_;
            box.metrics.reconstruction_accuracy_2 = non_observed_acc_;
            box.metrics.reconstruction_accuracy_3 = non_observed_acc_;
            box.metrics.reconstruction_accuracy_4 = non_observed_acc_;
            box.metrics.ahd_1 = box_size_*sqrt(3)*img_res_;
            box.metrics.ahd_3 = box_size_*sqrt(3)*img_res_;
            box.metrics.kappa_1 = -1;
            box.metrics.kappa_3 = -1;

        }
    }
    if (is_gt_empty) {
        if (is_observed) {
            box.metrics.emptyness_l1 = ComputeMetrics::getEmptynessL1(rec_vector);
        } else {
            box.metrics.emptyness_l1 = non_observed_l1_;
        }
    }
    if (is_dkl_metrics_) {
        if (is_observed) {
            box.metrics.dkl = ComputeMetrics::getKLDivergence(rec_vector, gt_vector);
        } else {
            box.metrics.dkl = non_observed_dkl_;
        }
    }
    if (is_ot_metrics_) {
        if (!is_gt_empty) {
            if (is_observed) {
                ComputeOTMetrics::OccFreeWasserstein wasserstein;
                wasserstein = ot_metrics_.normalizeAndGetSinkhornDistanceSigned(rec_vector, gt_vector);
                box.metrics.wasserstein_occ_remapped = wasserstein.occ;
                box.metrics.wasserstein_free_remapped = wasserstein.free;
            } else {
                box.metrics.wasserstein_occ_remapped = non_observed_wd_;
                box.metrics.wasserstein_free_remapped = non_observed_wd_;

            }
        }
    }
}

void ComparatorRecGt::processSliceInListOfBoxesSample(int start, int end) {
    assert(sampleOfBoxes_.begin() + end <=sampleOfBoxes_.end());
    unsigned int counter{0};
    for (auto it = sampleOfBoxes_.begin()+start; it!=sampleOfBoxes_.begin()+end; ++it){
        //we want to calculate the limit only
        //and we save to disk only at the end
        calcLimitMetricFromBox(*it);
    }
}

void ComparatorRecGt::processSliceInListOfBoxes(int start, int end) {
    assert(listOfBoxes_.begin() + end <= listOfBoxes_.end());
    if (unit_test_) {
        cout << "Flag unit_test_ set to: " << boolalpha << unit_test_ << endl;
    }
    for (auto it = listOfBoxes_.begin()+start; it!=listOfBoxes_.begin()+end; ++it){
        if (!unit_test_) {  
            calcMetricFromBox(*it);
        } else {
            unitTestBox(*it);
        }
    }
}
//We loop on GT dataset and get N samples of a mimimum level occupancy <level>
void ComparatorRecGt::addNSamplesWithAMinLevelOfOcc(std::vector<ComparatorDatatypes::BoxToProcessLimitOnly>& samples, size_t n_target, double level) const {
    size_t init_n = samples.size();
    for (auto box : listOfBoxesLimitOnly_) {
        if(isBoxNotEmptyWithLevel(box, gtMatVect_.at(box.index_in_vect), level)) {
            samples.push_back(box);
            if ((samples.size() - init_n) >= n_target) {
                return;
            }
        }
    }
    cout << "WARNING, not enough boxes found with the desired level of occupancy" << endl;
    assert(samples.size()>0);//if it is zeros, we just stop
    return;
}
void ComparatorRecGt::addNSamplesEmpty(std::vector<ComparatorDatatypes::BoxToProcessLimitOnly>& samples, size_t n_target) const {
    size_t init_n = samples.size();
    for (auto box : listOfBoxesLimitOnly_) {
        if(isBoxOnlyZeros(box, gtMatVect_.at(box.index_in_vect))) {
            samples.push_back(box);
            if ((samples.size() - init_n) >= n_target) {
                return;
            }
        }
    }
    cout << "WARNING, not enough empty boxes found" << endl;
    assert(samples.size()>0);//if it is zeros, we just stop
    return;
}
void ComparatorRecGt::saveCubeToImg(const ComparatorDatatypes::BoxToProcessLimitOnly& box, const cv::Mat_<uint8_t>& refMat, const std::string& fname, const std::string& dir) const {
    BoxTools::saveSampleCubeToImg(box.ranges, refMat, fname, dir);
}
void ComparatorRecGt::saveCubeToImg(const ComparatorDatatypes::BoxToProcess& box, const cv::Mat_<uint8_t>& refMat, const std::string& fname, const std::string& dir) const {
    BoxTools::saveSampleCubeToImg(box.ranges, refMat, fname, dir);
}

void ComparatorRecGt::calcLimitWDDataset(int xp_number, bool sample_with_occ, bool sample_half_half, bool sample_empty_only, bool sample_with_ratio, size_t n_samples, double occ_level, double dataset_ratio){
    //check that we can write results to disk at the end:
    std::string dir = output_dir_ + "/"; 
    char suffix[40];
    
    //load gt dataset
    getListOfBoxesLimitOnly();
    loadGTOnlyDataset();

    std::random_device rd;
    static std::default_random_engine engine(rd());
    if (sample_with_occ) {
        
        sprintf(suffix, "_sample_with_occ_%.02f", occ_level);
        res_base_filename_ = dir + "limit_dataset_" + to_string(xp_number) + suffix;
        resf_short_ = ofstream{res_base_filename_ +".csv",ios::out};
        checkFile(resf_short_);
        resf_short_.close();
        
        cout << "will sample at max " << n_samples << " with an occupancy ratio of " << occ_level << endl;
        std::shuffle(listOfBoxesLimitOnly_.begin(), listOfBoxesLimitOnly_.end(), engine); 
        sampleOfBoxes_.clear();
        addNSamplesWithAMinLevelOfOcc(sampleOfBoxes_, n_samples, occ_level); 
    } else if (sample_half_half) {
        
        sprintf(suffix, "_sample_half_half_occ_%.02f", occ_level);
        res_base_filename_ = dir + "limit_dataset_" + to_string(xp_number) + suffix;
        resf_short_ = ofstream{res_base_filename_ + ".csv",ios::out};
        checkFile(resf_short_);
        resf_short_.close();
        
        cout << "will sample at max " << n_samples/2 << " with an occupancy ratio of " << occ_level << "and as many empty." << endl;
        std::shuffle(listOfBoxesLimitOnly_.begin(), listOfBoxesLimitOnly_.end(), engine); 
        sampleOfBoxes_.clear();
        addNSamplesWithAMinLevelOfOcc(sampleOfBoxes_, n_samples/2, occ_level); 
        addNSamplesEmpty(sampleOfBoxes_, sampleOfBoxes_.size()); 
    } else if (sample_empty_only) {
        sprintf(suffix, "_sample_empty_only");
        res_base_filename_ = dir + "limit_dataset_" + to_string(xp_number) + suffix;
        resf_short_ = ofstream{res_base_filename_ + ".csv",ios::out};
        checkFile(resf_short_);
        resf_short_.close();
        cout << "will sample at max " << n_samples << " empty samples." << endl;
        std::shuffle(listOfBoxesLimitOnly_.begin(), listOfBoxesLimitOnly_.end(), engine); 
        //and then we sample
        sampleOfBoxes_.clear();
        addNSamplesEmpty(sampleOfBoxes_, n_samples); 
    } else if (sample_with_ratio) {
        cout << "will sample a ratio of the gt dataset " << dataset_ratio << endl;
    //in the following case, we sample a ratio of the dataset, without looking if the boxes are occupied or free:
    //select a random sample on which to calculate the limit (10%)
        sprintf(suffix, "_sample_with_ratio_%.02f", dataset_ratio);
        res_base_filename_ = dir + "limit_dataset_" + to_string(xp_number) + suffix;
        resf_short_ = ofstream{res_base_filename_ + ".csv",ios::out};
        checkFile(resf_short_);
        resf_short_.close();
        
        assert(dataset_ratio > 0);
        assert(dataset_ratio < 1);
        size_t n_samples = std::floor(listOfBoxesLimitOnly_.size() * dataset_ratio);
        assert(n_samples > 0);
        
        sampleOfBoxes_.resize(n_samples);
        std::sample(listOfBoxesLimitOnly_.begin(), listOfBoxesLimitOnly_.end(), sampleOfBoxes_.begin(), n_samples, engine);

    } else {
        cout << "not doing anything" << endl;
    }
    
    int tot = sampleOfBoxes_.size();
    int k = tot / n_threads_;
    cout << "sample size is: " << tot << endl;
    
    std::vector<std::shared_ptr<std::thread> > threads;
    threads.resize(n_threads_);
    for (int i{0}; i < n_threads_; ++i) {
        if (i == n_threads_ - 1) {
            threads[i].reset(new std::thread(&ComparatorRecGt::processSliceInListOfBoxesSample,this,i*k, tot));

        } else {
            threads[i].reset(new std::thread(&ComparatorRecGt::processSliceInListOfBoxesSample,this,i*k, (i+1)*k));
        }
    }
    for (int i{0}; i < n_threads_; ++i) {
        threads[i]->join();
    }
    saveResultsLimitToDisk(tot, sampleOfBoxes_); 
    cout << "csv results are stored in: " << dir << endl;
}

void ComparatorRecGt::compare(int xp_number){
    //check that we can write results to disk at the end:
    std::string dir = output_dir_ + "/"; 
    res_base_filename_ = dir + "comparaison_new_" + to_string(xp_number);
    resf_short_ = ofstream{res_base_filename_ + ".csv",ios::out};
    std::cout << "Test writing on: " << res_base_filename_ << std::endl;
    checkFile(resf_short_);
    resf_short_.close();

    getListOfBoxes();
    loadCompleteDataset();
    
    //and now we multithread:
    int tot = listOfBoxes_.size();
    int k = tot/n_threads_;
    std::cout << " number of boxes to process: " << tot << std::endl;
    std::vector<std::shared_ptr<std::thread> > threads;
    threads.resize(n_threads_);
    for (int i{0}; i < n_threads_; ++i) {
        if (i == n_threads_ - 1) {
            threads[i].reset(new std::thread(&ComparatorRecGt::processSliceInListOfBoxes,this,i*k, tot));

        } else {
            threads[i].reset(new std::thread(&ComparatorRecGt::processSliceInListOfBoxes,this,i*k, (i+1)*k));
        }
    }
    for (int i{0}; i < n_threads_; ++i) {
        threads[i]->join();
    }

    saveResultsToDisk(tot); 
}

void ComparatorRecGt::saveResultsLimitToDisk(int count, const std::vector<ComparatorDatatypes::BoxToProcessLimitOnly>& samples) {
    //write the results to disk
    char end_name[20];
    sprintf(end_name, "_%09d.csv", count);
    resf_short_ = ofstream{res_base_filename_ + end_name ,ios::out};
    resf_short_ << "x y z"
                << " random_emptyness_l1"
                << " random_wd_occ random_wd_free" 
                << " ideal_conv1_wd_occ ideal_conv1_wd_free" 
                << " ideal_conv2_wd_occ ideal_conv2_wd_free" 
                << " ideal_conv3_wd_occ ideal_conv3_wd_free" 
                << " n_gt_points"
                << endl;

    for (auto box : samples) {
        resf_short_ << fixed << setprecision(2) 
            << box.mbbox.xmin << " "
            << box.mbbox.ymin << " "
            << box.mbbox.zmin
            << setprecision(8)
            << " " << box.random.emptyness_l1
            << " " << box.random.wasserstein_occ_remapped << " " << box.random.wasserstein_free_remapped
            << " " << box.ideal_conv1.wasserstein_occ_remapped << " " << box.ideal_conv1.wasserstein_free_remapped
            << " " << box.ideal_conv2.wasserstein_occ_remapped << " " << box.ideal_conv2.wasserstein_free_remapped
            << " " << box.ideal_conv3.wasserstein_occ_remapped << " " << box.ideal_conv3.wasserstein_free_remapped
            << " " << box.random.n_gt_points
            << endl;
    }
    resf_short_.close();
}

void ComparatorRecGt::saveResultsToDisk() {
    saveResultsToDisk(0);
}
void ComparatorRecGt::saveResultsToDisk(int count) {
    char end_name[20];
    sprintf(end_name, "%09d.csv", count);
    resf_short_ = ofstream{res_base_filename_ + end_name ,ios::out};
    std::cout <<"Writing results to: " << res_base_filename_ + end_name << std::endl;
    resf_short_ << "x y z";
    if (is_dkl_metrics_) {
        resf_short_ << " dkl";
    }
    if (is_ot_metrics_) {
        resf_short_ << " wasserstein_occ_remapped wasserstein_free_remapped wasserstein_direct";
    }
    if (is_dist_metrics_) {
        resf_short_ << " surface_coverage_1 surface_coverage_2 surface_coverage_3 surface_coverage_4"
                    << " reconstruction_accuracy_1 reconstruction_accuracy_2 reconstruction_accuracy_3 reconstruction_accuracy_4" 
                    << " volumetric_information"
                    << " n_match_1 n_match_2 n_match_3 n_match_4"
                    << " ahd_1 ahd_3 kappa_1 kappa_3";
    }
    resf_short_ << " emptyness_l2 emptyness_l1 n_rec_points n_gt_points";
    resf_short_ << endl;

    for (auto box : listOfBoxes_) {
        resf_short_ << fixed << setprecision(2) 
            << box.mbbox.xmin << " "
            << box.mbbox.ymin << " "
            << box.mbbox.zmin
            << setprecision(8);
        if (is_dkl_metrics_) {
            resf_short_ << " " << box.metrics.dkl;
        }
        if (is_ot_metrics_) {
            resf_short_ << " " << box.metrics.wasserstein_occ_remapped << " " << box.metrics.wasserstein_free_remapped << " " << box.metrics.wasserstein_direct;
        }
        if (is_dist_metrics_) {
            resf_short_ << " " << box.metrics.surface_coverage_1 << " " << box.metrics.surface_coverage_2 << " " << box.metrics.surface_coverage_3 << " " << box.metrics.surface_coverage_4
                        << " " << box.metrics.reconstruction_accuracy_1  << " " << box.metrics.reconstruction_accuracy_2 << " " << box.metrics.reconstruction_accuracy_3 << " " << box.metrics.reconstruction_accuracy_4 
                        << " " << box.metrics.volumetric_information;
            resf_short_ << " " << box.metrics.n_match_1 << " " << box.metrics.n_match_2 << " " << box.metrics.n_match_3 << " " << box.metrics.n_match_4 ;
            resf_short_ << " " << box.metrics.ahd_1 << " " << box.metrics.ahd_3 << " " << box.metrics.kappa_1 << " " << box.metrics.kappa_3;
        }
        resf_short_ << " " << box.metrics.emptyness_l2 << " " << box.metrics.emptyness_l1; 
        resf_short_ << " " << box.metrics.n_rec_points << " " << box.metrics.n_gt_points; 
        resf_short_ << endl;
    }
    resf_short_.close();
    cout << "we have processed at least: " << counter_target_ << " voxels." << endl;
    string dir = output_dir_+ "/";
    string out_info = dir + "compare_ot_gt_info.txt";
    std::cout <<"Writing info to: " << out_info << std::endl;
    ofstream out_file_info{out_info,ios::out};
    checkFile(out_file_info);
    out_file_info << "xmin " << constrained_bbox_.xmin << endl
                  << "xmax " << constrained_bbox_.xmax << endl
                  << "ymin " << constrained_bbox_.ymin << endl
                  << "ymax " << constrained_bbox_.ymax << endl
                  << "zmin " << constrained_bbox_.zmin << endl
                  << "zmax " << constrained_bbox_.zmax << endl
                  << "res " << target_res_ << endl;
    out_file_info.close();
     
    cout << "csv results are stored in: " << dir << endl;

}
void ComparatorRecGt::save_sample_debug(const ComparatorDatatypes::BoxToProcess& box, int debugcount) const {
    string dir = base_dir_ + "res/obs/" + xp_name_ + "/";
    ofstream debugf = ofstream(dir + "sample_debug.csv", ios::app);
    checkFile(debugf);
    debugf << debugcount << " "
           << fixed << setprecision(2)
           << box.mbbox.xmin << " "
           << box.mbbox.ymin << " "
           << box.mbbox.zmin << endl;
    
    char fname1[40];
    sprintf(fname1, "debug_gt_%02d_", debugcount);
    saveSampleCubeToImg(box.ranges, gtMatVect_.at(box.index_in_vect), fname1, dir);
    char fname2[40];
    sprintf(fname2, "debug_ot_%02d_", debugcount);
    saveSampleCubeToImg(box.ranges, recMatVect_.at(box.index_in_vect), fname2, dir);
}
