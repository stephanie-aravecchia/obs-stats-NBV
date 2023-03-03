#ifndef COMPARE_REC_GT_H
#define COMPARE_REC_GT_H
#include "comparator/datatypes.h"
#include "comparator/compute_ot_metrics.h"
#include "comparator/box_tools.h"


//This class compares 2 datasets ground-truth and reconstruction
//The datasets are composed of :images, image list, and the bounding-box and resolution of the images
//This comparator loads the same sample of ground-truth and reconstruction, compares it, and save the metrics of the comparaison to disk 
class ComparatorRecGt : protected BoxTools{

    protected:
    //public:
        std::string base_dir_;
        std::string xp_name_;
        std::string output_dir_;
        std::vector<std::string> gt_imglist_;
        std::vector<std::string> rec_imglist_;
        double target_res_;
        double img_res_;
        size_t tot_box_;
        size_t counter_{0};
        size_t counter_target_;
        int n_slices_;
        int istart_;
        int jstart_;
        int kstart_;
        int nx_;
        int ny_;
        int nz_;
        std::vector<cv::Mat_<uint8_t> > gtMatVect_;
        std::vector<cv::Mat_<uint8_t> > recMatVect_;
        std::vector<ComparatorDatatypes::BoxToProcess> listOfBoxes_;
        std::vector<ComparatorDatatypes::BoxToProcessLimitOnly> listOfBoxesLimitOnly_;
        std::vector<ComparatorDatatypes::BoxToProcessLimitOnly> sampleOfBoxes_;
        std::ofstream resf_;  
        std::ofstream resf_short_;  
        std::ofstream resf_vect_short_;  
        std::string res_base_filename_;
        bool keep_complete_file_{false}; 
        ComparatorDatatypes::PixBBox space_pix_bbox_; 
        ComparatorDatatypes::BBox space_bbox_; 
        ComparatorDatatypes::BBox constrained_bbox_; 
        ComparatorDatatypes::Offsets rec_offset_;
        ComparatorDatatypes::Offsets gt_offset_;
        double ground_thres_; //height starting from the bottom we don't use
        ComputeOTMetrics ot_metrics_;
        bool is_dist_metrics_; 
        bool is_dkl_metrics_; 
        bool is_ot_metrics_; 
        int n_threads_;
        double ot_reg_;
        int ot_maxiter_;
        double occupancy_thres_; //to consider a voxel occupied only to count n_gt points and n_rect points, the metrics use cov_thresholds 
        double reg_distance_; //to consider a match between rec and gt voxels
        double ot_stopthres_;

        //constants to use when a box is not observed and threshold to define it is or not observed
        double divergence_to_unknown_thres_;
        double non_observed_vi_;
        double non_observed_cov_;
        double non_observed_acc_;
        double non_observed_l1_;
        double non_observed_dkl_;
        double non_observed_wd_;
        
        //to calculate the limit of metrics on the dataset:
        //the values are (sigma, ksize, additional noise)
        std::array<ComparatorDatatypes::NoiseOnGtParams, 3> noise_on_gt_ = std::array{
            ComparatorDatatypes::NoiseOnGtParams(0.05,7,0),
            ComparatorDatatypes::NoiseOnGtParams(0.08,7,0.05),
            ComparatorDatatypes::NoiseOnGtParams(0.2,11,0.1)
        };
        //we define the couple of thresholds used for distance metric such as coverage.
        //the values are (occupancy likelihood, registration distance)
        std::array<ComparatorDatatypes::CovThresholds,4> cov_thresholds = std::array{
            ComparatorDatatypes::CovThresholds(0.8, 0.05),
            ComparatorDatatypes::CovThresholds(0.8, 0.1),
            ComparatorDatatypes::CovThresholds(0.7, 0.1),
            ComparatorDatatypes::CovThresholds(0.7, 0.15),
        };
        //to test a selection of boxes only
        bool unit_test_;

        //get the list of img path we will load from the img list on disk and the directory
        //we start loading with an offset
        void getImgList(std::vector<std::string >&, const std::string&, const std::string&, int) const;
        
        //get starting index to start the comparaison at a coordinate which is a multiple of res
        //i.e if xmin is -5.28 and res .1, we start the comparaison a -5.20
        //this is mandatory for the comparaison to the observation in next stages
        int getStartingIndex(double) const;
        //get starting distance after the initial bounding box to match the starting index
        double getStartingDistance(double) const;

        //load a succession of box_size images after start, into a 3D Mat, 
        //the index of the pixels from the image that will be loaded in the 3D mat depend on the offsets and the width and height of the 3dmat
        void load3DMat(cv::Mat_<uint8_t>& mat, const std::vector<std::string >& imglist, const ComparatorDatatypes::Offsets& offset, int start) const;

        void loadCompleteDataset();
        void loadGTOnlyDataset();
        void reshapeMatVect();
        void getListOfBoxes();
        void getListOfBoxesLimitOnly();
        void saveResultsToDisk(int);
        void saveResultsLimitToDisk(int, const std::vector<ComparatorDatatypes::BoxToProcessLimitOnly>&);

        bool isQueryBoxIncludedInRef(const ComparatorDatatypes::BBox&, const ComparatorDatatypes::BBox&) const; 
        bool isQueryBoxIncludedInRef(const ComparatorDatatypes::PixBBox&, const ComparatorDatatypes::PixBBox&) const; 
        
        //Determine which metric will be used, call it appropriately, and store the results
        void calcMetricFromBox(ComparatorDatatypes::BoxToProcess& box);
        //calculate the limit of the Wasserstein distance, with a reconstruction drawned randomly from gt
        void calcLimitMetricFromBox(ComparatorDatatypes::BoxToProcessLimitOnly& box);
        void calcLimitMetricFromBoxTestOnly(ComparatorDatatypes::BoxToProcessLimitOnly& box, std::ofstream&);
        void computeLimitMetricsOnVect(ComparatorDatatypes::Metrics&, std::vector<double>&, std::vector<double>&);

        void unitTestBox(ComparatorDatatypes::BoxToProcess&);
        void unitTestBox(ComparatorDatatypes::BoxToProcessLimitOnly&);
        
        void getVectorFromBox(const ComparatorDatatypes::BoxToProcess&, const cv::Mat_<uint8_t>&, std::vector<double>&) const;
        void getVectorFromBox(const ComparatorDatatypes::BoxToProcessLimitOnly&, const cv::Mat_<uint8_t>&, std::vector<double>&) const;
        void getNoisyGtVector(const ComparatorDatatypes::BoxToProcessLimitOnly&, const cv::Mat_<uint8_t>&, std::vector<double>&, double uniform_noise_level=.1, bool use_fixed_sigma=false, double sigma=0.) const;
        
        bool isBoxOnlyZeros(const ComparatorDatatypes::BoxToProcessLimitOnly&, const cv::Mat_<uint8_t>&) const;        
        bool isBoxOnlyZeros(const ComparatorDatatypes::BoxToProcess&, const cv::Mat_<uint8_t>&) const;        
        double getBoxOccupancyRate(const ComparatorDatatypes::BoxToProcessLimitOnly&, const cv::Mat_<uint8_t>&) const;        
        double getBoxOccupancyRate(const ComparatorDatatypes::BoxToProcess&, const cv::Mat_<uint8_t>&) const;        
        //if the box contains at least <level> of not empty points, returns true
        bool isBoxNotEmptyWithLevel(const ComparatorDatatypes::BoxToProcessLimitOnly&, const cv::Mat_<uint8_t>&,double) const;        
        void addNSamplesWithAMinLevelOfOcc(std::vector<ComparatorDatatypes::BoxToProcessLimitOnly>&, size_t, double) const; 
        void addNSamplesEmpty(std::vector<ComparatorDatatypes::BoxToProcessLimitOnly>&, size_t) const; 


        void processSliceInListOfBoxes(int, int); 
        void processSliceInListOfBoxesSample(int, int); 
        
        //utility functions to convert coordinates of the BBOX from pixels to meters
        void bboxPixToMeters(const ComparatorDatatypes::PixBBox&, ComparatorDatatypes::BBox&) const;
        
        
        //debug
        void saveCubeToImg(const ComparatorDatatypes::BoxToProcess&, const cv::Mat_<uint8_t>&, const std::string&, const std::string&) const;
        void saveCubeToImg(const ComparatorDatatypes::BoxToProcessLimitOnly&, const cv::Mat_<uint8_t>&, const std::string&, const std::string&) const;
        bool isSampleDebug(const ComparatorDatatypes::BBox& mbox, double x, double y, double z) const;
        void save_sample_debug(const ComparatorDatatypes::BoxToProcess& box, int i) const;

    public:
        ComparatorRecGt();
        ComparatorRecGt(const std::string& base_dir, const std::string& xp, const std::string& output_dir, double target_res, double img_res, 
                        ComparatorDatatypes::PixBBox pixbox, ComparatorDatatypes::BBox metricsbox,
                        ComparatorDatatypes::Offsets gt_offset, ComparatorDatatypes::Offsets rec_offset, double ground_thres,
                        bool dist_metrics, bool dkl_metrics, bool ot_metrics, int nthreads, 
                        double ot_reg, int ot_maxiter, double occ_thres, double ot_stop_thres=1e-9,
                        double divergence_to_unknown=0.1, double non_observed_vi=0, double non_observed_cov=0, double non_observed_acc=0,
                        double non_observed_l1=500, double non_observed_dkl=6600, double non_observed_wd=100, bool unit_test=false, const std::string& test_filename="");
        //determine on which images we can compare, and compare them
        void compare(int xp_number);
        //calculate the limit of the wasserstein distance based on a sample of GT
        void calcLimitWDDataset(int xp_number, bool sample_with_occ, bool sample_half_half, bool sample_empty_only, bool sample_with_ratio, size_t n_samples, double occ_level, double dataset_ratio);
        void saveResultsToDisk();
};

#endif
