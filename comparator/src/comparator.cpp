#include "comparator/compare_rec_gt.h"
#include "comparator/dataset_info.h"
#include <ros/ros.h>
#include <assert.h>
#include <thread>
#include <algorithm>
#include <signal.h>
#include <unistd.h>

using namespace std;

static ComparatorRecGt* global_it;

//Save when ctrl-C
static void sig_handler(int signo) {
    if (signo == SIGINT) { 
        global_it->saveResultsToDisk(); 
    }
    exit(0);
};

//This code is a ROS wrapper on ComparatorRecGT, mostly to ease the use with launchfiles
class Comparator {
    protected:
    public : 
    Comparator() {
        
        ros::NodeHandle nh_("~");
        std::string base_dir;
        std::string xp_name;
        std::string output_dir;
        std::string test_file_name;
        double ot_reg;
        double ot_stop_thres;
        int ot_max_iter;
        int nthreads;
        int xpnum;
        double output_res;
        double ground_thres;
        bool limit_only;
        bool sample_with_occ;
        bool sample_half_half;
        bool sample_with_ratio;
        bool sample_empty_only;
        size_t nsamples;
        double occ_level_selection;
        double occ_thres; 
        double dataset_ratio; 

        bool do_ot_metrics;
        bool do_hausdorff_metrics;
        bool do_dkl_metrics;
        double divergence_to_unknown_thres;
        double non_observed_vi;
        double non_observed_cov;
        double non_observed_acc;
        double non_observed_l1;
        double non_observed_dkl;
        double non_observed_wd;
        bool unit_test;

        nh_.param<std::string>("base_dir",base_dir,"");
        nh_.param<std::string>("xp_name",xp_name,"");
        nh_.param<std::string>("output_dir",output_dir,"");
        nh_.param<bool>("limit_only",limit_only,false);
        nh_.param<bool>("do_ot_metrics",do_ot_metrics,true);
        nh_.param<bool>("do_hausdorff_metrics",do_hausdorff_metrics,true);
        nh_.param<bool>("do_dkl_metrics",do_dkl_metrics,false);
        nh_.param<bool>("limit_only",limit_only,false);
        nh_.param<bool>("sample_with_occ",sample_with_occ,false);
        nh_.param<bool>("sample_half_half",sample_half_half,false);
        nh_.param<bool>("sample_with_ratio",sample_with_ratio,false);
        nh_.param<bool>("sample_empty_only",sample_empty_only,false);
        nh_.param<double>("ot_reg",ot_reg,1.);
        nh_.param<int>("ot_maxiter",ot_max_iter,1000);
        nh_.param<double>("ot_stop_thres",ot_stop_thres,1e-9);
        nh_.param<double>("output_res",output_res,0.2);
        nh_.param<double>("ground_thres",ground_thres,0);
        nh_.param<int>("n_threads",nthreads,1);
        nh_.param<int>("xp_num",xpnum, 0);
        int n;
        nh_.param<int>("n_samples",n,5000);
        nsamples = static_cast<size_t>(n);
        nh_.param<double>("occ_level_selection",occ_level_selection,0.001);//ratio of occupied vox in cube
        nh_.param<double>("occ_thres_hausdorff",occ_thres,0.55);
        nh_.param<double>("dataset_ratio",dataset_ratio,0.1);
        nh_.param<double>("divergence_to_unknown_thres",divergence_to_unknown_thres,0.1);
        nh_.param<double>("non_observed_vi",non_observed_vi,0.0);
        nh_.param<double>("non_observed_cov",non_observed_cov,0.0);
        nh_.param<double>("non_observed_acc",non_observed_acc,0.0);
        nh_.param<double>("non_observed_l1",non_observed_l1,500.0);
        nh_.param<double>("non_observed_dkl",non_observed_dkl,6600.0);
        nh_.param<double>("non_observed_wd",non_observed_wd,100.0);
        nh_.param<bool>("unit_test",unit_test,false);
        nh_.param<std::string>("unit_test_file",test_file_name,"");

        DatasetInfo dataset_info(base_dir, xp_name);
        ComparatorDatatypes::PixBBox pixbox = dataset_info.getDatasetPixBbox();
        ComparatorDatatypes::BBox metricsbox = dataset_info.getDatasetMetricsBbox();
        ComparatorDatatypes::Offsets gt_offsets = dataset_info.getDatasetGTOffsets();
        ComparatorDatatypes::Offsets ot_offsets = dataset_info.getDatasetOTOffsets();
        double img_res = dataset_info.getDatasetResolution();
        cout << "dataset information loaded : "
            << "pixbox is : "<< pixbox << endl 
            << "gt offsets: " << gt_offsets << endl
            << "ot offsets: " << ot_offsets << endl
            << "img resolution: " << img_res << endl 
            << "Now starting the comparison." << endl; 
        
        //Regular version of the code, to compare the reconstruction to the GT
        ComparatorRecGt comparator(base_dir, xp_name, output_dir, output_res, img_res, pixbox, metricsbox, 
                gt_offsets, ot_offsets, ground_thres, do_hausdorff_metrics, do_dkl_metrics, do_ot_metrics, nthreads, 
                ot_reg, ot_max_iter, occ_thres, ot_stop_thres,
                divergence_to_unknown_thres, non_observed_vi, non_observed_cov, non_observed_acc, non_observed_l1,
                non_observed_dkl, non_observed_wd, unit_test, test_file_name);
        cout << "Comparator instance created" << endl;
        global_it = &comparator;
        signal(SIGINT, sig_handler);

        //If we want to compute the limit only of the dataset
        if (limit_only) {
            cout << "We will calculate the limit of the Wasserstein Distance on the dataset only" << endl;
            if (sample_with_occ) {
                cout << "Starting sample_with_occ" << endl;
                comparator.calcLimitWDDataset(xpnum, sample_with_occ, false, false, false, nsamples, occ_level_selection, dataset_ratio);
            }
            if (sample_half_half) {
                cout << "Starting sample_half_half" << endl;
                comparator.calcLimitWDDataset(xpnum, false, sample_half_half, false, false, nsamples, occ_level_selection, dataset_ratio);
            } 
            if (sample_empty_only) {
                cout << "Starting sample_empty_only" << endl;
                comparator.calcLimitWDDataset(xpnum, false, false, sample_empty_only, false, nsamples, occ_level_selection, dataset_ratio);
            }
            if (sample_with_ratio) {
                cout << "Starting sample_with_ratio" << endl;
                comparator.calcLimitWDDataset(xpnum, false, false, false, sample_with_ratio, nsamples, occ_level_selection, dataset_ratio);
            }
        } else {
            cout << "We will proceed to a complete comparaison" << endl;
            comparator.compare(xpnum);

        }
	cout << "Comparaison complete." << endl;
	exit(0);
    }
};

int main(int argc, char** argv) {
    ros::init(argc, argv, "comparator");
    Comparator comparator;
    return 0;
};
