#ifndef COMPUTE_METRICS_H
#define COMPUTE_METRICS_H
#include <vector>
#include <Eigen/Eigen>
#include "comparator/datatypes.h"


//Implementation of the metrics computation
namespace ComputeMetrics {
    
        //Occupied space metrics

        //Get the sum of the KLD between every element of rec and gt
        double getKLDivergence(const std::vector<double>&, const std::vector<double>&); 
        
        //Get the sum of the entropy of the elements 
        double getVolumetricInformation(const std::vector<double>&); 
        
        //Metrics below are based on Euclidean distance
        //In surface coverage and in accuracy, we compare 3D grids GT and REC.
        //If the distance between a point in GT and its closest point in REC is less than a registration distance,
        //the GT point is considered as observed
        //The Surface Coverage is then the percentage of GT points considered observed on to the total number of GT points
        //The reconstruction Accuracy is the symetric, where we inverse GT and REC.
        double getSurfaceCoverage(const ComparatorDatatypes::BoxToProcess& box, const cv::Mat_<uint8_t>& gtMat, const cv::Mat_<uint8_t>& recMat, 
                        double reg_dist, double occ_thres, int img_width, int img_height, int box_size, double img_res);
        double getReconstructionAccuracy(const ComparatorDatatypes::BoxToProcess& box, const cv::Mat_<uint8_t>& gtMat, const cv::Mat_<uint8_t>& recMat, 
                        double reg_dist, double occ_thres, int img_width, int img_height, int box_size, double img_res);
        double getBoxDistMetrics(const ComparatorDatatypes::BoxToProcess& box, const cv::Mat_<uint8_t>& refMat, const cv::Mat_<uint8_t>& queryMat, 
                        double reg_dist, double occ_thres, int img_width, int img_height, int box_size, double img_res);
        int getBoxNMatch(const ComparatorDatatypes::BoxToProcess& box, const cv::Mat_<uint8_t>& refMat, const cv::Mat_<uint8_t>& queryMat, 
                        double reg_dist, double occ_thres, int img_width, int img_height, int box_size, double img_res);
        
        //The Regular Hausdorff distance
        double getHausdorffDistance(const ComparatorDatatypes::BoxToProcess& box, const cv::Mat_<uint8_t>& gtMat, const cv::Mat_<uint8_t>& recMat, 
                        double occ_thres, int img_width, int img_height, int box_size, double img_res);
        //The Average Hausdorff Distance: symmetric HD between query and rec, and rec and query
        double getAverageHausdorffDistance(const ComparatorDatatypes::BoxToProcess& box, const cv::Mat_<uint8_t>& gtMat, const cv::Mat_<uint8_t>& recMat, 
                        double occ_thres, int img_width, int img_height, int box_size, double img_res);

        struct Confusion {
            int tp = 0;
            int fp = 0;
            int tn = 0;
            int fn = 0;
            Confusion() {}
            Confusion(int tp, int fp, int tn, int fn):
                tp(tp), fp(fp), tn(tn), fn(fn) {}
        };
        //Computes the confusion matrix of the reconstruction vs ground_truth
        Confusion getConfusionMatrix(const ComparatorDatatypes::BoxToProcess& box, const cv::Mat_<uint8_t>& gtMat, const cv::Mat_<uint8_t>& recMat, 
                        double occ_thres, int img_width, int img_height, int box_size, double img_res);

        //Compute the Cohen's Kappa Metric
        //ref_mat is the GT
        double getKappa(const ComparatorDatatypes::BoxToProcess& box, const cv::Mat_<uint8_t>& gtMat, const cv::Mat_<uint8_t>& recMat, 
                        double occ_thres, int img_width, int img_height, int box_size, double img_res);
        
        
        double getSurfaceCoverage(const ComparatorDatatypes::BoxToProcess& box, const cv::Mat_<uint8_t>& gtMat, const cv::Mat_<uint8_t>& recMat, 
                        double reg_dist, double occ_thres, int img_width, int img_height, int box_size, double img_res);
        //L_1 and L_2 norm, metrics for empty space
        double getEmptynessL1(const std::vector<double>&);
        double getEmptynessL2(const std::vector<double>&);
        
        //Several utility functions
        //returns the number of points in the box above the occupancy threshold
        int getNpoints(const ComparatorDatatypes::BoxToProcess& box, const cv::Mat_<uint8_t>& refMat, double occ_thres, 
                        int img_width, int img_height, int box_size);
        int getNpoints(const std::vector<double>& vect, double occ_thres);
        
        void normalize(const std::vector<double>&, std::vector<double>&);
        
        //utility function to normalize double to color and the inverse
        uint8_t getNormalizedPixValue(double v);
        double pix2proba(uint8_t);
        
        //debug functions to make sure GT is composed only of 0 and 1
        bool isGTVectorBinary(const std::vector<double>&);
        bool isGTImageBinary(const cv::Mat_<uint8_t>&);
        bool isValue0or1(double);
        bool isValue0or1(uint8_t);
};

#endif