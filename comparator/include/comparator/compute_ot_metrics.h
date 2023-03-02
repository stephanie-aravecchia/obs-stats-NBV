#ifndef COMPUTE_OT_METRICS_H
#define COMPUTE_OT_METRICS_H
#include <vector>
#include <array>
#include <random>
#include <Eigen/Eigen>
#include "functions/Function.h"

//This class is implementing some optimal transport algorithms, 
//It follows the python implementation in Python Optimal Transport
//https://pythonot.github.io/

//The core of the algorithms follows:
//Cuturi, Marco. "Sinkhorn distances: Lightspeed computation of optimal transport." Advances in neural information processing systems 26 (2013).

class ComputeOTMetrics {
    
    protected:
        unsigned int cube_nrows_; 
        unsigned int cube_ncols_; 
        unsigned int cube_nslices_; 
        Eigen::MatrixXd M_;
        Eigen::MatrixXd K_;
        Eigen::VectorXd uniform_;
        double reg_;
        int numItermax_;
        double stopThres_;
    
    public:
        ComputeOTMetrics() {};
        //1d constructor cols , reg, numItermax, stopThres
        ComputeOTMetrics(unsigned int, double, int, double);
        //2d constructor rows, cols, reg, numItermax, stopThres
        ComputeOTMetrics(unsigned int, unsigned int, double, int, double);
        //3d constructor slices, rows, cols, reg, numItermax, stopThres
        ComputeOTMetrics(unsigned int, unsigned int, unsigned int, double, int, double);

        struct OccFreeWasserstein {
            double occ;
            double free;
            OccFreeWasserstein() {}
            OccFreeWasserstein(double occ, double free):
                occ(occ), free(free) {}
        };

        //takes std::vectors as input, transform to Eigen::Vector, remap to [-1,1], split to positive and negative
        //distribution, normalize (add noise and sum to one) each distribution
        //and then perform the optimization on positive (occ) and negative (free), returns the wasserstein distance
        //on occupied spacea and on free space
        OccFreeWasserstein normalizeAndGetSinkhornDistanceSigned(std::vector<double>& a, std::vector<double>& b)const;
        //same, but normalize the wasserstein distance
        OccFreeWasserstein normalizeAndGetSinkhornDistanceSignedNormalized(std::vector<double>& a, std::vector<double>& b) const;
        
        //takes std::vectors as input, transform to Eigen::Vector and call the overloaded function
        double normalizeAndGetSinkhornDistance(std::vector<double>& a, std::vector<double>& b) const;
        
        //normalize the data (add noise and sum to one) and call getSinkhornDistance
        double normalizeAndGetSinkhornDistance(const Eigen::VectorXd& a, const Eigen::VectorXd& b) const;
        
        //Sinkhorn algorithm, as in POT librairy, with small adjustments
        double getSinkhornDistance(const Eigen::VectorXd& a, const Eigen::VectorXd& b) const;
        
        void getOccAndFreeDistri(const Eigen::VectorXd&, Eigen::VectorXd&, Eigen::VectorXd&) const;

        //add a small constant and sum to one, when the input vector is a const, the 2nd passed by ref is normalized
        void addConstAndNormalize(const Eigen::VectorXd&, Eigen::VectorXd&) const;
        //add a small consant and sum to one
        void addConstAndNormalize(Eigen::VectorXd&) const;
        void getCommonCDF(Eigen::VectorXd&, const Eigen::VectorXd&, const Eigen::VectorXd&) const;       
        void getCDF(Eigen::VectorXd&, const Eigen::VectorXd&) const;       
        void getCDF(std::vector<double>&, const std::vector<double>&) const;       
        void getCDF(Function&, const Eigen::VectorXd&) const;       
        void getCDF(Function&, const std::vector<double>&) const;       
        void drawRandomDistriFromGT(const std::vector<double>&, std::vector<double>&) const;       
        void setCostMatrixSquareDist();
};

#endif



