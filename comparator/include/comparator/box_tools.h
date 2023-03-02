#ifndef BOX_TOOLS_H
#define BOX_TOOLS_H
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <vector>
#include <array>
#include <fstream>
#include <iostream>
#include <random>


//This class compares 2 datasets ground-truth and reconstruction
//The datasets are composed of :images, image list, and the bounding-box and resolution of the images
//This comparator loads the same sample of ground-truth and reconstruction, compares it, and save the metrics of the comparaison to disk 
class BoxTools {

    //protected:
    public:
        int box_size_;
        int img_width_;
        int img_height_;

        //load a succession of box_size images into a 3D Mat
        //this is a test function only, the images must all have the same size
        void load3DMat(cv::Mat_<uint8_t>& mat, const std::vector<std::string >& imglist) const;
        

        void getVectorFromBox(const std::vector<cv::Range>&, const cv::Mat_<uint8_t>&, std::vector<double>&) const;
        void getVectorFromSingleBox(const cv::Mat_<double>&, std::vector<double>&) const;
        void getVectorFromSingleBox(const cv::Mat_<uint8_t>&, std::vector<double>&) const;
        void get3DProbaImageFromRange(const std::vector<cv::Range>&, const cv::Mat_<uint8_t>&, cv::Mat_<double>&) const;
        void get3DProbaImageFromSingleCube(const cv::Mat_<uint8_t>&, cv::Mat_<double>&) const;
        void getVectorOfImagesFromCube(const cv::Mat_<double>&, std::vector<cv::Mat_<double>>&, int) const;
        void getCubeFromVectorOfImages(cv::Mat_<double>&, const std::vector<cv::Mat_<double>>&, int) const;
        
        
        cv::Size getRandomKSize() const;
        double getRandomKSigma() const;
        void addRandomNoise(std::vector<double>&, double) const;
        void normMinMax(std::vector<double>&) const;
        void drawRandomDistri(size_t, std::vector<double>&) const;
        
        
        //returns true if thw box contains only zeros 
        bool isBoxOnlyZeros(const std::vector<cv::Range>&, const cv::Mat_<uint8_t>&) const;        
        //if the box contains at least <level> of not empty points, returns true
        bool isBoxNotEmptyWithLevel(const std::vector<cv::Range>&, const cv::Mat_<uint8_t>&, double) const;        
        //args are the probability vector, divergence from the unknown to consider a box observed
        bool isBoxObserved(const std::vector<double>&, double) const;        
        double getBoxOccupancyRate(const std::vector<cv::Range>&, const cv::Mat_<uint8_t>&) const;        
        
        
        //add noise to a ground truth box and retuns a vector, to use to calculate the wd limit of the dataset
        //if the last three arguments are not provided, we pick a random sigma and add uniform noise
        void getNoisyGtVector(const std::vector<cv::Range>&, const cv::Mat_<uint8_t>&, std::vector<double>&, double sigma=0.5, int ksize=7, double additionnal_uniform_noise=0) const;
        void getNoisyGtVectorDebug(const std::vector<cv::Range>&, const cv::Mat_<uint8_t>&, std::vector<double>&, double uniform_noise_level=.1, bool use_fixed_sigma=false, double sigma=0.) const;
        
        //debug
        //last args are: fname, dir
        void saveSampleCubeToImg(const std::vector<cv::Range>&, const cv::Mat_<uint8_t>&, const std::string&, const std::string&) const;
        void saveSingleCubeToImg(const cv::Mat_<uint8_t>&, const std::string&, const std::string&) const;
        void saveSingleCubeToImg(const cv::Mat_<double>&, const std::string&, const std::string&) const;
        void saveBoxVectorToImg(const std::vector<double>&, const std::string&, const std::string&) const;
        
        
        //utility function to normalize double to color and the inverse
        uint8_t getNormalizedPixValue(double v) const;
        double pix2proba(uint8_t) const;
        
        //utility funciton to check we can write the output file
        void checkFile(std::ofstream&) const;

    public:
        BoxTools(){};
        //Constructors for tests only:
        //Instantiate an object with a box_size only, useful to tests functions using vectors only
        BoxTools(int box_size);
        //Instantiate an object with box_size, img_width and img_height, useful to load images and do tests on the 3dMat
        //All the images must have the same shape, and only box_size images will be processed
        BoxTools(int box_size, int img_width, int img_height);
};

#endif
