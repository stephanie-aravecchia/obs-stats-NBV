#include "comparator/box_tools.h"
#include <ros/ros.h>
#include <assert.h>
#include <algorithm>

using namespace std;

BoxTools::BoxTools(int box_size) : box_size_(box_size) {}
BoxTools::BoxTools(int box_size, int img_width, int img_height) : 
    box_size_(box_size), img_width_(img_width), img_height_(img_height) {}

//This is just for tests
void BoxTools::load3DMat(cv::Mat_<uint8_t>& mat, const std::vector<std::string >& imglist) const {
    for (int k{0}; k < box_size_; k++) {
        assert(k>=0);
        assert(k < imglist.size());
        cv::Mat_<uint8_t> tmp = cv::imread(imglist[k], cv::IMREAD_UNCHANGED);
        assert(tmp.size().width == img_width_);
        assert(tmp.size().height == img_height_);
        for (int i{0}; i < tmp.size().width; i++) {
            for (int j{0}; j < tmp.size().height; j++) {
                //remainder: opencv: tmp() refers to tmp(row, col)
                mat(i, j, k) = tmp(j,i);
                //This is the same that before all the changes
            }
        }
    }
}
void BoxTools::getVectorFromBox(const std::vector<cv::Range>& ranges, 
            const cv::Mat_<uint8_t>& refMat, vector<double>& vect) const {
    vect.clear();
    for (int k{ranges[2].start}; k < ranges[2].end; k++) {//img first
            for (int j{ranges[1].start}; j < ranges[1].end; j++) {//then row
                for (int i{ranges[0].start}; i < ranges[0].end; i++) {//then col
                assert(i>=0);
                assert(j>=0);
                assert(k>=0);
                assert(i<img_width_);
                assert(j<img_height_);
                assert(k<box_size_);
                vect.push_back(pix2proba(refMat(i,j,k)));//this is (col, row, img)
            }
        }
    }
}
void BoxTools::getVectorFromSingleBox(const cv::Mat_<double>& unitBox, vector<double>& vect) const {
    vect.clear();
    for (int k{0}; k < box_size_; k++) {//img first
            for (int j{0}; j < box_size_; j++) {//then row
                for (int i{0}; i < box_size_; i++) {//then col
                 vect.push_back(unitBox(i,j,k));//this is (col, row, img)
            }
        }
    }
}
void BoxTools::getVectorFromSingleBox(const cv::Mat_<uint8_t>& unitBox, vector<double>& vect) const {
    vect.clear();
    cout << "int box" << endl;
    for (int k{0}; k < box_size_; k++) {//img first
            for (int j{0}; j < box_size_; j++) {//then row
                for (int i{0}; i < box_size_; i++) {//then col
                vect.push_back(unitBox(i,j,k));//this is (col, row, img)
            }
        }
    }
}
void BoxTools::get3DProbaImageFromSingleCube(const cv::Mat_<uint8_t>& refMat, cv::Mat_<double>& out) const {
    assert(out.total() == box_size_ * box_size_ * box_size_);
    for (int k{0}; k < box_size_; k++) {//img first
            for (int j{0}; j < box_size_; j++) {//then row
                for (int i{0}; i < box_size_; i++) {//then col
                out(i,j,k)=(pix2proba(refMat(i,j,k)));//this is (col, row, img)
            }
        }
    }
}
void BoxTools::get3DProbaImageFromRange(const std::vector<cv::Range>& ranges, const cv::Mat_<uint8_t>& refMat, cv::Mat_<double>& out) const {
    assert(out.total() == box_size_ * box_size_ * box_size_);
    assert((ranges[2].end - ranges[2].start) == box_size_);
    assert((ranges[1].end - ranges[1].start) == box_size_);
    assert((ranges[0].end - ranges[0].start) == box_size_);
    for (int k{ranges[2].start}; k < ranges[2].end; k++) {//img first
            for (int j{ranges[1].start}; j < ranges[1].end; j++) {//then row
                for (int i{ranges[0].start}; i < ranges[0].end; i++) {//then col
                assert(i>=0);
                assert(j>=0);
                assert(k>=0);
                assert(i<img_width_);
                assert(j<img_height_);
                assert(k<box_size_);
                out(i-ranges[0].start,j-ranges[1].start,k-ranges[2].start)=(pix2proba(refMat(i,j,k)));//this is (col, row, img)
            }
        }
    }
}

void BoxTools::getVectorOfImagesFromCube(const cv::Mat_<double>& pcube, std::vector<cv::Mat_<double>>& vect, int order) const {
    assert(vect.size() == box_size_);
    const int mat_sz[2] = {box_size_,box_size_};
    switch(order) {
        case 0://k slicing
            for (int k{0}; k < box_size_; k++) {
                vect[k] = cv::Mat_<double>(2, mat_sz);
                for (int j{0}; j < box_size_; j++) {
                    for (int i{0}; i < box_size_; i++) {
                        vect[k](i,j) = pcube(i,j,k);//here
                    }
                }
            }
            break;
        case 1://j slicing
            for (int k{0}; k < box_size_; k++) {
                vect[k] = cv::Mat_<double>(2, mat_sz);
                for (int j{0}; j < box_size_; j++) {
                    for (int i{0}; i < box_size_; i++) {
                        vect[k](i,j) = pcube(i,k,j);//here
                    }
                }
            }
            break;
        case 2://i slicing
            for (int k{0}; k < box_size_; k++) {
                vect[k] = cv::Mat_<double>(2, mat_sz);
                for (int j{0}; j < box_size_; j++) {
                    for (int i{0}; i < box_size_; i++) {
                        vect[k](i,j) = pcube(k,i,j);//here
                    }
                }
            }
            break;
        default:
            cerr << "Switch case not implemented" << endl;
            break;
    }
}
void BoxTools::getCubeFromVectorOfImages(cv::Mat_<double>& pcube, const std::vector<cv::Mat_<double>>& vect, int order) const {
    assert(vect.size() == box_size_);
    assert(pcube.total() == (box_size_*box_size_*box_size_));
    switch(order) {
        case 0://k slicing
            for (int k{0}; k < box_size_; k++) {
                for (int j{0}; j < box_size_; j++) {
                    for (int i{0}; i < box_size_; i++) {
                        pcube(i,j,k) = vect[k](i,j);//here
                    }
                }
            }
            break;
        case 1://j slicing
            for (int k{0}; k < box_size_; k++) {
                for (int j{0}; j < box_size_; j++) {
                    for (int i{0}; i < box_size_; i++) {
                        pcube(i,k,j) = vect[k](i,j);//here
                    }
                }
            }
            break;
        case 2://i slicing
            for (int k{0}; k < box_size_; k++) {
                for (int j{0}; j < box_size_; j++) {
                    for (int i{0}; i < box_size_; i++) {
                        pcube(k,i,j) = vect[k](i,j);//here
                    }
                }
            }
            break;
        default:
            cerr << "Switch case not implemented" << endl;
            break;
    }
}

//Get a random kernel size, with a user chosen distribution
cv::Size BoxTools::getRandomKSize() const {
    static std::random_device rd;
    static std::default_random_engine engine(rd());
    std::uniform_int_distribution<int> distri(1, 100);
    int m = distri(engine);
    //TODO define that elsewhere:
    static int n3 = 55; //ksize = (3,3)
    static int n1 = 20; //ksize = (1,1)
    static int n5 = 20; //ksize = (5,5)
    static int n7 = 5; //ksize = (7,7)
    assert(n1 + n3 + n5 + n7 == 100);
    if (m <= n3) {
        return cv::Size(3,3);
    } else if (m <= n3+n1) {
        return cv::Size(1,1);
    } else if (m <= n3+n1+n5) {
        return cv::Size(5,5);
    } else {
        return cv::Size(7,7);
    }
}

//get a random ksigma, with a uniform distribution 
double BoxTools::getRandomKSigma() const {
    //TODO define that elsewhere:
    static double ksigma = .1;
    static double msigma = .5;
    //in the process we are applying twice the convolution, to have at the end the correct sigma, we divide it by 2
    static std::random_device rd;
    static std::default_random_engine engine(rd());
    static std::normal_distribution<double> distribution(msigma, ksigma/2);
    return distribution(engine);
}

void BoxTools::addRandomNoise(std::vector<double>& vect, double noise_level) const {
    assert(vect.size()>0);
    static std::random_device rd;
    static std::default_random_engine engine(rd());
    std::uniform_real_distribution<double> random_dist(0.0,noise_level);
    auto add_noise = [&engine, &random_dist](double& a){
        double noise = random_dist(engine);
        a += std::max(0.0, std::min(noise, 1.0));
        };
    std::for_each(vect.begin(), vect.end(),add_noise);
}

void BoxTools::normMinMax(std::vector<double> &vect) const {
    //we normalize vect (min max)
    assert(vect.size() > 0);
    double min = *std::min_element(vect.begin(), vect.end());
    double max = *std::max_element(vect.begin(), vect.end());
    assert((max - min) > 1e-6);
    auto norm = [min, max](double &x)
        { x = (x - min) / (max - min); };
    std::for_each(vect.begin(), vect.end(), norm);
}

void BoxTools::drawRandomDistri(size_t size, std::vector<double>& p) const {
    assert(size>0);
    p.resize(size);
    static std::random_device rd;
    static std::default_random_engine engine(rd());
    static std::uniform_real_distribution<double> random_dist(0.0,1.0);
    for (size_t i{0}; i<size; ++i){
        p[i] = random_dist(engine);
    }
} 

bool BoxTools::isBoxOnlyZeros(const std::vector<cv::Range>& ranges, const cv::Mat_<uint8_t>& refMat) const {
    assert((ranges[2].end - ranges[2].start) == box_size_);
    assert((ranges[1].end - ranges[1].start) == box_size_);
    assert((ranges[0].end - ranges[0].start) == box_size_);
    for (int k{ranges[2].start}; k < ranges[2].end; k++) { //img first
        for (int j{ranges[1].start}; j < ranges[1].end; j++) { //then row
            for (int i{ranges[0].start}; i < ranges[0].end; i++) { //then col
                assert(i >= 0);
                assert(j >= 0);
                assert(k >= 0);
                assert(i < img_width_);
                assert(j < img_height_);
                assert(k < box_size_);
                if (refMat(i, j, k) > 0) {
                    return false;
                }
            }
        }
    }
    return true;
}
bool BoxTools::isBoxNotEmptyWithLevel(const std::vector<cv::Range>& ranges, const cv::Mat_<uint8_t>& refMat, double level) const {
    assert((ranges[2].end - ranges[2].start) == box_size_);
    assert((ranges[1].end - ranges[1].start) == box_size_);
    assert((ranges[0].end - ranges[0].start) == box_size_);
    int thres{floor(box_size_*box_size_*box_size_*level)};
    int counter{0};
    for (int k{ranges[2].start}; k < ranges[2].end; k++) { //img first
        for (int j{ranges[1].start}; j < ranges[1].end; j++) { //then row
            for (int i{ranges[0].start}; i < ranges[0].end; i++) { //then col
                assert(i >= 0);
                assert(j >= 0);
                assert(k >= 0);
                assert(i < img_width_);
                assert(j < img_height_);
                assert(k < box_size_);
                if (refMat(i, j, k) > 0) {
                    ++counter;
                    if (counter>thres) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}
//if the max divergence is .1 to the unknown, we don't process the box
bool BoxTools::isBoxObserved(const std::vector<double>& rec_vector, double divergence) const {
    double minv = *std::min_element(rec_vector.begin(), rec_vector.end());
    double maxv = *std::max_element(rec_vector.begin(), rec_vector.end());
    if ((minv > (0.5-divergence)) && (maxv < (0.5+divergence))) {
        return false;
    } else {
        return true;
    }
} 
double BoxTools::getBoxOccupancyRate(const std::vector<cv::Range>& ranges, const cv::Mat_<uint8_t>& refMat) const {
    assert((ranges[2].end - ranges[2].start) == box_size_);
    assert((ranges[1].end - ranges[1].start) == box_size_);
    assert((ranges[0].end - ranges[0].start) == box_size_);
    int counter{0};
    for (int k{ranges[2].start}; k < ranges[2].end; k++) { //img first
        for (int j{ranges[1].start}; j < ranges[1].end; j++) { //then row
            for (int i{ranges[0].start}; i < ranges[0].end; i++) { //then col
                assert(i >= 0);
                assert(j >= 0);
                assert(k >= 0);
                assert(i < img_width_);
                assert(j < img_height_);
                assert(k < box_size_);
                if (refMat(i, j, k) == 255) {
                    ++counter;
                }
            }
        }
    }
    return static_cast<double>(counter)/(box_size_*box_size_*box_size_);
}

void BoxTools::getNoisyGtVector(const std::vector<cv::Range>& ranges, 
            const cv::Mat_<uint8_t>& refMat, vector<double>& vect, double sigma, int ksize, double additionnal_uniform_noise) const {
    
    //we create a 3d mat noisyGT, we initialize it with a copy of the gt selected cube
    const int mat_sz[3] = {box_size_,box_size_,box_size_};
    cv::Mat_<double> noisygt = cv::Mat_<double>(3, mat_sz);
    get3DProbaImageFromRange(ranges, refMat, noisygt);
    //we apply noise (3 times 2d_filter, in all directions)
    //first, create two vectors of 10 elements containing 10x10 images
    //imgvect will contain the images before we apply the filter
    //output vect will be the blurred images
    //both are updated at each step
    std::vector<cv::Mat_<double>> imgvect;
    imgvect.resize(box_size_);
    std::vector<cv::Mat_<double>> outputvect;
    outputvect.resize(box_size_);
    const int im_sz[2] = {box_size_,box_size_};
    for (int k{0}; k < box_size_; k++) {
        outputvect[k] = cv::Mat_<double>(2, im_sz);
        imgvect[k] = cv::Mat_<double>(2, im_sz);
    }
    getVectorOfImagesFromCube(noisygt, imgvect,0);
    cv::Size kernelsize = cv::Size(7,7);
    //then we slice the cube in the appropriate direction into a vector of image, we store taht in imgvect
    //do the convolution on all the images (temporary result saved in outputvect)
    //we now transform the convoluted images back to a 3dmat noisygt
    //we repeat on the other dimensions
    
    //the "slicing" will be z
    // done before the if getVectorOfImagesFromCube(noisygt, imgvect,0);
    for (int k{0}; k < box_size_; k++) {
        cv::GaussianBlur(imgvect[k], outputvect[k],kernelsize,sigma, 0.0, cv::BORDER_CONSTANT);
    }
    getCubeFromVectorOfImages(noisygt, outputvect,0);    
    
    //the "slicing" will be j
    getVectorOfImagesFromCube(noisygt, imgvect,1);
    for (int k{0}; k < box_size_; k++) {
        cv::GaussianBlur(imgvect[k], outputvect[k],kernelsize,sigma, 0.0, cv::BORDER_CONSTANT);
    }
    getCubeFromVectorOfImages(noisygt, outputvect,1);    
    
    //the "slicing" will be i
    getVectorOfImagesFromCube(noisygt, imgvect,2);
    for (int k{0}; k < box_size_; k++) {
        cv::GaussianBlur(imgvect[k], outputvect[k],kernelsize,sigma, 0.0, cv::BORDER_CONSTANT);
    }
    getCubeFromVectorOfImages(noisygt, outputvect,2);    
    //we get the vector back from this noisy cube
    getVectorFromSingleBox(noisygt, vect);
    normMinMax(vect);
    //and finally, we add a uniform noise on the complete vector, if uniform_noise_level is provided
    if (fabs(additionnal_uniform_noise)>1e-6) {
        addRandomNoise(vect, additionnal_uniform_noise);
    } 
} 

void BoxTools::getNoisyGtVectorDebug(const std::vector<cv::Range>& ranges, 
            const cv::Mat_<uint8_t>& refMat, vector<double>& vect, double uniform_noise_level, bool use_fixed_sigma, double sigma) const {
    bool debug = false;
    std::string dir;
    char suffix[40];
    if (debug) {
        dir = "/home/saravecchia/data/bags/husky/xp_metrics/";
        sprintf(suffix,"%.02f_%.02f", sigma, uniform_noise_level);
    }
    //If the box contains only 0, or sigma = 0, we skip this step:
    if (!(isBoxOnlyZeros(ranges, refMat))) {
        //we create a 3d mat noisyGT, we initialize it with a copy of the gt selected cube
        const int mat_sz[3] = {box_size_,box_size_,box_size_};
        cv::Mat_<double> noisygt = cv::Mat_<double>(3, mat_sz);
        get3DProbaImageFromRange(ranges, refMat, noisygt);
        //We save the original GT, with the two methods to make sure they give the same results:
        if (debug) {
            char fname[60];
            sprintf(fname, "cube2img_orig_gt_%s", suffix);
            saveSampleCubeToImg(ranges, refMat, fname, dir);
            getVectorFromSingleBox(noisygt, vect);
            char fname1[60];
            sprintf(fname1, "vect2img_orig_gt_%s", suffix);
            saveBoxVectorToImg(vect, fname1, dir);
        } 
        //we apply noise (3 times 2d_filter, in all directions)
        //first, create two vectors of 10 elements containing 10x10 images
        //imgvect will contain the images before we apply the filter
        //output vect will be the blurred images
        //both are updated at each step
        std::vector<cv::Mat_<double>> imgvect;
        imgvect.resize(box_size_);
        std::vector<cv::Mat_<double>> outputvect;
        outputvect.resize(box_size_);
        const int im_sz[2] = {box_size_,box_size_};
        for (int k{0}; k < box_size_; k++) {
            outputvect[k] = cv::Mat_<double>(2, im_sz);
            imgvect[k] = cv::Mat_<double>(2, im_sz);
        }
        //if no sigma is provided (zero)
        if (!use_fixed_sigma) {
            sigma = getRandomKSigma();
        }
        getVectorOfImagesFromCube(noisygt, imgvect,0);
        if (sigma>1e-6){ 
            cv::Size kzise = cv::Size(7,7);
            //then we slice the cube in the appropriate direction into a vector of image, we store taht in imgvect
            //do the convolution on all the images (temporary result saved in outputvect)
            //we now transform the convoluted images back to a 3dmat noisygt
            //we repeat on the other dimensions
            
            //the "slicing" will be z
            // done before the if getVectorOfImagesFromCube(noisygt, imgvect,0);
            for (int k{0}; k < box_size_; k++) {
                cv::GaussianBlur(imgvect[k], outputvect[k],kzise,sigma, 0.0, cv::BORDER_CONSTANT);
            }
            getCubeFromVectorOfImages(noisygt, outputvect,0);    
            
            //the "slicing" will be j
            getVectorOfImagesFromCube(noisygt, imgvect,1);
            for (int k{0}; k < box_size_; k++) {
                cv::GaussianBlur(imgvect[k], outputvect[k],kzise,sigma, 0.0, cv::BORDER_CONSTANT);
            }
            getCubeFromVectorOfImages(noisygt, outputvect,1);    
            
            //the "slicing" will be i
            getVectorOfImagesFromCube(noisygt, imgvect,2);
            for (int k{0}; k < box_size_; k++) {
                cv::GaussianBlur(imgvect[k], outputvect[k],kzise,sigma, 0.0, cv::BORDER_CONSTANT);
            }
            getCubeFromVectorOfImages(noisygt, outputvect,2);    
        }    
        //we get the vector back from this noisy cube
        getVectorFromSingleBox(noisygt, vect);
        normMinMax(vect);
    }
    else {
        //if the box contains only zeros, we initialize the vector with zeros
        vect = std::vector<double>(box_size_*box_size_*box_size_,0);
    }
    //and finally, we add a uniform noise on the complete vector, if uniform_noise_level is provided
    if (fabs(uniform_noise_level)>1e-6) {
        addRandomNoise(vect, uniform_noise_level);
    } else {
    }
    //We save the noisy gt, with the two methods to make sure they give the same results:
    if (debug) {
        char fname2[60];
        sprintf(fname2, "noisy_gt_%s", suffix);
        saveBoxVectorToImg(vect, fname2, dir);
    }
}
void BoxTools::saveSingleCubeToImg(const cv::Mat_<uint8_t>& refMat, const std::string& fname, const std::string& dir) const {
    cv::Mat_<uint8_t> slice_img = cv::Mat_<uint8_t>(box_size_, box_size_);
    for (int k{0}; k < box_size_; k++) {
        for (int i{0}; i<box_size_; i++) {
            for (int j{0}; j<box_size_; j++) {                    
                //(row, col)            //(col, row, img)
                slice_img(j,i) = refMat(i,j,k);
            }
        }
        char suffix[10];
        sprintf(suffix,"_%02d.png", k);
        cv::imwrite(dir + fname + suffix, slice_img);
    }
}
void BoxTools::saveSingleCubeToImg(const cv::Mat_<double>& refMat, const std::string& fname, const std::string& dir) const {
    cv::Mat_<uint8_t> slice_img = cv::Mat_<uint8_t>(box_size_, box_size_);
    for (int k{0}; k < box_size_; k++) {
        for (int i{0}; i<box_size_; i++) {
            for (int j{0}; j<box_size_; j++) {                    
                //(row, col)                             //(col, row, img)
                slice_img(j,i) = getNormalizedPixValue(refMat(i,j,k));
            }
        }
        char suffix[10];
        sprintf(suffix,"_%02d.png", k);
        cv::imwrite(dir + fname + suffix, slice_img);
    }
}


void BoxTools::saveSampleCubeToImg(const std::vector<cv::Range>& ranges, const cv::Mat_<uint8_t>& refMat, const std::string& fname, const std::string& dir) const {
    int rows = ranges[1].end-ranges[1].start;
    int cols = ranges[0].end-ranges[0].start;
    cv::Mat_<uint8_t> slice_img = cv::Mat_<uint8_t>(rows, cols);
    for (int k{ranges[2].start}; k < ranges[2].end; k++) {
        for (int i{ranges[0].start}; i<ranges[0].end; i++) {
            for (int j{ranges[1].start}; j<ranges[1].end; j++) {                    
                assert(i>=0);
                assert(j>=0);
                assert(k>=0);
                assert(i<img_width_);
                assert(j<img_height_);
                assert(k<box_size_);
                assert(j-ranges[1].start>=0);                
                assert(i-ranges[0].start>=0);
                assert(j-ranges[1].start<rows);
                assert(i-ranges[0].start<cols);
                //(row, col)                                                    //(col, row, img)
                slice_img(j-ranges[1].start,i-ranges[0].start) = refMat(i,j,k);
            }
        }
        char suffix[10];
        sprintf(suffix,"_%02d.png", k);
        cv::imwrite(dir + fname + suffix, slice_img);
    }
}
void BoxTools::saveBoxVectorToImg(const std::vector<double>& imgVect, const std::string& fname, const std::string& dir) const {
    int cube_size = box_size_*box_size_*box_size_;
    int img_stride = box_size_*box_size_;
    int row_stride = box_size_;
    assert(cube_size == imgVect.size());
    cv::Mat_<uint8_t> slice_img = cv::Mat_<uint8_t>(box_size_, box_size_);
    for (int k{0}; k < box_size_; k++) {//img first
        for (int j{0}; j < box_size_; j++) {//then row
            for (int i{0}; i < box_size_; i++) {//then col
                //this is (row,col)
                slice_img(j,i) = getNormalizedPixValue(imgVect.at(k*img_stride+j*row_stride+i));
            }
        } 
        char suffix[10];
        sprintf(suffix,"_%02d.png", k);
        cv::imwrite(dir + fname + suffix, slice_img);

    }
}

uint8_t BoxTools::getNormalizedPixValue(double v) const {
    return static_cast<uint8_t>(std::max(0., std::min((v*255), 255.)));
}
double BoxTools::pix2proba(uint8_t v) const {
    return std::max(0., std::min(v/255.,1.0));
}

void BoxTools::checkFile(ofstream & f) const {
    if (!f) {
        cerr << "Output file cannot be opened" << endl;
        exit(EXIT_FAILURE);
    }
}