#ifndef POLICY_H
#define POLICY_H
#include <numeric>
#include <vector>
#include <math.h>
#include "explo_quality_aware/Constants.h"
#include "obs_stat/Tools.h"
#include "obs_stat/stat_point_type.h"

enum class PolicyChoice{NOBS, SIGMA, RANGE, ANGLES, ALL};

//Information Gain formulation between an initial and an updated PointStat (PointCloud of ObsStat)
class Policy {
    
    protected:
        PolicyChoice policy_choice;
        std::vector<double> n_obs_igs;
        std::vector<double> sigma_igs;
        std::vector<double> range_min_igs;
        std::vector<double> angles_igs;
        double n_obs_gain;
        double sigma_gain;
        double range_min_gain;
        double angles_gain;
        double gain;

        int getNAnglesFromV(uint32_t v) const {
            int n_angles = 0;
            for (uint8_t i=0; i<32; i++) {
                if (v & (1 << i)) {
                    ++n_angles;
                }
            }
            return n_angles;
        };
    //https://en.wikipedia.org/wiki/Central_limit_theorem
    double gain_nobs(int a, int b) const {
        if (a < 1) {
            return 0.;
        } 
        if (b <= 0) {
            return log(a/MIN_LOG_VAL);
        }
        if (a < b) {
            return 0;
        } 
        return log(static_cast<double>(a)/b);
    };
    //https://en.wikipedia.org/wiki/Central_limit_theorem
    //This is exactly the same than what we called quickmax
    double gain_sigma(double a, double b) const {
        //if the value is still 0, we have no gain
        if (a<MIN_LOG_VAL) {
            return 0.;
        } 
        //if not, but b is zero, clamp b
        if (b<MIN_LOG_VAL) {
            b = MIN_LOG_VAL;
        }
        return ((a>b) ? (log(a) - log(b)) : 0.);
    };
    //Could be updated here:
    //log((n+1)/n) - log(sum_(n+1)(range)/sum_(n)(range))
    double gain_range_min(double a, double b) const {
        //if a > b we have no gain
        if (a>b) {
            return 0;
        }
        return b - a;
    };

    double gain_angles(uint32_t a, uint32_t b) const {
        double na = static_cast<double>(getNAnglesFromV(a));
        double nb = static_cast<double>(getNAnglesFromV(b));
        //if the value is still 0, we have no gain
        if (na<MIN_LOG_VAL) {
            return 0.;
        } 
        //if not, but b is zero, clamp b
        if (nb<MIN_LOG_VAL) {
            nb = MIN_LOG_VAL;
        }
        return ((na>nb) ? (log(na) - log(nb)) : 0.);
    };
    
    protected:
    //And we just sum on all the elements of the grid per stat
    void update_igs(const PointStat& updated, const PointStat& initial) {
        n_obs_igs.push_back(gain_nobs(updated.n_obs, initial.n_obs));
        range_min_igs.push_back(gain_range_min(updated.r_min, initial.r_min));
        angles_igs.push_back(gain_angles(updated.angular_viewpoints_2d,initial.angular_viewpoints_2d));
        
        //For sigma, we need to deal with the particular case that it is not defined before the first observation.
        //Nonetheless, there is in reality an information gain to observe the voxel for the first time
        //If this is the first observation, initial.sigma_angle is -.5. 
        //For the log to be defined, we just add that +.5 to both values we compare
        if (fabs(initial.sigma_angle - INVALID_SIGMA) > MIN_POSITIVE_VALUE) {
            sigma_igs.push_back(gain_sigma(updated.sigma_angle, initial.sigma_angle));
        } else {
            sigma_igs.push_back(gain_sigma(updated.sigma_angle -INVALID_SIGMA, initial.sigma_angle -INVALID_SIGMA));
        }
    };

    public:
        Policy() {};
        Policy(PolicyChoice p): policy_choice(p) {};
        double getPolicyFromGrids(const grid_tools::StatGrid& updated_grid, const grid_tools::StatGrid& initial_grid, const grid_tools::GridSize& grid_size) {
            n_obs_igs.clear();
            sigma_igs.clear();
            range_min_igs.clear();
            angles_igs.clear();
            for (size_t ix{0}; ix < grid_size.nx; ix++) {
                for (size_t iy{0}; iy < grid_size.ny; iy++) {
                    for (size_t iz{0}; iz < grid_size.nz; iz++) {
                        if (updated_grid[ix][iy][iz].n_obs<MAX_NOBS_FOR_IG) {
                        update_igs(updated_grid[ix][iy][iz], initial_grid[ix][iy][iz]);
                        }
                    }
                }
            }
            n_obs_gain = std::accumulate(n_obs_igs.begin(), n_obs_igs.end(), 0.);
            sigma_gain = std::accumulate(sigma_igs.begin(), sigma_igs.end(), 0.);
            range_min_gain = std::accumulate(range_min_igs.begin(), range_min_igs.end(), 0.);
            angles_gain = std::accumulate(angles_igs.begin(), angles_igs.end(), 0.);
            switch (policy_choice) {
                case PolicyChoice::NOBS:
                    return  n_obs_gain;
                case PolicyChoice::SIGMA:
                    return  sigma_gain;
                case PolicyChoice::RANGE:
                    return  range_min_gain;
                case PolicyChoice::ANGLES:
                    return  angles_gain;
                case PolicyChoice::ALL:
                    std::cout << "Policy not correctly implemented, the sum does not make sense, not the same order of magnitude. Todo" << std::endl;
                    exit(EXIT_FAILURE);
                default:
                    std::cout << "Policy not implemented" << std::endl;
                    exit(EXIT_FAILURE);
            }
        }; 
};
 
#endif