
#include <vector>
#include <Eigen/Core>
#include "ceres/ceres.h"
#include "gflags/gflags.h"
#include "glog/logging.h"

struct PathSegment {
    Eigen::Vector3f A;
    Eigen::Vector3f B;
    Eigen::Vector3f AB,u;
    double length;
    
    PathSegment(const Eigen::Vector3f & A, const Eigen::Vector3f & B) :
        A(A), B(B) {
            update();
    }

    PathSegment() {}

    Eigen::Vector3f center() const {
        return (A+B)/2;
    }

    void update() {
        AB = B - A;
        length = AB.norm();
        u = AB / length;
    }

    template <class T>
        T px(const T & lambda) const {
            return T(A.x()) + lambda * T(u.x());
        }

    template <class T>
        T py(const T & lambda) const {
            return T(A.y()) + lambda * T(u.y());
        }

    template <class T>
        T pz(const T & lambda) const {
            return T(A.z()) + lambda * T(u.z());
        }



};


struct Path {
    Eigen::Vector3f start;
    std::vector<PathSegment> edges;
    Eigen::Vector3f end;

    Path() {}

    bool loadFromFile(const std::string & fname) {
        float x1,y1,z1,x2,y2,z2,x3,y3,z3;
        int i;
        FILE * fp = fopen(fname.c_str(),"r");
        char buffer[1024];
        // search the first line
        while (!feof(fp)) {
            if (fgets(buffer,1023,fp) == NULL) {
                return false;
            }
            if (sscanf(buffer," %d %e %e %e %e %e %e %e %e %e ",&i,
                        &x1,&y1,&z1, &x2,&y2,&z2, &x3,&y3,&z3) == 10) {
                break;
            }
        }
        start << x1,y1,z1;
        while (!feof(fp)) {
            if (fgets(buffer,1023,fp) == NULL) {
                break;
            }
            if (sscanf(buffer," %d %e %e %e %e %e %e %e %e %e ",&i,
                    &x1,&y1,&z1, &x2,&y2,&z2, &x3,&y3,&z3) != 10) {
                break;
            }
            if (i != 0) {
                break;
            }
            Eigen::Vector3f P,A,B; 
            P << x1,y1,z1;
            A << x2,y2,z2;
            B << x3,y3,z3;
            PathSegment ps(A,B);
            edges.push_back(ps);
        }
        if (edges.size()<2) {
            return false;
        }
        // The last segment is actually not a segment but the last point
        const PathSegment & pse = edges[edges.size()-1];
        assert(pse.length < 1e-3);
        end = pse.center();
        edges.resize(edges.size()-1);
        return true;
    }

};

template <class T> 
    T sqr(T x) {return x*x;}

struct PathCost {
    PathCost(const Path & path, unsigned int index, double weight=1)
        : path(path), index(index), weight(weight) { 
        }

    template <typename T>
        bool operator()(const T* const P1, 
                T* residuals) const {
            assert(index < path.edges.size());
            T x2 = path.edges[index].px(P1[0]); 
            T y2 = path.edges[index].py(P1[0]);
            T z2 = path.edges[index].pz(P1[0]);
            T x1,y1,z1;
            if (index==0) { 
                x1 = T(path.start.x()); 
                y1 = T(path.start.y()); 
                z1 = T(path.start.z());
            } else if (index == path.edges.size()-1) {
                x1 = T(path.end.x()); 
                y1 = T(path.end.y()); 
                z1 = T(path.end.z());
            } else {
                return false;
            }

            residuals[0] = sqr(x1-x2)+ sqr(y1-y2)+sqr(z1-z2);
            residuals[0] *= T(weight);

            return true;
        }

    template <typename T>
        bool operator()(const T* const P1, const T*P2, 
                T* residuals) const {
            assert(index < path.edges.size()-1);
            T x1 = path.edges[index].px(P1[0]); 
            T y1 = path.edges[index].py(P1[0]);
            T z1 = path.edges[index].pz(P1[0]);
            T x2 = path.edges[index+1].px(P2[0]); 
            T y2 = path.edges[index+1].py(P2[0]);
            T z2 = path.edges[index+1].pz(P2[0]);

            residuals[0] = sqr(x1-x2)+ sqr(y1-y2)+sqr(z1-z2);
            residuals[0] *= T(weight);

            return true;
        }


    const Path & path;
    unsigned int index;
    double weight;
};


struct SegmentCost {
    SegmentCost(const PathSegment & ps, double weight=1) : 
        ps(ps), weight() {}

    template <typename T>
        bool operator()(const T* const P, 
                T* residuals) const {
            if (P[0] < T(0)) {
                residuals[0] = T(weight)*sqr(P[0]);
            } else if (P[0] > T(ps.length)) {
                residuals[0] = T(weight)*sqr(P[0]-T(ps.length));
            } else {
                residuals[0] = T(0);
            }
            return true;
        }
    const PathSegment & ps;
    double weight;
};


using namespace ceres;




int main(int argc, char * argv[]) {

    if (argc < 2) {
        printf("Usage: %s <pathfile>\n",argv[0]);
        return -1;
    }
    Path path;
    if (!path.loadFromFile(argv[1])) {
        return -1;
    }


    Problem problem;
    Solver::Options options;
    options.function_tolerance = 1e-3;
    options.parameter_tolerance = 1e-4;
    options.max_num_iterations = 1000;
    options.minimizer_progress_to_stdout = false;
    options.num_threads = 1;
    options.eta = 1e-2;
    options.max_solver_time_in_seconds = 1e32;
    options.use_nonmonotonic_steps = false;
    // options->minimizer_type = ceres::LINE_SEARCH;

    CHECK(StringToTrustRegionStrategyType("levenberg_marquardt",
                &options.trust_region_strategy_type));
    CHECK(StringToDoglegType("traditional_dogleg", &options.dogleg_type));
    options.use_inner_iterations = false;

    double P[path.edges.size()];
    LossFunction *loss_function = NULL;
    CostFunction *cost_function = NULL;
    
    for (size_t i=0;i<path.edges.size();i++) {
        P[i] = path.edges[i].length/2;
        if (i == 0) {
            cost_function = new ceres::AutoDiffCostFunction<PathCost, 1, 1>(new PathCost(path, i, 1));
            problem.AddResidualBlock(cost_function,loss_function,P+i);
        } else if (i == path.edges.size()-1) {
            cost_function = new ceres::AutoDiffCostFunction<PathCost, 1, 1>(new PathCost(path, i, 1));
            problem.AddResidualBlock(cost_function,loss_function,P+i);
        } else {
            cost_function = new ceres::AutoDiffCostFunction<PathCost, 1, 1, 1>(new PathCost(path, i, 1));
            problem.AddResidualBlock(cost_function,loss_function,P+i,P+i+1);
        }
        cost_function = new ceres::AutoDiffCostFunction<SegmentCost,1,1>(new SegmentCost(path.edges[i],1e3));
        problem.AddResidualBlock(cost_function,loss_function,P+i);
    }
    Solver::Summary summary;
    Solve(options, &problem, &summary);
    std::cout << summary.FullReport() << "\n";

    FILE * fp=fopen("optpath","w");
    fprintf(fp,"0 %e %e %e\n",path.start.x(),path.start.y(),path.start.z());
    for (size_t i=0;i<path.edges.size();i++){
        fprintf(fp,"0 %e %e %e\n",
                path.edges[i].px(P[i]),
                path.edges[i].py(P[i]),
                path.edges[i].pz(P[i]));
    }
    fprintf(fp,"0 %e %e %e\n",path.end.x(),path.end.y(),path.end.z());
    fprintf(fp,"\n\n");
    fclose(fp);

    return 0;
}


