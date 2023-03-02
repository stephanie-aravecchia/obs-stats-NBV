#include <assert.h>
#include "functions/Function.h"

using namespace std;

void Function::print(FILE * fp) const {
    FDB::const_iterator it = vals.begin();
    for (;it != vals.end();it++) {
        fprintf(fp,"%f\t%f\n",it->first,it->second);
    }
}

void Function::print(const char * fname) const {
    FILE * fp = fopen(fname,"w");
    if (fp == NULL) {
        perror("Function::print");
        return;
    }
    print(fp);
    fclose(fp);
}

pair<double,double> Function::fmin() const 
{
	FDB::const_iterator it = vals.begin();
	if (it == vals.end()) return pair<double,double>(0.0,0.0);
	double m = it->second;
	double x = it->first;
	for(;it!=vals.end();it++) 
		if (it->second < m) {
			m = it->second;
			x = it->first;
		}
	return pair<double,double>(x,m);
}

pair<double,double> Function::fmax() const 
{
	FDB::const_iterator it = vals.begin();
	if (it == vals.end()) return pair<double,double>(0.0,0.0);
	double m = it->second;
	double x = it->first;
	for(;it!=vals.end();it++) 
		if (it->second > m) {
			m = it->second;
			x = it->first;
		}
	return pair<double,double>(x,m);
}


Function Function::map(const Function & f) const
{
	Function res;
	FDB::const_iterator it = vals.begin();
	for (;it != vals.end();it++)
		res.set(it->first,f(it->second));
	return res;
}

Function::Support Function::minima(double dead_zone) const 
{
    Function d = derivative();
    Support s = d.zeros(), res; 

    std::map<double,double> mins;
    // First, store all candidates in a map, ordered by value
    for (Support::const_iterator it = s.begin(); it != s.end(); it++) {
        double d2 = d.derivativeAt(*it);
        if (d2 > 0) {
            mins.insert(std::map<double,double>::value_type((*this)(*it),*it));
        }
    }
    // Then remove anything closer than deadzone
    for (std::map<double,double>::iterator it = mins.begin();
            it != mins.end(); it++) {
        res.insert(it->second);
        std::vector<std::map<double,double>::iterator> to_erase;
        for (std::map<double,double>::iterator jt = mins.begin();
                jt != mins.end(); jt++) {
            if (jt == it) continue;
            if (fabs(jt->second-it->second) < dead_zone) {
                to_erase.push_back(jt);
            }
        }
        for (size_t i=0;i<to_erase.size();i++) {
            mins.erase(to_erase[i]);
        }
    }
    return res;
}

Function::Support Function::maxima(double dead_zone) const 
{
    Function d = derivative();
    Support s = d.zeros(), res; 

    typedef std::map< double,double,std::greater<double> > Maxs;
    Maxs maxs;
    // First, store all candidates in a map, ordered by value
    for (Support::const_iterator it = s.begin(); it != s.end(); it++) {
        double d2 = d.derivativeAt(*it);
        if (d2 < 0) {
            maxs.insert(Maxs::value_type((*this)(*it),*it));
        }
    }
    // Then remove anything closer than deadzone
    for (Maxs::iterator it = maxs.begin();
            it != maxs.end(); it++) {
        res.insert(it->second);
        std::vector<Maxs::iterator> to_erase;
        for (Maxs::iterator jt = maxs.begin();
                jt != maxs.end(); jt++) {
            if (jt == it) continue;
            if (fabs(jt->second-it->second) < dead_zone) {
                to_erase.push_back(jt);
            }
        }
        for (size_t i=0;i<to_erase.size();i++) {
            maxs.erase(to_erase[i]);
        }
    }
    return res;
}


Function::Support Function::zeros() const
{
	Support res;
    if (vals.empty()) {
        return res;
    }
	FDB::const_iterator it = vals.begin(),itp = it;
    if (it->second==0.0) {
        res.insert(it->first);
    }
    it++;
	for (;it != vals.end();it++,itp++) {
        if (it->second==0.0) {
            res.insert(it->first);
        } else if ((itp->second < 0) && (it->second > 0)) {
            res.insert(itp->first - itp->second * (it->first - itp->first) / (it->second - itp->second));
        } else if ((itp->second > 0) && (it->second < 0)) {
            res.insert(itp->first + itp->second * (it->first - itp->first) / (it->second - itp->second));
        }
    }
	return res;
}


Function Function::primitive(double x0, double dx) const
{
	Function res;
    double sum = 0;
#if 0
    double fx = xmin();
    double lx = xmax();
    for (double x = fx;x <= lx; x += dx) {
        sum += (*this)(x) * dx;
        res.set(x,sum);
    }
#else
	FDB::const_iterator it = vals.begin(), itp = it;
    res.set(it->first,sum);
    it++;
	for (;it != vals.end();it++,itp++) {
        sum += 0.5*(it->second+itp->second) * (it->first-itp->first);
        res.set(it->first,sum);
    }
#endif
    res = res - (*this)(x0);
	return res;
}

Function Function::derivative(double h) const
{
	Function res;
	if (h <= 0) {
		h = mindx()/10;
	}
	FDB::const_iterator it = vals.begin();
	for (;it != vals.end();it++) {
        res.set(it->first,derivativeAt(it->first,h));
	}
	return res;
}

double Function::derivativeAt(double x, double h) const
{
    if (h <= 0) {
        h = mindx()/10;
    }
    return ((*this)(x+h/2) - (*this)(x-h/2))/h;
}


#ifdef USE_GSL

#include "gsl/gsl_multifit.h"
Function Function::spline_approximation(double node_spacing,double sampling) {
    double fm = xmin();
    double fM = xmax();
    // First compute the number of nodes and then approximate dx
    if (size() <= 10) {
        return *this;
    }
    int n_nodes = std::min<int>(floor((fM-fm)/node_spacing)+1,size()-4);
    double dx = (fM-fm)/n_nodes;
    // printf("size %d n_nodes %d dx %f f [%f,%f]\n",int(size()),n_nodes,dx,fm,fM);
    assert((signed)size() >= n_nodes + 4);
    gsl_vector *P = gsl_vector_alloc(n_nodes+4);
    gsl_vector *B = gsl_vector_alloc(size());
    gsl_matrix *A = gsl_matrix_alloc(size(),n_nodes+4);
    gsl_matrix_set_zero(A);
    size_t i=0;
    for (const_iterator it=begin();it!=end();it++,i++) {
        assert(i < size());
        int node = floor((it->first - fm)/dx);
        double u = (it->first - (fm + node*dx))/dx;
        if (fabs(it->first - fM) < 1e-9) {
            node = n_nodes-1;
            u = 1;
        }
        assert(node < n_nodes);
        double pu[4] = {1, u, u*u, u*u*u};
        // Compute 
        // M=(1/6) * [1 4 1 0; -3 0 3 0 ; 3 -6 3 0; -1 3 -3 1];
        // [1 u u^2 u^3] * M
        double m[4] = {
            (1*pu[0] + -3*pu[1] +  3*pu[2] + -1*pu[3]) / 6,
            (4*pu[0] +  0*pu[1] + -6*pu[2] +  3*pu[3]) / 6,
            (1*pu[0] +  3*pu[1] +  3*pu[2] + -3*pu[3]) / 6,
            (0*pu[0] +  0*pu[1] +  0*pu[2] +  1*pu[3]) / 6
        };
        // printf("i %d node %d x %f y %f u %f m %f %f %f %f\n",int(i),node,it->first,it->second,u,m[0],m[1],m[2],m[3]);

        for (size_t k=0;k<4;k++) {
            assert(i < size());
            assert(int(node+k) < n_nodes+4);
            gsl_matrix_set(A,i,node+k,m[k]);
        }
        gsl_vector_set(B,i,it->second);
    }
    double chisq = 0;
    gsl_matrix *cov = gsl_matrix_alloc(n_nodes+4,n_nodes+4);
    gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc(size(),n_nodes+4);
    gsl_multifit_linear(A,B,P,cov,&chisq,work);
#if 0
    printf("A:\n");
    for (size_t i=0;i<size();i++) {
        for (size_t j=0;j<n_nodes+4;j++) {
            printf("%8.4f ",gsl_matrix_get(A,i,j));
        }
        printf("\n");
    }
    printf("B:\n");
    for (size_t i=0;i<size();i++) {
        printf("%8.4f\n",gsl_vector_get(B,i));
    }
    printf("P:\n");
    for (size_t i=0;i<n_nodes+4;i++) {
        printf("%8.4f\n",gsl_vector_get(P,i));
    }
#endif
    Function output;
    for (double x=fm;x<fM;x+=sampling) {
        int node = floor((x - fm)/dx);
        assert(node < n_nodes);
        double u = (x - (fm+node*dx))/dx;
        double pu[4] = {1, u, u*u, u*u*u};
        // Compute 
        // M=(1/6) * [1 4 1 0; -3 0 3 0 ; 3 -6 3 0; -1 3 -3 1];
        // [1 u u^2 u^3] * M
        double m[4] = {
            (1*pu[0] + -3*pu[1] +  3*pu[2] + -1*pu[3]) / 6,
            (4*pu[0] +  0*pu[1] + -6*pu[2] +  3*pu[3]) / 6,
            (1*pu[0] +  3*pu[1] +  3*pu[2] + -3*pu[3]) / 6,
            (0*pu[0] +  0*pu[1] +  0*pu[2] +  1*pu[3]) / 6
        };
        double y = 0;
        for (size_t k=0;k<4;k++) {
            y += m[k] * gsl_vector_get(P,node+k);
        }
        output.set(x,y);
    }
    gsl_matrix_free(A);
    gsl_vector_free(B);
    gsl_matrix_free(cov);
    gsl_multifit_linear_free(work);
    gsl_vector_free(P);
    return output;
}

AngularFunction AngularFunction::spline_approximation(double node_spacing,double sampling) {
    // First create a function that does not jump around pi
    Function tmp;
    bool first = true;
    double prev = 0;
    for (const_iterator it=begin();it!=end();it++) {
        if (first) {
            tmp.set(it->first,it->second);
            prev = it->second;
            first = false;
        } else {
            double dtheta = remainder(it->second-prev,2*M_PI);
            tmp.set(it->first,prev+dtheta);
            prev = it->second;
        }
    }
    // Now we can do the approximation
    return AngularFunction(tmp.spline_approximation(node_spacing,sampling));
}

#endif
