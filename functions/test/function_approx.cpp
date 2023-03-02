#include <stdlib.h>
#include <time.h>
#include "functions/Function.h"


int main() {
    long int seed = time(NULL);
    // seed = 1533818610;
    printf("Seed: %ld\n",seed);
    srand48(seed);
	Function f;
	unsigned int i;

	for (i=0;i<50;i++) {
        double x = drand48()*10;
        double y = sin(x) + drand48()/3;
		f.set(x,y);
	}
    Function f2 = f;

    Function g = f.spline_approximation(1.0,0.01);
    Function g2 = f2.spline_approximation(1.0,0.01);

    f.print("data");
    g.print("approx");
    f2.print("data2");
    g2.print("approx2");

	AngularFunction af;

	for (i=0;i<50;i++) {
        double x = drand48()*10;
        double y = 2*M_PI*sin(x) + drand48()/2;
		af.set(x,y);
	}

    // af = af.select(af.support());

    AngularFunction ag = af.spline_approximation(1.0,0.01);

    af.print("adata");
    ag.print("aapprox");

	return 0;

}

