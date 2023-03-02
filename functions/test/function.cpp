#include <stdlib.h>
#include "functions/Function.h"


int main() {
	Function f(Function::INTER_NEAREST);
	unsigned int i;

	for (i=0;i<=1000;i++) {
		f.set((2*M_PI*i)/1000.,sin((2*M_PI*i)/1000.));
	}

    Function f2(f);
    f2.setInterpolationMode(Function::INTER_LINEAR);

    Function g = f.derivative();
    g.setInterpolationMode(Function::INTER_LINEAR);
    Function h = f.primitive(0,0.1)-1;
    f2.print("f2");
    f.print("f");
    g.print("g");
    h.print("h");

	printf("F(0.00)   : %f %f %f %f\n",f(0),f2(0), g(0), h(0));
	printf("F(0.45)   : %f %f %f %f\n",f(0.45),f2(0.45),g(0.45),h(0.45));
	printf("F(1.00): %f %f %f %f\n",f(1.0),f2(1.0),g(1.0),h(1.0));

    Function fsum = f + g;
    Function fmul = f * g;
    Function fsubd = f - 5.0;
    Function::Support z = fsubd.zeros(); 
    Function fsin = f.map(sin); // = sin o f
    Function fcompose = f.map(g); // =  g o f
    Function finv = f.inverse(); // f^-1

	return 0;

}

