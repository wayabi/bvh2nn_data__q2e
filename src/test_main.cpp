#include "qua.h"
#include "rot2.h"
#include <stdio.h>
#include <stdlib.h>

using namespace boost;
using namespace boost::math;
typedef quaternion<double> Q;

/*
int main(int argc, char** argv){
	double E = 0.000001;
	srand(time(NULL));

	double x = M_PI*((rand()%1000-500)/1000.0);
	double y = M_PI*((rand()%1000-500)/1000.0);
	double z = M_PI*((rand()%1000-500)/1000.0);
	printf("x = %f\ny = %f\nz = %f\n", x, y, z);
	for(int i=0;i<12;++i){
		qua::RotSeq rs = (qua::RotSeq)i;
		Q q = qua::e2q(x, y, z, rs);
		double xx = 0;
		double yy = 0;
		double zz = 0;
		qua::q2e(q, xx, yy, zz, rs);
		double xx2 = 0;
		double yy2 = 0;
		double zz2 = 0;
		Q q2 = conj(q);
		qua::q2e(q2, xx2, yy2, zz2, rs);
		if((x-xx)*(x-xx)+(y-yy)*(y-yy)+(z-zz)*(z-zz) > E &&
			(x-xx2)*(x-xx2)+(y-yy2)*(y-yy2)+(z-zz2)*(z-zz2) > E){
			printf("seq[%d] is invalid(x:%f, y:%f, z:%f)\n", i, xx, yy, zz);
			printf("\t\t(x:%f, y:%f, z:%f)\n", xx2, yy2, zz2);
		}else{
			printf("seq[%d] OK\n", i);
		}
	}
	return 0;
}
*/
