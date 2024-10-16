#include <stdio.h>
#include <math.h>
#include "nr.h"

#define CONV_DO(i) ((double)(i))

void get_eps(float *eps){
	*eps = 0.f;
    float tmp = 1.0f;
    while ((1.0f + tmp) != 1.0f) {
        tmp /= 2.0f;
    }
    *eps = tmp * 2.0f; //마지막에서 하나 앞의 n을 반환
}
void get_eps_do(double *epsd){
    *epsd = 0.0;
    double tmp = 1.0;
    while ((1.0 + tmp) != 1.0) {
        tmp /= 2.0;
    }
    *epsd = tmp * 2.0; //마지막에서 하나 앞의 n을 반환
}
void machar_do(int *ibeta, int *it, int *irnd, int *ngrd, int *machep, int *negep, int *iexp, int *minexp, int *maxexp, double *eps, double *epsneg, double *xmin, double *xmax){
	int i,itemp,iz,j,k,mx,nxres;
	double a,b,beta,betah,betain,one,t,temp,temp1,tempa,two,y,z,zero;

	one=CONV_DO(1);
	two=one+one;
	zero=one-one;
	a=one;
	do {
		a += a;
		temp=a+one;
		temp1=temp-a;
	} while (temp1-one == zero);
	b=one;
	do {
		b += b;
		temp=a+b;
		itemp=(int)(temp-a);
	} while (itemp == 0);
	*ibeta=itemp;
	beta=CONV_DO(*ibeta);
	*it=0;
	b=one;
	do {
		++(*it);
		b *= beta;
		temp=b+one;
		temp1=temp-b;
	} while (temp1-one == zero);
	*irnd=0;
	betah=beta/two;
	temp=a+betah;
	if (temp-a != zero) *irnd=1;
	tempa=a+beta;
	temp=tempa+betah;
	if (*irnd == 0 && temp-tempa != zero) *irnd=2;
	*negep=(*it)+3;
	betain=one/beta;
	a=one;
	for (i=1;i<=(*negep);i++) a *= betain;
	b=a;
	for (;;) {
		temp=one-a;
		if (temp-one != zero) break;
		a *= beta;
		--(*negep);
	}
	*negep = -(*negep);
	*epsneg=a;
	*machep = -(*it)-3;
	a=b;
	for (;;) {
		temp=one+a;
		if (temp-one != zero) break;
		a *= beta;
		++(*machep);
	}
	*eps=a;
	*ngrd=0;
	temp=one+(*eps);
	if (*irnd == 0 && temp*one-one != zero) *ngrd=1;
	i=0;
	k=1;
	z=betain;
	t=one+(*eps);
	nxres=0;
	for (;;) {
		y=z;
		z=y*y;
		a=z*one;
		temp=z*t;
		if (a+a == zero || fabs(z) >= y) break;
		temp1=temp*betain;
		if (temp1*beta == z) break;
		++i;
		k += k;
	}
	if (*ibeta != 10) {
		*iexp=i+1;
		mx=k+k;
	} else {
		*iexp=2;
		iz=(*ibeta);
		while (k >= iz) {
			iz *= *ibeta;
			++(*iexp);
		}
		mx=iz+iz-1;
	}
	for (;;) {
		*xmin=y;
		y *= betain;
		a=y*one;
		temp=y*t;
		if (a+a != zero && fabs(y) < *xmin) {
			++k;
			temp1=temp*betain;
			if (temp1*beta == y && temp != y) {
				nxres=3;
				*xmin=y;
				break;
			}
		}
		else break;
	}
	*minexp = -k;
	if (mx <= k+k-3 && *ibeta != 10) {
		mx += mx;
		++(*iexp);
	}
	*maxexp=mx+(*minexp);
	*irnd += nxres;
	if (*irnd >= 2) *maxexp -= 2;
	i=(*maxexp)+(*minexp);
	if (*ibeta == 2 && !i) --(*maxexp);
	if (i > 20) --(*maxexp);
	if (a != y) *maxexp -= 2;
	*xmax=one-(*epsneg);
	if ((*xmax)*one != *xmax) *xmax=one-beta*(*epsneg);
	*xmax /= (*xmin*beta*beta*beta);
	i=(*maxexp)+(*minexp)+3;
	for (j=1;j<=i;j++) {
		if (*ibeta == 2) *xmax += *xmax;
		else *xmax *= beta;
	}
}
#undef CONV_DO



int main(){
	int ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp;
	float eps, epsneg, xmin, xmax;
    double epsd, epsnegd, xmind, xmaxd;

	machar(&ibeta, &it, &irnd, &ngrd, &machep, &negep, &iexp, &minexp, &maxexp,
			&eps, &epsneg, &xmin, &xmax);
    printf("***machar***\n");
	printf("Machine Accuracy (float): \t%0.20f\n", eps);
    machar_do(&ibeta, &it, &irnd, &ngrd, &machep, &negep, &iexp, &minexp, &maxexp,
			&epsd, &epsnegd, &xmind, &xmaxd);
    printf("Machine Accuracy (double): \t%0.20f\n", epsd);

    printf("***original get_eps***\n");
	get_eps(&eps);
	printf("Machine Accuracy (float): \t%0.20f\n", eps);
    get_eps_do(&epsd);
    printf("Machine Accuracy (double): \t%0.20f\n", epsd);

	return 0;
}