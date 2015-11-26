#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<unistd.h>
double f0(double t,double *x);
double f1(double t, double *x);
double dphi(double t,double *x);
double   dx(double t,double *x);
double   dy(double t,double *x);
double  dpy(double t,double *x);
double  dpx(double t,double *x);
double   da(double t,double *x);
int whose_close(double h,double t,unsigned int ptsc,double *pts);
double *vec_cpy(unsigned int dim,double *dest,double *sourse);
double *vector_add(unsigned int dim,double *dest,double *source1,double *source2);
double *vector_scalar_mult(unsigned int dim,double *dest,double *source,double lambda);
double vector_scalar_prod(unsigned int dim,double *source1,double *source2);
double * vector_function(unsigned int dim,double *dest, double (*f[])(double,double*),double*x,double t);
double norm_max(unsigned int dim,double *source);
double norm_evkl(unsigned int dim,double *source);
double distance(unsigned int dim,double *x1,double *x2,double (*norm)(unsigned int,double*));
double step_mult(double err,double eps);
double no_grows(double err,double eps);
double runge_step(double eps,double *h,unsigned int dim,double (*f[])(double,double*),double*x,double t,int ptsc,double *pts,double(*s_m)(double,double));
int osc_while(double *x,double t);
double harmonic_oscillator(double eps,double x0,double dx0,int(*condition)(double*,double));
double traj34(double eps,double *x);
double *X(double eps,double *u,double* beta);
double *test_X(double eps,double *u,double *beta);
double *test_X2(double eps,double *u,double *beta);
double *dF(double eps,unsigned int dim,unsigned int pdim,double *df,double*(*F)(double,double*,double*),double *pt);
int shot2(double eps,double *a,double*(*F) (double,double*,double*));
