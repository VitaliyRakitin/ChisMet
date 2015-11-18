#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#include "vect.h"
#include "RK_koeff.h"

/* Параметры метода РК */
#define RK_NUM 7
#define H_DEFAULT 1e-15
#define FAC 	0.9
#define FACMIN	0.5
#define FACMAX 	1.3
#define POW 	5.0

#define H_0 0.001

#include "vect.h"
#include "change_me.h"


double new_step(double EPS, double h, double t, struct vect v, struct vect k[DIM], struct vect a, struct vect b[DIM], struct vect p, struct vect pp,double (*f[])(struct  vect,double));
double max(double variation);
double min(double variation);
double max2(double a, double b);
double Variation(double h, double t, struct vect v, struct vect k[DIM], struct vect p, struct vect pp);



int runge(double EPS, const char*filename);
void new_parametrs(double h, double t, struct vect *v, struct vect k[DIM], struct vect p);

int initials(struct vect *v);




void RK_koeff(struct vect B[RK_NUM], struct vect *A, struct vect *P, struct  vect *PP);


void runge_step(double h,struct vect v, double t,double (*f[])(struct  vect,double), struct vect b[RK_NUM],struct vect a, struct vect K[DIM]);

