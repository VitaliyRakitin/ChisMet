//
//  functions.h
//  
//
//  Created by Виталий Ракитин on 26.09.15.
//
//


int runge(double *global_error, double EPS, const char*, double [CHECK_POSITIONS * VAL_NUM ]);

double nearest(double t);
double exp_int(double, double, double );

double  f1(double x, double y, double u, double v);
double  f2(double x, double y, double u, double v);
double  f3(double x, double y, double u, double v);
double  f4(double x, double y, double u, double v);

void new_k(double h, double [K_NUM][VAL_NUM], double x, double y, double u, double  v);

double max(double variation);
double min(double variation);
double max2(double, double);

double new_h(double EPS, double h, double [K_NUM][VAL_NUM], double x, double y, double u, double v);

double Variation(double h, double [K_NUM][VAL_NUM], double x, double y, double u, double v);

double new_parametr(double h, double old, double [K_NUM][VAL_NUM], double x, double y, double u, double v, int );

double new_parametr1(double h, double old, double [K_NUM][VAL_NUM], double x, double y, double u, double v,int );


