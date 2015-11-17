//
//  task1.h
//  
//
//  Created by Виталий Ракитин on 25.09.15.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define K_NUM 7 // количество чисел Рунге
#define VAL_NUM 4 // количество переменных
#define WORK_TIME 3.14159269 // конечная точка по времени
#define CHECK_POSITIONS 4 // количество точек для проверки локальной погрешности
#define CHECK_POINTS 3 // количество степеней для проверки локальной погрешности

//начальные условия
#define x0 0.0
#define y0 1.0
#define u0 0.0
#define v0 0.0

//начальный шаг и оценка пограешности
//#define EPS 1e-11
#define h0 	1e-3
#define H_DEFAULT 1e-15
#define FAC 	0.9
#define FACMIN	0.5
#define FACMAX 	2.0
#define POW 	5.0

//коэффициенты из метода Рунге-Кутта
#define a2 1.0/5.0
#define a3 3.0/10.0
#define a4 4.0/5.0
#define a5 8.0/9.0
#define a6 1.0
#define a7 1.0

#define b21 1.0/5.0

#define b31 3.0/40.0
#define b32 9.0/40.0

#define b41 44.0/45.0
#define b42 -56.0/15.0
#define b43 32.0/9.0

#define b51 19372.0/6561.0
#define b52 -25360.0/2187.
#define b53 64448.0/6561.0
#define b54 -212.0/729.0

#define b61 9017.0/3168.0
#define b62 -355./33.
#define b63 46732.0/5247.0
#define b64 49.0/176.0
#define b65 -5103.0/18656.0

#define b71 35.0/384.0
#define b72 0.
#define b73 500./1113.
#define b74 125.0/192.0
#define b75 -2187.0/6784.
#define b76 11./84.

#define p1 35.0/384.0
#define p2 0.
#define p3 500./1113.
#define p4 125.0/192.0
#define p5 -2187.0/6784.
#define p6 11./84.
#define p7 0.

#define pp1 5179.0/57600.0
#define pp2 0.
#define pp3 7571.0/16695.0
#define pp4 393.0/640.0
#define pp5 -92097.0/339200.0
#define pp6 187.0/2100.0
#define pp7 1./40.

#include "functions.h" 
