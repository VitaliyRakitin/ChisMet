/*
 * В данном файле содержатся все изначальные функции и начальные условия
 */

 #define DIM 2 // размерность пространства (количество уравнений)
 #define WORK_TIME 30000 // Время работы программы

/* Начальные условия */
 #define X_0 1
 #define Y_0 0


/* Основные функции */ 
double  f1(struct vect v, double t);
double  f2(struct vect v, double t);

void Create_Vect_func(double (*f[DIM]) (struct vect v, double t));