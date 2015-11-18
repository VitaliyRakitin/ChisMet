#include "my.h"
/*
 * Математический осцилятор
 */

/* Исходные функции системы */ 
double  f1(struct vect v, double t){
	//t=t;
	return v.array[1];
}
double  f2(struct vect v, double t){
	//t=t;
	return -v.array[0];
}

/* Создание вектор-функции из имеющихся начальных функций */
void Create_Vect_func(double (*f[DIM]) (struct vect v, double t)){
f[0] = f1;
f[1] = f2;
}


/* Начальные условия */
int initials(struct vect *vect){
	init(vect,DIM);
	if ((vect->length <= 0)||(vect->array == NULL)) return -2;

	vect->array[0] = X_0;
	vect->array[1] = Y_0;
	
	return 0;
}
