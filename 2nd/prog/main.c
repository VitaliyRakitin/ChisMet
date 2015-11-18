#include "my.h"


void ErrCheck(int errflag);

void func1(struct vect a);

int main(void){

	struct vect a,b;
	init(&a,4);
	for (int i =0;i<4;i++)
		a.array[i] = i;
	init_copy(&b,a);
	print(a);
/*	
init(b,4);
	
		b.array[i] = 2;
	}
	print(minus(a,b));
	printf("%e \n", norm(a));
*/


 	runge( 1e-7, "out.csv");
	return 0;
}

void ErrCheck(int errflag){
	switch(errflag){
	case -1: {printf("Ошибка выделения памяти.\n"); break;}
	case -2: {printf("Неинициализированный вектор.\n"); break;}
	default: {printf("Другая ошибка.\n");break;}
	}
}

void func1(struct vect a){
	for (int i = 0;i<a.length;i++)
		a.array[i] = a.array[i] * 2;
}