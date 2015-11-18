#include "my.h"

/* Векторизация параметров метода Рунге-Кутта */

void RK_koeff(struct vect B[RK_NUM], struct vect *A, struct vect *P, struct  vect *PP)
{	// Матрица B
	double b[RK_NUM][RK_NUM-1] = { 	{0,0,0,0,0,0},
									{b21,0,0,0,0,0},
									{b31,b32,0,0,0,0},
									{b41,b42,b43,0,0,0},
									{b51,b52,b53,b54,0,0},
									{b61,b62,b63,b64,b65,0},
									{b71,b72,b73,b74,b75,b76} };

	
	for (int i = 0; i < RK_NUM; i++){
		init(B+i ,i); 
		for (int j = 0; j< i; j++)
			B[i].array[j] = b[i][j];
	}


	//матрицы a,p,pp
	double p[RK_NUM]  = {p1,p2,p3,p4,p5,p6,p7};
	double pp[RK_NUM] = {pp1,pp2,pp3,pp4,pp5,pp6,pp7};
	double a[RK_NUM]  = {1,a2,a3,a4,a5,a6,a7};

	init(A , RK_NUM);
	init(P , RK_NUM); 
	init(PP, RK_NUM);


	for (int i=0;i<RK_NUM;i++)	A->array[i] = a[i];
	
	for (int i=0;i<RK_NUM;i++)	P->array[i] = p[i];
	
	for (int i=0;i<RK_NUM;i++)	PP->array[i] = pp[i];
	
}