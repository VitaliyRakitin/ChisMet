#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<unistd.h>
#include"prince.h"
/* Сия программа нацелена решать задачу 34,
 * которая сводится к краевой задаче с переменными
 * phi, x, y, py, px, a,
 * и краевыми условиями
 * phi(0)=0, phi(1)=1, x(0)=0, y(0)=1, y(1)=0, px(1)=0
 *
 * Для большей универсальности размерность векторного поля остаётся плавающей
 * */

double f0(double t,double *x){
  return x[1];
}

double f1(double t, double *x){
  return -x[0];
}
/* Векторная алгебра
 * функция копирования проверена
 */
double *vec_cpy(unsigned int dim,double *dest,double *sourse){
  unsigned int i;
  for (i=0;i<dim;i++) dest[i]=sourse[i];
  return dest;
}

double *vad(unsigned int dim,double *dest,double *source1,double *source2){
  unsigned int i;
  double cs1[dim]; vec_cpy(dim,cs1,source1);
  double cs2[dim]; vec_cpy(dim,cs2,source2);
  for (i=0;i<dim;i++) dest[i] = cs1[i] + cs2[i];
  return dest;
}
/* умножение вектора на скаляр работает*/
double *vsm(unsigned int dim,double *dest,double *source,double lambda){
  unsigned int i;
  double cs[dim]; vec_cpy(dim,cs,source); // а вдруг dest==source
  for (i=0;i<dim;i++) dest[i]=cs[i]*lambda;
  return dest;
}
/* Скалярное произведение годное*/
double vsp(unsigned int dim,double *source1,double *source2){
  unsigned int i; double sum=0;
  double cs1[dim]; vec_cpy(dim,cs1,source1);
  double cs2[dim]; vec_cpy(dim,cs2,source2);
  for (i=0;i<dim; i++) sum+=cs1[i] * cs2[i];
  return sum;
}
/* подсчёт значения вектор-функции векторного аргумента работает*/
double * vfn(unsigned int dim,double *dest, double (*f[])(double,double*),double*x,double t){
  unsigned int i;
  double cs[dim]; vec_cpy(dim,cs,x);// здесь это особенно важно
  for (i=0; i<dim; i++) dest[i] = f[i](t,cs);
  return dest;
}

double runge_step(double *h,unsigned int dim,double (*f[])(double,double*),double*x,double t){
/* Один шаг метода Рунге--Кутты
 * включая корректировку шага
 * в массив h сохраняется информация о шагах
 * h[0] сделанный шаг, h[1] рекомендация следующего
 */
  double k[7][dim]; // числа Рунге
  double zz[dim],zz1[dim];   // хранилище промежуточных расчётов
  double err = 0; //Ошибка на шаге
  double b[7][6] = { {0,0,0,0,0,0},{b21,0,0,0,0,0},{b31,b32,0,0,0,0},{b41,b42,b43,0,0,0},{b51,b52,b53,b54,0,0},{b61,b62,b63,b64,b65,0},{b71,b72,b73,b74,b75,b76} };
  unsigned int i,m,j;   // счётчик для тестов
  vsm(dim,k[0],vfn(dim,zz,f,x,t),h[0] ); // k1 = h * f(x) вроде правильно
 
  for(m=0;m<7;m++){
    vec_cpy(dim,zz,x);
    for(j=0;j<m;j++){
     vad(dim,zz,zz,vsm(dim,zz1,k[j],b[m][j]) );// добавить b_mj k_j  
    }
    vsm(dim,k[m],vfn(dim,zz,f,x,t),h[0]);
  }
 
  for(i=0;i<dim;i++) printf("%f\t",k[5][i]); printf("\n");
  return err;
}


int main(int argc,char**argv){
  double (*f[2]) (double,double*); double x,y; unsigned int n; double xx[2]={2.,3.};
  double yy[2]; double h[2]={0.01,0.01}; double t;
  char op[2][10]={"+","*"};
  if (argc<4) {printf("type number (1,2) of function and two operands\n"); return -1;}
  f[0]=f0;
  f[1]=f1;
  x = atof(argv[2]); y = atof(argv[3]); n= atoi(argv[1]) - 1;
  printf("%f %s %f=%f\n",x,op[n],y, runge_step(h,2,f,xx,t) );
  for (n=0;n<2;n++) printf("%f\t",xx[n]);printf("\n");
  vsm(2,yy,xx,4);
  for (n=0;n<2;n++) printf("%f\t",yy[n]);printf("\n");
  printf("%f\n",vsp(2,xx,yy));
  vfn(2,xx,f,yy,1);
  for (n=0;n<2;n++) printf("%f\t",xx[n]);printf("\n");
  return 0;
}
