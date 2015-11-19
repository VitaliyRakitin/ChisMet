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

#define alpha 0.1

double f0(double t,double *x){
  t= t;
  return x[1];
}

double f1(double t, double *x){
  t=t;
  return -x[0];
}

// Векторное поле для задачи 34
double dphi(double t,double *x){  t=t;  	  return  x[1];}
double   dx(double t,double *x){  t=t;  	  return  x[2];}
double   dy(double t,double *x){  t=t;  	  return  x[3]*(1 + alpha * t*t*x[1]*x[1]);}
double  dpy(double t,double *x){  t=t; 		  return -x[4];}
double  dpx(double t,double *x){  t=t;  	  return  x[5];}
double   da(double t,double *x){  t=t; x[0]=x[0]; return  0.;}

/*
 * Близко ли чисто t+h[0] к одному из pts[] (справа) и к какому
 * возвращает -1, если не к какому, номер точки иначе
 * */
int whose_close(double h,double t,unsigned int ptsc,double *pts){
  unsigned int i;
  for (i=0;i<ptsc;i++) if (t < pts[i] && t+h > pts[i]) return i;
  return -1;
}


/* Векторная алгебра
 * функция копирования проверена
 */
double *vec_cpy(unsigned int dim,double *dest,double *sourse){
  unsigned int i;
  for (i=0;i<dim;i++) dest[i]=sourse[i];
  return dest;
}

double *vector_add(unsigned int dim,double *dest,double *source1,double *source2){
  unsigned int i;
  double cs1[dim]; vec_cpy(dim,cs1,source1);
  double cs2[dim]; vec_cpy(dim,cs2,source2);
  for (i=0;i<dim;i++) dest[i] = cs1[i] + cs2[i];
  return dest;
}
/* умножение вектора на скаляр работает*/
double *vector_scalar_mult(unsigned int dim,double *dest,double *source,double lambda){
  unsigned int i;
  double cs[dim]; vec_cpy(dim,cs,source); // а вдруг dest==source
  for (i=0;i<dim;i++) dest[i]=cs[i]*lambda;
  return dest;
}
/* Скалярное произведение годное*/
double vector_scalar_prod(unsigned int dim,double *source1,double *source2){
  unsigned int i; double sum=0;
  double cs1[dim]; vec_cpy(dim,cs1,source1);
  double cs2[dim]; vec_cpy(dim,cs2,source2);
  for (i=0;i<dim; i++) sum+=cs1[i] * cs2[i];
  return sum;
}
/* подсчёт значения вектор-функции векторного аргумента работает*/
double * vector_function(unsigned int dim,double *dest, double (*f[])(double,double*),double*x,double t){
  unsigned int i;
  double cs[dim]; vec_cpy(dim,cs,x);// здесь это особенно важно
  for (i=0; i<dim; i++) dest[i] = f[i](t,cs);
  return dest;
}
/* Норма максимум модуля вроде работает*/
double norm_max(unsigned int dim,double *source){
   unsigned int i;
   double norm=fabs(source[0]);
   for(i=0;i<dim;i++) if (norm < fabs(source[i])) norm=fabs(source[i]);
   return norm;
}
/* Расстояние, задаваемое через норму, годно */
double distance(unsigned int dim,double *x1,double *x2,double (*norm)(unsigned int,double*)){
  double z[dim];
  vector_scalar_mult(dim,z,x2,-1);
  vector_add(dim,z,z,x1);
  return norm(dim,z);
}

double step_mult(double err,double eps){
  return fmin(2,fmax(0.3,0.9*pow(eps/err,1./7.) ) );
}
double no_grows(double err,double eps){
  return fmin(1,fmax(0.3,0.9*pow(eps/err,1./7.) ) );
}

double runge_step(double eps,double *h,unsigned int dim,double (*f[])(double,double*),double*x,double t,int ptsc,double *pts,double(*s_m)(double,double)){
/* Один шаг метода Рунге--Кутты
 * включая корректировку шага
 * в массив h сохраняется информация о шагах
 * h[0] сделанный шаг, h[1] рекомендация следующего
 * Сюда намертво вшиты:
 *   Метод Дормана—Принса 6-го порядка, 
 */
  double k[7][dim]; // числа Рунге
  double zz[dim],zz1[dim];   // хранилище промежуточных расчётов
  double x6[dim],x7[dim]; // точки, посчитанные 6-м и 7-м методом
  double err = 0; //Ошибка на шаге
  double a[7] = {1.,a2,a3,a5,a6,a7};
  double b[7][6] = { {0,0,0,0,0,0},{b21,0,0,0,0,0},{b31,b32,0,0,0,0},{b41,b42,b43,0,0,0},{b51,b52,b53,b54,0,0},{b61,b62,b63,b64,b65,0},{b71,b72,b73,b74,b75,b76} };
  double p[2][7] = { {p1,p2,p3,p4,p5,p6,0},{pp1,pp2,pp3,pp4,pp5,pp6,pp7} };
  unsigned int i,m,j;   // счётчик для тестов
  for(m=0;m<7;m++){  // подсчёт чисел Рунге K
    vec_cpy(dim,zz,x);
    for(j=0;j<m;j++){
      vector_add(dim,zz ,zz ,vector_scalar_mult(dim,zz1,k[j],b[m][j]*h[0]) );// добавить b_mj k_j 
    }
    vector_function(dim,k[m] ,f,zz ,t+h[0]   *a[m]); // km = f0(zz)
  }
  // подсчёт следующей точки траектории
  // похуже
  for(i=0,vec_cpy(dim,x6,x);i<6;i++){
    vector_add( dim , x6 , x6 , vector_scalar_mult(dim, zz, k[i] , p[0][i]*h[0]) );
  }
  // получше (7-й порядок)
  for(i=0,vec_cpy(dim,x7,x);i<7;i++) 
    vector_add( dim , x7 , x7 , vector_scalar_mult(dim, zz, k[i], p[1][i] *h[0]) );
  // Вот здесь вшита норма максимум модуля
  err = distance(dim, x6, x7, norm_max);
  // А тут вшита процедура выбора шага
  h[1] = h[0] * s_m(err,eps);
  if (err>eps) { h[0]=h[1]; return runge_step(eps,h,dim,f,x,t,ptsc,pts,s_m);}
  if( (i = 1 + whose_close(h[0],t,ptsc,pts))>0 ){
    printf("Точка\n");
    h[0] = pts[i-1]-t;
    return runge_step(eps,h,dim,f,x,t,0,NULL,no_grows);
  }
  // Тут уже происходят страшные вещи, меняются входные иксы
  vec_cpy(dim,x,x6);
  return err;
}

int osc_while(double *x,double t){
  return  (t<40);
}

double harmonic_oscillator(double eps,double x0,double dx0,int(*condition)(double*,double)){
  double (*f[2]) (double,double*) = {f0,f1};
  double h[2]={0.0001,0.0001}, t=0 , x[2]={x0,dx0}, err=0, pts[1]={M_PI_2};
  for (t=0; condition( x,t); t+=h[0] ){
    h[0] = h[1];
    err+=runge_step(eps,h,2,f,x,t,1,pts,step_mult);
  }
  printf("\tt\t\tx\t\tdx\tglobal_error\t|%.0fcos t %s %.0fsin t - x|\n%e\t%e\t%e\t%e\t%e\n",x0,(dx0>0?"+":""),dx0,t,x[0],x[1],err,fabs(x0*cos(t) + dx0*sin(t)-x[0]));
  return err;
}

int main(int argc,char**argv){
  /*double x[6] = {0,0,1,23,23,23};
  int i;*/
  if (argc<5) {printf("%s eps x0 y0 T\n",argv[0]); return -1;}
  double e1 = harmonic_oscillator(atof(argv[1]),atof(argv[2]),atof(argv[3]),osc_while);
  double e2 = harmonic_oscillator(atof(argv[1])/100.,atof(argv[2]),atof(argv[3]),osc_while);
  double e3 = harmonic_oscillator(atof(argv[1])/10000.,atof(argv[2]),atof(argv[3]),osc_while);
  printf("%e\n",fabs((e1-e2)/(e2-e3)));
  
  return 0;
}
