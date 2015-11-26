#include"fun.h"
#include"prince.h"
int whose_close(double h,double t,unsigned int ptsc,double *pts){
  unsigned int i;
  for (i=0;i<ptsc;i++) if (t < pts[i] && t+h > pts[i]) return i;
  return -1;
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
  err = distance(dim, x6, x7, norm_evkl);
  // А тут вшита процедура выбора шага
  h[1] = h[0] * s_m(err,eps);
  if (err>eps) { h[0]=h[1]; return runge_step(eps,h,dim,f,x,t,ptsc,pts,s_m);}
  if( (i = 1 + whose_close(h[0],t,ptsc,pts))>0 ){
    h[0] = pts[i-1]-t;
    return runge_step(eps,h,dim,f,x,t,0,NULL,no_grows);
  }
  // Тут уже происходят страшные вещи, меняются входные иксы
  vec_cpy(dim,x,x6);
  return err;
}
