#include"fun.h"
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

double norm_evkl(unsigned int dim,double *source){
  unsigned int i;
  double norm = 0;
  for(i=0;i<dim;i++) norm+=source[i]*source[i];
  return sqrt(norm);
}
/* Расстояние, задаваемое через норму, годно */
double distance(unsigned int dim,double *x1,double *x2,double (*norm)(unsigned int,double*)){
  double z[dim];
  vector_scalar_mult(dim,z,x2,-1);
  vector_add(dim,z,z,x1);
  return norm(dim,z);
}
