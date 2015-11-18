#include "my.h"

/* Подсчёт коеффициентов Рунге.
 */

void runge_step(double h,struct vect v, double t,double (*f[])(struct  vect,double), struct vect b[RK_NUM],struct vect a, struct vect K[DIM]){

  struct vect new_v;
  init(&new_v,DIM); 

  for (int i = 0; i < RK_NUM; i++){
    /* задаём new_v */
    for (int k = 0; k < DIM; k++)
       new_v.array[k] = scalar(K[k],b[i]);    
    
    for (int j = 0; j < DIM; j++)
      K[j].array[i] = h * f[j](sum(v,new_v),t + h * a.array[i] );
  }
}