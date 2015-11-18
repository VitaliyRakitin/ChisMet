/* Векторная алгебра.
 * Структура --- вектор.
 *
 * Errors:
 * -1 --- ошибка выделения памяти;
 * -2 --- неопределён вектор
 */

#include "my.h"

//инициализация
int init(struct vect *v, int length){
	if (length == 0){
		v->length =0;
		v->array = NULL;
		return 0;
	}

	v->length = length;
	v->array = (double *)malloc(length*sizeof(double));
	if (v->array == NULL) return -1;
	for (int i=0;i< length; i++)
		v->array[i] = 0;
	return 0;
}

/* Копируем vect1 в vect2.*/
int copy(struct vect vect1, struct vect *vect2){
	if (vect1.length ==0){
		if (vect2->array != NULL)
			free(vect2->array);
		vect2->array = NULL;
		vect2->length = 0;
		return 0;
	}
	if ((vect1.length < 0)||(vect1.array == NULL)) return -2;

	free(vect2->array);

	double *new_array;
	new_array = (double *)malloc(vect1.length*sizeof(double));
	if (new_array == NULL) return -1;

	for (int i = 0; i< vect1.length; i++)
		new_array[i] = vect1.array[i];
	
	vect2->array = new_array;
	vect2->length = vect1.length;
	return 0;
}

/* Инициализатор копирования (создаётся a как копия b) */
int init_copy(struct vect *a, struct vect b){
	int errflag;
	errflag = init(a,b.length);
	if (errflag == 0) errflag = copy(b,a);
	return errflag;
}

/* Скалярное произведение.
 * Замечание: если длинна одного больше, то второй подгоняется нулями.
 */
double scalar(struct vect vect1, struct vect vect2){
	if ((vect1.length <= 0)||(vect2.length <= 0)) {
		return 0;
	}

	double result = 0;

	if (vect1.length >= vect2.length)
		for (int i= 0; i< vect1.length;i++)
			result += vect1.array[i] * vect2.array[i];
	else
		for (int i= 0; i< vect1.length;i++)
			result += vect1.array[i] * vect2.array[i];

	return result;
}

/* подсчёт значения вектор-функции векторного аргумента 
 * размерность заведомо DIM, так как большего в нашей программе не надо
 */
struct vect vf(struct vect v,  double t, double (*f[])(struct vect,double)){
  struct vect new_vect;
  init_copy(&new_vect,v);// здесь это особенно важно

  for (int i=0; i<DIM; i++) new_vect.array[i] = f[i](v,t);
  return new_vect;
}

/* Сумма */
struct vect sum(struct vect a, struct vect b){
	struct vect result;
	init(&result, a.length);
	for (int i=0;i < a.length;i++)
		result.array[i] = a.array[i] + b.array[i];
	return result;
}

/* Разность */
struct vect minus(struct vect a, struct vect b){
	struct vect result;
	init(&result, a.length);
	for (int i=0;i < a.length;i++)
		result.array[i] = a.array[i] - b.array[i];
	return result;
}

/* Евклидова норма */
double norm(struct vect a){
	double sum = 0;
	for (int i=0; i < a.length; i++)
		sum += a.array[i] * a.array[i];
	return sum;
}


double norm_max(struct vect a){
	struct vect v;
	init_copy(&v,a);
	double max = 0;
	for (int i=0; i < v.length; i++){
		if (v.array[i] < 0) v.array[i] = v.array[i] * (-1);
		if (max < v.array[i]) max = v.array[i];
	}
	return max;
}
double norm_sum(struct vect a){
	struct vect v;
	init_copy(&v,a);
	double sum = 0;
	for (int i=0; i < v.length; i++){
		if (v.array[i] < 0) v.array[i] = v.array[i] * (-1);
		sum += v.array[i] ;
	}
	return sum;
}



/* стандартный вывод */
void print(struct vect v){
	if (v.length<=0) return;
	for (int i=0;i<v.length;i++)
		printf("%e ", v.array[i]);
	printf("\n");
}

/* вывод в открытй файл */
void fprint(struct vect v, FILE *f){
	if (v.length<=0) return;
	for (int i=0;i<v.length;i++)
		fprintf(f,"\t%e", v.array[i]);
}
