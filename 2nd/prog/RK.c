/* 
 * Основная функция метода Рунге-Кутта
 */

 #include "my.h"
 
 int runge(double EPS, const char*filename)
{
    struct vect v; // основной вектор координат
	struct vect k[DIM]; // числа Рунгe

	//*инициализация чисел Рунге */
	for (int i=0;i<DIM;i++)
    	init(k + i,RK_NUM);
    
    /* Вектор функций и расчётные коэффициенты */
    double (*f[DIM]) (struct vect v, double t);
	struct vect a, p, pp; 
	struct vect b[RK_NUM];

	/* инициализация вектора функций и коэффициентов */
	Create_Vect_func(f); //
	RK_koeff(b,&a,&p,&pp);


    double h = H_0; // шаг
    double t = 0; // время
	
    double errflag = 0;

    
    FILE *file = fopen(filename, "w");

    if (file == NULL) errflag = -1;
    else{
        fprintf(file, "t \t x \t y\n");

        /* вводим начальные условия */
        initials(&v);
        printf("%d\n",a.length);
        print(v);
        		//runge_step(&h,v,t, f, b, a, k);// вводим начальные значения чисел Рунге
	    
	    //h = new_h(EPS,h,k,v); // обновляем шаг
        
        //new_k(h,k,x,y,u,v); // ежели он изменился -- обновим числа Рунге

        for (t = 0;t<=WORK_TIME+h;t+=h){
            // если мы в нужной точке T, 3T/4, T/2, T/4, то flag == 1 

        	/* обновим координаты */

    
    /*    	

    		fprintf(file, "%e", t);
            fprint(v,file);
            fprintf(file, "\n");
    */
            h = new_step( EPS,  h,  t,  v,  k,  a,  b,  p,  pp, f);
			
			/* обновляем числа Рунге */
    

//   printf("(!)\n" );
            fprintf(file, "%e", t);
            fprint(v,file);
            fprintf(file, "\n");


            runge_step(h,v,t, f, b, a, k);

            new_parametrs(h,t,&v,k,p);

            if (h<H_DEFAULT) {errflag = -3; break;} // на случай особых точек
        }
    }
    printf("final step = %e;\n", h);
    fclose(file);
    return errflag;
}

/* функция обновления координат */
void new_parametrs(double h, double t, struct vect *v, struct vect k[DIM], struct vect p){
	for (int i = 0; i < DIM; i++){
		v->array[i] = v->array[i] + scalar(k[i],p);
	}
}