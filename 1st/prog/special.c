//
//  special.c
//  
//
//  Created by Виталий Ракитин on 26.09.15.
//
//

#include "task1.h"

// Основная функция. На вход: EPS --- порядок погрешности
// filename --- файл для вывода, local --- значения локальной погрешности в точках
int runge(double *global_error, double EPS, const char*filename, double local[CHECK_POSITIONS * VAL_NUM ])
{
    double x,y,u,v; // x y a b
    double h=0;

    double k[K_NUM][VAL_NUM]; // числа рунге
    int errflag = 0;
    int check_time = CHECK_POSITIONS; // количество оставшихся точек на временной шкале для проверки лок погрешности
    int flag = 0;
    double h_next = 0; // механическая переменная для хранения следующего шага

    *global_error = 0;

    FILE *file = fopen(filename, "w");

    if (file == NULL) errflag = -1;
    else{
        fprintf(file, "t,  x,  y\n");

        //вводим начальные условия
        x = x0;    y = y0;    u = u0;   v = v0;   h = EPS/h0;
        new_k(h,k,x,y,u,v); // начальные числа Рунге
        h = new_h(EPS,h,k,x,y,u,v); // обновляем шаг
        new_k(h,k,x,y,u,v); // ежели он изменился -- обновим числа Рунге

        for (double t = 0;t<=WORK_TIME+h;t+=h){
            // если мы в нужной точке T, 3T/4, T/2, T/4, то flag == 1 

            if (flag == 0)
                if (check_time > 0) // если меньше, то точек не осталось, увы
                        if ((nearest(t) - t <= h)&&(nearest(t) - t >= 0)){
                            // если мы сюда попали, то расстояние до нужной точки меньше h
                            // объявим следующий шаг этим расстоянием
                            h_next = nearest(t) - t;
                        }


            //подсчет глобальной погрешности
            (*global_error) = Variation(h,k,x,y,u,v) + (*global_error) * exp_int(h,x,y);

        // обновляем координаты
            x = new_parametr(h,x,k,x,y,u,v,0);
            y = new_parametr(h,y,k,x,y,u,v,1);
            u = new_parametr(h,u,k,x,y,u,v,2);
            v = new_parametr(h,v,k,x,y,u,v,3);
            fprintf(file, "%e,  %e,  %e\n", t,x,y);
            // если мы в нужной точке T, 3T/4, T/2, T/4, то flag == 1 
            if (flag == 1) {
                local[CHECK_POSITIONS*(CHECK_POSITIONS - check_time)] = x;
                local[CHECK_POSITIONS*(CHECK_POSITIONS - check_time) + 1] = y;
                local[CHECK_POSITIONS*(CHECK_POSITIONS - check_time) + 2] = u;
                local[CHECK_POSITIONS*(CHECK_POSITIONS - check_time) + 3] = v;
                h = h_next; // возвращаем исходное значение шага 
                h_next = 0; 
                flag = 0;
                check_time--; // мы же уже дошли до этой точке, осталось на одну меньше
            }

            h = new_h(EPS,h,k,x,y,u,v); //обновляем шаг

            if (h_next != 0){
                double p;
                p = h;
                h = h_next;
                h_next = p;//сохраняем предыдуший шаг на случай, если новый слишком маленький
                flag = 1; // для того, чтобы мы видели вывод
            }
            
            new_k(h,k,x,y,u,v); //обновляем числа Рунге
           
            if (h<H_DEFAULT) {errflag = -2; break;} // на случай особых точек
        }
    }
    printf("final step = %e;\n", h);
    fclose(file);
    return errflag;
}

double nearest(double t){
    double near = 0;
    double segment = WORK_TIME/4.; 
    if ((t >= 0)&&(t <= segment)) near = segment;
    if ((t >     segment)&&(t <= 2 * segment)) near = 2 * segment;
    if ((t > 2 * segment)&&(t <= 3 * segment)) near = 3 * segment;
    if ((t > 3 * segment)&&(t <= 4 * segment)) near = 4 * segment;
    return near; 
}

double exp_int(double h, double x, double y){
    return 1;
    //return p = h * exp(pow(x,2) + pow(y,2) + sqrt(5 - 8 * pow(x,2) + 4 * pow(x,4) - 8 * pow(y,2) + 8 * pow(x,2) * pow(y,2) + 4 * pow(y,4)));
}

// обновление шага
double new_h(double EPS, double h, double k[K_NUM][VAL_NUM], double x, double y, double u, double v){
    double var = Variation(h,k,x,y,u,v);
    double p = 1.0/(POW + 1.0);

    while (var > EPS){
        h = h*min(max(FAC * pow(EPS/var,p))); 
        new_k(h,k,x,y,u,v); 
        var = Variation(h,k,x,y,u,v); 
        if (h<H_DEFAULT) break;
    }
    h =  h*min(max(FAC * pow(EPS/var,p)));
    return h;
}

double max(double variation){
    if (FACMIN >  variation) return FACMIN;
    return variation;
}

double min(double variation){
    if (FACMAX >  variation) return variation;
    return FACMAX;
}

double max2(double a, double b){
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    if (a > b) return a;
    return b;
}

// погрешность
double Variation(double h, double k[K_NUM][VAL_NUM], double x, double y, double u, double v){ 
    double x1 = new_parametr(h,x,k,x,y,u,v,0);
    double y1 = new_parametr(h,y,k,x,y,u,v,1);
    double u1 = new_parametr(h,u,k,x,y,u,v,2);
    double v1 = new_parametr(h,v,k,x,y,u,v,3);
    
    double xx = new_parametr1(h,x,k,x,y,u,v,0);
    double yy = new_parametr1(h,y,k,x,y,u,v,1);
    double uu = new_parametr1(h,u,k,x,y,u,v,2);
    double vv = new_parametr1(h,v,k,x,y,u,v,3);

    //double a = x1 - xx;
    //double b = y1 - yy;
    //double c = u1 - uu;
    //double d = v1 - vv;

   // return max2(max2(a,b),max2(c,d));
    return sqrt((x1 - xx)*(x1 - xx) + (y1 - yy)*(y1 - yy) + (u1 - uu)*(u1 - uu) + (v1 - vv)*(v1 - vv)); // растояние от точки до точки

}

double new_parametr(double h, double old, double k[K_NUM][VAL_NUM], double x, double y, double u, double v, int pos){
  return old + p1 * k[0][pos] + p2 * k[1][pos] + p3 * k[2][pos] + p4 * k[3][pos] + p5 * k[4][pos] + p6 * k[5][pos];
}

double new_parametr1(double h, double old, double k[K_NUM][VAL_NUM], double x, double y, double u, double v, int pos){
   return old + pp1 * k[0][pos] + pp2 * k[1][pos] + pp3 * k[2][pos] + pp4 * k[3][pos] + pp5 * k[4][pos] + pp6 * k[5][pos] + pp7 * k[6][pos];
}