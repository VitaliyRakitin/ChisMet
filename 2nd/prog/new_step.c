#include "my.h"


double new_step(double EPS, double h, double t, struct vect v, struct vect k[DIM], struct vect a, struct vect b[DIM], 
                struct vect p, struct vect pp, double (*f[])(struct  vect,double)){


    double var = Variation(h, t, v, k, p, pp);
    double coef = 1.0/(POW + 1.0);
    
    while (var > EPS){

        struct vect k_new[DIM];
        for (int i=0;i<DIM;i++) init_copy(k_new + i, k[i]);

        h = h*min(max(FAC * pow(EPS/var,coef)));
 
        runge_step(h,v,t, f, b, a, k_new);
        var = Variation(h, t, v, k_new, p, pp);
        
        if (h<H_DEFAULT) break;
    }

    h = h*min(max(FAC * pow(EPS/var,coef)));
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
double Variation(double h, double t, struct vect v, struct vect k[DIM], struct vect p, struct vect pp){
    struct vect v_1, v_2;
    init_copy(&v_2,v);
    init_copy(&v_1,v);

    new_parametrs(h,t,&v_1,k,p);
    new_parametrs(h,t,&v_2,k,pp);

    v_1 = minus(v_1,v_2);

    return norm(v_1); // растояние от точки до точки

}