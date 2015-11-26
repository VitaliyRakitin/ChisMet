#include"fun.h"
/* Сия программа нацелена решать задачу 34,
 * которая сводится к краевой задаче с переменными
 * phi, x, y, py, px, a,
 * и краевыми условиями
 * phi(0)=0, phi(1)=1, x(0)=0, y(0)=1, y(1)=0, px(1)=0
 *
 * Для большей универсальности размерность векторного поля остаётся плавающей
 * */
#define alpha 0.
double f0(double t,double *x){ t= t;  return x[1];}
double f1(double t,double *x){  t=t;  return -x[0];}
// Векторное поле для задачи 34
double dphi(double t,double *x){  t=t;  	  return  x[1];}
double   dx(double t,double *x){  t=t;  	  return  x[2];}
double   dy(double t,double *x){  t=t;  	  return  x[3]*(1 + alpha * t*t*x[1]*x[1]);}
double  dpy(double t,double *x){  t=t; 		  return -x[4];}
double  dpx(double t,double *x){  t=t;  	  return  x[5];}
double   da(double t,double *x){  t=t; x[0]=x[0]; return  0.;}
double dB(double t,double *x){ return x[3]*x[3]*(1+alpha*t*t*x[1]*x[1]);}
/*
 * Близко ли чисто t+h[0] к одному из pts[] (справа) и к какому
 * возвращает -1, если не к какому, номер точки иначе
 * */

int osc_while(double *x,double t){
  return  (t<40) && x[0]<100000.;
}

double harmonic_oscillator(double eps,double x0,double dx0,int(*condition)(double*,double)){
  double (*f[2]) (double,double*) = {f0,f1};
  double h[2]={0.0001,0.0001}, t=0 , x[2]={x0,dx0}, err=0, pts[1]={40.};
  for (t=0; condition( x,t); t+=h[0] ){
    h[0] = h[1];
    err+=runge_step(eps,h,2,f,x,t,1,pts,step_mult);
  }
  printf("\t\t\tПогрешность на шаге %e\n",eps);
  printf("\tt\t\tx\t\tdx\tglobal_error\t|%.0fcos t %s %.0fsin t - x|\n%e\t%e\t%e\t%e\t%e\n",x0,(dx0>0?"+":""),dx0,t,x[0],x[1],err,fabs(x0*cos(t) + dx0*sin(t)-x[0]));
  return err;
}
/*
 * Функция для пристрелки для краевой задачи из 34. Кладёт в x конечные значения.
 * Здесь вшита alpha из define
 */
double traj34(double eps,double *x){
  double err=0,delta=0,mu,t,h[2]={0.0001,0.0001},pts[1]={1.},k,a,b;
  double (*f[6])(double,double*) = {dphi,dx,dy,dpy,dpx,da};
  for(t=0;t<1;t+=h[0]){
    h[0] = h[1];
    err=runge_step(eps,h,6,f,x,t,1,pts,step_mult);// это неправильный подсчёт шага
    a = (1 + alpha*(t+h[0])*(t+h[0])*x[1]*x[1])/2.;
    b = x[3]*alpha*(t+h[0])*(t+h[0])*x[1];
    k = b*b+b+1.+a*a;
    mu = fmax(3.*k,fmax(sqrt(3.*k/2. - 3./4.),pow(3.*a*a/16.,1./3.)));
    delta = err + exp(h[0]*sqrt(mu))*delta;
  }
  return delta;
}

double *X(double eps,double *u,double* beta){
  double x[6]={0,0,1,beta[0],beta[1],-beta[1]};
  traj34(eps,x);
  u[0]=x[0]-1; u[1] = x[2];
  return u;
}

double *test_X(double eps,double *u,double *beta){
  eps = eps;
  u[0]=beta[0]*beta[0]/2.;
  return u;
}
double *test_X2(double eps,double *u,double *beta){
  eps = eps;
  u[0]=sin(beta[0]);
  u[1]=exp(beta[1]);
  return u;
}


double *dF(double eps,unsigned int dim,unsigned int pdim,double *df,double*(*F)(double,double*,double*),double *pt){
  unsigned int i,j;
  double u1[dim],u[pdim][dim],dpt[pdim][pdim];// u[i][j] - i-я переменная выросла, j-я компонента вектора значений
  F(eps,u1,pt);
  for(j=0;j<pdim;j++) for(i=0;i<pdim;i++) dpt[i][j]=pt[j];
  for(j=0;j<pdim;j++){
    dpt[j][j] = pt[j] + sqrt(eps);
    F(eps,u[j],dpt[j]);
  }
  for(j=0;j<pdim;j++) for(i=0;i<dim;i++) df[pdim*i + j] = (u[j][i] - u1[i])/sqrt(eps);
  return df;
}

/*
 * Обратную матрицу брать дело не царское, когда размерность больше двух.
 * Поэтому универсальный метод подождёт
 * df[0] df[1]		 df[3] -df[1]|__________1_________
 * df[2] df[3] 		-df[2]  df[0]|df[0]df[3]-df[2]df[1]
 *
 */
int shot2(double eps,double *a,double*(*F) (double,double*,double*)){
  double z[2]; double df[4],det;
  F(eps,z,a);
  printf("%e\t&%e\t&",z[0],z[1]);
  if (norm_max(2,z)<eps) return 0;
  dF(eps,2,2,df,F,a);
  det = df[0]*df[3] - df[2]*df[1];
  a[0] -= (z[0]*df[3] - z[1]*df[1])/det;
  a[1] -= (z[1]*df[0] - z[0]*df[2])/det;
  return 1;
}




int main(int argc,char**argv){
  double u[2],v[2],a[2],df[4],b[2],eps;
  int ne_ok=1;
  unsigned int i;
  if (argc<5) {printf("%s eps x0(2x0) y0 T\n",argv[0]); return -1;}
  eps = atof(argv[1]);
  printf("Гармонический осциллятор\n");
  double e1 = harmonic_oscillator(eps,atof(argv[2]),atof(argv[3]),osc_while);
  double e2 = harmonic_oscillator(eps/100.,atof(argv[2]),atof(argv[3]),osc_while);
  double e3 = harmonic_oscillator(eps/10000.,atof(argv[2]),atof(argv[3]),osc_while);
  printf("%e\n",fabs((e1-e2)/(e2-e3)));
  a[0]=0; a[1]=0; b[0] = atof(argv[2]); b[1] = atof(argv[3]);
  X(atof(argv[1]),u,a);
  a[0]+=sqrt(atof(argv[1]));
  X(atof(argv[1]),v,a);
  printf("%e\n",(u[0]-v[0])/sqrt(eps));
  dF(eps,2,2,df,X,a);
  dF(eps,1,1,df,test_X,b);
  dF(eps,2,2,df,test_X2,b);
  printf("%e\t%e\n",df[0],fabs(df[0]-b[0]));
  printf("%e\t%e\t%e\t%e\t%e\t%e\n%e\t%e\t%e\t%e\t%e\t%e\n",
     df[0],df[1], cos(b[0]),0., fabs(df[0]-cos(b[0])),	 fabs(df[1]),
     df[2],df[3], 0.,exp(b[1]), fabs(df[2])		, fabs(df[3]-exp(b[1])));
  printf("Начало стрельбы\n");
  a[0]=0;a[1]=0;
  printf("$\\phi(1)-1$&\t $y(1)$ &$\\hat \\beta_1$ &\t$\\hat\\beta_2$ &\t $\\Delta\\beta_1$&\t$\\Delta\\beta_2$&\t ne ok \\\\\\hline\n");
  for (ne_ok=1,i=0; i<100 && ne_ok; i++){
    b[0]=a[0];b[1]=a[1];
    ne_ok = shot2(eps,a,X);
    printf("%e &\t%e &\t%e &\t%e &%d \\\\\\hline\n",b[0],b[1],fabs(b[0]-9.),fabs(b[1]-30.),ne_ok);
  }
  return 0;
}
