/* Векторная алгебра */

struct vect{
	int length;
	double * array;
};

int init(struct vect *v, int length);
int copy(struct vect vect1, struct vect *vect2);
int init_copy(struct vect *a, struct vect b); 
double scalar(struct vect vect1, struct vect vect2);
struct vect vf(struct vect v,  double t, double (*f[])(struct vect,double));
struct vect sum(struct vect a, struct vect b);
struct vect minus(struct vect a, struct vect b);
double norm(struct vect a);
double norm_max(struct vect a);
double norm_sum(struct vect a);
void print(struct vect v);
void fprint(struct vect v, FILE *f);