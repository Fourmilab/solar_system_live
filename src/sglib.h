/*

        Definitions for the simple graphics library

*/

#define PI 3.14159265358979323846

#define abs(x) (((x) < 0) ? (-(x)) : (x))
#define fzero(x) (abs((x)) < 0.01)
#define torad(x) ((x) * (PI / 180.0))

#define X  0                       /* Coordinate indices */
#define Y  1
#define Z  2
#define T  3

#define TRUE  1
#define FALSE 0

typedef double point[3];           /* Three dimensional point */
typedef double vector[4];          /* Homogeneous coordinate vector */
typedef double matrix[4][4];       /* Transformation matrix */

extern matrix cT;                  /* Current transformation matrix */

extern char *sgalloc(int x);       /* Memory allocator interface */

extern void vecget(vector v, double x, double y, double z);
extern void vecput(double *x, double *y, double *z, vector v);
extern void veccopy(vector vo, vector v);
extern void vecxmat(vector vo, vector v, matrix m);

extern void pointget(point p, double x, double y, double z);
extern void pointcopy(point po, point p);

extern double vecdot(point a, point b);
extern void veccross(point o, point a, point b);
extern void vecadd(point o, point a, point b);
extern void vecsub(point o, point a, point b);
extern void vecscal(point o, point a, double s);
extern double vecmag(point a);
extern void vecnorm(point o, point a);
extern void vecprint(vector v);

extern void matmul(matrix o, matrix a, matrix b);
extern void matident(matrix a);
extern void matcopy(matrix o, matrix a);
extern void matprint(matrix a);

extern void mattran(matrix m, double tx, double ty, double tz);
extern void matscal(matrix m, double sx, double sy, double sz);
extern void matrot(matrix m, double theta, int j);
extern void matpers(matrix m, double alpha, double zn, double zf);
extern void matorie(matrix m, double a, double b, double c,
                              double d, double e, double f,
                              double p, double q, double r);
extern void matshad(matrix m, double x, double y, double z, int w);

extern void tran(double tx, double ty, double tz);
extern void scal(double sx, double sy, double sz);
extern void rot(double theta, int j);
extern void pers(double alpha, double zn, double zf);
extern void orie(double a, double b, double c,
                 double d, double e, double f,
                 double p, double q, double r);
extern void shad(double x, double y, double z, int w);

extern void push(void);
extern void pop(void);
extern void then(void);
