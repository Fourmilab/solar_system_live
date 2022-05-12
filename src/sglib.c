/*

        The Simple Graphics Library

        Implemented by John Walker in March of 1988.

        This is a simple three dimensional transformation and modeling
        library based on Jim Blinn's modeling primitives, as documented
        in Jim Blinn's Corner, IEEE Computer Graphics and Applications,
        October 1987, and subsequent columns.

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "sglib.h"

/*  Coordinate system transforms  */

matrix cT = {                      /*  Current transformation matrix */
        {1.0, 0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0, 0.0},
        {0.0, 0.0, 1.0, 0.0},
        {0.0, 0.0, 0.0, 1.0}
};

struct ctstack {
        struct ctstack *ctlast;    /* Previous coordinate system */
        matrix ctsm;               /* Saved coordinate system */
};

static struct ctstack *cts = NULL; /* Coordinate system stack */
static int ctdepth = 0;            /* Coordinate system nesting depth */

/*  Utility routines  */

/*  SGALLOC  --  Allocate a buffer and check for out of memory  */

char *sgalloc(int x)
{
        char *b;

        b = malloc(x);
        if (b == NULL) {
           fprintf(stderr, "\nBoom!!!  Memory capacity exceeded!\n");
           abort();
        }
        return b;
}

/*  General vector routines  */

/*  VECGET  --  Set vector from X, Y, and Z coordinates  */

void vecget(vector v, double x, double y, double z)
{
        v[X] = x;
        v[Y] = y;
        v[Z] = z;
        v[T] = 1.0;
}

/*  VECPUT  --  Store vector into X, Y, and Z coordinates  */

void vecput(double *x, double *y, double *z, vector v)
{
        double w;

        w = v[T];
        *x = v[X] / w;
        *y = v[Y] / w;
        *z = v[Z] / w;
}

/*  VECCOPY  --  Copy vector to another  */

void veccopy(vector vo, vector v)
{
        register int i;

        for (i = X; i <= T; i++)
           vo[i] = v[i];
}

/*  VECXMAT  --  Multiply a vector by a matrix  */

void vecxmat(vector vo, vector v, matrix m)
{
        register int i, j;
        register double sum;

        for (i = 0; i < 4; i++) {
           sum = 0;
           for (j = 0; j < 4; j++) {
              sum += v[j] * m[j][i];
           }
           vo[i] = sum;
        }
}

/*  Vector algebra routines which operate on points  */

/*  POINTGET  --  Set point from X, Y, and Z coordinates  */

void pointget(point p, double x, double y, double z)
{
        p[X] = x;
        p[Y] = y;
        p[Z] = z;
}

/*  POINTCOPY  --  Copy point to another  */

void pointcopy(point po, point p)
{
        po[X] = p[X];
        po[Y] = p[Y];
        po[Z] = p[Z];
}

/*  VECDOT  --  Computes the dot (inner) product of two vectors and
                returns the result as a double.  Since this will frequently
                be used on points as well as vectors, only the first
                three terms are computed.  */

double vecdot(point a, point b)
{
        int i;
        double product;

        product = 0.0;
        for (i = 0; i < 3; i++) {
           product += a[i] * b[i];
        }

        return product;
}

/*  VECCROSS  --  Computes the cross product of two vectors and stores
                  the result in a third.  This actually works on points;
                  if a vector is passed, the fourth item is ignored.  */

void veccross(point o, point a, point b)
{
        point r;

        r[X] = a[Y] * b[Z] - a[Z] * b[Y];
        r[Y] = a[Z] * b[X] - a[X] * b[Z];
        r[Z] = a[X] * b[Y] - a[Y] * b[X];

        pointcopy(o, r);
}

/*  VECADD  --  Add two vectors and store the sum in a third.
                Operates on points.  */

void vecadd(point o, point a, point b)
{
        o[X] = a[X] + b[X];
        o[Y] = a[Y] + b[Y];
        o[Z] = a[Z] + b[Z];
}

/*  VECSUB  --  Subtracts vector b from vector a and stores the
                result in vector o.  Expects points as arguments. */

void vecsub(point o, point a, point b)
{
        o[X] = a[X] - b[X];
        o[Y] = a[Y] - b[Y];
        o[Z] = a[Z] - b[Z];
}

/*  VECSCAL  --  Multiply vector by a scalar and store the result
                 in a second vector.  Expects points.  */

void vecscal(point o, point a, double s)
{
        o[X] = a[X] * s;
        o[Y] = a[Y] * s;
        o[Z] = a[Z] * s;
}

/*  VECMAG  --  Returns magnitude of a vector.  This expects a point
                and uses only the first three terms.  */

double vecmag(point a)
{
        return sqrt(a[X] * a[X] + a[Y] * a[Y] + a[Z] * a[Z]);
}

/*  VECNORM  --  Normalise vector and store normalised result in
                 a second vector.  Works on points.  */

void vecnorm(point o, point a)
{
        vecscal(o, a, 1.0 / vecmag(a));
}

/*  VECPRINT  --  Print a vector  */

void vecprint(vector v)
{
        int j;

        fprintf(stderr, "+-----------------------------------------+\n");
        fprintf(stderr, "|");
        for (j = 0; j < 4; j++) {
           fprintf(stderr, " %9.4f", v[j]);
        }
        fprintf(stderr, " |\n");
        fprintf(stderr, "+-----------------------------------------+\n");
}

/*  General matrix routines  */

/*  MATMUL  --  Multiply two 4 X 4 matrices, storing copy in a third.  */

void matmul(matrix o, matrix a, matrix b)
{
        register int i, j, k;
        register double sum;

        for (i = 0; i < 4; i++) {
           for (k = 0; k < 4; k++) {
              sum = 0.0;
              for (j = 0; j < 4; j++) {
                 sum += a[i][j] * b[j][k];
              }
              o[i][k] = sum;
           }
        }
}

/*  MATIDENT  --  Set a matrix to the identity matrix  */

void matident(matrix a)
{
        register int i, j;

        for (i = 0; i < 4; i++) {
           for (j = 0; j < 4; j++) {
              a[i][j] = (i == j) ? 1.0 : 0.0;
           }
        }
}

/*  MATCOPY  --  Copy a matrix to another  */

void matcopy(matrix o, matrix a)
{
        register int i, j;

        for (i = 0; i < 4; i++) {
           for (j = 0; j < 4; j++) {
              o[i][j] = a[i][j];
           }
        }
}

/*  MATPRINT  --  Print a matrix  */

void matprint(matrix a)
{
        int i, j;

        fprintf(stderr, "+-----------------------------------------+\n");
        for (i = 0; i < 4; i++) {
           fprintf(stderr, "|");
           for (j = 0; j < 4; j++) {
              fprintf(stderr, " %9.4f", a[i][j]);
           }
           fprintf(stderr, " |\n");
        }
        fprintf(stderr, "+-----------------------------------------+\n");
}

/*  Transformation matrix construction routines  */

/*  MATTRAN  --  Build translation matrix  */

void mattran(matrix m, double tx, double ty, double tz)
{
        matident(m);
        m[T][X] = tx;
        m[T][Y] = ty;
        m[T][Z] = tz;
}

/*  MATSCAL  --  Build scaling matrix  */

void matscal(matrix m, double sx, double sy, double sz)
{
        matident(m);
        m[X][X] = sx;
        m[Y][Y] = sy;
        m[Z][Z] = sz;
}

/*  MATROT  --  Build rotation matrix.  THETA is the rotation
                angle, in radians, and J is the axis about which
                the rotation is to be performed, expressed as one
                 of the manifest constants X, Y, or Z. */

void matrot(matrix m, double theta, int j)
{
        double s, c;

        s = sin(theta);
        c = cos(theta);

        matident(m);
        switch (j) {

           case X:
              m[1][1] = m[2][2] = c;
              m[1][2] = -s;
              m[2][1] = s;
              break;

           case Y:
              m[0][0] = m[2][2] = c;
              m[0][2] = s;
              m[2][0] = -s;
              break;

           case Z:
              m[0][0] = m[1][1] = c;
              m[0][1] = -s;
              m[1][0] = s;
              break;

           default:
              fprintf(stderr, "\nInvalid axis (J) argument %d to matrot.\n",
                 j);
              abort();
        }
}

/*  MATPERS  --  Build perspective transformation matrix.  ALPHA is
                 the field of view, ZN is the near clipping plane,
                 and ZF is the far clipping plane.  */

void matpers(matrix m, double alpha, double zn, double zf)
{
        double s, c, q;

        s = sin(alpha / 2.0);
        c = cos(alpha / 2.0);
        q = s / (1.0 - zn / zf);
        matident(m);
        m[X][X] = m[Y][Y] = c;
        m[Z][Z] = q;
        m[T][Z] = - q * zn;
        m[Z][T] = s;
        m[T][T] = 0.0;
}

/*  MATORIE  --  Specify explicit orientation  */

void matorie(matrix m, double a, double b, double c,
                       double d, double e, double f,
                       double p, double q, double r)
{
        matident(m);
        m[0][0] = a;
        m[1][0] = b;
        m[2][0] = c;
        m[0][1] = d;
        m[1][1] = e;
        m[2][1] = f;
        m[0][2] = p;
        m[1][2] = q;
        m[2][2] = r;
}

/*  MATSHAD  --  Specify matrix for fake shadow generation.  The
                 light source is at X, Y, and Z, and W is FALSE
                 for a light source at infinity and TRUE for a
                 local light source. */

void matshad(matrix m, double x, double y, double z, int w)
{
        matident(m);
        m[0][0] = z;
        m[1][1] = z;
        m[2][0] = -x;
        m[2][1] = -y;
        m[2][2] = 0.0;
        m[2][3] = w ? -1.0 : 0.0;
        m[3][3] = z;
}

/*  Current coordinate system transformation composition routines  */

/*  TRAN  --  Compose translation matrix  */

void tran(double tx, double ty, double tz)
{
        matrix m, m1;

        mattran(m, tx, ty, tz);
        matmul(m1, m, cT);
        matcopy(cT, m1);
}

/*  SCAL  --  Build scaling matrix  */

void scal(double sx, double sy, double sz)
{
        matrix m, m1;

        matscal(m, sx, sy, sz);
        matmul(m1, m, cT);
        matcopy(cT, m1);
}

/*  ROT  --  Build rotation matrix.  THETA is the rotation
             angle, in radians, and J is the axis about which
             the rotation is to be performed, expressed as one
             of the manifest constants X, Y, or Z. */

void rot(double theta, int j)
{
        matrix m, m1;

        matrot(m, theta, j);
        matmul(m1, m, cT);
        matcopy(cT, m1);
}

/*  PERS  --  Build perspective transformation matrix.  ALPHA is
              the field of view, ZN is the near clipping plane,
              and ZF is the far clipping plane.  */

void pers(double alpha, double zn, double zf)
{
        matrix m, m1;

        matpers(m, alpha, zn, zf);
        matmul(m1, m, cT);
        matcopy(cT, m1);
}

/*  ORIE  --  Specify explicit orientation  */

void orie(double a, double b, double c,
          double d, double e, double f,
          double p, double q, double r)
{
        matrix m, m1;

        matorie(m, a, b, c, d, e, f, p, q, r);
        matmul(m1, m, cT);
        matcopy(cT, m1);
}

/*  SHAD  --  Compose matrix for fake shadow generation.  The
              light source is at X, Y, and Z, and W is FALSE
              for a light source at infinity and TRUE for a
              local light source. */

void shad(double x, double y, double z, int w)
{
        matrix m, m1;

        matshad(m, x, y, z, w);
        matmul(m1, m, cT);
        matcopy(cT, m1);
}


/*  Coordinate system push and pop routines  */

/*  PUSH  --  Save coordinate system  */

void push(void)
{
        struct ctstack *c;

        c = (struct ctstack *) sgalloc(sizeof(struct ctstack));
        c->ctlast = cts;
        matcopy(c->ctsm, cT);
        cts = c;
        ctdepth++;
}

/*  POP  -- Restore coordinate system  */

void pop(void)
{
        struct ctstack *c;

        if (ctdepth <= 0) {
           fprintf(stderr, "\nCoordinate system popped when none pushed.\n");
           abort();
        }
        c = cts;
        cts = c->ctlast;
        matcopy(cT, c->ctsm);
        free(c);
        ctdepth--;
}

/*  THEN  --  Pop old coordinate system, push new one  */

void then(void)
{
        pop();
        push();
}
