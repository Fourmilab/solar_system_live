/*

        Asteroid orbit calculation

*/

#include "vplanet.h"

/*  Elements of currently-tracked asteroid */

int aTracked = FALSE;           /* Tracking an asteroid ? */

struct asteroid_info ast_info = {0, "-none-"};  /* Information about currently-tracked asteroid */

/*  gKEPLER  --  Super-accurate solution to the equation of Kepler
                 for eccentricities from 0 through 1.  */

static int gKepler(const double e, const double t, const double q,
                   double *v, double *r)
{
    double f, x, d, m1;
    int i;
#define K   0.01720209895       /* Gaussian gravitational constant */

    if (e >= 1.0) {

#define LANDGRAF
#ifdef LANDGRAF

        /* Iterative solution for near-parabolic motion, based
           on the program by Werner Landgraf which implements the
           method of Karl Stumpff.  This program is presented in
           chapter 35 of Meeus, "Astronomical Algorithms", 2nd ed.

           This method is suitable for eccentricities in the
           range from 0.98 through 1.1.  We use it only for
           eccentricities greater than or equal to 1. */

        unsigned int z;
        const unsigned int maxiter = 4096;
        double q1, g, q2, q3, s, c = 1.0 / 3.0, d = 1.0E-12, d1 = 10000,
               f,
               l, s0, s1, z1, y, g1;

        if (t == 0) {
            *r = q;
            *v = 0;
            return TRUE;
        }

        q1 = (K / (2 * q)) * sqrt((1 + e) / q);
        g = (1 - e) / (1 + e);
        q2 = q1 * t;
        s = 2 / (3 * abs(q2));
        s = 2 / tan(2 * atan(pow(tan(atan(s) / 2), c)));
        if (t < 0) {
            s = -s;
        }
        if (e != 1.0) {
            l = 0;

            while (TRUE) {
                s0 = s;
                z = 1;
                y = s * s;
                g1 = (-y) * s;
                q3 = q2 + (2 * g * s * (y / 3));

                while (TRUE) {
                    z += 1;
                    g1 = (-g1) * g * y;
                    z1 = (z - (z + 1) * g) / ((2 * z) + 1);
                    f = z1 * g1;
                    q3 = q3 + f;
                    if ((z > maxiter) || (abs(f) > d1)) {
                        goto noconv;
                    }
                    if (abs(f) <= d) {
                        break;
                    }
                }
                l += 1;
                if (l > maxiter) {
                    goto noconv;
                }
                while (TRUE) {
                    s1 = s;
                    s = ((2 * s * s * (s / 3)) + q3) / ((s * s) + 1);
                    if (abs(s - s1) <= d) {
                        break;
                    }
                }
                if (abs(s - s0) <= d) {
                    break;
                }
            }
        }
        *v = 2 * atan(s);
        *r = (q * (1 + e)) / (1 + (e * cos(*v)));
        if (*v < 0) {
            *v += 2 * PI;
        }
        return TRUE;

noconv:;
#ifdef AST_DEBUG
fprintf(stderr, "gKepler failed to converge for e = %.4f t = %.4f q = %.4f\n", e, t, q);
#endif
        return FALSE;
#else

        /* Gauss method for high eccentricities. */

        double a, b, c, a2, apow, b0, b1, b2, b3, c1, d1, w1, t1, t2, s, g, v1;
        int s1;

        static double B[] = {
            1, 0, -0.017142857, -0.003809524,
            -0.001104267, -0.000367358, -0.000131675, -0.000049577
        };
        static double S[] = {
            1, -0.8, 0.04571429, 0.01523810, 0.00562820, 0.00218783,
            0.00087905, 0.00036155
        };

        a = sqrt((1 + 9 * e) / 10);
        b = 5 * (1 - e) / (1 + 9 * e);
        c = sqrt(5 * (1 + e) / (1 + 9 * e));

        b1 = 3 * a * K * t / sqrt(2 * q * q * q);
        b2 = 1;
        do {
           w1 = b2 * b1;
           b3 = atan2(2.0, w1);
           t1 = sin(b3 / 2) / cos(b3 / 2);
           s1 = sgn(t1);
           t1 = abs(t1);
           t2 = pow(t1, (1 / 3.0)) * s1;
           g = atan(t2);
           s = 2 * cos(2 * g) / sin(2 * g);
           a2 = b * s * s;
           b0 = b2;
           b2 = 0;

           if ((abs(a2) > 0.3) && (e >= 1)) {
               return FALSE;
           }

           apow = 1.0;
           for (i = 0; i < 8; i++) {
               b2 += B[i] * apow;
               apow *= a2;
           }
        } while (abs(b2 - b0) > 1E-8);

        c1 = 0;
        apow = 1.0;
        for (i = 0; i < 8; i++) {
            c1 += S[i] * apow;
            apow *= a2;
        }
        c1 = sqrt(1 / c1);
        v1 = c * c1 * s;
        d1 = 1 / (1 + a2 * c1 * c1);
        x = 2 * atan(v1);
        *r = q * d1 * (1 + v1 * v1);
#endif /* LANDGRAF */
    } else {

        /* Solution by binary search by Roger W. Sinnott, Sky and Telescope,
           Vol. 70, page 159 (August 1985). */

        double m, a1, ev;

        a1 = q / (1 - e);
        m = K * t * pow(a1, -1.5);

        f = sgn(m);
        m = abs(m) / (2 * PI);
        m = (m - ((int) m)) * 2 * PI * f;
        if (m < 0) {
            m += 2 * PI;
        }
        f = (m > PI) ? -1 : 1;
        if (m > PI) {
            m = (2 * PI) - m;
        }
        x = PI / 2;
        d = PI / 4;
        for (i = 0; i < 53; i++) {
            m1 = x - e * sin(x);
            x = x + sgn(m - m1) * d;
            d /= 2;
        }
        x *= f;
        ev = sqrt((1 + e) / (1 - e));
        *r = a1 * (1 - e * cos(x));
        x = 2 * atan2(ev * sin(x / 2), cos(x / 2));
    }
    if (x < 0) {
        x += 2 * PI;
    }
    *v = x;
    return TRUE;
}

static int cometary_orbit(const double peri_date, const double peri_dist,
                   const double eccentricity, const double arg_perihelion,
                   const double long_asc_node, const double inclination,
                   const double jd,
                   double *hx, double *hy, double *hz)
{

    double e = dtr(obliqeq(peri_date));     /* Obliquity of the ecliptic at perihelion date */
    double w = dtr(arg_perihelion);
    double n = dtr(long_asc_node);
    double i = dtr(inclination);

    double w1, w2, n1, n2, i1, i2, e1, e2, p7, p8, p9, q7, q8, q9, v, r;

    w1 = sin(w);
    w2 = cos(w);
    n1 = sin(n);
    n2 = cos(n);
    i1 = sin(i);
    i2 = cos(i);
    e1 = sin(e);
    e2 = cos(e);

    p7 = w2 * n2 - w1 * n1 * i2;
    p8 = (w2 * n1 + w1 * n2 * i2) * e2 - w1 * i1 * e1;
    p9 = (w2 * n1 + w1 * n2 * i2) * e1 + w1 * i1 * e2;
    q7 = -w1 * n2 - w2 * n1 * i2;
    q8 = (-w1 * n1 + w2 * n2 * i2) * e2 - w2 * i1 * e1;
    q9 = (-w1 * n1 + w2 * n2 * i2) * e1 + w2 * i1 * e2;

    if (gKepler(eccentricity, jd - peri_date, peri_dist,
                &v, &r)) {
        double x1, y1;

        x1 = r * cos(v);
        y1 = r * sin(v);


        *hx = p7 * x1 + q7 * y1;
        *hy = p8 * x1 + q8 * y1;
        *hz = p9 * x1 + q9 * y1;

    } else {
        return FALSE;
    }
    return TRUE;
}

static int asteroidal_orbit(const double epoch, const double mean_anomaly,
                     const double semimajor_axis,
                     const double eccentricity, const double arg_perihelion,
                     const double long_asc_node, const double inclination,
                     const double jd,
                     double *hx, double *hy, double *hz)
{
    return cometary_orbit(epoch -
        sqrt(semimajor_axis * semimajor_axis * semimajor_axis) * mean_anomaly *
            365.2422 / 360,
        semimajor_axis - semimajor_axis * eccentricity,
        eccentricity,
        arg_perihelion,
        long_asc_node,
        inclination,
        jd,
        hx, hy, hz);
}

/*  TRACKASTEROID  --  Determine location of asteroid at a given point.  */

void trackAsteroid(double jd, double *ra, double *dec, double *dist,
                   double *hlong, double *hlat, double *hrv, int quick)
{
    int ok;
    double x, y, z, sx, sy, sz, obelix, lon, hra, hdec;

    if (ast_info.cometary) {
        ok = cometary_orbit(ast_info.PeriDate, ast_info.PeriAU,
                            ast_info.Eccentricity, ast_info.ArgP,
                            ast_info.LANode, ast_info.Inclination,
                            jd, &x, &y, &z);
    } else {
        /* Fill in implicit fields for asteroids */
        ast_info.PeriDate = ast_info.Epoch -
                            sqrt(ast_info.SemiMajorAU * ast_info.SemiMajorAU *
                                 ast_info.SemiMajorAU) * ast_info.mAnomaly *
                                365.24219879 / 360;
        ast_info.PeriAU = ast_info.SemiMajorAU * (1 - ast_info.Eccentricity);
        ok = asteroidal_orbit(ast_info.Epoch, ast_info.mAnomaly, ast_info.SemiMajorAU,
                              ast_info.Eccentricity, ast_info.ArgP,
                              ast_info.LANode, ast_info.Inclination,
                              jd, &x, &y, &z);
    }

    if (ok) {
        obelix = dtr(obliqeq(jd));
        hra = atan2(y, x);
        hdec = atan2(z, sqrt(x * x + y * y));
        *hrv = sqrt(x * x + y * y + z * z);
        *hlong = rtd(atan2(sin(hra) * cos(obelix) + tan(hdec) * sin(obelix), cos(hra)));
        *hlat = rtd(asin(sin(hdec) * cos(obelix) - cos(hdec) * sin(obelix) * sin(hra)));

        if (!quick) {               /* Quick calculation doesn't need RA, Dec, distance */
            lon = dtr(planet_info[0].hlong);
            sx = planet_info[0].dist * cos(lon);
            sy = planet_info[0].dist * sin(lon) * cos(obelix);
            sz = planet_info[0].dist * sin(lon) * sin(obelix);

            x += sx;
            y += sy;
            z += sz;

            *ra = rtd(atan2(y, x));
            *dec = rtd(atan2(z, sqrt(x * x + y * y)));
            *dist = sqrt(x * x + y * y + z * z);
/*          *hlong = rtd(atan2(sin(dtr(*ra)) * cos(obelix) + tan(dtr(*dec)) * sin(obelix), cos(dtr(*ra)))); */
/*          *hlat = rtd(asin(sin(dtr(*dec)) * cos(obelix) - cos(dtr(*dec)) * sin(obelix) * sin(dtr(*ra)))); */
        }
    } else {
        /* Kepler equation solver didn't converge.  Signal is hrv = 0. */
        *hrv = *hlong = *hlat = *ra = *dec = *dist = 0.0;
    }
}
