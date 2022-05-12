/*

         Earth Viewer  --  Astronomical calculation routines

                          by John Walker
                       http://www.fourmilab.ch/

*/

#include "vplanet.h"

/*  Astronomical constants  */

#define epoch       2444238.5      /* 1980 January 0.0 */

/*  Constants defining the Sun's apparent orbit  */

#define elonge      278.833540     /* Ecliptic longitude of the Sun
                                      at epoch 1980.0 */
#define elongp      282.596403     /* Ecliptic longitude of the Sun at
                                      perigee */
#define eccent      0.016718       /* Eccentricity of Earth's orbit */
#define sunangsiz   0.533128       /* Sun's angular size, degrees, at
                                      semi-major axis distance */

/*  Elements of the Moon's orbit, epoch 1980.0  */

#define mmlong      64.975464      /* Moon's mean longitude at the epoch */
#define mmlongp     349.383063     /* Mean longitude of the perigee at the
                                      epoch */
#define mlnode      151.950429     /* Mean longitude of the node at the
                                      epoch */
#define minc        5.145396       /* Inclination of the Moon's orbit */
#define mecc        0.054900       /* Eccentricity of the Moon's orbit */
#define mangsiz     0.5181         /* Moon's angular size at distance a
                                      from Earth */
#define msmax       384401.0       /* Semi-major axis of Moon's orbit in km */
#define mparallax   0.9507         /* Parallax at distance a from Earth */

#define dsin(x) (sin(dtr((x))))     /* Sin from deg */
#define dcos(x) (cos(dtr((x))))     /* Cos from deg */

/*  UCTTOJ  --  Convert GMT date and time to astronomical
                Julian time (i.e. Julian date plus day fraction,
                expressed as a double).  */

double ucttoj(const long year, const int mon, const int mday,
              const int hour, const int min, const int sec)
{

    /* Algorithm as given in Meeus, Astronomical Algorithms, Chapter 7, page 61 */

    int a, b, m;
    long y;

    assert(mon  >= 0 && mon  < 12);
    assert(mday >  0 && mday < 32);
    assert(hour >= 0 && hour < 24);
    assert(min  >= 0 && min  < 60);
    assert(sec  >= 0 && sec  < 60);

    m = mon + 1;
    y = year;

    if (m <= 2) {
        y--;
        m += 12;
    }

    /* Determine whether date is in Julian or Gregorian calendar based on
       canonical date of calendar reform. */

    if ((year < 1582) || ((year == 1582) && ((mon < 9) || (mon == 9 && mday < 5)))) {
        b = 0;
    } else {
        a = ((int) (y / 100));
        b = 2 - a + (a / 4);
    }

    return (((long) (365.25 * (y + 4716))) + ((int) (30.6001 * (m + 1))) +
                mday + b - 1524.5) +
            ((sec + 60L * (min + 60L * hour)) / 86400.0);
}

/*  JTIME  --  Convert a Unix date and time (tm) structure to astronomical
               Julian time (i.e. Julian date plus day fraction,
               expressed as a double).  */

double jtime(const struct tm *t)
{
    return ucttoj(t->tm_year + 1900, t->tm_mon, t->tm_mday, t->tm_hour, t->tm_min, t->tm_sec);
}

/*  JYEAR  --  Convert  Julian  date  to  year,  month, day, which are
               returned via integer pointers to integers (note that year is a long).  */

void jyear(double td, long *yy, int *mm, int *dd)
{
    double z, f, a, alpha, b, c, d, e;

    td += 0.5;
    z = floor(td);
    f = td - z;

    if (z < 2299161.0) {
        a = z;
    } else {
        alpha = floor((z - 1867216.25) / 36524.25);
        a = z + 1 + alpha - floor(alpha / 4);
    }

    b = a + 1524;
    c = floor((b - 122.1) / 365.25);
    d = floor(365.25 * c);
    e = floor((b - d) / 30.6001);

    *dd = (int) (b - d - floor(30.6001 * e) + f);
    *mm = (int) ((e < 14) ? (e - 1) : (e - 13));
    *yy = (long) ((*mm > 2) ? (c - 4716) : (c - 4715));
}

/*  JHMS  --  Convert Julian time to hour, minutes, and seconds.  */

void jhms(double j, int *h, int *m, int *s)
{
    long ij;

    j += 0.5;                 /* Astronomical to civil */
    ij = (long) ((j - floor(j)) * 86400.0);
    *h = (int) (ij / 3600L);
    *m = (int) ((ij / 60L) % 60L);
    *s = (int) (ij % 60L);
}

/*  KEPLER  --  Solve the equation of Kepler.  */

static double kepler(double m, double ecc)
{
    double e, delta;
#define EPSILON 1E-6

    e = m = dtr(m);
    do {
        delta = e - ecc * sin(e) - m;
        e -= delta / (1 - ecc * cos(e));
    } while (abs(delta) > EPSILON);
    return e;
}

/*  SUNPOS  --  Calculate position of the Sun.  JD is the Julian  date
                of  the  instant for which the position is desired and
                APPARENT should be nonzero if  the  apparent  position
                (corrected  for  nutation  and aberration) is desired.
                The Sun's co-ordinates are returned  in  RA  and  DEC,
                both  specified  in degrees (divide RA by 15 to obtain
                hours).  The radius vector to the Sun in  astronomical
                units  is returned in RV and the Sun's longitude (true
                or apparent, as desired) is  returned  as  degrees  in
                SLONG.  */

void sunpos(const double jd, const int apparent,
            double *ra, double *dec, double *rv, double *slong)
{
    double t, t2, t3, l, m, e, ea, v, theta, omega, eps;

    /* Time, in Julian centuries of 36525 ephemeris days,
       measured from the epoch 1900 January 0.5 ET. */

    t = (jd - 2415020.0) / JulianCentury;
    t2 = t * t;
    t3 = t2 * t;

    /* Geometric mean longitude of the Sun, referred to the
       mean equinox of the date. */

    l = fixangle(279.69668 + 36000.76892 * t + 0.0003025 * t2);

    /* Sun's mean anomaly. */

    m = fixangle(358.47583 + 35999.04975 * t - 0.000150 * t2 - 0.0000033 * t3);

    /* Eccentricity of the Earth's orbit. */

    e = 0.01675104 - 0.0000418 * t - 0.000000126 * t2;

    /* Eccentric anomaly. */

    ea = kepler(m, e);

    /* True anomaly */

    v = fixangle(2 * rtd(atan(sqrt((1 + e) / (1 - e))  * tan(ea / 2))));

    /* Sun's true longitude. */

    theta = l + v - m;

    /* Obliquity of the ecliptic. */

    eps = obliqeq(jd);

    /* Corrections for Sun's apparent longitude, if desired. */

    if (apparent) {
       omega = fixangle(259.18 - 1934.142 * t);
       theta = theta - 0.00569 - 0.00479 * sin(dtr(omega));
       eps += 0.00256 * cos(dtr(omega));
    }

    /* Return Sun's longitude and radius vector */

    *slong = theta;
    *rv = (1.0000002 * (1 - e * e)) / (1 + e * cos(dtr(v)));

    /* Determine solar co-ordinates. */

    *ra =
    fixangle(rtd(atan2(cos(dtr(eps)) * sin(dtr(theta)), cos(dtr(theta)))));
    *dec = rtd(asin(sin(dtr(eps)) * sin(dtr(theta))));
}

/*  GMST  --  Calculate Greenwich Mean Siderial Time for a given
              instant expressed as a Julian date and fraction.  */

double gmst(const double jd)
{
    double t, theta0;

    /* Time, in Julian centuries of 36525 ephemeris days,
       measured from the epoch 1900 January 0.5 ET. */

    t = ((floor(jd + 0.5) - 0.5) - 2415020.0) / JulianCentury;

    theta0 = 6.6460656 + 2400.051262 * t + 0.00002581 * t * t;

    t = (jd + 0.5) - (floor(jd + 0.5));

    theta0 += (t * 24.0) * 1.002737908;

    theta0 = (theta0 - 24.0 * (floor(theta0 / 24.0)));

    return theta0;
}

/*  PHASE  --  Calculate phase of moon as a fraction:

    The  argument  is  the  time  for  which  the  phase is requested,
    expressed as a Julian date and fraction.  Returns  the  terminator
    phase  angle  as a percentage of a full circle (i.e., 0 to 1), and
    stores into pointer arguments  the  illuminated  fraction  of  the
    Moon's  disc, the Moon's age in days and fraction, the distance of
    the Moon from the centre of the Earth, and  the  angular  diameter
    subtended  by the Moon as seen by an observer at the centre of the
    Earth.
*/

double phase(
  const double  pdate,        /* Date for which to calculate phase */
  double  *pphase,            /* Illuminated fraction */
  double  *mage,              /* Age of moon in days */
  double  *dist,              /* Distance in kilometres */
  double  *angdia,            /* Angular diameter in degrees */
  double  *sudist,            /* Distance to Sun */
  double  *suangdia)                  /* Sun's angular diameter */
{

    double Day, N, M, Ec, Lambdasun, ml, MM, MN, Ev, Ae, A3, MmP,
       mEc, A4, lP, Varia, lPP, NP, y, x, Lambdamoon,
       MoonAge, MoonPhase,
       MoonDist, MoonDFrac, MoonAng,
       F, SunDist, SunAng;

    /* Calculation of the Sun's position */

    Day = pdate - epoch;                    /* Date within epoch */
    N = fixangle((360 / 365.2422) * Day);   /* Mean anomaly of the Sun */
    M = fixangle(N + elonge - elongp);      /* Convert from perigee
                                               co-ordinates to epoch 1980.0 */
    Ec = kepler(M, eccent);                 /* Solve equation of Kepler */
    Ec = sqrt((1 + eccent) / (1 - eccent)) * tan(Ec / 2);
    Ec = 2 * rtd(atan(Ec));                 /* True anomaly */
    Lambdasun = fixangle(Ec + elongp);      /* Sun's geocentric ecliptic
                                               longitude */
    /* Orbital distance factor */
    F = ((1 + eccent * cos(dtr(Ec))) / (1 - eccent * eccent));
    SunDist = SunSMAX / F;                  /* Distance to Sun in km */
    SunAng = F * sunangsiz;                 /* Sun's angular size in degrees */

    /* Calculation of the Moon's position */

    /* Moon's mean longitude */
    ml = fixangle(13.1763966 * Day + mmlong);

    /* Moon's mean anomaly */
    MM = fixangle(ml - 0.1114041 * Day - mmlongp);

    /* Moon's ascending node mean longitude */
    MN = fixangle(mlnode - 0.0529539 * Day);

    /* Evection */
    Ev = 1.2739 * sin(dtr(2 * (ml - Lambdasun) - MM));

    /* Annual equation */
    Ae = 0.1858 * sin(dtr(M));

    /* Correction term */
    A3 = 0.37 * sin(dtr(M));

    /* Corrected anomaly */
    MmP = MM + Ev - Ae - A3;

    /* Correction for the equation of the centre */
    mEc = 6.2886 * sin(dtr(MmP));

    /* Another correction term */
    A4 = 0.214 * sin(dtr(2 * MmP));

    /* Corrected longitude */
    lP = ml + Ev + mEc - Ae + A4;

    /* Variation */
    Varia = 0.6583 * sin(dtr(2 * (lP - Lambdasun)));

    /* True longitude */
    lPP = lP + Varia;

    /* Corrected longitude of the node */
    NP = MN - 0.16 * sin(dtr(M));

    /* Y inclination coordinate */
    y = sin(dtr(lPP - NP)) * cos(dtr(minc));

    /* X inclination coordinate */
    x = cos(dtr(lPP - NP));

    /* Ecliptic longitude */
    Lambdamoon = rtd(atan2(y, x));
    Lambdamoon += NP;

    /* Ecliptic latitude */
/*    BetaM = rtd(asin(sin(dtr(lPP - NP)) * sin(dtr(minc))));  */

    /* Calculation of the phase of the Moon */

    /* Age of the Moon in degrees */
    MoonAge = lPP - Lambdasun;

    /* Phase of the Moon */
    MoonPhase = (1 - cos(dtr(MoonAge))) / 2;

    /* Calculate distance of moon from the centre of the Earth */

    MoonDist = (msmax * (1 - mecc * mecc)) /
           (1 + mecc * cos(dtr(MmP + mEc)));

    /* Calculate Moon's angular diameter */

    MoonDFrac = MoonDist / msmax;
    MoonAng = mangsiz / MoonDFrac;

    /* Calculate Moon's parallax */

/*    MoonPar = mparallax / MoonDFrac;  */

    *pphase = MoonPhase;
    *mage = SynMonth * (fixangle(MoonAge) / 360.0);
    *dist = MoonDist;
    *angdia = MoonAng;
    *sudist = SunDist;
    *suangdia = SunAng;
    return fixangle(MoonAge) / 360.0;
}

/*  HIGHMOON  --  High precision calculation of the Moon's geocentric longitude,
                  latitude, and radius vector (distance) as given in Chapter 45
                  of Meeus' "Astronomical Algorithms".  */


#define NTERMS  60                    /* Number of terms for L, B, and R */

/*  Coefficients for L and R terms.  */

static const signed char lrCoeff[NTERMS][4] = {
    {0, 0, 1, 0},
    {2, 0,-1, 0},
    {2, 0, 0, 0},
    {0, 0, 2, 0},
    {0, 1, 0, 0},
    {0, 0, 0, 2},
    {2, 0,-2, 0},
    {2,-1,-1, 0},
    {2, 0, 1, 0},
    {2,-1, 0, 0},
    {0, 1,-1, 0},
    {1, 0, 0, 0},
    {0, 1, 1, 0},
    {2, 0, 0,-2},
    {0, 0, 1, 2},
    {0, 0, 1,-2},
    {4, 0,-1, 0},
    {0, 0, 3, 0},
    {4, 0,-2, 0},
    {2, 1,-1, 0},
    {2, 1, 0, 0},
    {1, 0,-1, 0},
    {1, 1, 0, 0},
    {2,-1, 1, 0},
    {2, 0, 2, 0},
    {4, 0, 0, 0},
    {2, 0,-3, 0},
    {0, 1,-2, 0},
    {2, 0,-1, 2},
    {2,-1,-2, 0},
    {1, 0, 1, 0},
    {2,-2, 0, 0},
    {0, 1, 2, 0},
    {0, 2, 0, 0},
    {2,-2,-1, 0},
    {2, 0, 1,-2},
    {2, 0, 0, 2},
    {4,-1,-1, 0},
    {0, 0, 2, 2},
    {3, 0,-1, 0},
    {2, 1, 1, 0},
    {4,-1,-2, 0},
    {0, 2,-1, 0},
    {2, 2,-1, 0},
    {2, 1,-2, 0},
    {2,-1, 0,-2},
    {4, 0, 1, 0},
    {0, 0, 4, 0},
    {4,-1, 0, 0},
    {1, 0,-2, 0},
    {2, 1, 0,-2},
    {0, 0, 2,-2},
    {1, 1, 1, 0},
    {3, 0,-2, 0},
    {4, 0,-3, 0},
    {2,-1, 2, 0},
    {0, 2, 1, 0},
    {1, 1,-1, 0},
    {2, 0, 3, 0},
    {2, 0,-1,-2}
};

/*  Coefficients for B terms.  */

static const signed char bCoeff[NTERMS][4] = {
    {0, 0, 0, 1},
    {0, 0, 1, 1},
    {0, 0, 1,-1},
    {2, 0, 0,-1},
    {2, 0,-1, 1},
    {2, 0,-1,-1},
    {2, 0, 0, 1},
    {0, 0, 2, 1},
    {2, 0, 1,-1},
    {0, 0, 2,-1},
    {2,-1, 0,-1},
    {2, 0,-2,-1},
    {2, 0, 1, 1},
    {2, 1, 0,-1},
    {2,-1,-1, 1},
    {2,-1, 0, 1},
    {2,-1,-1,-1},
    {0, 1,-1,-1},
    {4, 0,-1,-1},
    {0, 1, 0, 1},
    {0, 0, 0, 3},
    {0, 1,-1, 1},
    {1, 0, 0, 1},
    {0, 1, 1, 1},
    {0, 1, 1,-1},
    {0, 1, 0,-1},
    {1, 0, 0,-1},
    {0, 0, 3, 1},
    {4, 0, 0,-1},
    {4, 0,-1, 1},
    {0, 0, 1,-3},
    {4, 0,-2, 1},
    {2, 0, 0,-3},
    {2, 0, 2,-1},
    {2,-1, 1,-1},
    {2, 0,-2, 1},
    {0, 0, 3,-1},
    {2, 0, 2, 1},
    {2, 0,-3,-1},
    {2, 1,-1, 1},
    {2, 1, 0, 1},
    {4, 0, 0, 1},
    {2,-1, 1, 1},
    {2,-2, 0,-1},
    {0, 0, 1, 3},
    {2, 1, 1,-1},
    {1, 1, 0,-1},
    {1, 1, 0, 1},
    {0, 1,-2,-1},
    {2, 1,-1,-1},
    {1, 0, 1, 1},
    {2,-1,-2,-1},
    {0, 1, 2, 1},
    {4, 0,-2,-1},
    {4,-1,-1,-1},
    {1, 0, 1,-1},
    {4, 0, 1,-1},
    {1, 0,-1,-1},
    {4,-1, 0,-1},
    {2,-2, 0, 1}
};

static const long lTerms[NTERMS] = {
     6288774,
     1274027,
      658314,
      213618,
     -185116,
     -114332,
       58793,
       57066,
       53322,
       45758,
      -40923,
      -34720,
      -30383,
       15327,
      -12528,
       10980,
       10675,
       10034,
        8548,
       -7888,
       -6766,
       -5163,
        4987,
        4036,
        3994,
        3861,
        3665,
       -2689,
       -2602,
        2390,
       -2348,
        2236,
       -2120,
       -2069,
        2048,
       -1773,
       -1595,
        1215,
       -1110,
        -892,
        -810,
         759,
        -713,
        -700,
         691,
         596,
         549,
         537,
         520,
        -487,
        -399,
        -381,
         351,
        -340,
         330,
         327,
        -323,
         299,
         294,
           0
};

static const long rTerms[NTERMS] = {
   -20905355,
    -3699111,
    -2955968,
     -569925,
       48888,
       -3149,
      246158,
     -152138,
     -170733,
     -204586,
     -129620,
      108743,
      104755,
       10321,
           0,
       79661,
      -34782,
      -23210,
      -21636,
       24208,
       30824,
       -8379,
      -16675,
      -12831,
      -10445,
      -11650,
       14403,
       -7003,
           0,
       10056,
        6322,
       -9884,
        5751,
           0,
       -4950,
        4130,
           0,
       -3958,
           0,
        3258,
        2616,
       -1897,
       -2117,
        2354,
           0,
           0,
       -1423,
       -1117,
       -1571,
       -1739,
           0,
       -4421,
           0,
           0,
           0,
           0,
        1165,
           0,
           0,
        8752
};

static const long bTerms[NTERMS] = {
     5128122,
      280602,
      277693,
      173237,
       55413,
       46271,
       32573,
       17198,
        9266,
        8822,
        8216,
        4324,
        4200,
       -3359,
        2463,
        2211,
        2065,
       -1870,
        1828,
       -1794,
       -1749,
       -1565,
       -1491,
       -1475,
       -1410,
       -1344,
       -1335,
        1107,
        1021,
         833,
         777,
         671,
         607,
         596,
         491,
        -451,
         439,
         422,
         421,
        -366,
        -351,
         331,
         315,
         302,
        -283,
        -229,
         223,
         223,
        -220,
        -220,
        -185,
         181,
        -177,
         176,
         166,
        -164,
         132,
        -119,
         115,
         107
};

void highmoon(const double jd, double *l, double *b, double *r)
{
    double t, t2, t3, t4, lprime, d, m, mprime, f, a1, a2, a3, e[3],
           sigmaL, sigmaB, sigmaR, ang;
    int i;

    t = (jd - J2000) / JulianCentury;

    t4 = t * (t3 = t * (t2 = t * t));

    lprime = 218.3164591 + 481267.88134236 * t - 0.0013268 * t2 +
                t3 / 538841.0 - t4 / 65194000.0;

    d = 297.8502042 + 445267.1115168 * t - 0.00163 * t2 +
                t3 / 545868.0 - t4 / 113065000.0;

    m = 357.5291092 + 35999.0502909 * t - 0.0001536 * t2 +
                t3 / 24490000.0;

    mprime = 134.9634114 + 477198.8676313 * t + 0.008997 * t2 +
                t3 / 69699.0 - t4 / 14712000.0;

    f = 93.2720993 + 483202.0175273 * t - 0.0034029 * t2 -
                t3 / 3526000.0 + t4 / 863310000.0;

    a1 = 119.75 + 131.849 * t;

    a2 = 53.09 + 479264.29 * t;

    a3 = 313.45 + 481266.484 * t;

    e[0] = 1;
    e[1] = 1 - 0.002516 * t - 0.0000074 * t2;
    e[2] = e[1] * e[1];

    lprime = fixangr(dtr(lprime));
    d = fixangr(dtr(d));
    m = fixangr(dtr(m));
    mprime = fixangr(dtr(mprime));
    f = fixangr(dtr(f));
    a1 = fixangr(dtr(a1));
    a2 = fixangr(dtr(a2));
    a3 = fixangr(dtr(a3));

    sigmaL = sigmaB = sigmaR = 0;

    for (i = 0; i < NTERMS; i++) {
        ang = lrCoeff[i][0] * d + lrCoeff[i][1] * m +
              lrCoeff[i][2] * mprime + lrCoeff[i][3] * f;
        sigmaL += lTerms[i] * sin(ang) * e[abs(lrCoeff[i][1])];
        if (rTerms[i] != 0) {
            sigmaR += rTerms[i] * cos(ang) * e[abs(lrCoeff[i][1])];
        }

        ang = bCoeff[i][0] * d + bCoeff[i][1] * m +
              bCoeff[i][2] * mprime + bCoeff[i][3] * f;
        sigmaB += bTerms[i] * sin(ang) * e[abs(bCoeff[i][1])];
    }

    sigmaL += 3958.0 * sin(a1) + 1962.0 * sin(lprime - f) +
              318.0 * sin(a2);

    sigmaB += -2235.0 * sin(lprime) + 382.0 * sin(a3) + 175.0 * sin(a1 - f) +
              175.0 * sin(a1 + f) + 127.0 * sin(lprime - mprime) -
              115.0 * sin(lprime + mprime);

    *l = rtd(lprime) + sigmaL / 1000000.0;
    *b = sigmaB / 1000000.0;
    *r = 385000.56 + sigmaR / 1000.0;
}

/*  OBLIQEQ  --  Calculate the obliquity of the ecliptic for a given Julian
                 date.  This uses Laskar's tenth-degree polynomial fit
                 (J. Laskar, Astronomy and Astrophysics, Vol. 157, page 68 [1986])
                 which is accurate to within 0.01 arc second between AD 1000
                 and AD 3000, and within a few seconds of arc for +/-10000
                 years around AD 2000.  If we're outside the range in which
                 this fit is valid (deep time) we simply return the J2000 value
                 of the obliquity, which happens to be almost precisely the mean.  */

double obliqeq(const double jd)
{
#define Asec(x) ((x) / 3600.0)

    static const double oterms[10] = {
        Asec(-4680.93),
        Asec(   -1.55),
        Asec( 1999.25),
        Asec(  -51.38),
        Asec( -249.67),
        Asec(  -39.05),
        Asec(    7.12),
        Asec(   27.87),
        Asec(    5.79),
        Asec(    2.45)
    };

    double eps = 23 + (26 / 60.0) + (21.448 / 3600.0), u, v;
    int i;

    v = u = (jd - J2000) / (JulianCentury * 100);

    if (abs(u) < 1.0) {
        for (i = 0; i < 10; i++) {
            eps += oterms[i] * v;
            v *= u;
        }
    }
    return eps;
}

/* Periodic terms for nutation in longiude (delta \Psi) and
   obliquity (delta \Epsilon) as given in table 21.A of
   Meeus, "Astronomical Algorithms", first edition. */

#define NUTERMS 63

static const signed char nutArgMult[NUTERMS][5] = {
    {  0,  0,  0,  0,  1 },
    { -2,  0,  0,  2,  2 },
    {  0,  0,  0,  2,  2 },
    {  0,  0,  0,  0,  2 },
    {  0,  1,  0,  0,  0 },
    {  0,  0,  1,  0,  0 },
    { -2,  1,  0,  2,  2 },
    {  0,  0,  0,  2,  1 },
    {  0,  0,  1,  2,  2 },
    { -2, -1,  0,  2,  2 },
    { -2,  0,  1,  0,  0 },
    { -2,  0,  0,  2,  1 },
    {  0,  0, -1,  2,  2 },
    {  2,  0,  0,  0,  0 },
    {  0,  0,  1,  0,  1 },
    {  2,  0, -1,  2,  2 },
    {  0,  0, -1,  0,  1 },
    {  0,  0,  1,  2,  1 },
    { -2,  0,  2,  0,  0 },
    {  0,  0, -2,  2,  1 },
    {  2,  0,  0,  2,  2 },
    {  0,  0,  2,  2,  2 },
    {  0,  0,  2,  0,  0 },
    { -2,  0,  1,  2,  2 },
    {  0,  0,  0,  2,  0 },
    { -2,  0,  0,  2,  0 },
    {  0,  0, -1,  2,  1 },
    {  0,  2,  0,  0,  0 },
    {  2,  0, -1,  0,  1 },
    { -2,  2,  0,  2,  2 },
    {  0,  1,  0,  0,  1 },
    { -2,  0,  1,  0,  1 },
    {  0, -1,  0,  0,  1 },
    {  0,  0,  2, -2,  0 },
    {  2,  0, -1,  2,  1 },
    {  2,  0,  1,  2,  2 },
    {  0,  1,  0,  2,  2 },
    { -2,  1,  1,  0,  0 },
    {  0, -1,  0,  2,  2 },
    {  2,  0,  0,  2,  1 },
    {  2,  0,  1,  0,  0 },
    { -2,  0,  2,  2,  2 },
    { -2,  0,  1,  2,  1 },
    {  2,  0, -2,  0,  1 },
    {  2,  0,  0,  0,  1 },
    {  0, -1,  1,  0,  0 },
    { -2, -1,  0,  2,  1 },
    { -2,  0,  0,  0,  1 },
    {  0,  0,  2,  2,  1 },
    { -2,  0,  2,  0,  1 },
    { -2,  1,  0,  2,  1 },
    {  0,  0,  1, -2,  0 },
    { -1,  0,  1,  0,  0 },
    { -2,  1,  0,  0,  0 },
    {  1,  0,  0,  0,  0 },
    {  0,  0,  1,  2,  0 },
    { -1, -1,  1,  0,  0 },
    {  0,  1,  1,  0,  0 },
    {  0, -1,  1,  2,  2 },
    {  2, -1, -1,  2,  2 },
    {  0,  0, -2,  2,  2 },
    {  0,  0,  3,  2,  2 },
    {  2, -1,  0,  2,  2 },
};

static const long nutArgCoeff[NUTERMS][4] = {
    { -171996,   -1742,   92095,      89 },      /*  0,  0,  0,  0,  1 */
    {  -13187,     -16,    5736,     -31 },      /* -2,  0,  0,  2,  2 */
    {   -2274,      -2,     977,      -5 },      /*  0,  0,  0,  2,  2 */
    {    2062,       2,    -895,       5 },      /*  0,  0,  0,  0,  2 */
    {    1426,     -34,      54,      -1 },      /*  0,  1,  0,  0,  0 */
    {     712,       1,      -7,       0 },      /*  0,  0,  1,  0,  0 */
    {    -517,      12,     224,      -6 },      /* -2,  1,  0,  2,  2 */
    {    -386,      -4,     200,       0 },      /*  0,  0,  0,  2,  1 */
    {    -301,       0,     129,      -1 },      /*  0,  0,  1,  2,  2 */
    {     217,      -5,     -95,       3 },      /* -2, -1,  0,  2,  2 */
    {    -158,       0,       0,       0 },      /* -2,  0,  1,  0,  0 */
    {     129,       1,     -70,       0 },      /* -2,  0,  0,  2,  1 */
    {     123,       0,     -53,       0 },      /*  0,  0, -1,  2,  2 */
    {      63,       0,       0,       0 },      /*  2,  0,  0,  0,  0 */
    {      63,       1,     -33,       0 },      /*  0,  0,  1,  0,  1 */
    {     -59,       0,      26,       0 },      /*  2,  0, -1,  2,  2 */
    {     -58,      -1,      32,       0 },      /*  0,  0, -1,  0,  1 */
    {     -51,       0,      27,       0 },      /*  0,  0,  1,  2,  1 */
    {      48,       0,       0,       0 },      /* -2,  0,  2,  0,  0 */
    {      46,       0,     -24,       0 },      /*  0,  0, -2,  2,  1 */
    {     -38,       0,      16,       0 },      /*  2,  0,  0,  2,  2 */
    {     -31,       0,      13,       0 },      /*  0,  0,  2,  2,  2 */
    {      29,       0,       0,       0 },      /*  0,  0,  2,  0,  0 */
    {      29,       0,     -12,       0 },      /* -2,  0,  1,  2,  2 */
    {      26,       0,       0,       0 },      /*  0,  0,  0,  2,  0 */
    {     -22,       0,       0,       0 },      /* -2,  0,  0,  2,  0 */
    {      21,       0,     -10,       0 },      /*  0,  0, -1,  2,  1 */
    {      17,      -1,       0,       0 },      /*  0,  2,  0,  0,  0 */
    {      16,       0,      -8,       0 },      /*  2,  0, -1,  0,  1 */
    {     -16,       1,       7,       0 },      /* -2,  2,  0,  2,  2 */
    {     -15,       0,       9,       0 },      /*  0,  1,  0,  0,  1 */
    {     -13,       0,       7,       0 },      /* -2,  0,  1,  0,  1 */
    {     -12,       0,       6,       0 },      /*  0, -1,  0,  0,  1 */
    {      11,       0,       0,       0 },      /*  0,  0,  2, -2,  0 */
    {     -10,       0,       5,       0 },      /*  2,  0, -1,  2,  1 */
    {      -8,       0,       3,       0 },      /*  2,  0,  1,  2,  2 */
    {       7,       0,      -3,       0 },      /*  0,  1,  0,  2,  2 */
    {      -7,       0,       0,       0 },      /* -2,  1,  1,  0,  0 */
    {      -7,       0,       3,       0 },      /*  0, -1,  0,  2,  2 */
    {      -7,       0,       3,       0 },      /*  2,  0,  0,  2,  1 */
    {       6,       0,       0,       0 },      /*  2,  0,  1,  0,  0 */
    {       6,       0,      -3,       0 },      /* -2,  0,  2,  2,  2 */
    {       6,       0,      -3,       0 },      /* -2,  0,  1,  2,  1 */
    {      -6,       0,       3,       0 },      /*  2,  0, -2,  0,  1 */
    {      -6,       0,       3,       0 },      /*  2,  0,  0,  0,  1 */
    {       5,       0,       0,       0 },      /*  0, -1,  1,  0,  0 */
    {      -5,       0,       3,       0 },      /* -2, -1,  0,  2,  1 */
    {      -5,       0,       3,       0 },      /* -2,  0,  0,  0,  1 */
    {      -5,       0,       3,       0 },      /*  0,  0,  2,  2,  1 */
    {       4,       0,       0,       0 },      /* -2,  0,  2,  0,  1 */
    {       4,       0,       0,       0 },      /* -2,  1,  0,  2,  1 */
    {       4,       0,       0,       0 },      /*  0,  0,  1, -2,  0 */
    {      -4,       0,       0,       0 },      /* -1,  0,  1,  0,  0 */
    {      -4,       0,       0,       0 },      /* -2,  1,  0,  0,  0 */
    {      -4,       0,       0,       0 },      /*  1,  0,  0,  0,  0 */
    {       3,       0,       0,       0 },      /*  0,  0,  1,  2,  0 */
    {      -3,       0,       0,       0 },      /* -1, -1,  1,  0,  0 */
    {      -3,       0,       0,       0 },      /*  0,  1,  1,  0,  0 */
    {      -3,       0,       0,       0 },      /*  0, -1,  1,  2,  2 */
    {      -3,       0,       0,       0 },      /*  2, -1, -1,  2,  2 */
    {      -3,       0,       0,       0 },      /*  0,  0, -2,  2,  2 */
    {      -3,       0,       0,       0 },      /*  0,  0,  3,  2,  2 */
    {      -3,       0,       0,       0 },      /*  2, -1,  0,  2,  2 */
};

/*  NUTATION  --  Calculate the nutation in longitude, deltaPsi, and
                  obliquity, deltaEpsilon for a given Julian date
                  jd.  */

void nutation(const double jd, double *deltaPsi, double *deltaEpsilon)
{
    int i, j;
    double t = (jd - 2451545.0) / 36525.0, t2, t3, to10;
    double ta[5];
    double dp = 0, de = 0, ang;

    t3 = t * (t2 = t * t);

    /* Calculate angles.  The correspondence between the elements
       of our array and the terms cited in Meeus are:

       ta[0] = D  ta[0] = M  ta[2] = M'  ta[3] = F  ta[4] = \Omega

    */

    ta[0] = dtr(297.850363 + 445267.11148 * t - 0.0019142 * t2 +
                t3 / 189474.0);
    ta[1] = dtr(357.52772 + 35999.05034 * t - 0.0001603 * t2 -
                t3 / 300000.0);
    ta[2] = dtr(134.96298 + 477198.867398 * t + 0.0086972 * t2 +
                t3 / 56250.0);
    ta[3] = dtr(93.27191 + 483202.017538 * t - 0.0036825 * t2 +
                t3 / 327270);
    ta[4] = dtr(125.04452 - 1934.136261 * t + 0.0020708 * t2 +
                t3 / 450000.0);

    /* Range reduce the angles in case the sine and cosine functions
       don't do it as accurately or quickly. */

    for (i = 0; i < 5; i++) {
        ta[i] = fixangr(ta[i]);
    }

    to10 = t / 10.0;
    for (i = 0; i < NUTERMS; i++) {
        ang = 0;
        for (j = 0; j < 5; j++) {
            if (nutArgMult[i][j] != 0) {
                ang += nutArgMult[i][j] * ta[j];
            }
        }
        dp += (nutArgCoeff[i][0] + nutArgCoeff[i][1] * to10) * sin(ang);
        de += (nutArgCoeff[i][2] + nutArgCoeff[i][3] * to10) * cos(ang);
    }

    /* Return the result, converting from ten thousandths of arc
       seconds to radians in the process. */

    *deltaPsi = dtr(dp / (3600.0 * 10000.0));
    *deltaEpsilon = dtr(de / (3600.0 * 10000.0));
}

/*  ECLIPTOEQ  --  Convert celestial (ecliptical) longitude and latitude into
                   right ascension (in degrees) and declination.  We must supply
                   the time of the conversion in order to compensate correctly
                   for the varying obliquity of the ecliptic over time.  */

void ecliptoeq(const double jd, const double Lambda, const double Beta,
               double *Ra, double *Dec)
{
    double eps;

    /* Obliquity of the ecliptic. */

    eps = dtr(obliqeq(jd));

    *Ra = fixangle(rtd(atan2((cos(eps) * sin(dtr(Lambda)) -
                         (tan(dtr(Beta)) * sin(eps))), cos(dtr(Lambda)))));
    *Dec = rtd(asin((sin(eps) * sin(dtr(Lambda)) * cos(dtr(Beta))) +
                 (sin(dtr(Beta)) * cos(eps))));

}
