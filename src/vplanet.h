/*

                  Solar System Live  --  Definitions

                           by John Walker
                       http://www.fourmilab.ch/

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <limits.h>
#include <assert.h>
#include <stdint.h>

#include "bitmap.h"
#include "sglib.h"

#define FALSE   0
#define TRUE    1


#define TERMINC  100                  /* Circle segments for terminator */

#define PI 3.14159265358979323846  /* Assume not near black hole nor in
                                      Tennessee */
#define V       (void)

#define EOS     '\0'

/* Determine number of elements in an array */

#define ELEMENTS(array) (sizeof(array)/sizeof((array)[0]))

/*  Frequently used astronomical constants */

#define J2000               2451545.0       /* Julian day of J2000 epoch */
#define ModifiedJulian      2400000.5       /* Offset for modified Julian Day numbers */
#define JulianCentury       36525.0         /* Days in Julian century */
#define AstronomicalUnit    149597870.0     /* Astronomical unit in kilometres */
#define SunSMAX  (AstronomicalUnit * 1.000001018) /* Semi-major axis of Earth's orbit */
#define EarthRad            6378.14         /* Earth's equatorial radius, km (IAU 1976) */
#define LunatBase           2423436.0       /* Base date for E. W. Brown's numbered
                                               series of lunations (1923 January 16) */
#define SynMonth            29.53058868     /* Synodic month (mean time from new Moon to new Moon) */

#define max(x, y) (((x) > (y)) ? (x) : (y))               /* Maximum */
#define min(x, y) (((x) < (y)) ? (x) : (y))               /* Maximum */
#define sgn(x) (((x) < 0) ? -1 : ((x) > 0 ? 1 : 0))       /* Extract sign */
/* #define abs(x) ((x) < 0 ? (-(x)) : (x))  */            /* Absolute val */
#define fixangle(a) ((a) - 360.0 * (floor((a) / 360.0)))  /* Fix angle */
#define fixangr(a)  ((a) - (PI*2) * (floor((a) / (PI*2))))  /* Fix angle in radians */
#define dtr(x) ((x) * (PI / 180.0))                       /* Degree->Radian */
#define rtd(x) ((x) / (PI / 180.0))                       /* Radian->Degree */

struct planet {                 /* Planet information entry */
    double hlong;               /* FK5 Heliocentric longitude */
    double hlat;                /* FK5 Heliocentric latitude */
    double hrv;                 /* Heliocentric radius vector */
    double dhlong;              /* Dynamical Heliocentric longitude */
    double dhlat;               /* Dynamical Heliocentric latitude */
    double ra;                  /* Apparent right ascension */
    double dec;                 /* Apparent declination */
    double dist;                /* True distance from the Earth */
    double mag;                 /* Approximate magnitude */
    double lha;                 /* Local hour angle */
    double alt;                 /* Altitude above (-below) horizon */
    double az;                  /* Azimuth from South: West positive, East negative */
};
extern struct planet planet_info[11]; /* Calculated planetary information */

extern int aTracked;            /* Tracking an asteroid ? */

struct asteroid_info {          /* Asteroid information */
    int cometary;               /* Nonzero if cometary element format */
    char Name[80];              /* Name and number */
    double MagH;                /* IAU Commission 20 magnitude argument H */
    double MagG;                /* IAU Commission 20 magnitude argument G */
    double SemiMajorAU;         /* Semimajor axis in AU */
    double Eccentricity;        /* Eccentricity of orbit */
    double Inclination;         /* Inclination of orbit to ecliptic */
    double ArgP;                /* Argument of perihelion */
    double LANode;              /* Longitude of the ascending node */
    double mAnomaly;            /* Mean anomaly at the epoch */
    double Epoch;               /* Epoch of elements */
    double PeriDate;            /* Time of perihelion passage */
    double PeriAU;              /* Perihelion distance, AU */
};

extern struct asteroid_info ast_info;   /* Information about currently-tracked asteroid */

struct canned_orbit {
    unsigned short ulong;       /* Longitude */
    short ulat;                 /* Latitude */
    unsigned short urv;         /* Radius vector */
};
extern struct canned_orbit currentOrbits[];

/*  Select colour based on scheme:

        Cscheme(full_colour, black_on_white, white_on_black, night_vision)
*/

#define SCS_COLOUR          0             /* Normal colour display */
#define SCS_COLOUR_ON_WHITE 1             /* Colour on white background */
#define SCS_BLACK_ON_WHITE  2             /* Black print on white background */
#define SCS_WHITE_ON_BLACK  3             /* White print on black background */
#define SCS_NIGHT_VISION    4             /* Red to preserve night vision */

#define Cscheme(a, b, c, d, e) (skyColourScheme == SCS_COLOUR ? a : \
                              (skyColourScheme == SCS_COLOUR_ON_WHITE ? b : \
                                (skyColourScheme == SCS_BLACK_ON_WHITE ? c : \
                                  (skyColourScheme == SCS_WHITE_ON_BLACK ? d : e \
                                  ) \
                                )   \
                              ) \
                            )

/*  Default colour scheme for material in map.  */

#define CschemeD(a) Cscheme(a, a, cBlack, cWhite, cRed)

/* From vplanet.c */

extern char *progpath;          /* Program invocation path (argv[0]) */
extern double siteLat, siteLon; /* Observing site latitude and longitude */
#define maxOrbitalElements 10
extern char orbitalElements[maxOrbitalElements][132]; /* Elements of body being tracked */
extern char savedOrbitalElements[maxOrbitalElements][132]; /* Elements submitted with stateless request */
extern char *elementSpecification;  /* Orbital elements supplied by the user from WWW_elements */

/* From astrelem.c */

extern int element_rs(char *s);
extern void astrelem(const int preloaded);
extern void print_elements(FILE *fp);

/* From astro.c */

extern double ucttoj(const long year, const int mon, const int mday,
                     const int hour, const int min, const int sec);
extern double jtime(const struct tm *t);
extern void jyear(double td, long *yy, int *mm, int *dd);
extern void jhms(double j, int *h, int *m, int *s);
extern double gmst(const double jd);
extern void sunpos(const double jd, const int apparent,
                   double *ra, double *dec, double *rv, double *slong);
extern double phase(
                    const double  pdate,
                    double  *pphase, double  *mage, double  *dist,
                    double  *angdia, double  *sudist, double  *suangdia);
extern void highmoon(const double jd, double *l, double *b, double *r);
extern void nutation(const double jd, double *deltaPsi, double *deltaEpsilon);
extern void ecliptoeq(const double jd, const double Lambda, const double Beta,
                      double *Ra, double *Dec);
extern double obliqeq(const double jd);

/* From cometel.c */

extern void cometel(char *firstline);

/* From csvelements.c */

extern int csvElements(char *s);

/* From gifout.c */

extern int gifout(int cols, int rows, int rowlen, RGBQUAD *colmap, int colors,
                  int transparent, unsigned char *bitmap, FILE *outfile);

/* From htmlutil.c */

extern void write_html_content(FILE *fp, const char *s);
extern char *escape_html_content(char *s);

/* From jplelements.c */

extern int jplElements(const char *s);

/* From vsop87.c */

#define FULL_PLANET_INFO    0x8000          /* Calculate complete high-precision planet info */
extern void planets(const double jd, const int which);  /* Update planetary positions */

/* From packargs.c */

extern char *packargs(int argc, char *argv[]);
extern char *unpackargs(const char *oblock);

/* From uncgi.c */

extern void uncgi(void);

/* From planetp.c */

extern void updatePlanet(double jd, int normal, FILE *ofile, char *obSite);
extern void calcPlanet(double jd);
extern void buildPlanets(const double jd);

/* From pluto.c */
extern int plutoPrecise;        /* Precise position for Pluto required */
extern int pluto(double jd, double *l, double *b, double *r); /* Special calculation for Pluto */

/* From asteroid.c */

extern void selectAsteroid(char *aname, double aelem[8]);
extern void selectComet(char *aname, double aelem[8]);
extern void trackAsteroid(double jd, double *ra, double *dec, double *dist,
                          double *hlong, double *hlat, double *hrv, int quick);

/* From strlcpy.c */

size_t strlcpy(char *dst, const char *src, size_t siz);
