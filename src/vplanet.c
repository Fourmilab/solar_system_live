/*

              Solar System Live --  Main server program

                           by John Walker
                      http://www.fourmilab.ch/

*/

#include "vplanet.h"

/* #define TESTMODE */      /* Direct subsequent requests to Solar.t script */

/*  If CACHE_WARNING is defined, generated HTML will contain a warning
    to users not to link to the ephemeral cache files.  */
/* #define CACHE_WARNING */ /**/

#define TimeLimit   25                /* CPU time limit in seconds */

#ifdef TimeLimit
#define __USE_XOPEN_EXTENDED          /* Required to define sigset() with GCC library */
#include <signal.h>
#include <sys/time.h>
#include <sys/resource.h>
#endif

double siteLat = 47.0, siteLon = -7.0;
char orbitalElements[maxOrbitalElements][132];        /* Elements of body being tracked */
char savedOrbitalElements[maxOrbitalElements][132];   /* Elements submitted with dynamic request */
char *progpath;                       /* Program invocation path (argv[0]) */
static int directCall = FALSE;        /* Were we invoked directly, without uncgi first ? */
static int embeddedImage = FALSE;     /* Is this an embedded image specified with a "?di=" block ? */
char *elementSpecification = NULL;    /* Contents of orbit elements textbox in the request */

static int dto = 0;                   /* Date/time option */

static BITMAPFILEHEADER bfh;
static struct bmih {
    BITMAPINFOHEADER bmiHeader;
    RGBQUAD bmiColors[256];
} bmi, *tmbmap;
static int bmRowLen;
static unsigned char *bmBits = NULL, *tmbits;
static unsigned char ct[4];

static double jt, sunra, sundec, sunlong, gt;
static struct tm ctim;
static char viewdesc[132], dtutc[80];
static int transparent = -1, cBlue, cGreen, cBorder = -1, cWhite,
    cGrey, cLtGrey, cDkRed, cRed, cBlack;
static int dumpelem = FALSE;          /* Dump orbital elements ? */
static int stereoView = FALSE;
static int skyColourScheme = SCS_COLOUR;    /* Colour scheme */
static double stereoOcular = 6;       /* Ocular separation - cross, + wall */
#define stereoSepPixels 10
static int useimage = FALSE;          /* Use images for planets ? */
#ifdef SAVESET
static int bookmark = FALSE;          /* Save parameters for bookmark ? */
static int gotten = FALSE;            /* Was "get" method used ? */
#endif

#ifdef DBLOG
static FILE *db;
#define DB(x)   fprintf(db, "%s\n", x)
#else
#define DB(x)
#endif

extern char **environ;                  /* Environment variables */

#ifdef TimeLimit

/*  DINGALING  --  Handle SIGALRM when time runs out.  */

static void dingaling(int sig/*, int code, struct sigcontext *scp, char *addr*/)
{
    struct rusage u;

    if (sig == SIGALRM) {
        getrusage(RUSAGE_SELF, &u);
        if ((u.ru_utime.tv_sec + u.ru_stime.tv_sec) > TimeLimit) {
            printf("<h1><font color=\"#FF0000\">Time limit (%d seconds) exceeded.</font></h1>\n", TimeLimit);
            exit(0);
        }
        alarm(TimeLimit);             /* Wind the cat */
    }
}
#endif

/*  LOAD_BITMAP  --  Load a Microsoft Windows .BMP file into memory.  */

#define rB(x)   (x) = getc(fp)
#define rS(x)   dstat = fread(ct, 2, 1, fp); (x) = ((ct[1] << 8) | ct[0])
#define rL(x)   dstat = fread(ct, 4, 1, fp); (x) = (((long) ct[3] << 24) | ((long) ct[2] << 16) | (ct[1] << 8) | ct[0])

static void load_bitmap(const char *bmpfile)
{
    FILE *fp = fopen(bmpfile, "r");
    int dstat;

    if (bmBits != NULL) {
        free(bmBits);
        bmBits = NULL;
    }
    if (fp != NULL) {
        rS(bfh.bfType);
        rL(bfh.bfSize);
        rS(bfh.bfReserved1);
        rS(bfh.bfReserved2);
        rL(bfh.bfOffBits);

        rL(bmi.bmiHeader.biSize);
        rL(bmi.bmiHeader.biWidth);
        rL(bmi.bmiHeader.biHeight);
        rS(bmi.bmiHeader.biPlanes);
        rS(bmi.bmiHeader.biBitCount);
        rL(bmi.bmiHeader.biCompression);
        rL(bmi.bmiHeader.biSizeImage);
        rL(bmi.bmiHeader.biXPelsPerMeter);
        rL(bmi.bmiHeader.biYPelsPerMeter);
        rL(bmi.bmiHeader.biClrUsed);
        rL(bmi.bmiHeader.biClrImportant);

        assert(bmi.bmiHeader.biBitCount == 8);
        assert(bmi.bmiHeader.biCompression == BI_RGB);

        if (bmi.bmiHeader.biClrUsed == 0) {
            bmi.bmiHeader.biClrUsed = 1 << bmi.bmiHeader.biBitCount;
        }

        /* Read the colour palette table. */

        dstat = fread(&bmi.bmiColors[0], sizeof(RGBQUAD), bmi.bmiHeader.biClrUsed, fp);

        /* Compute the row length and read the bitmap. */

        bmRowLen = (bmi.bmiHeader.biWidth + 3) & ~3;

        fseek(fp, bfh.bfOffBits, SEEK_SET);
        bmBits = (unsigned char *) malloc(bmRowLen * bmi.bmiHeader.biHeight);
        dstat = fread(bmBits, bmRowLen * bmi.bmiHeader.biHeight, 1, fp);

        fclose(fp);
        dstat++;                /* Purely to get rid of "set but not used" warning */
    }
}

static void dump_ppm(const char *fname, const struct bmih *bi, const unsigned char *bits)
{
    FILE *fp;
    int i, j;

    if (fname[0] == '-') {
        fp = stdout;
    } else {
        fp = fopen(fname, "w");
    }
    assert(fp != NULL);
    fprintf(fp, "P6 %ld %ld %d\n",
        (long)( bi->bmiHeader.biWidth), (long) (bi->bmiHeader.biHeight), 255);
    for (i = 0; i < bi->bmiHeader.biHeight; i++) {
        unsigned const char *px = bits +
            (((bi->bmiHeader.biHeight - 1) - i) *
            ((bi->bmiHeader.biWidth + 3) & ~3));

        for (j = 0; j < bi->bmiHeader.biWidth; j++) {
            int p = *px++;

            putc(bi->bmiColors[p].rgbRed, fp);
            putc(bi->bmiColors[p].rgbGreen, fp);
            putc(bi->bmiColors[p].rgbBlue, fp);
        }
    }
    if (fp != stdout) {
        fclose(fp);
    }
}

/*  PARSEDEG  --  Parse latitude and longitude degree specifications.  */

static double parseDeg(const char *dec)
{
    if ((strchr(dec, 'd') != NULL) ||
        (strchr(dec, 'D') != NULL) ||
        (strchr(dec, ' ') != NULL) ||
        (strchr(dec, '°') != NULL)) {
        double dd = 0, mm = 0, ss = 0;
        char c1, c2;

        sscanf(dec, "%lf%c%lf%c%lf", &dd, &c1, &mm, &c2, &ss);
        return sgn(dd) * (abs(dd) + (mm / 60) + (ss / 3600));
    } else {
        return atof(dec);
    }
}

/*  LOCFILE  --  Expand local file name to full directory name from
                 which we're executing this program.  */

static char *locfile(const char *basename)
{
    static char pn[PATH_MAX];

    strcpy(pn, progpath);
    if (strrchr(pn, '/') != NULL) {
        *(strrchr(pn, '/') + 1) = EOS;
    } else {
        strcpy(pn, "./");
    }
    if (directCall || embeddedImage) {
        strcat(pn, "../cgi-executables/");
    }
    strcat(pn, basename);
    return pn;
}

/*  EDDEG  --  Edit degrees and minutes.  */

static char *eddeg(const double ds)
{
    static char buf[2][20];
    static int n = 0;
    long a = (long) ((abs(ds) + 0.00001) * 3600);
    int d, m, s;

    n = (n + 1) % 2;
    d = a / 3600;
    m = (a / 60) % 60;
    s = (a % 60);
    sprintf(buf[n], "%s%d\260", ds < 0 ? "-" : "", d);
    if (m > 0 || s > 0) {
        sprintf(buf[n] + strlen(buf[n]), "%d'", m);
        if (s > 0) {
            sprintf(buf[n] + strlen(buf[n]), "%d&quot;", s);
        }
    }
    return buf[n];
}

/*  EDLAT  --  Edit latitude.  */

static char *edlat(const double lat)
{
    static char slat[30];
    double ulat = fixangle(rtd(lat));

    if (ulat >= 180) {
        ulat = -(360 - ulat);
    }
    sprintf(slat, "%s%s", eddeg(abs(ulat)), ulat < 0 ? "S" : "N");
    return slat;
}

/*  EDLON  --  Edit longitude.  */

static char *edlon(const double lon)
{
    static char slon[30];
    double ulon = fixangle(rtd(lon));

    if (ulon >= 180) {
        ulon = -(360 - ulon);
    }
    sprintf(slon, "%s%s", eddeg(abs(ulon)), ulon < 0 ? "E" : "W");
    return slon;
}

/*  EDVPOS  --  Edit viewing position and return string.  */

static char *edvpos(const double tlat, const double tlon)
{
    static char vpos[80];

    sprintf(vpos, "%s %s",
        edlat(tlat), edlon(tlon));

    return vpos;
}

/*  PSRVECTOR  --  Draw a vector into the current image array.  */


#define PixA(x, y) ((tpixel + ((x) + (((DWORD) (tmhm1 - (y))) * (tmrowlen)))))
#define Opix(x, y) ((opixel + ((x) + (((DWORD) ((adiheight - 1) - (y))) * (rowlen)))))

static int rowlen, tmrowlen, adiheight;
static unsigned char *opixel;
static unsigned char *tpixel;

void psrvector(const unsigned int fx, const unsigned int fy,
               const unsigned int tx, const unsigned int ty, const unsigned int colour)
{
    unsigned int cfx, cfy, ctx, cty;
    int x, y, m, f, xinc, loopcnt, dx, dy;

#ifdef FlipY
    fy = (iarray[1] - 1) - fy;
    ty = (iarray[1] - 1) - ty;
#endif

    if (fy < ty) {
        cfx = fx;
        cfy = fy;
        ctx = tx;
        cty = ty;
    } else {
        cfx = tx;
        cfy = ty;
        ctx = fx;
        cty = fy;
    }
    dy = cty - cfy;
    dx = ctx - cfx;

    if (dx < 0) {
        xinc = -1;
        dx = -dx;
    } else  {
        xinc = 1;
    }

    x = cfx;
    y = cfy;

    if (dx > dy) {
        f = dy - dx / 2;
        m = dx - dy;
        loopcnt = dx + 1;

        while (loopcnt--) {
            *Opix(x, y) = colour;
            x += xinc;
            if (f > 0) {
                f -= m;
                y++;
            } else {
                f += dy;
            }
        }
    } else {
        f = -(dy / 2 - dx);
        m = dy - dx;
        loopcnt = dy + 1;

        while (loopcnt--) {
            *Opix(x, y) = colour;
            y++;
            if (f > 0) {
                x += xinc;
                f -= m;
            } else {
                f += dx;
            }
        }
    }
}

static double orrLat = 90;              /* Orrery latitude */
static double orrLong = 0;              /* Orrery longitude */
static int orrInner = TRUE;             /* Orrery inner solar system only ? */
static int orrScale = 2;                /* Distance scale */

static double earthlong, earthrv;       /* Heliocentric coordinates of the Earth */
static double minperh = 0.28;           /* Minimum perihelion (usually Mercury, but possible ast./comet) */

/* The fastest-changing aspect of Solar system orbits is the perihelion shift in
   the orbit of Mercury, and that amounts to only about 1.5 degrees per century (
   0.559974 +/= 0.000041 arc seconds per century, to be precise).  So, we only
   bother to recalculate the orbits when they're a century out of date. */

/* Sidereal rotation period of planets in days */

static double rotationPeriod[] = {
    87.969,
    224.701,
    365.256,
    686.980,
    4332.589,
    10759.22,
    30685.4,
    60189.0,
    90465.0,
    0.0                                 /* Asteroid period calculated from orbital elements */
};

/*  SETPROJECT  --  Define projection matrix for a given viewpoint.  */

static void setproject(const double orrlat, const double orrlong, const double stereo)
{
    matident(cT);
    if (stereo != 0.0) {
        rot(dtr(stereo), Y);                 /* Stereo ocular offset */
    }
    rot(dtr(180.0 - (orrlat - 90.0)), X);    /* Tilt to latitide */
    rot(dtr(90.0 - orrlong), Z);             /* Spin to longitude */
}

/*  PROJECT  --  Project a longitude, latitude, and radius vector into the
                 display space.  */

static void project(const double along, const double alat, const double arv, const int cr,
                    const int plno, const int nplan, const double diau,
                    int *u, int *v, int *z)
{
    double lrv, acoslat, vx, vy, vz;
    vector p, uvw;
#define MercOffset  16

    if (orrScale == 0) {
        /* Log projection */
        lrv = log10(1.0 + 9.0 * ((arv - minperh) / (diau)));
        if (lrv < 0) {
            lrv = 0;
        } else {
            lrv = (int) (MercOffset + (cr - MercOffset) * lrv);
        }
    } else if (orrScale == 1) {
        /* Real projection */
        lrv = (int) (MercOffset + (cr - MercOffset) * (arv / diau));
    } else {/* orrScale == 2 */
        /* Equal projection */
        lrv = (int) (MercOffset + (cr - MercOffset) * ((plno == 10 ? 4.5 : ((double) plno)) / ((double) nplan)));
    }

    /* Get rectangular co-ordinates of point. */
    acoslat = abs(cos(dtr(alat)));
    vecget(p, lrv * sin(dtr(along)) * acoslat,
              lrv * cos(dtr(along)) * acoslat,
              -lrv * sin(dtr(alat)));
    /* Transform by viewing matrix */
    vecxmat(uvw, p, cT);
    /* Store projected co-ordinates. */
    vecput(&vx, &vy, &vz, uvw);

    /* Return screen co-ordinates (relative to the centre of
       the screen.  The Z co-ordinate represents the distance
       of the point behind (if negative) or above (if positive)
       the screen.  This is used to sort the order in which the
       planet icons are drawn to that any obscurations are done
       in the correct order. */

    *u = (int) vx;
    *v = (int) vy;
    *z = (int) (vz * 2000);
}

/*  ORRERY  --  Generate Orrery view.  */

static void orrery(const int isize, const int tmwidth, const int tmheight,
                   const int diwidth, const int diheight, const double faketime,
                   const int stereo)
{
    int i, x, y, ix, iy, o = -1, t, tmhm1, cr, cpx = 0, cpy = 0,
        nplan = aTracked ? 10 : (orrInner ? 4 : 9),
        siwidth = diwidth, sbias, pass,
        planX[11], planY[11], planZ[11], planord[11];

    double diau = orrInner ? (aTracked ?
       max(1.7, min(5.5, ((ast_info.SemiMajorAU * (1.0 + ast_info.Eccentricity)) * 1.1)))
            : 1.7) :
           (aTracked ?
                max(50, min(150, ((ast_info.SemiMajorAU * (1.0 + ast_info.Eccentricity)) * 1.1)))
            : 50),
           eplan = aTracked ? (orrInner ? 5 : 9) : (orrInner ? 4 : 9);
    double dialim = diau * 2;
    BITMAPINFOHEADER *bh, *obmap;
    WORD ncol, offbits;
    double orrlat = orrLat, orrlong = orrLong;
#ifdef SaveOrbits
    FILE *fp = fopen("/tmp/orbits.c", "w");
#endif

/*
if (aTracked) {
fprintf(stderr, "Diau = %.2f, acomp = %.2f\n", diau,
    (ast_info.SemiMajorAU * (1.0 + ast_info.Eccentricity)) * 1.1);
}
*/
    if (stereo) {
        siwidth += diwidth + stereoSepPixels;
    }
    adiheight = diheight;
    bh = &bmi.bmiHeader;
    ncol = (WORD) (bh->biClrUsed == 0 ? (1L << bh->biBitCount) : bh->biClrUsed);
    offbits = (WORD) (bh->biSize + ncol * sizeof(RGBQUAD));
    tpixel = bmBits;
    tmrowlen = ((tmwidth + (sizeof(LONG) - 1)) / sizeof(LONG)) * sizeof(LONG);
    tmhm1 = tmheight - 1;

    if (tmbmap != NULL) {
        free(tmbmap);
        tmbmap = NULL;
    }
    rowlen = ((siwidth + (sizeof(LONG) - 1)) / sizeof(LONG)) * sizeof(LONG);
    obmap = (BITMAPINFOHEADER *) malloc(offbits + diheight * rowlen);
    memset((char *) obmap, 0, offbits + diheight * rowlen);
    memcpy((char *) obmap, (char *) bh, offbits);
    obmap->biWidth = siwidth;
    obmap->biHeight = diheight;
    obmap->biSizeImage = ((DWORD) diheight) * rowlen;
    opixel = ((unsigned char *) (obmap)) + offbits;
    tmbmap = (struct bmih *) obmap;
    tmbits = (unsigned char *) opixel;

    if ((skyColourScheme == SCS_BLACK_ON_WHITE) ||
        (skyColourScheme == SCS_COLOUR_ON_WHITE)) {
        for (y = 0; y < diheight; y++) {
            memset(Opix(0, y), cWhite, stereo ? siwidth : diwidth);
        }
    }

    /* If we're generating a stereo pair, draw the border between
       the two images.  This provides a reference which makes them
       easier to fuse. */

    if (stereo) {
        for (y = 0; y < diheight; y++) {
            for (x = diwidth; x < diwidth + stereoSepPixels; x++) {
                *Opix(x, y) = cBorder;
            }
        }
    }

    for (pass = 0; pass < (stereo ? 2 : 1); pass++) {
        sbias = (pass == 0) ? 0 : (diwidth + stereoSepPixels);

        setproject(orrlat, orrlong,
                   stereo ? (stereoOcular * (pass == 0 ? 1 : -1)) : 0.0);

        {
            double plong, plat, prv, sjd, timestep;
            int cx, cy, vx, vy, vz, fvx = 0, fvy = 0, ftf = TRUE, nsteps, colour = 0;
            char spi[sizeof(planet_info)];

#define inLim(x)        ((x) >= 0 && (x) < diheight)
#define MoveTo(x, y)    cpx = (x), cpy = (y)
#define LineTo(x, y)    if (inLim(cpx) && inLim(cpy) && inLim(x) && inLim(y)) { \
                            psrvector(sbias + cpx, cpy, \
                                      sbias + (x), (y), colour); \
                        } cpx = (x), cpy = (y)
#define nSteps      64

            /* Constrain use of canned orbits to decent interval around J2000.0 */

            if (abs((faketime - J2000) / JulianCentury) <= 5) {
                o = 0;
            }

            cx = diwidth / 2;
            cy = diheight / 2;
            cr = (min(diwidth, diheight) / 2) - (MercOffset / 2);

            for (i = 1; i <= min(nplan, 9); i++) {
                ftf = TRUE;

#ifdef SaveOrbits
                fprintf(fp, "\n    /* Planet %d */ \n\n", i);
#endif

                /* Test for zero radius vector which indicates planetary theory not
                   known for this date (Pluto only).  If that's the case, skip
                   drawing the planet. */

                if ((orrInner && (i > 4)) ||
                    ((i != 3) && (planet_info[i].hrv == 0.0))) {
                    continue;
                }
                plutoPrecise = FALSE;           /* Allow suspect values to close Pluto orbit */
                memcpy(spi, (char *) planet_info, sizeof spi);
                nsteps = nSteps;

                /* Draw the orbital path for the planet by stepping forward in
                   time one revolution period.  (N.b. for the planets whose
                   orbits are almost perfectly circular, it would be *enormously*
                   faster to just draw circles of the mean radius.  But for the
                   moment we do it the hard way for all of them.  This may be
                   slow, but it's *right*.)  */

                timestep = rotationPeriod[i - 1] / nsteps;
                for (sjd = faketime; sjd <= faketime + rotationPeriod[i - 1]; sjd += timestep) {

                    if (o >= 0 && i <= 9) {
                        plong = (((double) currentOrbits[o].ulong) * 360.0) / 65535.0;
                        plat = (((double) currentOrbits[o].ulat) * 180.0) / 32767.0;
                        prv = (((double) currentOrbits[o++].urv) * 50.0) / 65535.0;
                    } else {
                        if (i == 3) {
                            double sunra, sundec;

                            sunpos(sjd, TRUE, &sunra, &sundec, &prv, &plong);
                            plong += 180;
                            plat = 0;
                        } else {
                            planets(sjd, 1 << i);
                            plong = planet_info[i].hlong;
                            plat = planet_info[i].hlat;
                            prv = planet_info[i].hrv;
                        }
                    }
#ifdef SaveOrbits
                    if (fp != NULL && i <= 9) {
                        unsigned int ulong, urv;
                        int ulat;

                        ulong = (unsigned int) ((fixangle(plong) * 65535.0) / 360.0);
                        ulat = (int) ((plat * 32767.0) / 180.0);
                        urv = (unsigned int) ((prv * 65535.0) / 50.0);
                        fprintf(fp, "    { %6u, %6d, %6u },\n", ulong, ulat, urv);
                    }
#endif
                    if (prv > 0) {
                        project(plong, plat, prv, cr, i, eplan, diau,
                                &vx, &vy, &vz);
                        vx += cx;
                        vy = cy - vy;
                        if (ftf) {
                            MoveTo(fvx = vx, fvy = vy);
                            ftf = FALSE;
                        } else {
                            if (plat < 0) {
                                colour = Cscheme(cGreen, cGreen, cLtGrey, cGrey, cDkRed);
                            } else {
                                colour = Cscheme(cBlue, cBlue, cBlack, cWhite, cRed);
                            }
                            LineTo(vx, vy);
                        }
                    }
                }
                if (!ftf) {
                    /* Draw line to guarantee orbit is closed */
                    LineTo(fvx, fvy);
                }
                plutoPrecise = TRUE;
                memcpy((char *) planet_info, spi, sizeof spi);
            }

            /* If an asteroid or comet is being tracked, we have to handle it
               specially.  First of all, we can never use the canned planetary
               orbits to avoid calculation.  Second, unlike planetary orbits,
               asteroidal and cometary orbits can have large eccentricity,
               which means that we can't use a naive approximation to the
               number of steps needed to draw the orbit smoothly, but must
               instead take into account the orbital velocity at various
               locations on the orbit and adjust the time step accordingly.
               Finally, the orbit may not be closed either because it is
               genuinely parabolic or hyperbolic, or merely such an eccentric
               ellipse that the curve is truncated at the edge of the display. */

            if (aTracked) {
                double eSmA, lprv, aprv = 0;
                int arm, avx = 0, avy = 0, lvx = 0, lvy = 0;
                /* int aclosed = FALSE; */

                i = 10;
                memcpy(spi, (char *) planet_info, sizeof spi);

                /* Because the orbit may not be closed, we start at the
                   perihelion date and plot into the future until we
                   either arrive at the aphelion or we run into the
                   cutoff based on the scale of the display.  Then we go
                   back to the perihelion and plot backward until we
                   encounter the equivalent event on the other branch. */

                eSmA = ast_info.Eccentricity < 1.0 ? ast_info.SemiMajorAU :
                        40;

                for (arm = -1; arm <= 1; arm += 2) {
                    ftf = TRUE;                 /* No vector on first point */
                    sjd = ast_info.PeriDate;
                    lprv = 0;

                    while (TRUE) {
                        double ov = 0.5;

                        planets(sjd, 1 << 10);
                        plong = planet_info[i].hlong;
                        plat = planet_info[i].hlat;
                        prv = planet_info[i].hrv;

                        if (prv > 0) {
                            lvx = vx;
                            lvy = vy;
                            project(plong, plat, prv, cr, 10, eplan, diau,
                                    &vx, &vy, &vz);
                            vx += cx;
                            vy = cy - vy;
                            if (ftf) {
                                MoveTo(vx, vy);
                                ftf = FALSE;
                            } else {
                                if (plat < 0) {
                                    colour = cGreen;
                                } else {
                                    colour = cBlue;
                                }
                                LineTo(vx, vy);
                            }
                        }
#ifdef AST_DEBUG
else {
  fprintf(stderr, "Prv = %.4f  sjd = %.4f\n", prv, sjd);
  break;
}

if (!finite(prv)) {
  fprintf(stderr, "Not finite Prv = %.4f  sjd = %.4f\n", prv, sjd);
  break;
}
#endif

                        if ((prv == 0) || (prv > dialim) || (prv < lprv)) {
                            if (prv > 0) {
                                if (arm > 0) {
                                    /* If orbit closed, connect ends of orbit */
                                    if (aprv <= lprv && prv <= lprv) {
                                        LineTo(avx, avy);
                                        LineTo(avx, avy);           /* Set last pixel */
                                    }
                                } else {
                                    /* aclosed = prv < lprv; */
                                    aprv = prv;
                                    avx = lvx;
                                    avy = lvy;
                                }
                            }
                            break;
                        }
                        lprv = prv;

                        /* Orbital velocity.  Note that since we plot both branches of
                           the orbit outward from the perihelion, we always overestimate
                           the average velocity for the segment, guaranteeing that the
                           curve will be smooth at the expense of a few extra segments.
                           If the orbit of the body extends into the outer solar system,
                           the calculation of its orbital velocity may result in taking
                           the square root of a negative number due to the hokey way we
                           we estimate it.  If this situation arises, we simply use the
                           last valid orbital velocity computed (which will result in
                           more segments than needed) for the remaining part of the orbit. */

#ifdef AST_DEBUG
if (prv == 0) {
  fprintf(stderr, "Doooh!  OV calc Prv = %.4f  sjd = %.4f\n", prv, sjd);
  break;
}
#endif

                        if ((1 / (2 * eSmA)) < (1 / prv)) {
/*  fprintf(stderr, "Limiting OV, was %.4f at prv %.4f\n", ov, prv); */
                            ov = 42.1219 * sqrt((1 / prv) - (1 / (2 * eSmA)));
                        }

                        /* Time step to move 0.1 AU in inner solar system, correspondingly
                           more in full system. */

                        sjd += arm * (orrInner ? 1 : (prv < 1.0 ? 1 : max(1, (eSmA / 5)))) *
                                (((AstronomicalUnit / ov) / (24.0 * 60 * 60)) / 10);
#ifdef AST_DEBUG
if (!finite(sjd)) {
  fprintf(stderr, "Doooh!  SJD increment arm = %d  prv = %.4f  ov = %.4f\n", arm, prv, ov);
  break;
}
#endif
                    }
                }
                memcpy((char *) planet_info, spi, sizeof spi);
            }
        }

        /* Plug in well-known values for the heliocentric
           position of the Sun.  We only do this for the
           first image of a stereo pair, since the existing
           value will be re-used for the second. */

        if (pass == 0) {
            planord[0] = planX[0] = planY[0] = planZ[0] = 0;
        }

        /* Calculate transformed plot positions for the planet
           icons. */

        for (i = 1; i <= 10; i++) {
            int o;

            /* On the first image of a stereo pair, fill the table
               of planets in linear order.  For the second image,
               where we want to re-use the sort order of the
               original image (see comment below for why), find the
               planet in the order table and use that slot.  Note
               that we have to include a bail-out in case this
               planet wasn't rendered in the selected view. */

            if (pass == 0) {
                o = i;
            } else {
                for (o = 0; o <= 10; o++) {
                    if (planord[o] == i) {
                        break;
                    }
                }
                if (o > 10) {
                    continue;
                }
            }
            planord[o] = i;
            if (i <= nplan) {
                double plong, plat, prv;

                if (i == 3) {
                    plong = earthlong;
                    plat = 0;
                    prv = earthrv;
                } else {
                    plong = planet_info[i].hlong;
                    plat = planet_info[i].hlat;
                    prv = planet_info[i].hrv;
                }
                if ((orrInner && (i > 4 && i <= 9)) || (prv == 0.0)) {
                    planord[o] = planX[o] = planY[o] = planZ[o] = -1;
                    continue;           /* Planetary theory for Pluto invalid at this time */
                }
                project(plong, plat, prv, cr, i, eplan, diau,
                        &planX[o], &planY[o], &planZ[o]);

            } else {
                planord[o] = planX[o] = planY[o] = planZ[o] = -1;
            }
        }

        /* Sort planets by projected distance.  We do this only
           for the first image of a stereo pair.  The second stereo
           image is always drawn with the same order of planets
           as the first to avoid peekaboo planets when the
           slight difference in viewpoint changes the Z-depth. */

        if (pass == 0) {
            for (i = 0; i <= 10; i++) {
                int j;

                for (j = i + 1; j <= 10; j++) {
                    if (planZ[j] < planZ[i]) {
                        int t;

#define Swop(p) t = p[i]; p[i] = p[j]; p[j] = t
                        Swop(planX);
                        Swop(planY);
                        Swop(planZ);
                        Swop(planord);
#undef Swop
                    }
                }
            }
        }

        /* Plot the planet icons in order of projected distance
           from the view: most distant to nearest.  This performs
           a "painter's algorithm" rendering which ensures all
           obscurations are done in the right order. */

        for (i = 0; i <= 10; i++) {
            if (planord[i] >= 0) {
                ix = diwidth / 2 + planX[i];
                iy = diheight / 2 - planY[i];
                t = planord[i] * 32;  /* Select proper icon */
                if (planord[i] == 10 && ast_info.cometary) {
                    t += 32;
                }
                for (x = 0; x < 32; x++) {
                    int ax = (ix - 16) + x;

                    if (ax >= 0 && ax < diwidth) {
                        for (y = 0; y < 32; y++) {
                            int ay = (iy - 16) + y;

                            if (ay >= 0 && ay < diheight) {
                                unsigned char p = *PixA(x + t, y);

                                if (p != transparent) {
                                    *Opix(sbias + ax, ay) = p;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
#ifdef SaveOrbits
    fclose(fp);
#endif
}

/*  WRITE_CACHE_IMAGE_WARNING  --  Include warning in HTML source to deter
                                   users from linking to images in the cache.  */

#ifdef CACHE_WARNING
static void write_cache_image_warning(FILE *fp)
{
#define W(x)    fprintf(fp, "%s\n", x)
    W("<!--");
    W("                      WARNING!!!");
    W("");
    W("You may be looking at the source for this Solar System Live result");
    W("document in order to figure out where the image file comes from");
    W("so you can link to it from your site.  DON'T DO THAT!");
    W("");
    W("The image files generated by Solar System Live are stored in");
    W("a temporary \"cache\" directory, and are kept only long enough");
    W("to guarantee that users will be able to download them to their");
    W("browsers.  After 10 minutes or so, images are automatically deleted.");
    W("So, if you link to one of these images, it may work if you try it");
    W("right away, but the link is sure to break before long, leaving an");
    W("ugly broken image on your page.");
    W("");
    W("Rather than linking to the image, load it in your browser and then");
    W("save a local copy on your site and reference it from there.");
    W("If you want to link to an image which will be updated every time");
    W("users request it, create a custom request with:");
    W("  http://www.fourmilab.ch/solar/custom.html");
    W("then copy the URL it generates into a link on your page.");
    W("");
    W("Thank you for your understanding and cooperation.");
    W("-->");
#undef W
}
#endif

/*  WRITEPOST  --  Write balance of HTML file, including controls
                   which permit regenerating the view.  */

static void writepost(FILE *fp, const double vlat, const double vlon, const double alt,
    const char *bmpname, const int innersystem, const int orbittype, const int tmsize, const double jt)
{
    int i, j;
    double lat, lon;
    static const char checked[] = " checked=\"checked\"";

    lon = fixangle(rtd(vlon));
    if (lon >= 180) {
        lon = -(360 - lon);
    }
    lat = fixangle(rtd(vlat));
    if (lat >= 180) {
        lat = -(360 - lat);
    }

#define W(x)    fprintf(fp, "%s\n", x)
#define L(x, y)    fprintf(fp, "<a href=\"/solar/help/%s.html\"><b>%s:</b></a>", y, x);
#ifdef TESTMODE
    /* We always use the "get" method in TESTMODE so that the arguments are visible
       and can be edited for debugging. */
    W("<form method=\"get\" name=\"request\" action=\"/cgi-bin/Solar.t\">");
#else
#ifdef SAVESET
    fprintf(fp, "<form method=\"%s\" name=\"request\" action=\"/cgi-bin/Solar\">\n", bookmark ? "get" : "post");
    if (bookmark) {
        W("<input name=\"got\" type=\"hidden\" value=\"-xh\" />\n");
    }
#else
    fprintf(fp, "<form method=\"post\" name=\"request\" action=\"/cgi-bin/Solar\">\n");
#endif
#endif
    W("<center>");
    W("<p />");

                              L("Time", "timedate");

    fprintf(fp, "<input type=\"radio\" name=\"date\" onclick=\"0\" value=\"0\"%s />&nbsp;<a href=\"/solar/help/timedate.html#Now\">Now</a>\n",
        dto == 0 ? checked : "");
    fprintf(fp, "<input type=\"radio\" name=\"date\" onclick=\"0\" value=\"1\"%s />&nbsp;<a href=\"/solar/help/timedate.html#UTC\">UTC:</a>&nbsp;",
        dto == 1 ? checked : "");
    fprintf(fp, "<input type=\"text\" name=\"utc\" value=\"%s\" size=\"20\" onchange=\"document.request.date[1].checked=true;\" />\n",
        dtutc);
    fprintf(fp, "<input type=\"radio\" name=\"date\" onclick=\"0\" value=\"2\"%s />&nbsp;<a href=\"/solar/help/timedate.html#Julian\">Julian:</a>&nbsp;",
        dto == 2 ? checked : "");
    fprintf(fp, "<input type=\"text\" name=\"jd\" value=\"%.5f\" size=\"15\" onchange=\"document.request.date[2].checked=true;\" />\n",
        jt);
    W("<br />");

#ifdef SAVESET
    fprintf(fp, "  <input type=\"checkbox\" name=\"bmark\" value=\"-xs\"%s /> ",
        (bookmark && !gotten) ? "checked" : "");
    fprintf(fp, "<a href=\"/solar/help/saveset.html\"><b>Save settings</b></a>");
#endif
    W( "<input type=\"submit\" value=\"Update\" />");
    W("<br />");

                              L("Show", "icons");
    fprintf(fp, "<label><input type=\"radio\" name=\"img\" value=\"-k0\"%s />&nbsp;Icons</label>\n",
        !useimage ? checked : "");
    fprintf(fp, "<label><input type=\"radio\" name=\"img\" value=\"-k1\"%s />&nbsp;Images</label>\n",
        useimage ? checked : "");
    W("<br />");

                            L("Display", "in-out");
    fprintf(fp, "<label><input type=\"radio\" name=\"sys\" value=\"-Sf\"%s />&nbsp;Full&nbsp;system</label>\n",
        !innersystem ? checked : "");
    fprintf(fp, "<label><input type=\"radio\" name=\"sys\" value=\"-Si\"%s />&nbsp;Inner&nbsp;system</label>\n",
        innersystem ? checked : "");

                              L("Size", "size");
    fprintf(fp, " <input type=\"text\" name=\"imgsize\" value=\"%d\" size=\"6\" />\n",
        tmsize);

    fprintf(fp, "<input type=\"checkbox\" name=\"stereo\" value=\"-v\"%s /> \n",
        stereoView ? checked : "");

                             L("Stereo", "stereo");
    fprintf(fp, "<label><input type=\"radio\" name=\"eyes\" value=\"0\"%s />&nbsp;Cross</label>\n",
        stereoOcular >= 0 ? checked : "");
    fprintf(fp, "<label><input type=\"radio\" name=\"eyes\" value=\"1\"%s />&nbsp;Wall</label>\n",
        stereoOcular < 0 ? checked : "");

    W("<br />");
                             L("Orbits", "orbits");
    fprintf(fp, "<label><input type=\"radio\" name=\"orb\" value=\"-b1\"%s />&nbsp;Real</label>\n",
        orbittype == 1 ? checked : "");
    fprintf(fp, "<label><input type=\"radio\" name=\"orb\" value=\"-b0\"%s />&nbsp;Logarithmic</label>\n",
        orbittype == 0 ? checked : "");
    fprintf(fp, "<label><input type=\"radio\" name=\"orb\" value=\"-b2\"%s />&nbsp;Equal</label>\n",
        orbittype == 2 ? checked : "");

    W("<br />");

                         L("Observing site", "site");
    W(" Lat.");
    fprintf(fp, "<input type=\"text\" name=\"lat\" value=\"%s\" size=\"10\" />&nbsp;",
        eddeg(abs(lat)));
    fprintf(fp, "<label><input type=\"radio\" name=\"ns\" value=\"North\"%s />&nbsp;N</label>&nbsp;",
        lat > 0 ? checked : "");
    fprintf(fp, "<label><input type=\"radio\" name=\"ns\" value=\"South\"%s />&nbsp;S</label>\n",
        lat < 0 ? checked : "");

    W("Long.");
    fprintf(fp, "<input type=\"text\" name=\"lon\" value=\"%s\" size=\"10\" />&nbsp;",
        eddeg(abs(lon)));
    fprintf(fp, "<label><input type=\"radio\" name=\"ew\" value=\"East\"%s />&nbsp;E</label>&nbsp;",
        lon < 0 ? checked : "");
    fprintf(fp, "<label><input type=\"radio\" name=\"ew\" value=\"West\"%s />&nbsp;W</label>\n",
        lon > 0 ? checked : "");
    W("<br />");

                          L("Heliocentric", "heliocentric");
    fprintf(fp, " Lat. <input type=\"text\" name=\"hlat\" value=\"%s\" size=\"10\" />&nbsp;",
        eddeg(abs(orrLat)));
    fprintf(fp, "<label><input type=\"radio\" name=\"hns\" value=\"North\"%s />&nbsp;N</label>&nbsp;",
        orrLat > 0 ? checked : "");
    fprintf(fp, "<label><input type=\"radio\" name=\"hns\" value=\"South\"%s />&nbsp;S</label>\n",
        orrLat < 0 ? checked : "");

    W("Long.");
    fprintf(fp, "<input type=\"text\" name=\"hlon\" value=\"%s\" size=\"10\" />\n",
        eddeg(orrLong));
    W("<br />");

                            L("Colour scheme", "scheme");
#define Sel(x)  ((x) ? " selected=\"selected\"" : "")
    fprintf(fp, "\
&nbsp;<select name=\"scheme\" size=\"1\">\n\
<option value=\"0\"%s>Colour</option>\n\
<option value=\"1\"%s>Colour on white background</option>\n\
<option value=\"2\"%s>Black on white background</option>\n\
<option value=\"3\"%s>White on black background</option>\n\
<option value=\"4\"%s>Night vision (red)</option>\n\
</select>\n", Sel(skyColourScheme == 0), Sel(skyColourScheme == 1),
              Sel(skyColourScheme == 2), Sel(skyColourScheme == 3),
              Sel(skyColourScheme == 4));



    W("</center>");

    /* Append ephemeris after the control panel. */

    W("<p />");
                   L("Ephemeris", "ephemeris");
    updatePlanet(jt, TRUE, stdout, edvpos(vlat, vlon));

    /* Include elements for object being tracked, if any */

    if (aTracked) {
        L("Orbital elements for asteroid or comet", "elements");
        for (j = 7; j >= 0; j--) {
            if (strlen(orbitalElements[j]) > 0) {
                break;
            }
        }
        j = max(1, j + 1);
    } else {
        W("<b>To track an asteroid or comet, paste </b>");
        L("orbital elements below", "elements");
        j = 1;
    }
    fprintf(fp, "  <label><input type=\"checkbox\" name=\"edump\" value=\"-xe\"%s /> Echo elements</label>\n",
        dumpelem ? checked : "");
    W("<br />");
    fprintf(fp, "<textarea name=\"elements\" rows=\"%d\" cols=\"70\">", j);
    for (i = 0; i < j; i++) {
        write_html_content(fp, orbitalElements[i]);
        putc('\n', fp);
    }
    W("</textarea><br />");
    W("</form>");

    if (dumpelem && aTracked) {
        print_elements(stdout);
    }
#undef W
#undef L
}

/* Main program. */

int main(int argc, char **argv)
{
    double vlat, vlon, valt,
           tlat = 0, tlon = 0, talt = 35785;
    int i, f = 0,
        tmsize = 320, wppm = FALSE, trackelem = FALSE, html = FALSE,
        innersystem = FALSE, orbittype = 0, endopts = FALSE,
        dynimg = FALSE, stateless = FALSE;
    char *di, *cp, opt, *cachep;
    char cacheName[PATH_MAX];
    char tbuf[132];
    long cctime;
    FILE *ofile = stdout;
    static char *imagesource = "earth-topo.bmp";
    static char *mname[12] = {
        "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct",
        "Nov", "Dec" };

#ifdef DBLOG
    db = fopen("/tmp/dblog", "w");
#endif
    progpath = argv[0];               /* Save path used to call us */

    /*  Clear orbitalElements array.  */
    for (i = 0; i < 8; i++) {
        orbitalElements[i][0] = 0;
    }

    /*  If the WWW_di environment variable is defined, we've been invoked
        from a "stateless" mode document to generate a dynamic image with the
        arguments packed into the argument block.  Unpack the block into a
        synthesised argc/argv and transcribe the orbital elements (if any) to
        the orbitalElements array.  Existence of a WWW_di block overrides any
        arguments specified on the command line.  */

    di = getenv("QUERY_STRING");
    if ((di != NULL) && (strncmp(di, "di=", 3) == 0)) {
        unsigned char *ar = (unsigned char *) unpackargs(di + 3);
        int i, ne;
        char *cp;

        di += 3;
#ifdef DUMP_DI_ARGUMENT_BLOCK
        int na;

        fprintf(stderr, "IV = %02X  Version: %d\n", ar[0], ar[1]);
        na = ar[2];
        fprintf(stderr, "Arguments: %d\n", na);
        cp = (char *) (ar + 3);
        for (i = 1; i < na; i++) {
            fprintf(stderr, "   %d:  \"%s\"\n", i, cp);
            cp += strlen(cp) + 1;
        }

        ne = *((unsigned char *) cp);
        cp++;

        fprintf(stderr, "Elements: %d\n", ne);
        for (i = 0; i < ne; i++) {
            fprintf(stderr, "   %d:  \"%s\"\n", i, cp);
            cp += strlen(cp) + 1;
        }
return 0;
#endif

        if (ar[1] != 1) {
            fprintf(stderr, "%s: Invalid version %d in WWW_di option block.\n", progpath, ar[1]);
            abort();
        }
        argc = ar[2];
        argv = malloc(sizeof(char *) * argc);
        if (argv == NULL) {
            fprintf(stderr, "%s: Unable to allocate %d item argument vector for WWW_di option block.\n",
                progpath, argc);
            abort();
        }
        argv[0] = progpath;

        /*  Link arguments in block to new argv pointer table.  */

        cp = (char *) (ar + 3);
        for (i = 1; i < argc; i++) {
            argv[i] = cp;
            cp += strlen(cp) + 1;
        }

        /*  Transcribe elements (if any) to orbitalElements array.  */

        ne = *((unsigned char *) cp);
        if (ne > maxOrbitalElements) {
            ne = maxOrbitalElements;
        }
        cp++;
        for (i = 0; i < ne; i++) {
            strcpy(savedOrbitalElements[i], cp);
            cp += strlen(cp) + 1;
        }
        embeddedImage = TRUE;
    } else {
        di = NULL;
    }

    /*  If we are called with no command line arguments and no
        dynamic image parameter block was supplied, this is a
        direct invocation which expects us to decode the CGI
        arguments ourself.  We call the embedded uncgi()
        function to unpack them into environment variables
        and then synthesise a command line with options set from them.  */

    if ((di == NULL) && (argc == 1)) {
        char *arg, *arg2;
#define MaxArgs 32
        static char *orgv[MaxArgs + 1];
        int orgc = 1;
        int needuncgi = TRUE, i;

        /* We may have been invoked with a URL which embeds a call
           on uncgi, since people have copied these URLs created
           by the custom request page into their own pages over the
           years.  If so, we don't want to call uncgi() again and mess
           up the already-built environment.  We detect this case by
           scanning the environment to see if it already contains a
           variable which begins with "WWW_".  If so, we skip our
           own call on the uncgi() function. */

        for (i = 0; environ[i] != NULL; i++) {
            if (strncmp(environ[i], "WWW_", 4) == 0) {
                needuncgi = FALSE;
                break;
            }
        }

        if (needuncgi) {
            uncgi();
        }

        directCall = TRUE;

        orgv[0] = argv[0];
        orgv[orgc++] = "-f";
        orgv[orgc++] = "-t";

        /* Dynamic image request */
        if (getenv("WWW_dynimg") != NULL) {
            orgv[orgc++] = "-j";
        } else {
            orgv[orgc++] = "-h";
        }

        /* Set image size */
        if ((arg = getenv("WWW_imgsize")) != NULL) {
            char *s = malloc(strlen(arg) + 3);
            strcpy(s, "-i");
            strcpy(s + 2, arg);
            orgv[orgc++] = s;
        }

        /* Date and time */
        if ((arg = getenv("WWW_date")) != NULL) {
            switch (arg[0]) {
                case '0':               /* Now */
                default:
                    orgv[orgc++] = "-e0";
                    break;

                case '1':               /* UTC date and time */
                    if ((arg2 = getenv("WWW_utc")) != NULL) {
                        char *s = malloc(strlen(arg2) + 4);
                        strcpy(s, "-e1");
                        strcpy(s + 3, arg2);
                        orgv[orgc++] = s;
                    } else {
                        orgv[orgc++] = "-e0";
                    }
                    break;

                case '2':               /* Julian day */
                    if ((arg2 = getenv("WWW_jd")) != NULL) {
                        char *s = malloc(strlen(arg2) + 4);
                        strcpy(s, "-e2");
                        strcpy(s + 3, arg2);
                        orgv[orgc++] = s;
                    } else {
                        orgv[orgc++] = "-e0";
                    }
                    break;
            }
        } else {
            orgv[orgc++] = "-e0";       /* Now */
        }


        /* Heliocentric Latitude */
        if ((arg = getenv("WWW_hlat")) != NULL) {
            char *s = malloc(strlen(arg) + 5);
            strcpy(s, "-da");
            if ((arg2 = getenv("WWW_hns")) != NULL) {
                if (arg2[0] == 'l') {
                    arg2++;
                }
                if (arg2[0] == 'S') {
                    strcat(s, "-");
                }
            }
            strcat(s, arg);
            orgv[orgc++] = s;
        } else {
            orgv[orgc++] = "-da90";
        }

        /* Heliocentric Longitude */
        if ((arg = getenv("WWW_hlon")) != NULL) {
            char *s = malloc(strlen(arg) + 5);
            strcpy(s, "-do");
            strcat(s, arg);
            orgv[orgc++] = s;
        } else {
            orgv[orgc++] = "-do0";
        }

        /* Inner / full solar system selection */
        if (((arg = getenv("WWW_sys")) != NULL) &&
            (strncmp(arg, "-S", 2) == 0)) {
            orgv[orgc++] = strdup(arg);
        } else {
            orgv[orgc++] = "-Sf";
        }

        /* Orbit plotting selection */
        if (((arg = getenv("WWW_orb")) != NULL) &&
            (strncmp(arg, "-b", 2) == 0)) {
            char *s = malloc(4);
            strcpy(s, "-b");
            s[2] = arg[2];
            s[3] = EOS;
            orgv[orgc++] = s;
        } else {
            orgv[orgc++] = "-b0";
        }

        /* Stereo option and cross/wall eye selection */
        if (((arg = getenv("WWW_stereo")) != NULL) &&
            (strcmp(arg, "-v") == 0)) {
            int eyel = 1;
            char *s;
            if ((arg2 = getenv("WWW_eyes")) != NULL) {
                eyel += strlen(arg2);
            }
            s = malloc(strlen(arg) + eyel);
            strcpy(s, arg);
            if (arg2 != NULL) {
                strcat(s, arg2);
            }
            orgv[orgc++] = s;
        }

        /* Dump orbital elements option */
        if (((arg = getenv("WWW_edump")) != NULL) &&
            (strcmp(arg, "-xe") == 0)) {
            orgv[orgc++] = "-xe";
        }

        /* Image or icon selection */
        if (((arg = getenv("WWW_img")) != NULL) &&
            ((strcmp(arg, "-k0") == 0) ||
             (strcmp(arg, "-k1") == 0))) {
            orgv[orgc++] = strdup(arg);
        }

        /* Colour scheme */
        if (((arg = getenv("WWW_scheme")) != NULL) &&
            ((arg[0] >= '0') && (arg[0] <= '4'))) {
            char *s = malloc(4);
            strcpy(s, "-g");
            s[2] = arg[0];
            s[3] = EOS;
            orgv[orgc++] = s;
        }


        /* End of options.  Add terminating "-" and observer's
           latitude and longitude if specified. */

        orgv[orgc++] = "-";

        /* Observer Latitude */
        if ((arg = getenv("WWW_lat")) != NULL) {
            char *s = malloc(strlen(arg) + 2);
            s[0] = EOS;
            if ((arg2 = getenv("WWW_ns")) != NULL) {
                if (arg2[0] == 'l') {
                    arg2++;
                }
                if (arg2[0] == 'S') {
                    strcpy(s, "-");
                }
            }
            strcat(s, arg);
            orgv[orgc++] = s;
        }

        /* Observer Longitude */
        if ((arg = getenv("WWW_lon")) != NULL) {
            char *s = malloc(strlen(arg) + 2);
            s[0] = 0;
            if (((arg2 = getenv("WWW_ew")) != NULL) &&
                (arg2[0] != 'W')) {
                strcpy(s, "-");
            }
            strcat(s, arg);
            orgv[orgc++] = s;
        }

        /* Altitude */
        if ((arg = getenv("WWW_alt")) != NULL) {
            orgv[orgc++] = strdup(arg);
        }

        orgv[orgc] = NULL;

        if (orgc >= MaxArgs) {
            fprintf(stderr, "%s: Too many command line arguments (%d generated, %d max)\n",
                progpath, orgc, MaxArgs);
            abort();
        }
#undef MaxArgs

        for (i = 1; i < orgc; i++) {
            if (orgv[i] == NULL) {
                fprintf(stderr, "%s: Allocation of command line argument %d failed\n",
                    progpath, i);
                abort();
            }
        }

        argv = orgv;
        argc = orgc;

        /* If comet or asteroid orbital elements were specified, save
           a pointer to them for subsequent parsing. */
        if ((arg = getenv("WWW_elements")) != NULL) {
            /* If the last line of the elements is not terminated by
               a line feed, add one so alement_rs() in astrelem.c
               doesn't get confused. */
            if (arg[strlen(arg) - 1] != '\n') {
                char *s = malloc(strlen(arg) + 2);
                strcpy(s, arg);
                strcat(s, "\n");
                arg = s;
            }
            elementSpecification = arg;
        }

    }

    /*  Process command line options.  */

DB("Args");
    for (i = 1; i < argc; i++) {
        cp = argv[i];
DB(cp);
        if (strcmp(cp, "-") == 0) {   /* Bare "-" marks end of options, start of arguments */
            endopts = TRUE;
            continue;
        }
        if ((!endopts) && (*cp == '-' && !isdigit(cp[1]))) {
            opt = *(++cp);
            if (islower(opt))
                opt = toupper(opt);
            switch (opt) {

                case 'B':             /* Specify orbit type */
                    orbittype = atoi(cp + 1);
                    if (orbittype < 0 || orbittype > 2) {
                        orbittype = 0;
                    }
                    break;

                case 'D':             /* Set heliocentric viewpoint */
                    opt = cp[1];
                    if (islower(opt)) {
                        opt = toupper(opt);
                    }
                    if (opt == 'A') {
                        orrLat = parseDeg(cp + 2);
                    } else if (opt == 'O') {
                        orrLong = parseDeg(cp + 2);
                    }
                    break;

                case 'E':             /* Epoch */
                    /* Dates, however specified, are limited to the valid
                       domain for VSOP87, -4712 through +8000.  Dates
                       which are outside that range, or malformed,
                       flip back to the current date and time. */
                    switch (cp[1]) {
                        case '0':
                            dto = 0;
                            break;

                        case '1':
                            {
                                char dlm;

                                dto = 1;
                                ctim.tm_year = -5000;
                                ctim.tm_mon = ctim.tm_hour =
                                    ctim.tm_min = ctim.tm_sec = 0;
                                ctim.tm_mday = 1;
                                sscanf(cp + 2, "%d%c%d%c%d %d%c%d%c%d",
                                    &ctim.tm_year, &dlm, &ctim.tm_mon, &dlm, &ctim.tm_mday,
                                    &ctim.tm_hour, &dlm, &ctim.tm_min, &dlm, &ctim.tm_sec);
                                if (ctim.tm_year < -4712 || ctim.tm_year > 8000 ||
                                    ctim.tm_mon < 1 || ctim.tm_mon > 12 ||
                                    ctim.tm_mday < 1 || ctim.tm_mday > 31 ||
                                    ctim.tm_hour < 0 || ctim.tm_hour >= 24 ||
                                    ctim.tm_min < 0 || ctim.tm_min >= 60 ||
                                    ctim.tm_sec < 0 || ctim.tm_sec >= 60) {
                                    dto = 0;
                                } else {
                                    ctim.tm_year -= 1900;
                                    ctim.tm_mon--;
#ifndef Solaris
                                    ctim.tm_gmtoff = 0;
                                    ctim.tm_zone = "UTC";
#endif
                                    ctim.tm_isdst = 0;
                                    jt = jtime(&ctim);
                                    ctim.tm_wday = fmod(jt + 1.5, 7);
                                }
                                break;
                            }

                        case '2':
                            {   long a;
                                double ejt;

                                dto = 2;
                                jt = atof(cp + 2);
                                if (jt < 0.5 || jt > 4346655.5) {
                                    dto = 0;
                                } else {
                                    ejt = jt - 0.5;
                                    jyear(jt, &a, &ctim.tm_mon, &ctim.tm_mday);
                                    ctim.tm_wday = fmod(jt + 1.5, 7);
                                    ctim.tm_year = a - 1900;
                                    ctim.tm_mon--;
#ifndef Solaris
                                    ctim.tm_gmtoff = 0;
                                    ctim.tm_zone = "UTC";
#endif
                                    ctim.tm_isdst = 0;
                                    a = (long) (((ejt - ((long) ejt)) * 86400.0) + 0.5);
                                    ctim.tm_hour = a / 3600;
                                    ctim.tm_min = (a / 60) % 60;
                                    ctim.tm_sec = (a % 60);
                                }
                            }
                            break;
                    }
                    break;

                case 'F':
                    stateless = TRUE; /* Use stateless image protocol */
                    break;

                case 'G':             /* -G[0-4]  --  Colour scheme  */
                    skyColourScheme = atoi(cp + 1);
                    if (skyColourScheme < 0 || skyColourScheme > 4) {
                        skyColourScheme = 0;
                    }
                    break;

                case 'H':             /* Generate HTML output */
                    html = TRUE;
                    break;

                case 'I':             /* Image size */
                    tmsize = atoi(cp + 1);
                    tmsize = max(100, min(tmsize, 1024));
                    break;

                case 'J':             /* Generate dynamic image result */
                    dynimg = TRUE;
                    break;

                case 'K':             /* Specify icon file */
                    useimage = atoi(cp + 1);
                    break;

                case 'O':             /* Output .GIF file */
                    cachep = strstr(cp + 1, "/cache");
                    assert(cachep != NULL);
                    strcpy(cacheName, cachep);
                    ofile = fopen(cp + 1, "w");
                    if (ofile == NULL) {
                        ofile = stdout;
                    }
                    break;

                case 'S':             /* Set Inner or Outer system */
                    opt = cp[1];
                    if (islower(opt)) {
                        opt = toupper(opt);
                    }
                    if (opt == 'I') {
                        innersystem = TRUE;
                    } else if (opt == 'O') {
                        innersystem = FALSE;
                    }
                    break;

                case 'T':             /* Track asteroid or comet */
                    trackelem = TRUE;
                    break;

                case 'V':
                    switch (cp[1]) {
                        case '0':
                            stereoView = TRUE;
                            stereoOcular = abs(stereoOcular);
                            break;

                        case '1':
                            stereoView = TRUE;
                            stereoOcular = -abs(stereoOcular);
                            break;
                    }
                    break;

                case 'W':             /* Write PPM to standard output */
                    wppm = TRUE;
                    break;

                case 'X':             /* Diagnostic or extended output */
                    opt = cp[1];
                    if (islower(opt)) {
                        opt = toupper(opt);
                    }
                    switch (opt) {
                        case 'E':
                            dumpelem = TRUE;
                            break;

#ifdef SAVESET
                        case 'H':
                            gotten = TRUE;
                            break;

                        case 'S':
                            bookmark = TRUE;
                            break;
#endif
                    }
                    break;

                case '?':
                case 'U':
    fprintf(stderr,"\nSOLAR.  Call");
    fprintf(stderr,
       "\n             with solar [options] [latitude [[longitude] [altitude]]]");
    fprintf(stderr,"\n");
    fprintf(stderr,"\n        Options:");
    fprintf(stderr,"\n              -       Mark end of options");
    fprintf(stderr,"\n              -bn     Orbits: 0 = log, 1 = real, 2 = equal");
    fprintf(stderr,"\n              -dAnnn  Set heliocentric latitude");
    fprintf(stderr,"\n              -dOnnn  Set heliocentric longitude");
    fprintf(stderr,"\n              -e0     Epoch = Now");
    fprintf(stderr,"\n              -e1dt   Epoch = Date/time in UTC");
    fprintf(stderr,"\n              -e2jd   Epoch = Julian date jd");
    fprintf(stderr,"\n              -f      Use stateless image protocol");
    fprintf(stderr,"\n              -gN     Colour scheme N 0=colour, 1=black/white 2=white/black 3=night vision");
    fprintf(stderr,"\n              -h      Create HTML file");
    fprintf(stderr,"\n              -isize  Image size in pixels");
    fprintf(stderr,"\n              -j      Generate dynamic image");
    fprintf(stderr,"\n              -kN     Select  n = 0 icon, n = 1 image");
    fprintf(stderr,"\n              -ofname GIF output to fname");
    fprintf(stderr,"\n              -si/o   Show Inner/Outer system");
    fprintf(stderr,"\n              -t      Track asteroid or comet");
    fprintf(stderr,"\n              -u      Print this message");
    fprintf(stderr,"\n              -vn     Stereo n: 0 = cross, 1 = wall");
    fprintf(stderr,"\n              -w      Write PPM, not GIF output");
    fprintf(stderr,"\n              -xe     Echo asteroid/comet elements");
#ifdef SAVESET
    fprintf(stderr,"\n              -xh     Use action=get for form to save parameters");
    fprintf(stderr,"\n              -xs     Save settings for bookmark");
#endif
    fprintf(stderr,"\n");
    fprintf(stderr, "        Version %s\n", Version);
                    return 0;
            }
        } else {
            if (isdigit(cp[0]) || cp[0] == '+' || cp[0] == '-') {
                switch (f) {
                    case 0:
                        siteLat = parseDeg(cp);
                        f++;
                        break;

                    case 1:
                        siteLon = parseDeg(cp);
                        f++;
                        break;

                    case 2:
                        talt = atof(cp);
                        f++;
                        break;


                    default:
                        fprintf(stderr, "Extra coordinates ignored.\n");
                }
            }
        }
    }

    /*  If we received a WWW_di argument block, adjust the options accordingly
        to generate the dynamic image.  This is a big kludge and should be looked
        at closely when options are added or modified.   */

    if (di) {
        stateless = FALSE;
        html = FALSE;
        dynimg = TRUE;
    }

    orrInner = innersystem;
    orrScale = orbittype;

    /* Limit image size to 512 for stereo views (saves time and
       bandwidth, and 512 is really too big for most people to
       fuse anyway on typical screens. */

    if (stereoView) {
        tmsize = min(tmsize, 512);
    }

    if (dto == 0) {
        time(&cctime);
        ctim = *gmtime(&cctime);
        jt = jtime(&ctim);
    }
    sunpos(jt, TRUE, &sunra, &sundec, &earthrv, &sunlong);
    earthlong = sunlong + 180;
    gt = gmst(jt);
    sprintf(dtutc, "%d-%02d-%02d %d:%02d:%02d",
        ctim.tm_year + 1900, ctim.tm_mon + 1, ctim.tm_mday,
        ctim.tm_hour, ctim.tm_min, ctim.tm_sec);

    tlat = dtr(siteLat);
    tlon = dtr(siteLon);
    vlat = tlat;
    vlon = tlon;
    valt = talt;
    sprintf(viewdesc, "View from %s", edvpos(vlat, vlon));

    /* Load icon bitmap. */

    load_bitmap(locfile(useimage ? "solar_images.bmp" :
        Cscheme("yourtel-icons.bmp", "yourtel-icons.bmp", "yourtel-icons-b.bmp", "yourtel-icons-w.bmp", "yourtel-icons-r.bmp")));
    assert(bmBits != NULL);

    /* Find the colour table indices of colours we're going to use
       in drawing in the bitmap.  Note that we can reserve any needed
       colours which aren't used in the icons by including a dummy icon
       at the right end of the icon bitmap which references all the
       additional required colours. */

#define findColour(cindex, r, g, b)  if (bmi.bmiColors[i].rgbRed == r &&  \
                                         bmi.bmiColors[i].rgbGreen == g && \
                                          bmi.bmiColors[i].rgbBlue == b)    \
                                          cindex = i
    for (i = 0; i < bmi.bmiHeader.biClrUsed; i++) {
        findColour(transparent, 255, 0, 255);
        findColour(cBlack, 0, 0, 0);
        findColour(cBlue, 0, 0, 255);
        findColour(cGreen, 0, 255, 0);
        findColour(cBorder, 182, 182, 182);
        findColour(cWhite, 255, 255, 255);
        findColour(cGrey, 128, 128, 128);
        findColour(cLtGrey, 192, 192, 192);
        findColour(cDkRed, 128, 0, 0);
        findColour(cRed, 255, 0, 0);
    }

    /*  Set time limit to bail us out of possible loops due to
        errors in rendering code.  */

#ifdef TimeLimit
    sigset(SIGALRM, (void (*)()) dingaling);
    alarm(TimeLimit);
#endif

    /*  If we're tracking an asteroid or comet, read its orbital
        elements.  */

    if (trackelem) {
        astrelem(di != NULL);
    }

#ifdef SDEBUG
    /*  If this is a stateless request, assemble the argument block for the image request.  */

    if (stateless) {
        unsigned char *ha = packargs(argc, argv);
        int i, na, ne;
        char *cp;
        unsigned char *ar = unpackargs(ha);

        fprintf(stderr, "IV = %02X  Version: %d\n", ar[0], ar[1]);
        na = ar[2];
        fprintf(stderr, "Arguments: %d\n", na);
        cp = (char *) (ar + 3);
        for (i = 1; i < na; i++) {
            fprintf(stderr, "   %d:  \"%s\"\n", i, cp);
            cp += strlen(cp) + 1;
        }

        ne = *((unsigned char *) cp);
        cp++;

        fprintf(stderr, "Elements: %d\n", ne);
        for (i = 0; i < ne; i++) {
            fprintf(stderr, "   %d:  \"%s\"\n", i, cp);
            cp += strlen(cp) + 1;
        }
    }
#endif

    /* Calculate positions of objects at the requested epoch. */

    buildPlanets(jt);

    /*  Generate HTTP response header.  */

    if (html) {
        printf("Content-type: text/html\r\n");
        printf("\r\n");
    } else {
        if (dynimg) {
            char *ruri = getenv("REQUEST_URI");

            if (wppm) {
                printf("Content-type: image/x-portable-pixmap\r\n");
            } else {
                printf("Content-type: image/gif\r\n");
            }
            if (ruri != NULL) {
                printf("Content-Location: %s\r\n", ruri);
            }
            printf("Pragma: no-cache\r\n");
            /* printf("Cache-Control: no-cache\r\n");
               *** "Cache-Control: pricate" added by server for all cgi-bin
                   requests at Fourmilab. *** */
            printf("\r\n");
        }
        /*  If dynimg is not set, output is simply the image, written
            to standard output, with no HTTP response header.  This is
            used primarily for low-level debugging of image generation
            without the need to get the whole Web pipeline in the act. */
    }

    /*  If this is not a stateless request, generate the Orrery image.  */

    if (dynimg || (!stateless)) {

        /* Generate the Orrery image. */

        orrery(tmsize, bmi.bmiHeader.biWidth, bmi.bmiHeader.biHeight,
               tmsize, tmsize, jt, stereoView);
        if (wppm) {
            dump_ppm("-", tmbmap, tmbits);
        } else {
            int tmrowlen = (tmbmap->bmiHeader.biWidth + 3) & ~3;

            gifout(tmbmap->bmiHeader.biWidth, tmbmap->bmiHeader.biHeight, tmrowlen,
                   tmbmap->bmiColors, tmbmap->bmiHeader.biClrUsed,
                   cBorder,
                   tmbits, ofile);
        }

        if (ofile != stdout) {
            fclose(ofile);
        }
    }

    /* Generate HTML reply page referencing the image. */

    if (html) {
        static char *weekdays[] = { "Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat" };

        sprintf(tbuf, "%s %d %s %d %d:%02d", weekdays[ctim.tm_wday],
            ctim.tm_year + 1900, mname[ctim.tm_mon], ctim.tm_mday,
            ctim.tm_hour, ctim.tm_min);

        printf("\
<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>\n\
<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\"\n\
    \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n\
<html xmlns=\"http://www.w3.org/1999/xhtml\" xml:lang=\"en\" lang=\"en\">\n");

        printf("<head>\n<title>\nSolar System Live\n</title>\n</head>\n");
        printf("<body %s>\n", (skyColourScheme == SCS_NIGHT_VISION) ?
            "bgcolor=\"#000000\" text=\"#FF0000\"" :
            "bgcolor=\"#FFFFFF\"");

#ifdef UNDERC
        printf("<h1><img src=\"/images/travaux.gif\" alt=\"\" align=\"middle\" /> Under Construction</h1>\n");
#endif

        if (dto == 0) {
            printf("<h1>Solar System Live");
        } else {
            printf("<h1>Solar System: %s", tbuf);
        }
#ifdef SAVESET
        printf("%s</h1>\n", gotten ? " (Custom)" : "");
#else
        printf("</h1>\n");
#endif
        printf("<center>\n");
#ifdef CACHE_WARNING
        if (!stateless) {
            write_cache_image_warning(stdout);
        }
#endif
#ifdef ETEST
        printf("<img src=\"test.gif\" width=\"%d\" height=\"%d\" /><br />\n",
            (int) tmsize, (int) tmsize);
#else
        if (stateless) {
#ifdef OLDWAY
            printf("<img src=\"/cgi-bin/SolarImage%s?di=%s\" width=\"%d\" height=\"%d\" alt=\"\" /><br />\n",
#else
            printf("<img src=\"/cgi-bin/Solar%s?di=%s\" width=\"%d\" height=\"%d\" alt=\"\" /><br />\n",
#endif
#ifdef TESTMODE
                ".t",
#else
                "",
#endif
                packargs(argc, argv),
                (int) (stereoView ? (tmsize * 2 + stereoSepPixels) : tmsize), (int) tmsize);
        } else {
            printf("<img src=\"%s\" width=\"%d\" height=\"%d\" alt=\"\" /><br />\n",
                cacheName, (int) (stereoView ? (tmsize * 2 + stereoSepPixels) : tmsize), (int) tmsize);
        }
#endif
        printf("</center>\n");
        writepost(stdout, vlat, vlon, valt, imagesource, innersystem,
            orbittype,
            tmsize, jt);
    }

#ifdef ShowCPUtime
    {
        struct rusage u;

        getrusage(RUSAGE_SELF, &u);
        printf("CPU time: %ld seconds.<br />\n", u.ru_utime.tv_sec + u.ru_stime.tv_sec);
    }
#endif

    if (directCall) {
        FILE *postamble = fopen(locfile("solar-post.html"), "r");

#ifdef SHOW_ARGUMENTS
{
int i;

printf("<pre>\nEnvironment Variables:\n\n");
for (i = 0; environ[i] != NULL; i++) {
    printf("    %s\n", environ[i]);
}
printf("\nCommand Line Arguments:\n\n");
for (i = 1; i < argc; i++) {
    printf("    %d: %s\n", i, argv[i]);
}
printf("</pre>\n");
}
#endif

        if (postamble != NULL) {
            char s[1024];
            int l;

            while ((l = fread(s, 1, sizeof s, postamble)) > 0) {
                fwrite(s, 1, l, stdout);
            }
            fclose(postamble);
        }
    }

#ifdef DBLOG
    fclose(db);
#endif
    return 0;
}
