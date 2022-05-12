/*

                       Postage Stamp Rasteriser

        by John Walker, Autodesk Neuchatel, 25e septembre 1991

    Ever  needed  a  little  rasteriser?  A little one--small, simple,
    suitable, for example, for creating  preview  bitmaps  that  allow
    users  to  select images graphically rather than typing in stoopid
    little file names that mean nothing to nobody.

    PSTAMPR.C is a software component that performs vector and polygon
    scan  conversion   into   monochrome   bitmaps   of   user-defined
    resolution.   It  handles  binary  (light and dark) bitmaps of any
    addressable resolution and scan converts both vectors and polygons
    (concave  or  convex).   Its handling of concave polygons uses the
    latest, most accurate algorithm known  to  the  author  (patiently
    evolved within the belly of AutoCAD for over eight years).

    This  component  is  intended  to  be  clean, self-sufficient, and
    independent.  It isn't particularly efficient, but neither  is  it
    profligate  in  its  use  of  your CPU resources.  It's called the
    "Postage Stamp Rasteriser" because it's intended for making small,
    simple,   monochrome   images,  however,  there  are  no  inherent
    limitations  that  prevent  your  using  it  in   more   ambitious
    applications.  If performance is critical, however, you'll want to
    optimise some of the algorithms: AutoCAD's implementation of these
    functions is far more CPU efficient, but at the cost of enormously
    increased complexity of the code.  PSTAMPR should, in  all  cases,
    produce   the  same  bitmaps  as  AutoCAD's  rasteriser  would  if
    presented with the  same  vectors  and  polygons.   Built-in  code
    allows  exporting the bitmaps created with PSTAMPR as uncompressed
    TIFF files.

    You access PSTAMPR through the following functions.   For  clarity
    of  documentation,  the  functions  are  described  in  prototyped
    declaration form.  For compatibility  with  vintage  C  compilers,
    however, the actual declarations use the old-fashioned form.

    INITIALISATION
    ==============

    To use PSTAMPR to scan convert a vector image, first initialise it
    by calling:

        char *psrinit(unsigned int xsize, unsigned int ysize,
                      unsigned int background)

    where  xsize and ysize give, respectively, the width and height of
    the bitmap in  pixels  and  background  specifies  the  background
    colour  (zero or nonzero) to which the bitmap will be initialised.
    psrinit() returns a character pointer to the  bitmap,  dynamically
    allocated  with  malloc(),  or  NULL  if  the  bitmap could not be
    allocated.  The bitmap is stored as 8 bits per char.  Your program
    can  set  and  read the bits without knowing the precise format of
    the  bitmap  by  calling  the  psrsetpixel()   and   psrgetpixel()
    functions  (see  below).   If,  for efficiency's sake, you want to
    directly access the bitmap, refer to the code in  those  functions
    for  an example of how the bitmap is addressed.  The bitmap can be
    internally organised with pixel 0,0 at the bottom left or the  top
    left.   If  the  compile-time variable FlipY is defined, 0,0 is at
    the bottom left; if it isn't defined, then pixel 0,0 is at the top
    left.  Lines in the bitmap are padded to the next  byte  boundary;
    you  can  determine the number of bytes in each line of the bitmap
    by calling psrlinelen().

    The pointer returned by psrinit() is a "handle" you  pass  to  all
    the  other  functions of PSTAMPR.  The fact that all state is kept
    in the bitmap itself permits your program to use PSTAMPR to create
    any  number of bitmaps simultaneously, even in parallel threads of
    a multitasking program.  IMPORTANT--to dispose of the  bitmap  you
    *must*  call  psrterm()  with the handle of the bitmap rather than
    attempting to free the handle  directly  with  free().   psrinit()
    allocates  some local storage at the start of the bitmap and hence
    the handle it returns doesn't point to the start of the  allocated
    buffer.  Consequently, passing that pointer to free() will lead to
    disaster.

    PIXEL SCAN CONVERSION
    =====================

    To set an individual pixel in the bitmap, call:

        void psrsetpixel(char *bhandle, unsigned int x, unsigned int y,
                                        unsigned int colour)

    where bhandle is the handle to the bitmap returned by psrinit(), x
    and  y specify the co-ordinates of the pixel to be set, and colour
    specifies whether the pixel is to be set to zero (colour == 0)  or
    one (colour != 0).

    VECTOR SCAN CONVERSION
    ======================

    To draw a vector into a bitmap, call:

        void psrvector(char *bhandle, unsigned int fx, unsigned int fy,
                                      unsigned int tx, unsigned int ty,
                                      unsigned int colour)

    bhandle  is  the  handle  to  the  bitmap,  fx  and fy specify one
    endpoint of the vector while tx and ty give the other.  The colour
    of  the  vector  (whether  it  is  drawn  in  zero or one bits) is
    determined by whether colour is zero or nonzero.

    POLYGON SCAN CONVERSION
    =======================

    To draw a solid-filled polygon (which may be  concave  or  convex)
    into the bitmap, call:

        void psrpoly(char *bhandle, struct spolygon *poly,
                                    unsigned int colour)

    Once  again, bhandle is the handle to the bitmap.  The vertices of
    the polygon are given in the spolygon  structure  poly,  which  is
    declared as follows:

        struct spoint {
            unsigned int x, y;
        };

        struct spolygon {
            int npoints;                     Number of points in polygon
            struct spoint pt[MAXVERT + 1];   Actual points
        };

    The compile-time variable MAXVERT specifies the maximum vertices a
    polygon may contain.  This number must be known  when  PSTAMPR  is
    compiled  in  order  to  correctly  allocate  stack space used for
    tables of edges and intersections within psrpoly().  Finally,  the
    colour of the polygon (whether it is scan converted to zero or one
    bits) is specified by whether  the  argument  colour  is  zero  or
    nonzero.

    RETRIEVING INDIVIDUAL PIXELS
    ============================

    You  may  read  the  value  of an individual pixel from the bitmap
    with:

        int psrgetpixel(char *bhandle, unsigned int x, unsigned int y)

    where bhandle is the handle to the bitmap, and x and y specify the
    column  and row address of the pixel to be read.  The value of the
    pixel is returned by psrgetpixel():  0  or  1.   Reading  out  the
    entire bitmap with psrgetpixel() is staggeringly inefficient; it's
    far better to write code that accesses the bitmap directly.

    DIRECT BITMAP ACCESS
    ====================

    To access a row of pixels in the bitmap, you need  only  know  the
    row address and the number of bytes each row of pixels occupies in
    memory.  You can obtain the memory length of pixel rows with:

        int psrlinelen(char *bhandle)

    where bhandle is the handle to the bitmap.  The  function  returns
    the row length in bytes.  For example, the pixels in row 12 of the
    image myImage are stored at the address:

        char *row12, *myImage;

        row12 = myImage + (12 * psrlinelen(myImage));

    Pixels are stored 8 per byte, with the leftmost pixel in the  high
    order  (0x80)  bit.   Thus,  the pixel in column X of row12 may be
    extracted with the expression:

        ((*(row12 + (X >> 3)) & (0x80 >> (X & 7))) ? 1 : 0)

    Obviously, for efficiency you  should  call  psrlinelen()  only  a
    single  time  after creating the bitmap with psrinit(), then store
    its value in a local  variable  for  subsequent  accesses  to  the
    bitmap.

    EXPORTING THE IMAGE IN TIFF FORMAT
    ==================================

    You can write your bitmap as an uncompressed monochrome TIFF image
    file by calling:

        void psrtiff(char *bhandle, FILE *fd)

    where bhandle is the handle to the bitmap and fd is the handle  of
    a  binary  file  already  opened  for  writing.  The TIFF image is
    written starting at the current file position of the file fd,  and
    after  the image has been written, the current file position of fd
    will be left at the first byte following the TIFF image.  Since no
    seeks  are  performed  on  fd,  it may be any type of file, even a
    device file or pipe which cannot be positioned.

    RELEASING A BITMAP
    ==================

    To release the storage allocated for a bitmap, call:

        void psrterm(char *bhandle)

    where bhandle is the image handle.  The storage allocated for  the
    bitmap  in  psrinit() will be freed.  No references may be made to
    the bitmap using the handle after its storage is relinquished with
    psrterm().

    BUILT-IN TEST PROGRAM
    =====================

    If  you  compile  PSTAMPR.C  with  TestProgram  defined,  a main()
    function is included that performs a fairly  demanding  regression
    test on the PSTAMPR functions.  It creates a bitmap using a mix of
    vectors, polygons, and pixels, then exports the bitmap as  a  TIFF
    file.   The TIFF file is then read back and compared byte-for-byte
    against a known-correct canned TIFF file.  Extraneous  garbage  at
    the  end  of  the TIFF file is also reported.  If you compile with
    PosRast defined,  the  test  program  will  write  the  bitmap  to
    standard  output  in  the  PBM format used by Jef Poskanzer raster
    toolkit (normally ASCII mode, but binary if you also define Binary
    at  compile time).  This can be handy if you run into problems and
    want to use the PBMPLUS tools to analyse the generated bitmap (for
    example,  subtracting  what you got from a known correct bitmap to
    see which pixels differ).

    CONFIGURING FOR PERFORMANCE
    ===========================

    PSTAMPR   uses   the  <assert.h>  package  extensively  to  verify
    arguments to functions and to make  internal  consistency  checks.
    Once  you're  happy  that  PSTAMPR  is  working  properly, you may
    compile it with the symbol NDEBUG defined at compile  time,  which
    will  remove  the  assertions and the substantial compute overhead
    they create.

*/

#include    <stdio.h>
#include    <assert.h>
#ifdef lint
#include    <memory.h>
#endif

/*  If FlipY is defined the screen is addressed with pixel 0,0 at the
    bottom left.  If it isn't defined, pixel 0,0 is at the top left. */

#define FlipY

#define MAXVERT 10                    /* Maximum vertices in polygon */

#define HIGHSCRC 32767                /* Largest screen co-ordinate */

/* Screen point */

struct spoint {
    unsigned int x, y;
};

/* Screen polygon */

struct spolygon {
    int npoints;                      /* Number of points in polygon */
    struct spoint pt[MAXVERT + 1];    /* Actual points */
};

#define FALSE    0
#define TRUE     1

#define V        (void)

/* LINTLIBRARY */

#ifndef min
#define min(a, b)   ((a) < (b) ? (a) : (b))
#endif

#ifndef max
#define max(a, b)   ((a) > (b) ? (a) : (b))
#endif

extern char *malloc();

/*  PSRINIT  --  Initialise  the rasteriser.  Passed the image size in
                 pixels.  Allocates and clears the image buffer to the
                 specified  background  colour  and initialises static
                 variables used by the  draw  routines.   Returns  the
                 address of the image array or zero if the array could
                 not be allocated. */

char *psrinit(xsize, ysize, background)
  unsigned int xsize, ysize, background;
{
    int xbytes = (xsize + 7) / 8;     /* Calculate row width in bytes */
    char *pparray;
    int *iarray;

    /* Allocate image array */

    assert(xsize > 0);
    assert(ysize > 0);
    iarray = (int *) malloc((3 * sizeof(int)) + (ysize * xbytes));
    if (iarray == NULL) {
        return NULL;
    }
    pparray = (char *) (&iarray[3]);
    iarray[0] = xsize;
    iarray[1] = ysize;
    iarray[2] = xbytes;
    V memset(pparray, background ? 0xFF : 0, (int) (ysize * xbytes));
    return pparray;
}

/*  PSRLINELEN  -- Return the line length, in bytes, of the array. */

int psrlinelen(image)
  char *image;
{
    int *iarray = (int *) (image - (3 * sizeof(int)));

    return iarray[2];
}

/*  PSRVECTOR  --  Draw a vector into the current image array.  */

void psrvector(pparray, fx, fy, tx, ty, colour)
  char *pparray;
  unsigned int fx, fy, tx, ty, colour;
{
    unsigned int cfx, cfy, ctx, cty;
    int *iarray = (int *) (pparray - (3 * sizeof(int)));
    int x, y, m, f, xinc, loopcnt, dx, dy, xbytes = iarray[2];

    assert(fx >= 0 && fx < iarray[0]);
    assert(fy >= 0 && fy < iarray[1]);
    assert(tx >= 0 && tx < iarray[0]);
    assert(ty >= 0 && ty < iarray[1]);
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
            char *where = pparray + (x >> 3) + (y * xbytes);
            char mask = 0x80 >> (x & 0x7);

            *where = (*where & (~mask)) | (colour ? mask : 0);
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
            char *where = pparray + (x >> 3) + (y * xbytes);
            char mask = 0x80 >> (x & 0x7);

            *where = (*where & (~mask)) | (colour ? mask : 0);
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

/*  PSRSETPIXEL  --  Set a pixel in the image to a given value. */

static void psrsetpixel(image, x, y, value)
  char *image;
  unsigned int x, y, value;
{
    int *iarray = (int *) (image - (3 * sizeof(int)));
    int bit = 0x80 >> (x & 0x7);
    char *padr;

    assert(x >= 0 && x < iarray[0]);
    assert(y >= 0 && y < iarray[1]);

#ifdef FlipY
    y = (iarray[1] - 1) - y;          /* Flip picture to right way up */
#endif
    padr = image + (x >> 3) + (y * iarray[2]);
    *padr = (*padr & (~bit)) | (value ? bit : 0);
}

/*  PSRPOLY  --  Rasterise  a  polygon  of  a  given number of points. */

void psrpoly(image, poly, colour)
  char *image;
  struct spolygon *poly;
  unsigned int colour;
{

    int minint, maxint, thisint;
    struct edg {
            int ymax, ymin, xpos, xend;
            double slope;
    };
    struct edg edge[MAXVERT];
    unsigned int inter[MAXVERT];
    int i, j, k, l, m, n, nvert, nedge, nint, horedge, iterations, count;
    unsigned int s, t, u, y, bymax, bymin;

#define vertex poly->pt

    nvert = poly->npoints;

    bymin = HIGHSCRC;
    bymax = 0;
    for (i = 0; i < poly->npoints; i++) {
        bymin = min(bymin, (y = poly->pt[i].y));
        bymax = max(bymax, y);
    }

    horedge = FALSE;
    assert(nvert <= MAXVERT);

    if (nvert == 1) {
        psrsetpixel(image, poly->pt[0].x, poly->pt[0].y, colour);
        return;
    };

    for (i = j = 0; i < nvert; i++) {
        k = (i + 1) % nvert;
        if (vertex[i].y > vertex[k].y) {
            l = i;
            m = k;
        } else if (vertex[i].y == vertex[k].y) {
            horedge = TRUE;
            continue;
        } else {
            l = k;
            m = i;
        }
        edge[j].ymax = vertex[l].y;
        edge[j].ymin = vertex[m].y;
        edge[j].xpos = vertex[l].x;
        edge[j].xend = vertex[m].x;

        edge[j].slope =
                (double) ((long) edge[j].xend - (long) edge[j].xpos) /
                (double) ((long) edge[j].ymin - (long) edge[j].ymax);
        j++;
    }
    nedge = j;
    iterations = bymax - bymin;

    for (y = bymax, count = 0; count <= iterations; y--, count++) {
        nint = 0;
        minint = HIGHSCRC;
        maxint = 0;

        for (j = 0; j < nedge; j++) {
            if (y <= edge[j].ymax && y >= edge[j].ymin) {
                if (y == edge[j].ymin) {
                    thisint = inter[nint] = edge[j].xend;
                } else {
                    thisint = inter[nint] =
                       (edge[j].xpos - edge[j].slope * (edge[j].ymax - y)) +
                       0.4999;
                }
                if (thisint <= minint) {
                     minint = thisint;
                }
                if (thisint >= maxint) {
                    maxint = thisint;
                }
                nint++;
            }
        }

        if (nint > 1) {
            for (i = 0; i < (nint - 1); i++) {
                for (j = i + 1; j < nint; j++) {
                    if (inter[i] >= inter[j]) {
                        s = inter[i];
                        inter[i] = inter[j];
                        inter[j] = s;
                   }
                }
            }
            if (nint > 2) {
                for (i = 0; i < (nint - 1); i++) {
                    if (inter[i] == inter[i + 1]) {
                        for (j = 0; j < nvert; j++) {
                            if (vertex[j].y == y && vertex[j].x == inter[i]) {
                                break;
                            }
                        }
                        if (j < nvert) {
                            s = vertex[(j == 0) ? nvert - 1 : j - 1].y;
                            t = vertex[(j + 1) % nvert].y;
                            u = vertex[j].y;
                            if (!((u > s && u > t) || (u < s && u < t))) {
                                for (j = i + 1; j < (nint - 1); j++) {
                                    inter[j] = inter[j + 1];
                                }
                                nint--;
                            }
                        }
                    }
                }
            }
            if (horedge && (nint > 2)) {
                for (i = 0; i < nvert; i++) {
                    if ((y == vertex[i].y) &&
                        (y == vertex[(i + 1) % nvert].y)) {
                        l = vertex[i].x;
                        m = vertex[(i + 1) % nvert].x;
                        for (k = 0; k < nint - 1; k++) {
                            if (((l == inter[k]) && (m == inter[k + 1])) ||
                                ((m == inter[k]) && (l == inter[k + 1]))) {

                                u = 0;
                                s = vertex[((i + 2) % nvert)].y;
                                t = vertex[((i + nvert) - 1) % nvert].y;
                                if ((k & 1) == 0) {
                                    if (((s <= y) && (t >= y)) ||
                                        ((s >= y) && (t <= y))) {
                                        u = 1;
                                    }
                                } else { /* ((k & 1) == 1) */
                                    assert(l == inter[k] || l == inter[k + 1]);
                                    u = 1;
                                    if (((s <= y) && (t <= y)) ||
                                        ((s >= y) && (t >= y))) {
                                        u++;
                                    }
                                }
                                while (u-- > 0) {
                                    for (n = k; n < (nint - 1); n++)
                                        inter[n] = inter[n + 1];
                                    nint--;
                                }
                            }
                        }
                    }
                }
            }
            for (i = 0; i < (nint - 1); i += 2) {
                psrvector(image, inter[i], y, inter[i + 1], y, colour);
            }
        }
    }
}

/*  PSRGETPIXEL  --  Get a pixel from the image.  */

int psrgetpixel(image, x, y)
  char *image;
  unsigned int x, y;
{
    int *iarray = (int *) (image - (3 * sizeof(int)));
    char *padr;

    assert(x >= 0 && x < iarray[0]);
    assert(y >= 0 && y < iarray[1]);

#ifdef FlipY
    y = (iarray[1] - 1) - y;          /* Flip picture to right way up */
#endif
    padr = image + (x >> 3) + (y * iarray[2]);
    return (*padr & (0x80 >> (x & 0x7))) ? 1 : 0;
}

/*  PSRTERM  --  Terminate rasteriser and release storage.  */

void psrterm(image)
  char *image;
{
    int *iarray = (int *) (image - (3 * sizeof(int)));

    free((char *) iarray);
}

/*  PSRTIFF  --  Export image as a TIFF file.  */

/*  Write a short in LSB-first order. */

static void Putshort(w, fd)
  int w;
  FILE *fd;
{
    V putc((char) (w & 0xFF), fd);
    V putc((char) ((w >> 8) & 0xFF), fd);
}

/*  Write a long in LSB-first order. */

static void Putlong(l, fd)
  long l;
  FILE *fd;
{
    V putc((char) (l         & 0xFF), fd);
    V putc((char) ((l >>  8) & 0xFF), fd);
    V putc((char) ((l >> 16) & 0xFF), fd);
    V putc((char) ((l >> 24) & 0xFF), fd);
}

#define TIFF_VERSION    42            /* TIFF version indicator */

#define TIFF_LITTLEENDIAN       0x4949  /* Intel order.  Honest endian! */

#define TIFFTAG_OSUBFILETYPE            255     /* kind of data in subfile */
#define     FILETYPE_REDUCEDIMAGE       0x1     /* reduced resolution version */
#define TIFFTAG_IMAGEWIDTH              256     /* image width in pixels */
#define TIFFTAG_IMAGELENGTH             257     /* image height in pixels */
#define TIFFTAG_BITSPERSAMPLE           258     /* bits per channel (sample) */
#define TIFFTAG_COMPRESSION             259     /* data compression technique */
#define     COMPRESSION_NONE            1       /* dump mode */
#define TIFFTAG_PHOTOMETRIC             262     /* photometric interpretation */
#define     PHOTOMETRIC_MINISBLACK      1       /* min value is black */
#define TIFFTAG_STRIPOFFSETS            273     /* offsets to data strips */
#define TIFFTAG_SAMPLESPERPIXEL         277     /* samples per pixel */
#define TIFFTAG_ROWSPERSTRIP            278     /* rows per strip of data */
#define TIFFTAG_STRIPBYTECOUNTS         279     /* bytes counts for strips */
#define TIFFTAG_XRESOLUTION             282     /* pixels/resolution in x */
#define TIFFTAG_YRESOLUTION             283     /* pixels/resolution in y */
#define TIFFTAG_PLANARCONFIG            284     /* storage organization */
#define     PLANARCONFIG_CONTIG         1       /* single image plane */

typedef enum {
/*      BYTE = 1,                  8-bit unsigned integer */
/*      ASCII = 2,                 8-bit bytes w/ last byte null */
        SHORT = 3,              /* 16-bit unsigned integer */
        LONG = 4,               /* 32-bit unsigned integer */
        RATIONAL = 5            /* 64-bit fractional (numerator+denominator) */
} TIFFDataType;

static void TiffDir(fd, tag, type, count, offset)
  FILE *fd;
  int tag;
  TIFFDataType type;
  long count, offset;
{
    Putshort(tag, fd);
    Putshort((short) type, fd);
    Putlong(count, fd);
    Putlong(offset, fd);
}

void psrtiff(image, fd)
  char *image;
  FILE *fd;
{
    int *iarray = (int *) (image - (3 * sizeof(int))),
        width = iarray[0],
        height = iarray[1];
    int y;

    /* Painfully elaborated structure of the TIFF file.  Since we have
       to  write  the file in one pass without any opportunity to seek
       back to earlier addresses, we've precalculated the  offsets  of
       each  component in the file.  Note that this precludes any form
       of compression unless we wanted to  make  two  passes,  one  to
       determine  the  lengths  of the lines, and a second to actually
       emit the file. */

    int numdirent = 13;               /* Number of directory entries */
    long dirbase, axres, ayres, astrip;

    dirbase = 8;                           /* Directory base address */
    axres = dirbase + 2 + numdirent * 12L; /* Address of X resolution */
    ayres = axres + 8;                     /* Address of Y resolution */
    astrip = ayres + 8;                    /* Address of first data strip */

    Putshort(TIFF_LITTLEENDIAN, fd);  /* File type sentinel / byte order */
    Putshort(TIFF_VERSION, fd);       /* TIFF version of this file */
    Putlong(dirbase, fd);             /* Offset of first directory */

    /* Build directory items. */

#define Stl(x)  ((long) (x))

    Putshort(numdirent, fd);          /* Number of directory entries */
    TiffDir(fd, TIFFTAG_OSUBFILETYPE,    SHORT, 1L, Stl(FILETYPE_REDUCEDIMAGE));
    TiffDir(fd, TIFFTAG_IMAGEWIDTH,      SHORT, 1L, Stl(width));
    TiffDir(fd, TIFFTAG_IMAGELENGTH,     SHORT, 1L, Stl(height));
    TiffDir(fd, TIFFTAG_BITSPERSAMPLE,   SHORT, 1L, Stl(1));
    TiffDir(fd, TIFFTAG_COMPRESSION,     SHORT, 1L, Stl(COMPRESSION_NONE));
    TiffDir(fd, TIFFTAG_PHOTOMETRIC,     SHORT, 1L, Stl(PHOTOMETRIC_MINISBLACK));
    TiffDir(fd, TIFFTAG_STRIPOFFSETS,    LONG,  1L, astrip);
    TiffDir(fd, TIFFTAG_SAMPLESPERPIXEL, SHORT, 1L, Stl(1));
    TiffDir(fd, TIFFTAG_ROWSPERSTRIP,    LONG,  1L, Stl(height));
    TiffDir(fd, TIFFTAG_STRIPBYTECOUNTS, LONG,  1L, Stl(iarray[2] * height));
    TiffDir(fd, TIFFTAG_XRESOLUTION, RATIONAL,  1L, axres);
    TiffDir(fd, TIFFTAG_YRESOLUTION, RATIONAL,  1L, ayres);
    TiffDir(fd, TIFFTAG_PLANARCONFIG,    SHORT, 1L, Stl(PLANARCONFIG_CONTIG));

    /* Output the X and Y resolution as rational numbers. */

    for (y = 0; y < 2; y++) {
        Putlong(30L, fd);
        Putlong(1L, fd);
    }

    /* And finally, output the actual pixels. */

    for (y = 0; y < height; y++) {
        V fwrite(image + y * iarray[2], iarray[2], 1, fd);
    }
}

/*  Built-in test program.  */

#ifdef TestProgram
int main()
{
#define Xsize   132
#define Ysize   157         /* Odd Y size is very demanding on psrpoly() */
    char *pixels = psrinit(Xsize, Ysize, 0);
    struct spolygon sp;
    FILE *fp;
    int i;

#ifndef FlipY
    Hey!! FlipY must be defined when you build the TestProgram!!
#endif

/* Mode for opening TIFF binary file.  You may have to change this for
   your compiler. */

#ifdef unix
#define Wmode   "w+"
#else
#define Wmode   "wb+"
#endif

    /* Known-correct TIFF file output.  Embedded in C with XD.  The XD
       utility is listed on the American  Information  Exchange  under
       the symbol XD.JWALKER. */

    static unsigned char tiffcomp[2851] = {
        73,73,42,0,8,0,0,0,13,0,255,0,3,0,1,0,0,0,1,0,0,0,0,1,3,0,1,0,0,0,132,
        0,0,0,1,1,3,0,1,0,0,0,157,0,0,0,2,1,3,0,1,0,0,0,1,0,0,0,3,1,3,0,1,0,0,
        0,1,0,0,0,6,1,3,0,1,0,0,0,1,0,0,0,17,1,4,0,1,0,0,0,182,0,0,0,21,1,3,0,
        1,0,0,0,1,0,0,0,22,1,4,0,1,0,0,0,157,0,0,0,23,1,4,0,1,0,0,0,109,10,0,
        0,26,1,5,0,1,0,0,0,166,0,0,0,27,1,5,0,1,0,0,0,174,0,0,0,28,1,3,0,1,0,
        0,0,1,0,0,0,30,0,0,0,1,0,0,0,30,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,32,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,112,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,112,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,168,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,168,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,40,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,36,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,36,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,34,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,34,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,33,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,33,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,33,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,32,
        128,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,32,128,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        16,32,64,0,0,0,0,0,0,0,0,0,0,0,0,0,0,16,32,64,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,32,32,64,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,0,32,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,32,0,32,0,0,0,0,0,0,0,0,0,0,0,0,0,0,64,0,16,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,64,0,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,128,0,8,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,128,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,8,0,0,0,0,0,0,0,0,
        0,0,0,0,0,1,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,4,0,0,0,0,0,0,0,0,0,
        0,0,0,0,2,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,2,0,0,0,0,0,0,0,0,0,0,
        0,0,0,4,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
        0,0,8,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,128,0,0,0,0,0,0,0,0,0,0,
        0,0,8,0,0,0,128,0,0,0,0,0,0,0,0,0,0,0,0,16,0,0,0,64,0,0,0,0,0,0,0,0,0,
        0,0,0,16,0,0,0,64,0,0,0,0,0,0,0,0,0,0,0,0,32,0,0,0,64,0,0,0,0,0,0,0,0,
        0,0,0,0,32,0,32,0,32,0,0,0,0,0,0,0,0,0,0,64,0,32,0,0,0,32,0,16,0,0,0,
        0,0,0,0,0,0,0,64,0,32,0,16,0,0,0,0,0,0,0,0,0,0,0,0,64,0,32,0,16,0,0,0,
        0,0,0,0,0,0,0,0,0,128,0,32,0,16,0,0,0,0,0,0,0,0,0,0,0,0,128,0,32,0,8,
        0,0,0,0,0,0,0,0,0,0,0,1,0,0,32,0,8,0,0,0,0,0,0,0,0,0,0,0,1,0,0,32,0,4,
        0,0,0,0,0,0,0,0,0,0,0,1,0,0,32,0,4,0,0,0,0,0,0,0,0,0,0,0,2,0,0,32,0,2,
        0,0,0,0,0,0,0,0,0,0,0,2,0,0,32,0,2,0,0,0,0,0,0,0,0,0,0,0,4,0,0,32,0,2,
        0,0,0,0,0,0,0,0,0,0,0,4,0,0,32,0,1,0,0,0,0,0,0,0,0,0,0,0,8,0,0,32,0,1,
        0,0,0,0,0,0,0,0,0,0,0,8,0,0,32,0,0,128,0,0,0,0,0,0,0,0,0,0,8,0,0,32,0,
        0,128,0,0,0,0,0,0,0,0,0,0,16,0,0,32,0,0,128,0,0,0,0,0,0,0,0,0,0,16,0,
        0,32,0,0,64,0,0,0,0,0,0,0,0,0,0,32,0,0,32,0,0,64,0,0,0,0,0,0,0,0,0,0,
        32,0,0,32,0,0,32,0,0,0,0,0,0,0,0,0,0,64,0,0,32,0,0,32,0,0,0,0,0,0,0,0,
        0,0,64,0,0,32,0,0,16,0,0,0,0,0,0,0,0,0,0,64,0,0,32,0,0,16,0,0,0,0,0,0,
        0,0,0,0,128,0,0,32,0,0,16,0,0,0,0,0,0,0,0,0,0,128,0,0,32,0,0,8,0,0,0,
        0,0,0,0,0,0,1,0,0,0,32,0,0,8,0,0,0,0,0,0,0,0,0,1,0,0,0,32,0,0,4,0,0,0,
        0,0,0,0,0,0,1,0,0,0,32,0,0,4,0,0,0,0,0,0,0,0,0,2,0,0,0,32,0,0,4,0,0,0,
        0,0,0,0,0,0,2,0,0,0,32,0,0,2,0,0,0,0,0,0,0,0,0,4,0,0,0,32,0,0,2,0,0,0,
        0,0,0,0,0,0,4,0,0,0,32,0,0,1,0,0,0,0,0,0,0,0,0,8,0,0,0,32,0,0,1,0,0,0,
        0,0,0,0,0,0,8,0,0,0,32,0,0,0,128,0,0,0,0,0,0,0,0,8,0,0,0,32,0,0,0,128,
        0,0,0,0,0,0,0,0,16,0,0,0,32,0,0,0,128,0,0,0,0,0,0,0,0,16,0,0,0,32,0,0,
        0,64,0,0,0,0,0,0,0,0,32,0,0,0,32,0,0,0,64,0,0,0,0,0,0,0,0,32,0,0,0,32,
        0,0,0,32,0,0,0,0,0,0,0,0,64,0,0,0,32,0,0,0,32,0,0,0,0,0,0,0,0,64,8,0,
        0,32,0,0,128,32,0,0,0,0,0,0,0,0,64,12,0,0,32,0,1,128,16,0,0,0,0,0,0,0,
        0,128,6,0,0,32,0,3,0,16,0,0,0,0,0,0,0,0,128,7,0,0,32,0,7,0,8,0,0,0,0,
        0,0,0,1,0,7,0,0,32,0,7,0,8,0,0,0,0,0,0,0,1,0,7,128,0,32,0,15,0,4,0,0,
        0,0,0,0,0,2,0,3,192,0,32,0,30,0,4,0,0,0,0,0,0,0,2,0,3,224,0,32,0,62,0,
        4,0,0,0,0,0,0,0,2,0,3,240,0,32,0,126,0,2,0,0,0,0,0,0,0,4,0,1,248,0,32,
        0,252,0,2,0,0,0,0,0,0,0,4,0,1,248,0,32,0,252,0,1,0,0,0,0,0,0,0,8,0,1,
        252,0,32,1,252,0,1,0,0,0,0,0,0,0,8,0,1,254,0,32,3,252,0,1,0,0,0,0,0,0,
        0,8,0,0,255,0,32,7,248,0,0,128,0,0,0,0,0,0,16,0,0,255,128,32,15,248,0,
        0,128,0,0,0,0,0,0,16,0,0,255,192,32,31,248,0,0,64,0,0,0,0,0,0,32,0,0,
        127,224,32,63,240,0,0,64,0,0,0,0,0,0,32,0,0,127,224,32,63,240,0,0,32,
        0,0,0,0,0,0,64,0,0,127,240,32,127,240,0,0,32,0,0,0,0,0,0,64,0,0,127,248,
        32,255,240,0,0,32,0,0,0,0,0,0,64,0,0,63,252,33,255,224,0,0,16,0,0,0,0,
        0,0,128,0,0,63,254,35,255,224,0,0,16,0,0,0,0,0,0,128,0,0,63,255,39,255,
        224,0,0,8,0,0,0,0,0,1,0,0,0,63,255,39,255,224,0,0,8,0,0,0,0,0,1,0,0,0,
        31,255,175,255,192,0,0,8,0,0,0,0,0,2,0,0,0,31,255,255,255,192,0,0,4,0,
        0,0,0,0,2,0,0,0,31,255,255,255,192,0,0,4,0,0,0,0,0,2,0,0,0,15,255,255,
        255,128,0,0,2,0,0,0,0,0,4,0,0,0,15,255,255,255,128,0,0,2,0,0,0,0,0,4,
        0,0,0,15,255,255,255,128,0,0,1,0,0,0,0,0,8,0,0,0,15,255,255,255,128,0,
        0,1,0,0,0,0,0,8,0,0,0,7,255,255,255,0,0,0,1,0,0,0,0,0,16,0,0,0,7,255,
        255,255,0,0,0,0,128,0,0,0,0,16,0,0,0,7,255,255,255,0,0,0,0,128,0,0,0,
        0,16,0,0,0,3,255,255,254,0,0,0,0,64,0,0,0,0,32,0,0,0,3,255,255,254,0,
        0,0,0,64,0,0,0,0,32,0,0,0,3,255,255,254,0,0,0,0,64,0,0,0,0,64,0,0,0,3,
        255,255,254,0,0,0,0,32,0,0,0,0,64,0,0,0,1,255,255,252,0,0,0,0,32,0,0,
        0,0,64,0,0,0,1,224,0,124,0,0,0,0,16,0,0,0,0,128,0,0,0,1,224,0,124,0,0,
        0,0,16,0,0,0,0,128,0,0,0,0,240,0,248,0,0,0,0,8,0,0,0,1,0,0,0,0,0,240,
        0,248,0,0,0,0,8,0,0,0,1,0,0,0,0,0,240,0,248,0,0,0,0,8,0,0,0,2,0,0,0,0,
        0,240,0,248,0,0,0,0,4,0,0,0,2,0,0,0,0,0,120,1,240,0,0,0,0,4,0,0,0,2,0,
        0,0,0,0,120,1,240,0,0,0,0,2,0,0,0,4,0,0,0,0,0,120,1,240,0,0,0,0,2,0,0,
        0,4,0,0,0,0,0,60,1,224,0,0,0,0,2,0,0,0,8,0,0,0,0,0,60,3,224,0,0,0,0,1,
        0,0,0,8,0,0,0,0,0,60,3,224,0,0,0,0,1,0,0,0,16,0,0,0,0,0,60,3,224,0,0,
        0,0,0,128,0,0,16,0,0,0,0,0,30,3,192,0,0,0,0,0,128,0,0,16,0,0,0,0,0,30,
        7,192,0,0,0,0,0,64,0,0,32,0,0,0,0,0,30,7,192,0,0,0,0,0,64,0,0,32,0,0,
        0,0,0,15,7,128,0,0,0,0,0,64,0,0,64,0,0,0,0,0,15,7,128,0,0,0,0,0,32,0,
        0,64,0,0,0,0,0,15,15,128,0,0,0,0,0,32,0,0,128,0,0,0,0,0,15,143,128,0,
        0,0,0,0,16,0,0,128,0,0,0,0,0,7,143,0,0,0,0,0,0,16,0,0,128,0,0,0,0,0,7,
        143,0,0,0,0,0,0,16,0,1,0,0,0,0,0,0,7,159,0,0,0,0,0,0,8,0,1,0,0,0,0,0,
        0,7,223,0,0,0,0,0,0,8,0,2,0,0,0,0,0,0,3,222,0,0,0,0,0,0,4,0,2,0,0,0,0,
        0,0,3,254,0,0,0,0,0,0,4,0,2,0,0,0,0,0,0,3,254,0,0,0,0,0,0,2,0,4,0,0,0,
        0,0,0,1,252,0,0,0,0,0,0,2,0,4,0,0,0,0,0,0,1,252,0,0,0,0,0,0,2,0,8,0,0,
        0,0,0,0,1,252,0,0,0,0,0,0,1,0,8,0,0,0,0,0,0,1,252,0,0,0,0,0,0,1,0,16,
        0,0,0,0,0,0,0,248,0,0,0,0,0,0,0,128,16,0,0,0,0,0,0,0,248,0,0,0,0,0,0,
        0,128,16,0,0,0,0,0,0,0,248,0,0,0,0,0,0,0,128,32,0,0,0,0,0,0,0,112,0,0,
        0,0,0,0,0,64,32,0,0,0,0,0,0,0,112,0,0,0,0,0,0,0,64,64,0,0,0,0,0,0,0,112,
        0,0,0,0,0,0,0,32,64,0,0,0,0,0,0,0,112,0,0,0,0,0,0,0,32,128,0,0,0,0,0,
        0,0,32,0,0,0,0,0,0,0,16,128,0,0,0,0,0,0,0,32,0,0,0,0,0,0,0,16
    };

    if (pixels == NULL) {
        fprintf(stderr, "Unable to allocate pixel array.\n");
        return 1;
    }

    assert(psrlinelen(pixels) == ((Xsize + 7) / 8));

    psrvector(pixels, Xsize / 2, 0, Xsize / 2, Ysize - 1, 1);
    psrvector(pixels, Xsize / 2, Ysize - 1, 0, 0, 1);
    psrvector(pixels, Xsize / 2, Ysize - 1, Xsize - 1, 0, 1);
    psrvector(pixels, Xsize / 2, (3 * Ysize) / 4 + 2,
                      Xsize / 2, (7 * Ysize) / 8, 0);

    assert(psrgetpixel(pixels, Xsize / 4, (3 * Ysize) / 4) == 0);
    psrsetpixel(pixels, Xsize / 4, (3 * Ysize) / 4, 1);
    assert(psrgetpixel(pixels, Xsize / 4, (3 * Ysize) / 4) == 1);
    assert(psrgetpixel(pixels, Xsize / 2, (3 * Ysize) / 4) == 1);
    psrsetpixel(pixels, Xsize / 2, (3 * Ysize) / 4, 0);
    assert(psrgetpixel(pixels, Xsize / 2, (3 * Ysize) / 4) == 0);
    psrsetpixel(pixels, (3 * Xsize) / 4, (3 * Ysize) / 4, 1);

    sp.npoints = 4;
    sp.pt[0].x = Xsize / 3;
    sp.pt[0].y = Ysize / 2;
    sp.pt[1].x = Xsize / 2;
    sp.pt[1].y = 0;
    sp.pt[2].x = (2 * Xsize) / 3;
    sp.pt[2].y = Ysize / 2;
    sp.pt[3].x = Xsize / 2;
    sp.pt[3].y = Ysize / 3;
    psrpoly(pixels, &sp, 1);

    sp.npoints = 3;
    sp.pt[0].x = (9 * Xsize) / 20;
    sp.pt[0].y = Ysize / 4;
    sp.pt[1].x = Xsize / 2;
    sp.pt[1].y = (Ysize * 1) / 10;
    sp.pt[2].x = (11 * Xsize) / 20;
    sp.pt[2].y = Ysize / 4;
    psrpoly(pixels, &sp, 0);

    /* Make a TIFF file. */

    fp = fopen("pstampr.tif", Wmode);
    assert(fp != NULL);
    psrtiff(pixels, fp);

    /* Now verify that it's the same as a known-correct file. */

    rewind(fp);
    for (i = 0; i < sizeof tiffcomp; i++) {
        int fchar = getc(fp);

        if (tiffcomp[i] != fchar) {
            V fprintf(stderr, "Compare error on byte %d of TIFF file.\n",
                i);
            V fprintf(stderr, "  Expected 0x%02X, received 0x%02X\n",
                tiffcomp[i], fchar);
            return 1;
        }
    }
    if (getc(fp) != EOF) {
        V fprintf(stderr, "Unexpected garbage at end of TIFF file.\n");
    }
    fclose(fp);

#ifdef PosRast

    /* Write the bitmap to standard output in the PBM format  used  by
       Jef  Poskanzer's  PBMPLUS  raster  toolkit.   If  you  run into
       trouble, the tools in  this  toolkit  may  help  you  find  the
       problem. */

#ifdef Binary
    printf("P4\n%d %d\n", Xsize, Ysize);
    V fwrite(pixels, psrlinelen(pixels) * Ysize, 1, stdout);
#else
    printf("P1\n%d %d\n", Xsize, Ysize);
    {  int x, y, l = 0;

#ifdef FlipY
        for (y = Ysize - 1; y >= 0; y--)
#else
        for (y = 0; y < Ysize; y++)
#endif
        {
            for (x = 0; x < Xsize; x++) {
                if (++l > 70) {
                    putchar('\n');
                    l = 1;
                }
                putchar('0' + psrgetpixel(pixels, x, y));
            }
            putchar('\n');
        }
    }
#endif
#endif /* PosRast */

    psrterm(pixels);

    return 0;
}
#endif /* TestProgram */
