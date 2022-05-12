/*

             Stateless Argument Packing and Unpacking

                           by John Walker
                      http://www.fourmilab.ch/

*/

#include "vplanet.h"

/* If FOURMILAB is defined, use the private Fourmilab host key to
   encrypt stateless arguments.  Otherwise, include the default
   public key provided with the source distribution.  */

#ifdef FOURMILAB
#include "hostkey_FOURMILAB.h"
#else
#include "hostkey_PUB.h"
#endif

/*  PACKARGS  --  Pack command line and orbital elements into monolithic
                  option block, quoted for inclusion in a URL.  */

char *packargs(int argc, char *argv[])
{
    int h, i, j, l = 1, r;
    char *argblock, *cp, *hblock;
    unsigned char v;

    assert(argc < 256);

    /* Compute length of all arguments packed into string. */
    for (i = 1; i < argc; i++) {
        l += strlen(argv[i]) + 1;
    }

    /* Determine number of lines of orbital elements */
    for (j = (maxOrbitalElements - 1); j >= 0; j--) {
        if (strlen(orbitalElements[j]) > 0) {
            break;
        }
    }
    j++;

    /* Add length to encode orbital elements. */
    l += 1;
    for (i = 0; i < j; i++) {
        l += strlen(orbitalElements[i]) + 1;
    }

    l += 3;
    argblock = malloc(l);
    if (argblock == NULL) {
        fprintf(stderr, "%s: Unable to allocate %d byte argument block in packargs.\n", progpath, l);
        abort();
    }

    cp = argblock;

    r = (getppid() >> 3) ^ getpid() ^ time(NULL);
    h = 3 + (r & 0xF);
    for (i = 0; i < h; i++) {
        r ^= hostkey[17 + (r & 0x1F)];
    }
    *cp++ = r;                              /* Initial vector */
    *cp++ = 1;                              /* Format code */


    *cp++ = argc;                           /* Number of arguments */
    for (i = 1; i < argc; i++) {
        strcpy(cp, argv[i]);                /* Arguments as null terminated strings */
        cp += strlen(argv[i]) + 1;
    }

    *cp++ = j;                              /* Number of lines of orbital elements */
    for (i = 0; i < j; i++) {
        strcpy(cp, orbitalElements[i]);
        cp += strlen(orbitalElements[i]) + 1;
    }

    /* Compute and append argument block checksum. */
    v = 0;
    for (i = 0; i < (l - 1); i++) {
        v = ((v << 1) | (v >> 7)) ^ argblock[i];
    }
    *cp++ = v;
/* fprintf(stderr, "Checksum(%d) = %02X\n", l, v & 0xFF); */

    assert((cp - argblock) == l);           /* Does length match predicted length ? */

    hblock = malloc((l * 2) + 1);
    if (hblock == NULL) {
        fprintf(stderr, "%s: Unable to allocate %d byte encoded argument block in packargs.\n", progpath, (l * 2) + 1);
        abort();
    }

    h = 0;

    cp = hblock;
    sprintf(cp, "%02X", argblock[0] & 0xFF);
    cp += 2;
    v = argblock[0];
    for (i = 1; i < l; i++) {
        v = v ^ hostkey[h++] ^ argblock[i];
        h = h % (sizeof hostkey);
        sprintf(cp, "%02X", v & 0xFF);
        cp += 2;
    }
    cp++;
    assert((cp - hblock) == ((l * 2) + 1));

/* fprintf(stderr, "%s\n", hblock); */

    return hblock;
}

/*  UNPACKARGS  --  Unpack command line arguments and orbital elements
                    from monolithic option block.  */

char *unpackargs(const char *oblock)
{
    char *argblock, *cp;
    char v;
    unsigned char uv;
    int h, i, nb;

    assert(strlen(oblock) >= 9);
    argblock = malloc(nb = (strlen(oblock) / 2));
    if (argblock == NULL) {
        fprintf(stderr, "%s: Unable to allocate %d byte unpacked argument buffer.\n", progpath, nb);
        abort();
    }

    cp = argblock;
    for (i = 0; i < nb; i++) {
        int h;
        sscanf(oblock + (i * 2), "%2X", &h);
        *cp++ = h;
/* fprintf(stderr, "%02X", h); */
    }
/* fprintf(stderr, "\n"); */

    assert((cp - argblock) == nb);

    /* Decode option block with host key. */
    v = argblock[0];
    h = 0;
    for (i = 1; i < nb; i++) {
        char nv = argblock[i];
        argblock[i] = v ^ hostkey[h++] ^ argblock[i];
        h = h % (sizeof hostkey);
        v = nv;
    }

    /* Verify argument block checksum. */
    uv = 0;
    for (i = 0; i < (nb - 1); i++) {
        uv = ((uv << 1) | (uv >> 7)) ^ argblock[i];
    }
    if ((uv & 0xFF) != (argblock[nb - 1] & 0xFF)) {
        fprintf(stderr, "%s: Checksum error (%02X : %02X) in argument buffer.\n",
            progpath, uv & 0xFF, argblock[nb - 1] & 0xFF);
        abort();
    }

    return argblock;
}
