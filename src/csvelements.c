/*

    Scan asteroid and comet orbital elements from the compressed
    CSV format used in our object catalogues.  The formats are
    as follows (comma-separated fields broken out onto separate
    lines):

        Asteroids:

            A,
            Name/Number,
            Magnitude H,
            Magnitude G,
            Mean anomaly,
            Arg. perihelion,
            Long. node,
            Inclination,
            Eccentricity,
            Semimajor axis,
            Epoch (MJD)

        Comets:
            C,
            Name,
            Perihelion time,
            Perihelion AU,
            Eccentricity,
            Long. perihelion,
            Long. node,
            Inclination

*/

#include "vplanet.h"

static char *cp = NULL;

/*  CSVNEXT  --  Return next field from current CSV record.  Note
                 that these object catalogue records do not use
                 CSV quoting.  */

static char *csvnext(void)
{
    char *ef, *sr;

    if (cp == NULL || (*cp == 0)) {
        cp = NULL;
        return NULL;
    }

    ef = strchr(cp, ',');
    if (ef == NULL) {
        if (strlen(cp) > 0) {
            sr = cp;
            cp = NULL;
            return sr;
        }
        cp = NULL;
        return NULL;
    }
    *ef = 0;
    sr = cp;
    cp = ef + 1;
    return sr;
}

/*  CSVELEMENTS  --  Scan CSV elements from string and place in
                     ast_info.  Returns TRUE if valid elements
                     were stored and FALSE if an error was detected
                     in the elements.  aTracked is set to the
                     return value of the function.  */

int csvElements(char *s)
{
    char *f;
    char ac, ditch;
    int yy, mm;
    double epoch, dd;
/*
#define nextF() if ((f = csvnext()) == NULL) { fprintf(stderr, "CSV bailout\n"); return FALSE; } else fprintf(stderr, ">%s<\n", f)
*/
#define nextF() if ((f = csvnext()) == NULL) return FALSE

    cp = s;
    aTracked = FALSE;

    cp = s;
    nextF();
    ac = f[0];                        /* Save asteroid/comet switch */
    if (islower(ac)) {
        ac = toupper(ac);
    }

    nextF();                          /* Transcribe name */
    strlcpy(ast_info.Name, f, sizeof ast_info.Name);

    switch (ac) {

        /*  Asteroids  */

        case 'A':
            ast_info.cometary = FALSE;
            nextF();
            ast_info.MagH = atof(f);
            nextF();
            ast_info.MagG = atof(f);
            nextF();
            ast_info.mAnomaly = atof(f);
            nextF();
            ast_info.ArgP = atof(f);
            nextF();
            ast_info.LANode = atof(f);
            nextF();
            ast_info.Inclination = atof(f);
            nextF();
            ast_info.Eccentricity = atof(f);
            nextF();
            ast_info.SemiMajorAU = atof(f);
            nextF();
            ast_info.Epoch = atof(f) + ModifiedJulian;
            return aTracked = TRUE;

        /*  Comets  */

        case 'C':
            ast_info.cometary = TRUE;
            ast_info.MagH = 0;
            ast_info.MagG = 0;
            ast_info.mAnomaly = 0;
            nextF();
            if (sscanf(f, "%d%c%d%c%lf", &yy, &ditch, &mm, &ditch, &dd) < 5) {
                return FALSE;
            }
            if (yy < 100) {
                if (yy < 90) {
                    yy += 2000;
                } else {
                    yy += 1900;
                }
            }
            epoch = ucttoj(yy, mm - 1, (int) dd, 0, 0, 0);
            epoch += dd - ((int) dd);
            ast_info.PeriDate = ast_info.Epoch = epoch;
            nextF();
            ast_info.PeriAU = atof(f);
            nextF();
            ast_info.Eccentricity = atof(f);
            nextF();
            ast_info.ArgP = atof(f);
            nextF();
            ast_info.LANode = atof(f);
            nextF();
            ast_info.Inclination = atof(f);
            if (ast_info.Eccentricity < 1) {
                ast_info.SemiMajorAU = ast_info.PeriAU / (1.0 - ast_info.Eccentricity);
            } else {
                ast_info.SemiMajorAU = 0;    /* Parabolic */
            }
            return aTracked = TRUE;

        default:
            break;
    }
    return FALSE;
}
