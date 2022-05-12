/*

    Scan asteroid and comet orbital elements from the format
    used by the JPL Small-Body Orbital Elements databases:
        http://ssd.jpl.nasa.gov/sb_elem.html
    We automatically figure out whether the body is an asteroid
    or comet from the format of the item.  Data appear in columns,
    but column alignment is not perfect for some comets with
    very long names, and we are forgiving in general of poor
    alignment.

        Asteroids:
68858 2002 JW           53600  3.0045080 0.11827070  11.52479 240.11931  83.87259 146.4065300 14.60 0.15 MPO49525

            Number      (always purely numeric)
            Name        (may include blanks)
            Epoch       (Modified Julian Date)
            a           Semi-major axis of the orbit, AU
            e           Eccentricity
            i           Inclination in degrees
            w           Argument of the perihelion
            Node        Longitude of the ascending node
            M           Mean anomaly
            H           Absolute magnitude
            G           Magnitude slope parameter
            Ref         Reference (may include blanks)

        Comets:
73P-B/Schwassmann-Wachmann 3          52120 0.93750307 0.69401060  11.40605 198.78365  69.92112 20010127.97736 JPL 9

            Number      (\d*[CDP](\-[A-Z])?)?
            /           Slash always precedes name
            Name        (may include blanks, and/or be followed by a parenthesised designation)
            Epoch       (Modified Julian Date)
            q           Perihelion distance, AU
            e           Eccentricity
            i           Inclination in degrees
            w           Argument of the perihelion
            Node        Longitude of the ascending node
            Tp          Time of perihelion passage (YYYYMMDD.DDDDD)
            Ref         Reference (may include blanks)

*/

#include "vplanet.h"

/*  JPLELEMENTS  --  Scan JPL elements from string and place in
                     ast_info.  Returns TRUE if valid elements
                     were stored and FALSE if an error was detected
                     in the elements.  aTracked is set to the
                     return value of the function.  */

int jplElements(const char *s)
{
    int yy, mm;
    double dd;
    int n;
    const char *c, *f, *oname;

    aTracked = FALSE;

#ifdef ELDEBUG
    fprintf(stderr, "Parsing JPL elements \"%s\"\n", s);
#endif

    /* Skip leading blanks */

    oname = s;
    while ((*oname) && isspace(*oname)) {
        oname++;
    }

    /* Due to sloppy formatting of long comet names and inconsistent column
       layout between asteroid and comet tables, we find the epoch and hence
       the end of the object name by scanning along the line looking for the
       "a" (for asteroids) or "q" (for comets) field, which is always a positive
       decimal number.  Having found it, we can then work backwards to the epoch
       and from there to the end of the object name. */

    n = 0;
    c = oname + 1;
#ifdef ELDEBUG
        fprintf(stderr, "Looking forward from %s\n", c);
#endif
    while ((n >= 0) && (*c++)) {
#ifdef ELDEBUG_STATE_TRACE
        fprintf(stderr, "State %d >>%s\n", n, c);
#endif

        switch (n) {
            case 0:
                if (*c == ' ') {
                    n = 1;
                }
                break;

            case 1:
                if (isdigit(*c)) {
                    n = 2;
                } else {
                    if (*c != ' ') {
                        n = 0;
                    }
                }
                break;

            case 2:
                if (*c == '.') {
                    n = 3;
                } else if (isdigit(*c)) {
                    n = 2;
                } else if (*c == ' ') {
                    n = 1;
                } else {
                    n = 0;
                }
                break;

            case 3:
                if (isdigit(*c)) {
                    n = -1;
                } else {
                    n = 0;
                }
                break;
        }
    }

    if (n >= 0) {
#ifdef ELDEBUG
        fprintf(stderr, "Unable to find epoch; ended in state %d at %s\n", n, c);
#endif
        return FALSE;
    }
#ifdef ELDEBUG
        fprintf(stderr, "Located a/q field at %s\n", c);
#endif

    /* Back-track to start of epoch. */

    while ((c > s) && !isspace(*c)) {
        c--;
    }
    while ((c > s) && isspace(*c)) {
        c--;
    }
    while ((c > s) && (isdigit(*c) || (*c == '-'))) {
        c--;
    }
    if (c == s) {
#ifdef ELDEBUG
        fprintf(stderr, "Unable to find epoch; failed to back over epoch field.\n");
#endif
        return FALSE;
    }
    f = c + 1;

#ifdef ELDEBUG
        fprintf(stderr, "Ended backtracking at %s\n", c);
#endif


    /* Delete trailing white space from object name and save it. */

    while (isspace(*c)) {
        c--;
    }
    memcpy(ast_info.Name, oname, (c - oname) + 1);
    ast_info.Name[(c - oname) + 1] = 0;

#ifdef ELDEBUG
        fprintf(stderr, "Object name [%s], elements: %s\n", ast_info.Name, f);
#endif

    if (((c = strchr(s, '/')) != NULL) &&
        (c > s) &&
        isalpha(c[-1])) {
            char periDate[132];

            ast_info.cometary = TRUE;
            ast_info.MagH = 0;
            ast_info.MagG = 0;
            ast_info.mAnomaly = 0;

            if (sscanf(f, "%lf %lf %lf %lf %lf %lf %s",
                    &ast_info.Epoch, &ast_info.PeriAU, &ast_info.Eccentricity,
                    &ast_info.Inclination, &ast_info.ArgP, &ast_info.LANode, periDate) == 7) {
                ast_info.Epoch += ModifiedJulian;
                sscanf((periDate + strlen(periDate)) - 10, "%02d%lf", &mm, &dd);
                periDate[strlen(periDate) - 10] = 0;
                yy = atoi(periDate);
                ast_info.PeriDate = ucttoj(yy, mm - 1, (int) dd, 0, 0, 0) + (dd - ((int) dd));
                if (ast_info.Eccentricity < 1) {
                    ast_info.SemiMajorAU = ast_info.PeriAU / (1.0 - ast_info.Eccentricity);
                } else {
                    ast_info.SemiMajorAU = 0;    /* Parabolic */
                }
                return aTracked = TRUE;
            } else {
#ifdef ELDEBUG
                fprintf(stderr, "Failed to parse JPL cometary elements\n");
#endif
            }
    } else {
            ast_info.cometary = FALSE;
            if (sscanf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
                    &ast_info.Epoch, &ast_info.SemiMajorAU, &ast_info.Eccentricity,
                    &ast_info.Inclination, &ast_info.ArgP, &ast_info.LANode, &ast_info.mAnomaly,
                    &ast_info.MagH, &ast_info.MagG) == 9) {
                ast_info.Epoch += ModifiedJulian;
                return aTracked = TRUE;
            } else {
#ifdef ELDEBUG
                fprintf(stderr, "Failed to parse JPL asteroidal elements\n");
#endif
            }
    }

    return FALSE;
}
