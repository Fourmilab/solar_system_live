/*

           Planet information panel

*/

#include "vplanet.h"

/*  CALCPLANET  --  Calculate planetary positions and altitude and azimuth from
                    viewer's position.  */

void calcPlanet(double jd)
{
    int i;
    double igmst, latsin, latcos;

    planets(jd, 0xFFFF);
    igmst = gmst(jd);
    latsin = sin(dtr(siteLat));
    latcos = cos(dtr(siteLat));
    for (i = 0; i <= (aTracked ? 10 : 9); i++) {
        planet_info[i].lha = dtr(fixangle((igmst * 15) - siteLon - planet_info[i].ra));
        planet_info[i].az = rtd(atan2(sin(planet_info[i].lha), cos(planet_info[i].lha) * latsin -
                                tan(dtr(planet_info[i].dec)) * latcos));
        planet_info[i].alt = rtd(asin(latsin * sin(dtr(planet_info[i].dec)) +
                                latcos * cos(dtr(planet_info[i].dec)) * cos(planet_info[i].lha)));
    }
}

/*  REFRESHLASTALT  --  Refresh last altitude to determine rising from setting.  */

static double lastalt[10] = {-9999};        /* Last calculated altitudes */
static double lastsaltime = -9999;          /* Time of last altitude calculation */

static void refreshLastSalt(double jd)
{
    int i;

    calcPlanet(lastsaltime = (jd - (1.0 / 48.0)));
    for (i = 0; i <= (aTracked ? 10 : 9); i++) {
        lastalt[i] = planet_info[i].alt;
    }
}

/*  BUILDPLANETS  --  Build table of planetary positions.  */

void buildPlanets(const double jd)
{
    if (lastsaltime > jd || ((jd - lastsaltime) > (1.0 / 24.0))) {
        refreshLastSalt(jd);
    }
    calcPlanet(jd);
    lastsaltime = jd;
}

/*  UPDATEPLANET  --  Generate table of planetary positions.  */

void updatePlanet(double jd, int normal, FILE *ofile, char *obSite)
{
    int i;
    char tbuf[80], obuf[90];
    char *plump;
    double d, md, mangdia, suangdia, setalt, upalt;
    static char *plnames[] = { "Sun", "Mercury", "Venus", "Moon",
                               "Mars", "Jupiter", "Saturn", "Uranus",
                               "Neptune", "Pluto", "?" };

#define Set(x) memcpy(obuf + (x) + (8 - strlen(tbuf)), tbuf, strlen(tbuf))

    plnames[10] = ast_info.Name;
    fprintf(ofile, "<pre>\n");
    fprintf(ofile,
    "<b>              Right                   Distance    From %s:\n", obSite);
    fprintf(ofile,
    "            Ascension    Declination      (AU)   Altitude Azimuth</b>\n");

    phase(jd, &d, &d, &md, &mangdia, &d, &suangdia);
    assert(lastsaltime != -9999);

    for (i = 0; i <= (aTracked ? 10 : 9); i++) {
        /* If outside the known calculable orbit of Pluto, hide the bogus numbers */
        if (i == 9 && planet_info[i].hrv <= 0) {
            continue;
        }

        memset(obuf, ' ', (sizeof obuf) - 1);
        obuf[(sizeof obuf) - 1] = 0;

        /* If we're about to process the name of a tracked asteroid or
           comet, escape any HTML meta-character to thwart attempt to inject
           HTML into the reply document.  If escape_html_content replaces the
           argument string with a quoted one, we'll orphan that buffer, but
           it isn't worth worrying about since this function is only called
           once in Solar System Live and the program exits shortly thereafer. */

        if (i == 10) {
            plnames[i] = escape_html_content(plnames[i]);
        }

        /* If the name of the object doesn't fit in the space
           available, try to break it onto the next line in a
           semi-intelligent fashion. */

        if (strlen(plnames[i]) > 10) {
            char *cp, c;

            c = plnames[i][10];
            plnames[i][10] = EOS;

            cp = NULL;
            if (strrchr(plnames[i], ' ')) {
                cp = strrchr(plnames[i], ' ');
            }
            if (strrchr(plnames[i], '/') && strrchr(plnames[i], '/') > cp) {
                cp = strrchr(plnames[i], '/');
            }
            if (strrchr(plnames[i], '.') && strrchr(plnames[i], '.') > cp) {
                cp = strrchr(plnames[i], '.');
            }
            if (strrchr(plnames[i], '-') && strrchr(plnames[i], '-') > cp) {
                cp = strrchr(plnames[i], '-');
            }
            if (strrchr(plnames[i], '(') && strrchr(plnames[i], '(') > cp) {
                cp = strrchr(plnames[i], '(');
            }
            if (strrchr(plnames[i], ')') && strrchr(plnames[i], ')') > cp) {
                cp = strrchr(plnames[i], ')');
            }

            if (cp != NULL) {
                char a = cp[1];

                plump = cp + 1;
                cp[1] = EOS;
                sprintf(tbuf, "<b>%s</b>", plnames[i]);
                cp[1] = a;
            } else {
                sprintf(tbuf, "<b>%s</b>", plnames[i]);
                plump = plnames[i] + 10;
            }

            plnames[i][10] = c;

            while (*plump && isspace(*plump)) {
                plump++;
            }
            if (strlen(plump) == 0) {
                plump = NULL;
            }
        } else {
            sprintf(tbuf, "<b>%s</b>", plnames[i]);
            plump = NULL;
        }
        memcpy(obuf, tbuf, strlen(tbuf));

        {   double m = fixangle(planet_info[i].ra) / 15;

            sprintf(tbuf, "%2dh %2dm %2ds", (int) m, ((int) (m * 60)) % 60,
                    (int) fmod(m * 3600, 60.0));
        }
        Set(22);
        {   double m = abs(planet_info[i].dec);

            sprintf(tbuf, "%c%d\260 %4.1f'", planet_info[i].dec < 0 ? '-' : '+',
                    (int) m, fmod(m * 60, 60));
        }
        Set(35);
        if (i == 3) {
            sprintf(tbuf, "%.1f ER", md / EarthRad);
        } else {
            sprintf(tbuf, "%.3f", planet_info[i].dist);
        }
        Set(45);
        sprintf(tbuf, "%.3f", planet_info[i].alt);
        Set(55);
        sprintf(tbuf, "%.3f", planet_info[i].az);
        Set(64);
        setalt = -0.5667;
        upalt = 0;

        /* (Approximately) compensate for the effect of refraction by
           the Earth's atmosphere so that "rising" and "setting" appears
           at the time the object actually becomes visible or disappears
           rather than considering pure geometry of touching the
           horizon.  These calculations follow the guidance of Meeus in
           chapter 14 of Astronomical Algorithms.

           Note that the Sunrise is considered to be the first appearance
           of the limb above the horizon.

           Moonrise is even more complicated.  First of all, in addition to
           refraction we must also compensate for the parallax of the
           Moon which, in turn, varies due to the eccentricity of the Moon's
           orbit.  Further (departing from normal astronomical convention),
           if we wish to also time moonrise from first appearance of the
           limb rather than the centre of the disc (we adopt this nonstandard
           criterion to prevent squawks from folks who see the Moon rising
           out their window before we've labeled it as "rising"), we must add
           the semidiameter of the Moon to our refraction correction.  That,
           of course, *also* changes from perigee to apogee, so there goes
           another correction into the stew.  Ain't this fun?

           Then there's the question of how long "rising" and "setting"
           lasts.  For a planet it's almost instantaneous, but we
           compromise and consider rising and setting to be the interval
           between refraction-corrected rise/set and geometric (about half
           a degree to the horizon.  For the Sun and Moon, we must adjust the
           end of rise/set times for the diameter of the disc.  */

        if (i == 0) {               /* Sun */
            setalt -= suangdia / 2;
            upalt = setalt + suangdia;
        } else if (i == 3) {        /* Moon */
#define msmax       384401.0                    /* Semi-major axis of Moon's orbit in km */
#define mparallax   0.9507                      /* Parallax at distance msmax from Earth */
            setalt += 0.7275 * (mparallax * (md / msmax));
            setalt -= mangdia / 2;
            upalt = setalt + mangdia;
        }
        if ((planet_info[i].alt > setalt) && (planet_info[i].alt < upalt)) {
            strcpy(tbuf, (planet_info[i].alt > lastalt[i]) ? "Rising" : "Setting");
        } else {
            if (planet_info[i].alt > upalt) {
                strcpy(tbuf, (abs(planet_info[i].lha) < dtr(0.5)) ?
                         "Transit" : "Up");
            } else {
                strcpy(tbuf, "Set");
            }
        }
        strcpy(obuf + 73, tbuf);                /* Also truncates to proper length */
        fprintf(ofile, "%s\n", obuf);
        if (plump != NULL) {
            fprintf(ofile, "<b>%s</b>\n", plump);
        }
        lastalt[i] = planet_info[i].alt;        /* Save last altitude */
    }
    fprintf(ofile, "</pre>\n");
}
