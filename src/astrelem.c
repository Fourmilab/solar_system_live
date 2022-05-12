/*

    Scan asteroid elements in the 8-line format used by the Minor
    Planet Center (and in the Minor Planet Electronic Circulars).

    Designed and implemented by John Walker in October of 1994.
    Modified for use in Solar System Live by the same perp in
    February of 1995.

     (4646) Kwee
    Epoch 1993 Jan. 13.0 TT = JDT 2449000.5  (M-P)          Oishi
    M 220.20030              (2000.0)            P               Q
    n   0.26003234     Peri.  350.28914     +0.79481320     +0.60665475
    a   2.4309973      Node   332.34432     -0.55450494     +0.71559595
    e   0.1921938      Incl.    1.92064     -0.24656897     +0.34625491
    P   3.79           H   14.0           G   0.15
    From observations at 4 oppositions, 1960-1990.  Ref. MPC 17195.
              1         2         3         4         5         6         7
    01234567890123456789012345678901234567890123456789012345678901234567890

    Elements for comets (identified by the content of the first line), are
    handed off to cometel.c, which parses all the myriad ways in which they
    are and have been published.

*/

#include "vplanet.h"

static int alreadyLoaded, aloadi;

/*  RS  --  Read orbital elements into the specified string buffer,
            which must be at least 132 characters in length.  If
            alreadyLoaded is set, the elements are retrieved from
            the savedOrbitalElements array; otherwise they are
            read from standard input or parsed from the data block
            saved from the WWW_elements environment variable set
            when uncgi() processed the query string.  */

int element_rs(char *s)
{
    if (alreadyLoaded) {
        if (savedOrbitalElements[aloadi][0] == 0) {
            return FALSE;
        }
        strlcpy(s, savedOrbitalElements[aloadi++], 132);
    } else {
        if (elementSpecification != NULL) {
            char *e = strchr(elementSpecification, '\n');
            char c;
            if (e == NULL) {
                return FALSE;
            }
            c = e[1];           /* Save character after new line */
            e[1] = EOS;         /* Set string terminator after new line */
            /* Protect against buffer overflow by lopping off any material
               in the string after character 131. */
            if (strlen(elementSpecification) > 131) {
                strcpy(elementSpecification + 130, "\n");
            }
            strlcpy(s, elementSpecification, 132); /* Return this line to caller */
            e[1] = c;           /* Restore start of next line */
            elementSpecification = e + 1;       /* Advance to next line */
        } else {
            if (fgets(s, 132, stdin) == NULL) {
                return FALSE;
            }
        }
    }
    while (strlen(s) > 0 && s[strlen(s) - 1] < ' ') {
        s[strlen(s) - 1] = 0;
    }
    return TRUE;
}

#define s   orbitalElements

/*  SC  --  Copy string, deleting all white space.  */

static void sc(char *b, const int line, const int start, const int end) {
    char *c = s[line];
    int i, n;

    for (i = 0, n = start; n <= end; n++) {
        if (c[n] == 0) {
            break;
        }
        if (!isspace(c[n])) {
            b[i++] = c[n];
        }
    }
    b[i] = 0;
}

/*  ZAPELEM  -- Delete invalid elements from text buffer.  */

static void zapelem(void)
{
    int i;

    for (i = 1; i < 8; i++) {
        s[i][0] = EOS;
    }
    strcpy(s[0], "* Specified elements invalid *");
}

void astrelem(const int preloaded)
{
    int i, p, recov, stealthComet = FALSE;
#define IBN 4
    char ib[IBN][132];
    char anum[82], aname[82], epoch[82], meanan[82], semimaj[82],
         eccen[82], argperi[82], longnode[82], inclin[82], magh[82], magg[82],
         saname[82];
    double cepoch = 0;
    char *u, *v;

    alreadyLoaded = preloaded;
    aloadi = 0;

        aTracked = FALSE;
        p = 0;
        while (TRUE) {
#define ibm(x)  ib[(p + (IBN - (x))) % IBN]
            if (!element_rs(ib[p])) {
                return;
            }

            /*  Check for elements in CSV format.  If so, invoke
                the special scanner for them.  */

            if (ib[p][0] == '=') {
                strlcpy(orbitalElements[0], ib[p], 132);
                if (!csvElements(ib[p] + 1)) {
                    zapelem();
                }
                return;
            }

            /*  Check for elements in JPL format.  Our test is somewhat tacky--if
                the line is more than 80 characters long, it's deemed JPL format.
                For the moment, I know of no potential ambiguity, but should problems
                occur, we can always handle the failure to parse in JPL format by
                trying other formats (since JPL elements are single-line and
                attempting to parse them won't "eat" additional lines from the
                elements input stream.  */

            if (strlen(ib[p]) > 80) {
                strlcpy(orbitalElements[0], ib[p], 132);
                if (!jplElements(ib[p])) {
                    zapelem();
                }
                return;
            }

            if ((strstr(ib[p], "COMET") != NULL) ||
                (strstr(ib[p], "P/") != NULL) ||            /* MPC 8 line periodic comet */
                (strstr(ib[p], "C/") != NULL)) {            /* MPC 8 line parabolic comet */
                cometel(ib[p]);
                return;
            }

#ifdef ELDEBUG
            fprintf(stderr, "Astrelem processing \"%s\"\n", ib[p]);
#endif

            if (strncmp(ib[p], "Epoch ", 6) == 0 &&
                strstr(ib[p], "TT") != NULL &&
                strchr(ib[p], '=') != NULL &&
                strstr(ib[p], "JDT") != NULL) {

                /* MPECs for recovered objects contain an "Id." line after the name
                   and before the Epoch, instead of a "From" line at the end of the
                   elements.  If this is the case with these elements, move the ID
                   to the end of the elements and shift everything else up one line. */

                recov = strncmp(ibm(1), "Id", 2) == 0;

                strlcpy(s[0], ibm(recov ? 2 : 1), 132);
                strlcpy(s[1], ib[p], 132);
                for (i = 0; i < 5; i++) {
                    if (!element_rs(s[i + 2])) {
                        zapelem();
                        return;
                    }
                    if (i == 1 && strncmp(s[2], "T ", 2) == 0 &&
                        strncmp(s[3], "q  ", 3) == 0) {
                        i--;
                        continue;
                    }
                }
                p = (p + 1) % IBN;
                break;
            }
            p = (p + 1) % IBN;
        }

        /*  Scan asteroid number, if one has been assigned. */

        strlcpy(saname, s[0], sizeof saname);
        /*  Some MPECs contain additional information right justified in
            the first line, for example "Earth MOID = 0.0003 AU" for the
            Earth-grazer 2009 DD45.  We test whether the name record contains
            five or more embedded spaces and if so truncate the record at the
            start of the string of spaces. */
        if ((u = strstr(saname, "     ")) != NULL) {
            *u = 0;
        }
        if ((u = strchr(saname, '(')) != NULL && (v = strchr(u + 1, ')')) != NULL) {
            u += 1;
            *v++ = 0;
            strcpy(anum, u);
        } else {
            v = saname;
            strlcpy(anum, "", sizeof anum);
        }

        /*  Scan name or other designation */

        while (isspace(*v)) {
            v++;
        }
        while (strlen(v) > 0 && isspace(v[strlen(v) - 1])) {
            v[strlen(v) - 1] = 0;
        }
        strlcpy(aname, v, sizeof aname);

        /* Epoch isn't always precisely aligned in MPC postings.
           Allow for some slack. */

        v = strstr(s[1], "JDT");
        if (v == NULL) {
            strlcpy(epoch, "0", sizeof epoch);  /* Zero should get somebody's attention */
        } else {
            strlcpy(epoch, v + 3, sizeof epoch);
        }

        /* Elements for some periodic comets are posted as
           if they were asteroids, distinguished only by the
           use of a "T " line for perihelion date instead of
           a "M " for mean anomaly; we call these "stealth
           comets".  Recognise them and set everything correctly. */

        if (strncmp(s[2], "T ", 2) == 0) {
            char y[20], m[20], daf[20];
            static char mname[] = "janfebmaraprmayjunjulaugsepoctnovdec";
            int dint;
            double dfrac;

            stealthComet = TRUE;

            /* Undo the wreckage done by the asteroid name
               and number scanner above. */

            u = s[0];
            while (*u && isspace(*u)) {
                u++;
            }
            strcpy(aname, u);

            sscanf(u = (s[2] + 2), "%s %s", y, m);
            for (i = 0; i < 12; i++) {
                if (strncasecmp(m, mname + i * 3, 3) == 0) {
                    sprintf(m, "%d", i + 1);
                    break;
                }
            }
            if (i >= 12) {
                return;               /* Bogus month name */
            }

            /* Tiptoe up to the start of the year, since in the
               case of "Sept.199x" the month runs into the first
               digit of the year.  This code is resilient in case
               a space is added later. */

            while (*u && isspace(*u)) {
                u++;
            }
            while (*u && !isspace(*u)) {
                u++;
            }
            while (*u && isspace(*u)) {
                u++;
            }
            while (*u && !isdigit(*u)) {
                u++;
            }
            sscanf(u, "%s", daf);
            cepoch = ucttoj(atoi(y), atoi(m) - 1, (dint = atoi(daf)), 0, 0, 0);
            if ((dfrac = atof(daf)) != ((double) dint)) {
                cepoch += dfrac - dint;
            }
        } else {
            sc(meanan, 2, 2, 15);
        }
        sc(semimaj, 4, 2, 15);
        sc(eccen, 5, 2, 15);
        sc(argperi, 3, 26, 35);
        sc(longnode, 4, 26, 35);
        sc(inclin, 5, 26, 35);
        sc(magg, 6, 39, 48);
        sc(magh, 6, 20, 30);

        /* Mean anomaly is blank in some MPECs, meaning zero.
           Home Planet correctly interprets the blank field as zero,
           but it might confuse other programs.  Replace it with 0.0 if
           blank. */

        if (strlen(meanan) == 0) {
            strlcpy(meanan, "0.0", sizeof meanan);
        }

        /* Store information in ast_info structure. */

        ast_info.cometary = stealthComet;
        if (stealthComet) {
            strlcpy(ast_info.Name, aname, sizeof ast_info.Name);
        } else {
            /* Paranoia to truncate cleanly rather than buffer overflow
               in the case of absurdly long name and/or number */
            char namebuffy[(sizeof aname) + (sizeof anum) + 2];
            sprintf(namebuffy, "%s %s", aname, anum);
            strlcpy(ast_info.Name, namebuffy, sizeof ast_info.Name);
        }
        if (stealthComet) {
            ast_info.MagH = ast_info.MagG = 0.0;
        } else {
            ast_info.MagH = atof(magh);
            ast_info.MagG = atof(magg);
        }
        ast_info.ArgP = atof(argperi);
        ast_info.LANode = atof(longnode);
        ast_info.Inclination = atof(inclin);
        ast_info.Eccentricity = atof(eccen);
        ast_info.SemiMajorAU = atof(semimaj);
        if (stealthComet) {
            ast_info.PeriDate = cepoch;
            ast_info.PeriAU = ast_info.SemiMajorAU * (1 - ast_info.Eccentricity);
        } else {
            ast_info.mAnomaly = atof(meanan);
        }
        ast_info.Epoch = atof(epoch);
        aTracked = TRUE;
        return;
}

/*  EJD  --  Edit Julian Date.  */

static char *ejd(const double j)
{
    long yy;
    int mm, dd, h, m, s;
    static char t[80];

    jhms(j, &h, &m, &s);
    jyear(j, &yy, &mm, &dd);
    sprintf(t, "%.5f (%ld/%d/%d %d:%02d:%02d)", j, yy, mm, dd, h, m, s);
    return t;
}

/*  PRINT_ELEMENTS  --  Print asteroid or comet elements.  */

void print_elements(FILE *fp)
{
    fprintf(fp, "<h3>Elements for %s (%s):</h3>\n<blockquote>\n",
        escape_html_content(ast_info.Name),
        ast_info.cometary ? "cometary" : "asteroidal");
    fprintf(fp, "Epoch: %s<br />\n", ejd(ast_info.Epoch));
    fprintf(fp, "o (argument of perihelion): %.5f<br />\n", ast_info.ArgP);
    fprintf(fp, "O (longitude of ascending node): %.5f<br />\n", ast_info.LANode);
    fprintf(fp, "i (inclination): %.5f<br />\n", ast_info.Inclination);
    fprintf(fp, "e (eccentricity): %.7f<br />\n", ast_info.Eccentricity);

    if (ast_info.Eccentricity < 1.0) {
        fprintf(fp, "a  (semimajor axis): %.7f AU<br />\n", ast_info.SemiMajorAU);
    }

    if (ast_info.cometary) {
        fprintf(fp, "T  (perihelion date): %s<br />\n", ejd(ast_info.PeriDate));
        fprintf(fp, "q  (perihelion distance): %.5f AU<br />\n", ast_info.PeriAU);
    } else {
        fprintf(fp, "M  (mean anomaly): %.5f<br />\n", ast_info.mAnomaly);
        fprintf(fp, "T  (perihelion date): %s<br />\n",
                ejd(ast_info.Epoch - sqrt(ast_info.SemiMajorAU * ast_info.SemiMajorAU *
                    ast_info.SemiMajorAU) * ast_info.mAnomaly * 365.24219879 / 360));
        fprintf(fp, "q  (perihelion distance): %.5f AU<br />\n",
                ast_info.SemiMajorAU * (1 - ast_info.Eccentricity));
        fprintf(fp, "H  (mean visual magnitude): %.2f<br />\n", ast_info.MagH);
        fprintf(fp, "G  (magnitude slope factor): %.2f<br />\n", ast_info.MagG);
    }
    fprintf(fp, "</blockquote>\n");
}
