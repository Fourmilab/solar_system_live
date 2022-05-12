/*

    Scan comet elements in the formats used by the Minor Planet Center
    (and in the Minor Planet Electronic Circulars).

    Designed and implemented by John Walker in November of 1994.
    Modified for use in Solar System Live by the same perp in
    February of 1995.  Fixed on 4 April 1997 to handle new
    style elements in which periodic comets are distinguished
    entirely by the presence of an eccentricity.

    Modified on 19 April 1997 to make the unused last line of
    elements for periodic comets (the "a =" line which gives
    the semimajor axis and period) optional.  Apparently
    elements missing that line squeaked through before the
    4 April extension to parse new style periodic elements.
    The line will, if present, be echoed to the elements box
    but no harm is done if it is not supplied.  The old code
    is kept around for reference, disabled by OBSOLETE_1997_04_19.

    PERIODIC COMET SHOEMAKER 4 (1994k)

                . . .

         T = 1994 Oct. 31.231 TT          Peri. = 196.521
         e = 0.52851                      Node  =  92.603   2000.0
         q = 2.92503 AU                   Incl. =  25.345
           a =  6.20375 AU     n = 0.063786     P =  15.45 years
              1         2         3         4         5         6         7
    01234567890123456789012345678901234567890123456789012345678901234567890

    COMET MUELLER (1994c)

                    . . .

         T = 1993 Dec.  3.990 TT          Peri. = 102.512
                                          Node  =   4.933   2000.0
         q = 1.81118 AU                   Incl. = 145.454
              1         2         3         4         5         6         7
    01234567890123456789012345678901234567890123456789012345678901234567890

*/

#include "vplanet.h"

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

/*  SCT  --  Copy string, trimming white space from start and end.  */

static void sct(char *b, const int line, const int start, const int end)
{
    char *c = s[line], sb[132];

    memcpy(sb, c + start, (end - start) + 1);
    sb[(end - start) + 1] = 0;
    while (strlen(sb) > 0 && isspace(sb[strlen(sb) - 1])) {
        sb[strlen(sb) - 1] = 0;
    }
    c = sb;
    while (*c && isspace(*c)) {
        c++;
    }
    strcpy(b, c);
}

/*  SCR  --  Copy string, replacing all sequences of white space
             with a single space character.  */

static void scr(char *b, const int line, const int start, const int end)
{
    char *c = s[line], sb[132];

    memcpy(sb, c + start, (end - start) + 1);
    sb[(end - start) + 1] = 0;

    while (*c != 0) {
        if (isspace(*c)) {
            *b++ = ' ';
            do {
                c++;
            } while (*c && isspace(*c));
        } else {
            *b++ = *c++;
        }
    }
    *b++ = 0;
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

void cometel(char *firstline)
{
    int i, iline = 0;
    static char cname[82], periau[82], perilong[82],
                eccen[82], longnode[82], inclin[82], semimaj[82], period[82];
    double xepoch = 0, epoch, dfrac;
    int periodic = -1, dint;
    char *u, *v;
    int mpc8line = FALSE;

#ifdef RSTRACE
    Trace = fopen("/t/solar.t", "w");
    fprintf(Trace, "%s\n", firstline);
#endif
#ifdef ELDEBUG
    fprintf(stderr, "Cometel processing \"%s\"\n", firstline);
#endif
    aTracked = FALSE;
    strcpy(s[0], firstline);
    if (strncmp(s[0], "COMET ", 6) == 0) {
        sct(cname, 0, 6, 80);
        periodic = FALSE;
    } else if ((u = strstr(s[0], "C/")) != NULL) {
        strlcpy(cname, u, sizeof cname);
        periodic = FALSE;
    } else if ((u = strstr(s[0], "P/")) != NULL) {
        strlcpy(cname, u, sizeof cname);
        periodic = TRUE;
    } else if (strncmp(s[0], "PERIODIC COMET ", 15) == 0) {
        sct(cname, 0, 15, 80);
        periodic = TRUE;
    } else if ((u = strstr(s[0], "COMET P/")) != NULL) {
        strlcpy(cname, u + 6, sizeof cname);
        periodic = TRUE;
     }


#ifdef ELDEBUG
    fprintf(stderr, "  Comet name \"%s\"\n", cname);
#endif

    /* Scan until we either find a line which begins with the "T ="
       specification or run into another comet header (indicating
       this circular didn't contain elements. */

    iline = 1;
    while (TRUE) {
        char t[132];

        if (!element_rs(s[iline])) {
            zapelem();
            return;
        }
        if (strncmp(s[iline], "COMET ", 6) == 0) {
            sct(cname, iline, 6, 80);
            periodic = FALSE;
            strlcpy(s[0], s[iline], 132);
            continue;
        } else if (strncmp(s[iline], "PERIODIC COMET ", 15) == 0) {
            sct(cname, iline, 15, 80);
            periodic = TRUE;
            strlcpy(s[0], s[iline], 132);
            continue;
        }
        sc(t, iline, 0, 80);
#ifdef ELDEBUG
            fprintf(stderr, "  Cometel scanning %s\n", t);
#endif

        /* Check for explicit epoch specification in MPC 8 line
           format.  If present, save in xepoch, which will override the
           default of epoch equal to the perihelion time. */

        if ((strncmp(s[iline], "Epoch ", 6) == 0) &&
            ((u = strstr(t, "TT=JDT")) != NULL) &&
            isdigit(u[6])) {
            xepoch = atof(u + 6);
#ifdef ELDEBUG
            fprintf(stderr, " Explicit epoch %.4f\n", xepoch);
#endif
            iline++;
            continue;
        }

        /* Check for perihelion time in MPC 8 line format.  */

        if ((strncmp(s[iline], "T ", 2) == 0) && isdigit(t[1]) &&
            (strstr(s[iline], " TT") != NULL)) {
            s[iline][1] = '=';
            sc(t, iline, 0, 80);
            mpc8line = TRUE;
#ifdef ELDEBUG
            fprintf(stderr, "  Ahhhh.  This appears to be MPC 8-line format.\n");
            fprintf(stderr, "    Rescanning fudged date %s\n", t);
#endif
        }

        if (strncmp(t, "T=", 2) == 0 && isdigit(t[2])) {
            char y[20], m[20], daf[20];
            static char mname[] = "janfebmaraprmayjunjulaugsepoctnovdec";

            sscanf(u = strchr(s[iline], '=') + 1, "%s %s", y, m);
            for (i = 0; i < 12; i++) {
                if (strncasecmp(m, mname + i * 3, 3) == 0) {
                    sprintf(m, "%d", i + 1);
                    break;
                }
            }
            if (i >= 12) {
                continue;             /* Bogus month name */
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
            epoch = ucttoj(atoi(y), atoi(m) - 1, (dint = atoi(daf)), 0, 0, 0);
            if ((dfrac = atof(daf)) != ((double) dint)) {
                epoch += dfrac - dint;
            }

#ifdef ELDEBUG
            fprintf(stderr, "  Epoch: %.4f\n", epoch);
#endif

            if (mpc8line) {

                s[iline][1] = ' ';      /* Restore original syntax for next time */

                /*  Current MPC 8 line format.

                    Parabolic format:

    C/2005 M1 (Christensen)
    T 2005 Jan. 24.728 TT                                   MPC
    q   2.94607              (2000.0)            P               Q
                       Peri.  140.756       +0.258376       +0.961898
                       Node   143.962       -0.920172       +0.273233
    e   1.0            Incl.    8.741       -0.294152       -0.009823
    From 14 observations 2005 June 17-19.
              1         2         3         4         5         6         7
    01234567890123456789012345678901234567890123456789012345678901234567890

                    Periodic format:

    P/2005 L4 (Christensen)
    T 2005 Aug. 24.033 TT                                   MPC
    q   2.36711              (2000.0)            P               Q
    n   0.118160       Peri.   24.531       +0.606258       +0.742788
    a   4.11298        Node   284.072       -0.769436       +0.457542
    e   0.42448        Incl.   17.032       -0.201048       +0.488795
    P   8.34
    From 34 observations 2005 June 3-18.
              1         2         3         4         5         6         7
    01234567890123456789012345678901234567890123456789012345678901234567890

                */

                /*  Line 2 (q)  */
                iline++;
                if (!element_rs(s[iline])) {
                    zapelem();
                    return;
                }
                scr(t, iline, 0, 80);
#ifdef ELDEBUG
                fprintf(stderr, "  MPC 8/2: %s\n", t);
#endif

                v = strstr(t, "q ");
                if (v == NULL) {
                    zapelem();
                    return;
                }
                sscanf(v + 2, "%s", periau);

                /* Line 3 ([n] Peri.)  */

                iline++;
                if (!element_rs(s[iline])) {
                    zapelem();
                    return;
                }
                scr(t, iline, 0, 80);
#ifdef ELDEBUG
                fprintf(stderr, "  MPC 8/3: %s\n", t);
#endif
                v = strstr(t, "Peri. ");
                if (v == NULL) {
                    zapelem();
                    return;
                }
                sscanf(v + 5, "%s", perilong);

                /* Line 4 ([a] Node)  */

                iline++;
                if (!element_rs(s[iline])) {
                    zapelem();
                    return;
                }
                scr(t, iline, 0, 80);
#ifdef ELDEBUG
                fprintf(stderr, "  MPC 8/4: %s\n", t);
#endif
                v = strstr(t, "Node ");
                if (v == NULL) {
                    zapelem();
                    return;
                }
                sscanf(v + 5, "%s", longnode);

                /* Line 5 (e Incl.)  */

                iline++;
                if (!element_rs(s[iline])) {
                    zapelem();
                    return;
                }
                scr(t, iline, 0, 80);
#ifdef ELDEBUG
                fprintf(stderr, "  MPC 8/5: %s\n", t);
#endif
                v = strstr(t, "e ");
                if (v == NULL) {
                    zapelem();
                    return;
                }
                sscanf(v + 2, "%s", eccen);

                v = strstr(t, "Incl. ");
                if (v == NULL) {
                    zapelem();
                    return;
                }
                sscanf(v + 5, "%s", inclin);

                /* Read up rest of elements, up to 8 lines, so they're
                   echoed back to the requester. */

                for (i = iline + 1; i < 8; i++) {
                    if (!element_rs(s[i])) {
                        break;
                    }
#ifdef ELDEBUG
                    fprintf(stderr, "  MPC 8/%d: %s\n", i, s[i]);
#endif
                }

            } else {

                /* Old informal MPC format. */

                v = strchr(u, '=');
                if (v != NULL) {
                    sc(perilong, 1, (v + 1) - s[1], 80);
                } else {
                    continue;
                }

                /* Process second line. */

                if (!element_rs(s[2])) {
                    zapelem();
                    return;
                }
                sc(t, 2, 0, 80);
                if (strncmp(t, "e=", 2) == 0) {
                    sscanf(u = strchr(s[2], '=') + 1, "%s", eccen);
                    v = strchr(u, '=');
                    periodic = TRUE;
                } else {
                    strcpy(eccen, "1");         /* Parabolic */
                    v = strchr(s[2], '=');
                }
                if (v != NULL) {
                    sscanf(v + 1, "%s", longnode);
                } else {
                    continue;
                }

                /* Process third line. */

                if (!element_rs(s[3])) {
                    zapelem();
                    return;
                }
                u = strchr(s[3], '=');
                if (u != NULL) {
                    v = strchr(u + 1, '=');
                }
                if (u != NULL && v != NULL) {
                    sscanf(u + 1, "%s", periau);
                    sscanf(v + 1, "%s", inclin);
                } else {
                    continue;
                }

                /* Process fourth line if comet is periodic. */

                strcpy(semimaj, "");
                strcpy(period, "Parabolic");
#ifdef OBSOLETE_1997_04_19
                if (periodic) {
                    if (!element_rs(s[4])) {
                        zapelem();
                        return;
                    }
                    u = strchr(s[4], '=');
                    if (u != NULL) {
                        v = strchr(u + 1, '=');
                        if (v != NULL) {
                            v = strchr(v + 1, '=');
                        }
                    }
                    if (u != NULL && v != NULL) {
                        char t1[30], t2[30];

                        sscanf(u + 1, "%s %s", t1, t2);
                        sprintf(semimaj, "%s %s", t1, t2);
                        sscanf(v + 1, "%s %s", t1, t2);
                        sprintf(period, "%s %s", t1, t2);
                    } else {
                        continue;
                    }
                }
#else
                if (periodic) {
                    /* Read unused semimajor axis and period line, if
                       present, so it doesn't disappear from the elements
                       box. */
                    (void) element_rs(s[4]);
                }
#endif
            }
            break;
        }
    }

    /* Store information in ast_info structure. */

    ast_info.cometary = TRUE;
    strlcpy(ast_info.Name, cname, sizeof ast_info.Name);
    ast_info.MagH = 0;
    ast_info.MagG = 0;
    ast_info.mAnomaly = 0;
    ast_info.ArgP = atof(perilong);
    ast_info.LANode = atof(longnode);
    ast_info.Inclination = atof(inclin);
    ast_info.Eccentricity = atof(eccen);
    ast_info.PeriAU = atof(periau);
    ast_info.PeriDate = ast_info.Epoch = epoch;
    if (xepoch != 0) {
        ast_info.Epoch = xepoch;
    }
    if (ast_info.Eccentricity < 1.0) {
        ast_info.SemiMajorAU = ast_info.PeriAU / (1.0 - ast_info.Eccentricity);
    } else {
        ast_info.SemiMajorAU = 6.5;     /* Limit parabolics to Jupiter's orbit */
    }

    aTracked = TRUE;
    return;
}
