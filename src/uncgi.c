/*
 * @(#)uncgi.c  1.37 11/01/03
 *
 * Unescape all the fields in a form and stick them in the environment
 * so they can be used without awful machinations.
 *
 * Call with an ACTION such as:
 *      http://foo.bar.com/cgi-bin/uncgi/myscript/extra/path/stuff
 *
 * Uncgi will run "myscript" from the SCRIPT_BIN directory (configured at
 * compile time) and set PATH_INFO to "/extra/path/stuff".
 *
 * Environment variable names are "WWW_" plus the field name.
 *
 * Copyright 1994, Steven Grimm <koreth@midwinter.com>.
 *
 * Permission is granted to redistribute freely and use for any purpose,
 * commercial or private, so long as this copyright notice is retained
 * and the source code is included free of charge with any binary
 * distributions.
 */

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>

#define LIBRARY

#define VERSION "1.10"

#define NOPARAMS void

#define PREFIX  "WWW_"
#define ishex(x) (((x) >= '0' && (x) <= '9') || ((x) >= 'a' && (x) <= 'f') || \
                  ((x) >= 'A' && (x) <= 'F'))


/* Define some directory macros for systems that don't have them. */
#ifndef S_IFDIR
# define S_IFDIR 040000
#endif

#ifndef S_ISDIR
# define S_ISDIR(x) ((x & 0170000) == S_IFDIR)
#endif

/*
 * Convert two hex digits to a value.
 */
static int htoi(const unsigned char *s)
{
        int     value;
        char    c;

        c = s[0];
        if (isupper(c))
                c = tolower(c);
        value = (c >= '0' && c <= '9' ? c - '0' : c - 'a' + 10) * 16;

        c = s[1];
        if (isupper(c))
                c = tolower(c);
        value += c >= '0' && c <= '9' ? c - '0' : c - 'a' + 10;

        return (value);
}

/*
 * Get rid of all the URL escaping in a string.  Modify it in place, since
 * the result will always be equal in length or smaller.
 */
static void url_unescape(unsigned char *str)
{
        unsigned char   *dest = str;

        while (str[0])
        {
                if (str[0] == '+')
                        dest[0] = ' ';
                else if (str[0] == '%' && ishex(str[1]) && ishex(str[2]))
                {
                        dest[0] = (unsigned char) htoi(str + 1);
                        str += 2;
                }
                else
                        dest[0] = str[0];

                str++;
                dest++;
        }

        dest[0] = '\0';
}

/*
 * Print the start of an error message.
 */
static void http_head(NOPARAMS)
{
        puts("Content-type: text/html\n");
        puts("<html><head><title>Error!</title></head><body>");
        puts("<h1>An error has occurred while processing your request.</h1>");
}

/*
 * Print a standard tagline for Un-CGI error messages.
 */
static void uncgi_tag(NOPARAMS)
{
        printf("<pre>\n\n\n\n\n\n\n</pre><cite><a href=\"http://www.");
        printf("midwinter.com/~koreth/uncgi.html\">Un-CGI %s", VERSION);
        puts("</a></cite>\n</body></html>");
}

/*
 * Print an HTML error string with the right HTTP header and exit.
 */
static void html_perror(const char *str)
{
        http_head();

        puts("<p>The following error was encountered while processing");
        puts("your query.");
        puts("<blockquote>");
        printf("%s: %s\n", str, strerror(errno));
        puts("</blockquote>");
        puts("<p>Please try again later.");
        uncgi_tag();

        exit(1);
}

#ifndef LIBRARY
static void ad_path_error(const char *msg)
{
        http_head();
        printf("The path to the script was invalid: %s\n", msg);

        exit(1);
}


static void no_suid_error(void)
{
        http_head();
        printf("Execution of setuid scripts isn't permitted.\n");

        exit(1);
}
#endif /* LIBRARY */


/*
 * Stuff a URL-unescaped variable, with the prefix on its name, into the
 * environment.  Uses the "=" from the CGI arguments.  Putting an "=" in
 * a field name is probably a bad idea.
 *
 * If the variable is already defined, append a '#' to it along with the
 * new value.
 *
 * If the variable name begins with an underline, strip whitespace from the
 * start and end and normalize end-of-line characters.
 */
static void stuffenv(char *var)
{
        char    *buf, *c, *s, *t, *oldval, *newval;
        int     despace = 0, got_cr = 0;

#ifdef DEBUG
        printf("Before unescape: %s\n", var);
#endif

        url_unescape((unsigned char *) var);

        /*
         * Allocate enough memory for the variable name and its value.
         */
        buf = malloc(strlen(var) + sizeof(PREFIX) + 2);
        if (buf == NULL)
                html_perror("stuffenv");

        strcpy(buf, PREFIX);
        if (var[0] == '_')
        {
                strcpy(buf + sizeof(PREFIX) - 1, var + 1);
                despace = 1;
        }
        else
                strcpy(buf + sizeof(PREFIX) - 1, var);

        /*
         * If, for some reason, there wasn't an = in the query string,
         * add one so the environment will be valid.
         *
         * Also, change periods to underscores so folks can get at "image"
         * input fields from the shell, which has trouble with periods
         * in variable names.
         */
        for (c = buf; *c != '\0'; c++)
        {
                if (*c == '.')
                        *c = '_';
                if (*c == '=')
                        break;
        }
        if (*c == '\0')
                c[1] = '\0';
        *c = '\0';

        /*
         * Do whitespace stripping, if applicable.  Since this can only ever
         * shorten the value, it's safe to do in place.
         */
        if (despace && c[1])
        {
#ifdef DEBUG
                printf("  Stripping whitespace.\n");
#endif
                for (s = c + 1; *s && isspace(*s); s++)
                        ;
                t = c + 1;
                while (*s)
                {
                        if (*s == '\r')
                        {
                                got_cr = 1;
                                s++;
                                continue;
                        }
                        if (got_cr)
                        {
                                if (*s != '\n')
                                        *t++ = '\n';
                                got_cr = 0;
                        }
                        *t++ = *s++;
                }

                /* Strip trailing whitespace if we copied anything. */
                while (t > c && isspace(*--t))
                        ;
                t[1] = '\0';
        }

        /*
         * Check for the presence of the variable.
         */
        if ((oldval = getenv(buf)))
        {
#ifdef DEBUG
                printf("  Variable %s exists with value %s\n", buf, oldval);
#endif
                newval = malloc(strlen(oldval) + strlen(buf) + strlen(c+1) + 3);
                if (newval == NULL)
                        html_perror("stuffenv: append");

                *c = '=';
                sprintf(newval, "%s#%s", buf, oldval);
                *c = '\0';

                /*
                 * Set up to free the entire old environment variable -- there
                 * really ought to be a library function for this.  It's safe
                 * to free it since the only place these variables come from
                 * is a previous call to this function; we can never be
                 * freeing a system-supplied environment variable.
                 */
                oldval -= strlen(buf) + 1; /* skip past VAR= */
        }
        else
        {
#ifdef DEBUG
                printf("  Variable %s doesn't exist yet.\n", buf);
#endif
                *c = '=';
                newval = buf;
        }

#ifdef DEBUG
        printf("  putenv %s\n", newval);
#endif
        putenv(newval);

        if (oldval)
        {
                /*
                 * Do the actual freeing of the old value after it's not
                 * being referred to any more.
                 */
                free(oldval);
                free(buf);
        }
}

/*
 * Scan a query string, stuffing variables into the environment.  This
 * should ideally just use strtok(), but that's not available everywhere.
 */
static void scanquery(char *q)
{
        char    *next = q;

        do {
                next = strchr(next, '&');
                if (next)
                        *next = '\0';

                stuffenv(q);
                if (next)
                        *next++ = '&';
                q = next;
        } while (q != NULL);
}

/*
 * Read a POST query from standard input into a dynamic buffer.  Terminate
 * it with a null character.
 */
static char *postread(NOPARAMS)
{
        char    *buf = NULL;
        int     size = 0, sofar = 0, got;

        buf = getenv("CONTENT_TYPE");
        if (buf == NULL || strcmp(buf, "application/x-www-form-urlencoded"))
        {
                http_head();
                puts("<p>No content type was passed to uncgi.");
                uncgi_tag();
                exit(1);
        }

        buf = getenv("CONTENT_LENGTH");
        if (buf == NULL)
        {
                http_head();
                puts("<p>The server did not tell uncgi how long the request");
                puts("was.");
                uncgi_tag();
                exit(1);
        }

        size = atoi(buf);
        buf = malloc(size + 1);
        if (buf == NULL)
                html_perror("postread");
        do
        {
                got = fread(buf + sofar, 1, size - sofar, stdin);
                sofar += got;
        } while (got && sofar < size);

        buf[sofar] = '\0';

        return (buf);
}

#ifndef LIBRARY
/*
 * Run a shell script.  We use this instead of the OS's "#!" mechanism
 * because that mechanism doesn't work too well on SVR3-based systems.
 */
static void runscript(char *shell, char *script)
{
        char    *argvec[4], **ppArg = argvec, *pz;

#ifdef EXECUTABLES_ONLY
        if (access(script, X_OK))
                html_perror("Script isn't executable");
#endif

        /*
         *  "shell" really points to the character following the "#!",
         *  not to the name of the shell program to run.  Skip any
         *  leading white space.
         */
        while (isspace( *shell ))
                shell++;

        /*
         *  We only run the shell if there is a full path name.
         */
        if (*shell != '/')
                return;

        *ppArg++ = shell;

        pz = shell;
        while ((*pz != '\0') && (! isspace( *pz ))) {
                pz++;
        }

#ifdef DEBUG
        printf("Skipped whitespace\n");
        fflush(stdout);
#endif

        /*
         *  We have found the end of the command.
         *  Clip off the trailing white space.
         */
        if (*pz != '\0') {
                *pz++ = '\0';
                while (isspace( *pz ))
                        pz++;

#ifdef DEBUG
                printf("Stripped shell is '%s', args are '%s'\n", shell, pz);
                fflush(stdout);
#endif

                /*
                 * If anything follows the command *AND* what follows is
                 * not a comment, then insert it into the command argument
                 * vect.
                 */

                if ((*pz != '\0') && (*pz != '#')) {
                        *ppArg++ = pz;

                        /*
                         *  Trim off anything after the first white space char.
                         */
                        while ((*pz != '\0') && (! isspace( *pz ))) {
                                pz++;
                        }

                        *pz = '\0';
                }
        }

        *ppArg++ = script;
        *ppArg   = (char*)NULL;

#ifdef DEBUG
        if (argvec[2] == NULL)
                printf("Executing '%s %s'\n\n", argvec[0], argvec[1]);
        else
                printf("Executing '%s %s %s'\n\n", argvec[0], argvec[1],
                        argvec[2]);
        fflush(stdout);
#endif

        execv(shell, argvec);
        /* Fall back to main() on error. */
}

/*
 * Figure out which part of the path information refers to a backend program.
 * This might be more than one path element if the backend is in a
 * subdirectory of the main script directory.
 */
static int find_program(char *scriptdir, char *pathinfo)
{
        struct stat     st;
        int     sdlen = strlen(scriptdir);
        int     proglen = 0;
        char    path[1000];

        strncpy(path, scriptdir, sizeof(path));

        while (pathinfo[proglen])
        {
                if (pathinfo[proglen] == '/' &&
                    ! strncmp(pathinfo + proglen, "/../", 4))
                {
                        bad_path_error(".. not allowed");
                }

                path[proglen + sdlen] = pathinfo[proglen];
                proglen++;
                if (pathinfo[proglen] == '/' || pathinfo[proglen] == '\0')
                {
                        path[proglen + sdlen] = '\0';
                        if (stat(path, &st) || !S_ISDIR(st.st_mode))
                                return (proglen);
                }
        }

        return (proglen);       /* whole path, if it's all directories */
}
#endif

/*
 * Main program, optionally callable as a library function.
 */
void uncgi(NOPARAMS)
{
        char    *query, *dupquery, *method;

#ifdef DEBUG
        printf("Content-type: text/plain\n\nUn-CGI %s (%s), debug mode\n",
                VERSION, id);
        fflush(stdout);
#endif

        /*
         * First, get the query string, wherever it is, and stick its
         * component parts into the environment.  Allow combination
         * GET and POST queries, even though that's a bit strange.
         */
        query = getenv("QUERY_STRING");
        if (query != NULL && strlen(query))
        {
                /* Ultrix doesn't have strdup(), so we do this the long way. */
                dupquery = malloc(strlen(query) + 1);
                if (dupquery)
                {
                        strcpy(dupquery, query);
                        scanquery(dupquery);
                }

#ifdef DEBUG
                printf("Scanned query string.\n");
                fflush(stdout);
#endif
        }

        method = getenv("REQUEST_METHOD");
        if (method != NULL && ! strcmp(method, "POST"))
        {
                query = postread();
                if (query[0] != '\0')
                        scanquery(query);
#ifdef DEBUG
                printf("Read POST query.\n");
                fflush(stdout);
#endif
        }

#ifndef NO_QUERY_OK
        if (query == NULL)
        {
                http_head();
                puts("<p>I couldn't find a query to process.\n</body></html>");
                exit(1);
        }
#endif
}


#ifndef LIBRARY
main(int argc, char **argv)
{
        char    *program, *pathinfo, *newpathinfo;
        char    *scriptname, *newscriptname;
        char    *ptrans, *ptend;
        char    *argvec[2], shellname[200];
        int     proglen;
        FILE    *fp;

        uncgi();

#ifndef SCRIPT_BIN
#error SCRIPT_BIN not defined.  Please compile via the Makefile.
#endif

        /*
         * Now figure out which program the caller *really* wants, and adjust
         * PATH_INFO and PATH_TRANSLATED to look right to that program.
         * Also point SCRIPT_NAME to it, so it can figure out what it's called.
         */
        pathinfo = getenv("PATH_INFO");
        if (pathinfo != NULL && pathinfo[0])
        {
                proglen = find_program(SCRIPT_BIN, pathinfo);

#ifdef DEBUG
                printf("Found program, path info is '%s'\n", pathinfo+proglen);
                fflush(stdout);
#endif

                scriptname = getenv("SCRIPT_NAME");
                if (scriptname == NULL)
                        scriptname = "";
                newscriptname = malloc(proglen + strlen(scriptname) + 1 +
                                        sizeof("SCRIPT_NAME="));
                if (newscriptname == NULL)
                        html_perror("scriptname");
                sprintf(newscriptname, "SCRIPT_NAME=%s", scriptname);
                strncat(newscriptname, pathinfo, proglen);
                putenv(newscriptname);

#ifdef DEBUG
                printf("Script name is '%s'\n", newscriptname);
                fflush(stdout);
#endif

                /*
                 * Figure out the path to the backend program.
                 */
                program = malloc(proglen + sizeof(SCRIPT_BIN));
                if (program == NULL)
                        html_perror("program");
                strcpy(program, SCRIPT_BIN);
                strncat(program + sizeof(SCRIPT_BIN) - 1, pathinfo, proglen);

#ifdef DEBUG
                printf("Program path is '%s'\n", program);
                fflush(stdout);
#endif

                /*
                 * Strip "program" from PATH_TRANSLATED.
                 * XXX - this depends on strcpy() copying front to back.
                 */
                ptrans = getenv("PATH_TRANSLATED");
                if (ptrans != NULL && strlen(ptrans))
                {
                        ptend = ptrans + strlen(ptrans) - strlen(pathinfo);
                        strcpy(ptend, ptend + proglen);
                        ptend = malloc(strlen(ptrans) +
                                        sizeof("PATH_TRANSLATED="));
                        if (ptend == NULL)
                                html_perror("ptend");
                        sprintf(ptend, "PATH_TRANSLATED=%s", ptrans);
                        putenv(ptend);
                }
                else
                {
                        /*
                         * Netscape server is too lazy to fill in the variable
                         * properly.
                         */
                        ptend = NULL;
                }

#ifdef DEBUG
                if (ptend)
                        printf("Translated path is '%s'\n", ptend);
                else
                        printf("No translated path provided by server.\n");
                fflush(stdout);
#endif

                pathinfo += proglen;
                newpathinfo = malloc(strlen(pathinfo) + sizeof("PATH_INFO="));
                if (newpathinfo == NULL)
                        html_perror("newpathinfo");
                sprintf(newpathinfo, "PATH_INFO=%s", pathinfo);
                putenv(newpathinfo);

#ifdef DEBUG
                printf("New pathinfo is '%s'\n", newpathinfo);
                fflush(stdout);
#endif
        }
        else
        {
                /*
                 * No PATH_INFO means no program to run.
                 */
                http_head();
                puts("<p>Whoever wrote this form doesn't know how to use");
                puts("the 'uncgi' program, because they didn't tell it");
                puts("what to run.");
                puts("<h5>(Bummer.)</h5>");
                uncgi_tag();

                exit(0);
        }

        /*
         * SVR3-based systems seem to have trouble running shell scripts.
         * So if that's what this is, run its shell explicitly.  Don't do
         * it on Linux since some versions have a bug in fgets() that
         * causes them to hang if no newline is encountered.
         */
#ifndef __linux__
        fp = fopen(program, "r");
        if (fp != NULL)
        {
#ifdef DEBUG
                printf("Opened program\n");
                fflush(stdout);
#endif

                if (fgets(shellname, sizeof(shellname), fp) != NULL)
                {
                        fclose(fp);
                        if (shellname[0] == '#' && shellname[1] == '!')
                                runscript(shellname + 2, program);
                }
                fclose(fp);
        }
#endif /* __linux__ */

#ifdef DEBUG
        printf("Executing '%s'\n\n", program);
        fflush(stdout);
#endif

        /*
         * Now execute the program.
         */
        argvec[0] = program;
        argvec[1] = NULL;

        execv(program, argvec);

        /*
         * If we get here, the exec failed.
         */
        ptend = strrchr(program, '/');
        if (ptend == NULL)
                ptend = program;
        else
                ptend++;
        html_perror(ptend);
}
#endif
