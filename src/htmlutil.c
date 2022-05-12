/*

        Utilities for creating HTML

*/

#include "vplanet.h"


/*  WRITE_HTML_CONTENT  --  Write a string to an output stream which will
                            appear in HTML content.  Markup characters
                            are escaped as text entities.  */

void write_html_content(FILE *fp, const char *s)
{
    char c;

    while ((c = *s++) != EOS) {
        switch (c) {
            case '&':
                fputs("&amp;", fp);
                break;

            case '<':
                fputs("&lt;", fp);
                break;

            case '>':
                fputs("&gt;", fp);
                break;

            default:
                putc(c, fp);
                break;
        }
    }
}

/*  ESCAPE_HTML_CONTENT  --  Replace any HTML meta-characters in a
                             string which will appear in HTML
                             content.  If the string contains any
                             characters which require escaping, a new
                             string will be allocated and returned.
                             Otherwise, the argument string is simply
                             returned.  If you call this repeatedly, it's
                             up to you to test that the string has been
                             replaced and free the buffer.  */

char *escape_html_content(char *s)
{
    char c;
    const char *p = s;
    char *result = s, *o;
    size_t esclen = strlen(s) + 1;

    while ((c = *p++) != EOS) {
        switch (c) {
            case '&':
                esclen += 4;
                break;

            case '<':
                esclen += 3;
                break;

            case '>':
                esclen += 3;
                break;

            default:
                break;
        }
    }

    if (esclen > (strlen(s) + 1)) {
        result = malloc(esclen);
        if (result == NULL) {
            fprintf(stderr, "escape_html_content: cannot allocate %d byte escaped string\n",
                (int) esclen);
            abort();
        }
        o = result;
        while ((c = *s++) != EOS) {
            switch (c) {
                case '&':
                    strcpy(o, "&amp;");
                    o += 5;
                    break;

                case '<':
                    strcpy(o, "&lt;");
                    o += 4;
                    break;

                case '>':
                    strcpy(o, "&gt;");
                    o += 4;
                    break;

                default:
                    *o++ = c;
                    break;
            }
        }
        *o++ = EOS;

        assert(esclen == (strlen(result) + 1));
    }

    return result;
}
