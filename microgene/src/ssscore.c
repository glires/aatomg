/* ssscore.c     */
/* Kohji OKAMURA */

/* Oct 07, 1999 */
/* Jan 29, 2026 */

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>


#define  MAX_LEN    64
#define  MAX_LINE  128


int   length;
char  line[MAX_LINE];
char  name[MAX_LINE];
char  amino_acid[MAX_LINE];


int main(int argc, char *argv[])
{
    char  aa;
    char  pr;
    double  ih, ie;
    double  score_h;
    double  score_e;
    int   i, j, n;

    while (fgets(line, MAX_LINE, stdin))
      {
        if (line[0] == '>')
          {
            line[strlen(line)-1] = '\0';
            strcpy(name, &line[1]);
            sscanf(name, "%*d %*d %s", amino_acid);
            length = strlen(amino_acid);
            score_h = 0.0;
            score_e = 0.0;
            i = 0;
            continue;
          }
        else if (line[0]=='#' || line[0]=='\n')
                 continue;

        if (sscanf(line, "%d %c %c %lf %lf %*lf", &n, &aa, &pr, &ih, &ie) == 5)
          {
            if (++i == MAX_LEN)
              {
                fprintf(stderr, "Peptide should consist of less than "
                                           "%d amino acids.\n", MAX_LEN);
                break;
              }
            assert(i == n);

            switch (pr)
              {
                case 'H':  score_h += ih;
                           break;
                case 'E':  score_e += ie;
                           break;
                default :  break;
              }
            if (i == length*2)
                fprintf(stdout, "%s %lf %lf\n",
                                 name, score_h, score_e);
          }
        else
            assert(0);
      }

    return  0;
}


#if  0
int main(int argc, char *argv[])
{
    char  aa;
    char  cf;                /* Chou-Fasman */
    char  gof;               /* Garnier-Osguthorpe-Robson */
    unsigned  score_h=0;
    unsigned  score_b=0;
    int   i=0, n;

    while (fgets(line, MAX_LINE, stdin))
        if (!strncmp(line, "Pos ", 4))
            break;

    fgets(line, MAX_LINE, stdin);    /* skip one line */
    assert(line[0] == '\n');

    while (fgets(line, MAX_LINE, stdin))
      {
        sscanf(line, "%d %c %*s %*s %*s %*s %c %c", &n, &aa, &cf, &gof);
        if (++i == MAX_LEN)
          {
            fprintf(stderr, "Peptide should consist of less than "
                                       "%d amino acids.\n", MAX_LEN);
            break;
          }
        assert(i == n);
        amino_acid[i-1] = aa;
        switch (cf)
          {
            case 'h':  score_h++;
                       break;
            case 'H':  score_h += 2;
                       break;
            case 'b':  score_b++;
                       break;
            case 'B':  score_b += 2;
                       break;
            default :  break;
          }
        switch (gof)
          {
            case 'h':  score_h += 2;
                       break;
            case 'H':  score_h += 4;
                       break;
            case 'b':  score_b += 2;
                       break;
            case 'B':  score_b += 4;
                       break;
            default :  break;
          }
      }

    amino_acid[i] = '\0';
    fprintf(stdout, "%s %s %s %u %u\n",
            argv[1], argv[2], amino_acid, score_h, score_b);

    return  0;
}
#endif
