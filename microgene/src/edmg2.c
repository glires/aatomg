/* edmg2.c       */
/* Kohji OKAMURA */

/* Jul 08, 1999  */
/* Sep 17, 1999  */
/* Oct 07, 1999  */


#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define  MAX_LINE  128
#define  MAX_LEN    64    /* peptide for aatobase */


char  line[MAX_LINE];
char  peptide[MAX_LEN];
char  command[MAX_LINE];
FILE  *fp, *f1;


int main(int argc, char *argv)
{
    int  i, n_aa=0;
    unsigned long  n_mg;

    fp = fopen("peptides", "r");
    while (fgets(line, MAX_LINE, fp) != NULL)
      {
        if (!sscanf(line, "%lu %d %s", &n_mg, &i, peptide))
            exit(-1);

        fprintf(stdout, ">%06lu %d %s\n%s%s\n",
                          n_mg, i, peptide, peptide, peptide);

        n_aa++;
      }

    return  n_aa;
}


# if  0    /* for gcg */
int main(int argc, char *argv)
{
    int  i, n_aa=0;
    unsigned long  n_mg;

    fp = fopen("peptides.uniq", "r");
    while (fgets(line, MAX_LINE, fp) != NULL)
      {
        if (!sscanf(line, "%lu %d %s", &n_mg, &i, peptide))
            exit(-1);

        remove("seq");
        remove("seq.p2s");

        f1 = fopen("seq", "w");
     /* fprintf(f1, "%s%s%s\n", peptide, peptide, peptide); */
        fprintf(f1, "%s%s\n", peptide, peptide);
        fclose(f1);

        system("reformat seq");
        system("peptidestructure -IN=seq -D");
        sprintf(command, "blue_ssscore %08lu %d < seq.p2s >> ./scorelist",
                                       n_mg, i);
        system(command);

        n_aa++;
      }

    return  n_aa;
}
# endif    /* for gcg */
