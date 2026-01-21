/* edmg1.c       */
/* Kohji OKAMURA */

/* Jul. 08, 1999 */
/* Jan. 06, 2000 */
/* Jan. 07, 2000 */


#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define  MAX_LEN  64    /* peptide for aatobase */


char  line[MAX_LEN*3];


int main(int argc, char *argv[])
{
    int  i, j;
    int  n_seq=0;
    unsigned long  n_mg;    /* microgene number */

    /*
     * Writing the input amino sequence.
     */
    fprintf(stdout, "00000000 3    %s\n", argv[1]);

    while (fgets(line, MAX_LEN*3, stdin) != NULL)
      {
        if (strncmp(line, "microgene ", 10))
            continue;
        sscanf(line+10, "%lu", &n_mg);

        /* This version doesn't read the first frame (i=4). */
        for (i=1; i<3; i++)
          {
            fgets(line, MAX_LEN*3, stdin);
            fprintf(stdout, "%08lu %d    ", n_mg, i);
            j = 0;
            while (line[j] != '\0')
              {
                if (isalpha(line[j]))
                    fputc(toupper(line[j]), stdout);
                j++;
              }
            fputc('\n', stdout);
            n_seq++;
          }

        /*
         * read out three lines which are DNA sequences
         */
        fgets(line, MAX_LEN*3, stdin); /* first frame */
        fgets(line, MAX_LEN*3, stdin); /* sense */
        fgets(line, MAX_LEN*3, stdin); /* anti */

        if (!atoi(argv[2]))
            continue;

        for (i=4; i<7; i++)
          {
            fgets(line, MAX_LEN*3, stdin);
            fprintf(stdout, "%08lu% d    ", n_mg, i);
            for (j=strlen(line)-1; j>-1; j--)
                if (isalpha(line[j]))
                    fputc(toupper(line[j]), stdout);
            fputc('\n', stdout);
            n_seq++;
          }
      }

    return  n_seq;
}
