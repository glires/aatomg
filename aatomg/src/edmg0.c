/* edmg.c */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define  MAX_SEQ_LEN  0x100


char  header[MAX_SEQ_LEN];
char  sense3[MAX_SEQ_LEN];
char  sense2[MAX_SEQ_LEN];
char  sense1[MAX_SEQ_LEN];
char  sense0[MAX_SEQ_LEN];
char  antis0[MAX_SEQ_LEN];
char  antis1[MAX_SEQ_LEN];
char  antis2[MAX_SEQ_LEN];
char  antis3[MAX_SEQ_LEN];


int main(int argc, char *argv[])
{
    int  flag_anti_sense=0;

    if (argc > 1)
        flag_anti_sense = atoi(argv[1]);

    while (fgets(header, MAX_SEQ_LEN, stdin))
      {
        if (strncmp(header, "microgene ", 10))
            continue;
        else
          {
            fgets(sense3, MAX_SEQ_LEN, stdin);
            fgets(sense2, MAX_SEQ_LEN, stdin);
            fgets(sense1, MAX_SEQ_LEN, stdin);
            fgets(sense0, MAX_SEQ_LEN, stdin);
            fgets(antis0, MAX_SEQ_LEN, stdin);
            fgets(antis1, MAX_SEQ_LEN, stdin);
            fgets(antis2, MAX_SEQ_LEN, stdin);
            fgets(antis3, MAX_SEQ_LEN, stdin);
          }

        if (strchr(sense3, '*'))
            continue;
        if (strchr(sense2, '*'))
            continue;
        if (strchr(sense1, '*'))
            continue;

        if (flag_anti_sense)
          {
            if (strchr(antis1, '*'))
                continue;
            if (strchr(antis2, '*'))
                continue;
            if (strchr(antis3, '*'))
                continue;
          }

        fputc('\n',   stdout);
        fputs(header, stdout);
        fputs(sense3, stdout);
        fputs(sense2, stdout);
        fputs(sense1, stdout);
        fputs(sense0, stdout);
        fputs(antis0, stdout);
        fputs(antis1, stdout);
        fputs(antis2, stdout);
        fputs(antis3, stdout);
      }

    return  0;
}
