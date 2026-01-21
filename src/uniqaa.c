/* uniqaa.c      */
/* Kohji OKAMURA */

/* Sep 17, 1999  */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define  MAX_LEN  512


unsigned long  u1, u2;
char  line[MAX_LEN], aa[MAX_LEN/2], prev[MAX_LEN/2];


int main(int argc, char *argv[])
{
    while (fgets(line, MAX_LEN, stdin) != NULL)
      {
        sscanf(line, "%lu %lu %s", &u1, &u2, aa);
        if (strcmp(prev, aa))
          {
            fputs(line, stdout);
            strcpy(prev, aa);
          }
        else
            strcpy(prev, aa);
      }
    return  0;
}
