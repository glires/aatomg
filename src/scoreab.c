/* scoreab.c */

/* Nov. 01, 1999 */
/* Kohji OKAMURA */


#include <stdio.h>
#include <stdlib.h>


#define  LIMIT_ALPHA   5.0
#define  LIMIT_BETA    8.0


short  flag;
char   line[256];
char   str1[16], str2[2], str3[128];
float  score_a, score_b;


int main(int argc, char *argv[])
{
    if (argc != 2)
        goto  usage;

    switch (argv[1][0])
      {
        case  'a':  flag = 0;
                    break;
        case  'b':  flag = 1;
                    break;
        default  :  usage:  fprintf(stderr, "usage: scoreab [a,b]\n");
                    exit(-1);
      }

    while (fgets(line, 256, stdin) != NULL)
      {
        sscanf(line, "%s %s %s %f %f", str1, str2, str3, &score_a, &score_b);
        if (flag)
          {
            if (score_b > LIMIT_BETA)
              fputs(line, stdout);
          }
        else
            if (score_a > LIMIT_ALPHA)
              fputs(line, stdout);
      }
    return  0;
}
