/* aatobase.c    */
/* Kohji OKAMURA */

/* Jun 29, 1999 */
/* Aug 21, 1999 */
/* modernized by Kohji with Kiro's assistance : Jan. 30, 2026
   - Fixed compiler warnings for modern gcc and clang */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define  MAX_LEN  64


int  length;
char  amino_acid[MAX_LEN];
char  dna[(MAX_LEN-1)*3+1];
char  dna2[(MAX_LEN-1)*3*2+1];
short  counter[MAX_LEN];
short  sflag = -1;

char  F[2][4] = {"TTT", "TTC"};
char  L[6][4] = {"TTA", "TTG", "CTT", "CTC", "CTA", "CTG"};
char  I[3][4] = {"ATT", "ATC", "ATA"};
char  M[1][4] = {"ATG"};
char  V[4][4] = {"GTT", "GTC", "GTA", "GTG"};
char  S[6][4] = {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"};
char  P[4][4] = {"CCT", "CCC", "CCA", "CCG"};
char  T[4][4] = {"ACT", "ACC", "ACA", "ACG"};
char  A[4][4] = {"GCT", "GCC", "GCA", "GCG"};
char  Y[2][4] = {"TAT", "TAC"};
char  H[2][4] = {"CAT", "CAC"};
char  Q[2][4] = {"CAA", "CAG"};
char  N[2][4] = {"AAT", "AAC"};
char  K[2][4] = {"AAA", "AAG"};
char  D[2][4] = {"GAT", "GAC"};
char  E[2][4] = {"GAA", "GAG"};
char  C[2][4] = {"TGT", "TGC"};
char  W[1][4] = {"TGG"};
char  R[6][4] = {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"};
char  G[4][4] = {"GGT", "GGC", "GGA", "GGG"};


int codon_n(char aa)
{
    switch (aa)
      {
        case 'M':  case 'W':
                   return  1;
        case 'F':  case 'Y':  case 'H':  case 'Q':  case 'N':
        case 'K':  case 'D':  case 'E':  case 'C':
                   return  2;
        case 'I':  return  3;
        case 'V':  case 'P':  case 'T':  case 'A':  case 'G':
                   return  4;
        case 'L':  case 'S':  case 'R':
                   return  6;
        default :  fprintf(stdout, "Amino acid sequence error: \'%c\'.\a\n",
                                                                        aa);
                   exit(-4);
      }
}


int incr_counter(void)
{
    int  i=0;

    while (1)
      {
        if (i >= length)
            return  0;
        counter[i]++;
        if (counter[i] >= codon_n(amino_acid[i]))
          {
            counter[i] = 0;
            i++;
            continue;
          }
        else
            return  1;
      }
}


void print_tcag(void)
{
    int  i;

    for (i=0; i<length; i++)
      {
        switch (amino_acid[i])
          {
            case 'F':  strcat(dna, F[counter[i]]);
                       break;
            case 'L':  strcat(dna, L[counter[i]]);
                       break;
            case 'I':  strcat(dna, I[counter[i]]);
                       break;
            case 'M':  strcat(dna, M[counter[i]]);
                       break;
            case 'V':  strcat(dna, V[counter[i]]);
                       break;
            case 'S':  strcat(dna, S[counter[i]]);
                       break;
            case 'P':  strcat(dna, P[counter[i]]);
                       break;
            case 'T':  strcat(dna, T[counter[i]]);
                       break;
            case 'A':  strcat(dna, A[counter[i]]);
                       break;
            case 'Y':  strcat(dna, Y[counter[i]]);
                       break;
            case 'H':  strcat(dna, H[counter[i]]);
                       break;
            case 'Q':  strcat(dna, Q[counter[i]]);
                       break;
            case 'N':  strcat(dna, N[counter[i]]);
                       break;
            case 'K':  strcat(dna, K[counter[i]]);
                       break;
            case 'D':  strcat(dna, D[counter[i]]);
                       break;
            case 'E':  strcat(dna, E[counter[i]]);
                       break;
            case 'C':  strcat(dna, C[counter[i]]);
                       break;
            case 'W':  strcat(dna, W[counter[i]]);
                       break;
            case 'R':  strcat(dna, R[counter[i]]);
                       break;
            case 'G':  strcat(dna, G[counter[i]]);
                       break;
            default :  fprintf(stderr, "Unknown amino acid symbol: "
                                       "\'%c\'.\a\n", amino_acid[i]);
                       exit(-3);
          }
      }

    strcpy(dna2, dna);
    strcat(dna2, dna);

    if (sflag > -1)
      {
        if (strstr(dna2, "TAA") || strstr(dna2, "TAG") || strstr(dna2, "TGA"))
            return;
        else if (sflag == 1)
                 if (strstr(dna2, "TTA") || strstr(dna2, "CTA") 
                                         || strstr(dna2, "TCA"))
                     return;
      }

    fputs(dna, stdout);
    fputc('\n', stdout);
}


void reset_str(char *str)
{
    int  i, length;

    length = strlen(str);
    for (i=0; i<length; i++)
        str[i] = '\0';
}


void print_usage(void)
{
    fprintf(stdout, "usage: aatobase AminoAcidSequence [sense flag]\n");
    exit(-1);
}


int main(int argc, char *argv[])
{
    unsigned long  sequences=0;

    /*
     * argv[1]: amino acid sequece
     * argv[2]: sense and anti sense flag
     *              0: sense only
     *              1: sense and anti sense
     */
    if (argc!=2 && argc!=3)
        print_usage();

    if  ((length=strlen(argv[1])) >= MAX_LEN)
      {
        fprintf(stderr, "The sequence is too long.\a\n");
        exit(-2);
      }
    strcpy(amino_acid, argv[1]);

    if (argc == 3)
        sflag = (short)atoi(argv[2]);

    do
      {
        print_tcag();
        sequences++;
        reset_str(dna);
      } while (incr_counter());

    return  (int)sequences;
}
