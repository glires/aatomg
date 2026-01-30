/* microgene.c   */
/* Kohji OKAMURA */

/* Jul 06, 1998 */
/* Jul 13, 1998 */
/* modernized by Kohji with Kiro's assistance : Jan. 30, 2026
   - Fixed compiler warnings for modern gcc and clang */


#include <time.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>


#define  EOL  '\n'
#define  MAX_SEQ_LEN  0x100


#define  ENCODE_BASE(base)       \
           (                     \
             base == 'T' ? 1 : ( \
             base == 'C' ? 2 : ( \
             base == 'A' ? 3 : ( \
             base == 'G' ? 4 : 0 \
           ))))


typedef  int  amino_acid;


FILE  *in_file;
char  file_name[128];
char  sequence[MAX_SEQ_LEN];
int  sequence_number;
struct flag
{
  unsigned  file:1;
  unsigned  rand:1;
  unsigned  stop:1;
  unsigned  term:1;
} flag;


void complement_sequence(char *seq, int len, char *comp)
{
    int  i;

    comp[len] = '\0';
    for (i=0; i<len; i++)
      {
        switch (seq[len-i-1])
          {
            case 'T':  comp[i] = 'A';
                       break;
            case 'C':  comp[i] = 'G';
                       break;
            case 'A':  comp[i] = 'T';
                       break;
            case 'G':  comp[i] = 'C';
                       break;
          }
      }
}


void print_sequences(char *seq, char *comp, char *rframe[], int len)
{
    int  i, j, surplus, aa_len, how_many_space[3];
    (void)seq;  /* unused parameter */

    fprintf(stdout, "\nmicrogene %u (%dbp)\n", sequence_number, len);

    for (i=0; i<3; i++)
      {
        aa_len = strlen(rframe[2-i]);
        for (j=0; j<2-i; j++)
            fputc(' ', stdout);
        for (j=0; j<aa_len; j++)
          {
            fputc(rframe[2-i][j], stdout);
            if (j != aa_len-1)
                fprintf(stdout, "  ");
            else
                fputc(EOL, stdout);
          }
      }

    fprintf(stdout, "%s\n", sequence);

    for (i=0; i<len; i++)
        fputc(comp[len-i-1], stdout);
    fputc(EOL, stdout);

    surplus = len % 3;
    switch (surplus)
      {
        case  0:  how_many_space[0] = 2;
                  how_many_space[1] = 1;
                  how_many_space[2] = 0;
                  break;
        case  1:  how_many_space[0] = 0;
                  how_many_space[1] = 2;
                  how_many_space[2] = 1;
                  break;
        case  2:  how_many_space[0] = 1;
                  how_many_space[1] = 0;
                  how_many_space[2] = 2;
                  break;
      }
    for (i=3; i<6; i++)
      {
        aa_len = strlen(rframe[i]);
        for (j=0; j<how_many_space[i-3]; j++)
            fputc(' ', stdout);
        for (j=0; j<aa_len; j++)
          {
            fputc(rframe[i][aa_len-j-1], stdout);
            if (j != aa_len-1)
                fprintf(stdout, "  ");
            else
                fputc(EOL, stdout);
          }
      }
}


void make_random_sequence(int bp, int how_many)
{
    int  i, j;

    srand((unsigned)(time(NULL) & 0xffff));
    for (i=0; i<how_many; i++)
      {
        for (j=0; j<bp; j++)
          {
            switch ((rand()) % 4 + 1)
              {
                case  1: fputc('T', stdout);
                         break;
                case  2: fputc('C', stdout);
                         break;
                case  3: fputc('A', stdout);
                         break;
                case  4: fputc('G', stdout);
                         break;
              }
          }
        fputc(EOL, stdout);
      }
}


int get_sequence(int file_flag)
{
    int  i, j=0;

  retry:
    if (file_flag)
      {
        if (fgets(sequence, MAX_SEQ_LEN, in_file) == NULL)
            return  -1;         /* end of file */
      }
    else
      {
        fflush(stdout);
        fprintf(stderr, "\n>");
        fgets(sequence, MAX_SEQ_LEN, stdin);
      }
    for (i = 0;; i++)
      {
        if (ENCODE_BASE(toupper(*(sequence+i))))
            sequence[j++] = toupper(*(sequence+i));
        else
          {
            if (sequence[i] == ' ' || sequence[i] == '\t')
                /* skip space and horizontal tab */
                continue;
            if (sequence[i] == '#' ||
                sequence[i] == EOL || sequence[i] == '\0')
              { /* process comment */
                /* exclude '\n' */
                sequence[j] = '\0';
                if (j)
                    return  j;  /* length of the sequence */
                else
                    goto  retry;
              }
            else
                return  -2;    /* error */
          }
      }
    /*
     * return   0> : length of the sequence
     *          0  : null line or comment line
     *         -1  : end of file
     *         -2  : sequence error
     */
}


amino_acid translation(char *codon)
{
    int  code;

    code = ENCODE_BASE(*codon)     * 0x100 +
           ENCODE_BASE(*(codon+1)) * 0x10  +
           (*(codon+2) == '\0' ? 0 : ENCODE_BASE(*(codon+2)));

    if (code > 0x419 && code < 0x425)
        return  'A';                     /* Ala, alanine */
    if (code > 0x140 && code < 0x143)
        return  'C';                     /* Cys, cysteine */
    if (code > 0x430 && code < 0x433)
        return  'D';                     /* Asp, aspartic acid */
    if (code > 0x432 && code < 0x435)
        return  'E';                     /* Glu, glutamic acid */
    if (code > 0x110 && code < 0x113)
        return  'F';                     /* Phe, phenylalanine */
    if (code > 0x439 && code < 0x445)
        return  'G';                     /* Gly, glycine */
    if (code > 0x230 && code < 0x233)
        return  'H';                     /* His, histidine */
    if (code > 0x310 && code < 0x314)
        return  'I';                     /* Ile, isoleucine */
    if (code > 0x332 && code < 0x335)
        return  'K';                     /* Lys, lysine */
    if ((code > 0x112 && code < 0x115) || (code > 0x209 && code < 0x215))
        return  'L';                     /* Leu, leucine */
    if (code == 0x314)
        return  'M';                     /* Met, methionine */
    if (code > 0x330 && code < 0x333)
        return  'N';                     /* Asn, asparagine */
    if (code > 0x219 && code < 0x225)
        return  'P';                     /* Pro, proline */
    if (code > 0x232 && code < 0x235)
        return  'Q';                     /* Gln, glutamine */
    if ((code > 0x239 && code < 0x245) || (code > 0x342 && code < 0x345))
        return  'R';                     /* Arg, arginine */
    if ((code > 0x119 && code < 0x125) || (code > 0x340 && code < 0x343))
        return  'S';                     /* Ser, serine */
    if (code > 0x319 && code < 0x325)
        return  'T';                     /* Thr, threonine */
    if (code > 0x409 && code < 0x415)
        return  'V';                     /* Var, valine */
    if (code == 0x144)
        return  'W';                     /* Trp, tryptophan */
    if (code > 0x130 && code < 0x133)
        return  'Y';                     /* Tyr, tyrosine */
    if (code == 0x133 || code == 0x134 || code == 0x143)
        return  '*';                     /* stop codon */
    else
        return  -1;                     /* translation error */
}


int terminal_translation(char *seq, char *rframe,
                         int shift, char *trip, int *i)
    /*
     * translate one or two bases at 3' terminal if necessary
     */
{
    int  aa_len, surplus;
    char  term_seq[5];

    /*
     * decide whether translate or return
     */
    aa_len = strlen(seq) / 3;
    surplus = strlen(seq) % 3;
    switch (surplus)
      {
        case  0:  if (*i == aa_len)
                      return  0;
                  break;
        case  1:  switch (shift)
                    {
                      case  0:  break;
                      case  1:  return  0;
                      case  2:  if (*i == aa_len)
                                    return  0;
                    }
                  break;
        case  2:  switch (shift)
                    {
                      case  0:  if (*i > aa_len)
                                    return  0;
                      case  1:  break;
                      case  2:  return  0;
                    }
                  break;
      }

    /*
     * concatenate the sequence and trancelate it
     */
    strcpy(term_seq, trip);
    strncat(term_seq, seq, 2);
    return  ((rframe[(*i)++] = tolower((char)translation(term_seq))) == '*');
}


int seq_to_aa(char *seq, char *rframe, int shift)
{
    char  *trip;
    int  i=0, stop_codon=0;
    amino_acid  aa;

    trip = seq + shift;
    while (*trip != '\0' && *(trip+1) != '\0')
      {
        if ((aa=translation(trip)) == -1)
            break;
        if ((rframe[i++]=(char)aa) == '*')
            stop_codon |= 1;
        if (*(trip+2) == '\0')
            break;
        trip += 3;
      }
    if (strlen(seq) >2)
        stop_codon |= terminal_translation(seq, rframe, shift, trip, &i);
    rframe[i] = '\0';
    return  stop_codon;
}


void check_option(int argc, char *argv[])
{
    while (argc--)
      {
        if (*argv[argc] == '-')
            switch (argv[argc][1])
              {
                case  'f':  flag.file |= 1;
                            break;
                case  'r':  flag.rand |= 1;
                            break;
                case  's':  flag.stop |= 1;
                            break;
                default:    ;
              }
      }
}


void file_open_error(char *name)
{
    fprintf(stderr, "File open error: %s\a\n", name);
    exit(-1);
}


int main(int argc, char *argv[])
{
    int  i;
    char  *complementary, *reading_frame[6];
    int  length;

    check_option(argc, argv);

    if (flag.rand)
      {
        for (i=1; i<argc-2; i++)
          {
            if (strncmp(argv[i], "-r", 2))
                continue;
            else
              {
                make_random_sequence(atoi(argv[i+1]), atoi(argv[i+2]));
                return  0;
              }
          }
      }


    if (flag.file)
      {
        for (i=1; i<argc-1; i++)
          {
            if (strncmp(argv[i], "-f", 2))
                continue;
            else
              {
                if ((in_file=fopen(argv[i+1], "r")) == NULL)
                    file_open_error(argv[i+1]);
              }
          }
      }

    /*
     * read sequence 
     * and allocate memory region for the complementary and reading frames
     */
    while ((length=get_sequence(flag.file)) > 0)
      {
        sequence_number++;
        complementary = (char *)malloc(length + 1);
                            /* 1: for null character */ 
        complement_sequence(sequence, length, complementary);
        for (i=0; i<6; i++)
            reading_frame[i] = (char *)malloc(length/3 + 2);
                            /* 2: for the last two bases and null character */

        /*
         * reset stop codon flag
         */
        flag.term &= 0;

        /* 
         * translate six reading frame and check stop codon
         */
        for (i=0; i<3; i++)
          {
            flag.term |= seq_to_aa(sequence, reading_frame[i], i);
            flag.term |= seq_to_aa(complementary, reading_frame[i+3], i);
          }


        if (!(flag.term && flag.stop))
            print_sequences(sequence, complementary, reading_frame, length);

        free(complementary);
        for (i=0; i<6; i++)
            free(reading_frame[i]);
      }

    if (flag.file)
        fclose(in_file);
    return  0;
}
