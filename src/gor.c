/****
    
 <BW-MGORf.c>                             

 Program for predicting secondary structure
 using BW-MGOR method
 
 Coded by Takeshi Kawabata  
   (Email : takawaba@lab.nig.ac.jp)

 last modified : October 1, 1999 
 modified by kohji@ims.u-tokyo.ac.jp : Oct. 17, 1999                 \* KOHJI *\

 Reference: 

 T.Kawabata and J.Doi " Improvement of Protein Secondary Structure Prediction 
  Using Binary Word Encoding " 
 PROTEINS: Structure, Function, and Genetics 27:36-46 (1997) 

  
 We made following modification after the paper was published :  
 (a) Infomation parameters are recalculated using 680 chains,
     which is selelcted from  pdb_select.1997-Oct 25% list. 

 (b) Length of binary word is extetended to 9 residues 

 (c) The filtering procedure is performed after prediction, 
     which removes one-residue alpha-helices and beta-sheets.     

 ## PERFORMANCE
 When Jack-knife test was peformed with "pdb_select.1997-Oct 25% list",
 Q3 = 65.1 %
 (Oct97-680.list 680 0.45 1.20 
  Q3 65.108 QH 67.0 QE 43.5 QC 73.7 CH 0.477 CE 0.385 CE 0.446 CT 0.615 Ftype 3).

 When you can prepare multiple alignment data, prediction accuracy will  
 more high value (may be about 68 %).

 ****HOW TO COMPILE in UNIX system
 
 This program is written in C language (no ANSI).
 
 ( i)  You must prepare three files:
       BW-MGORf.c (this file) Isma.dat and IsA.dat.

 (ii) You must rewrite infodire[]  in BW-MGORf.c which shows the location        
      of "Isma.dat" and "IsA.dat".
      For example, when "Isma.dat" and "IsA.dat" are    
      in directory "/usr/people/takawaba/DATA/",
      you must modify as follows:   

      static char infodire[] = "/usr/people/takawaba/DATA/";

 (iii) Type as follows
 
    cc BW-MGORf.c -lm -o BW-MGORf


  ****HOW TO USE

   We prepare three type input file. 
  (i) Plain one-letter text file (fasta format) such as:

      > SEQUENCE NAME, COMMENT, ETC....   
      RPDFCLEPPYTGPCKARIIRYFYNAKAGLCQTFVYGGCRAKRNNFKSAEDCMRTCGGA

      Type as follows :  BW-MGORf P [filename]  

  (ii) DSSP format file

      Type as follows :  BW-MGORf D [filename]  

      We can see also the prediction accuracy. 

  (iii) Multiple Sequence Alignment in ClustalW format

      Type as follow :  BW-MGORf M [filename]  


  **** OUTPUT OPTION

 if you add '-O V' option to the command, 
 you can obtain following vertical output 

 (example)
 #    AA  Pr  IH     IE     IC     Re
 1     M   _  -5.990 -6.268  4.701   
 2     E   _  -2.444 -1.491  0.705   
 3     I   _  -1.735 -1.177  0.047   
 4     T   _  -1.638 -1.353  0.239   
 5     N   _  -1.315 -1.895  0.345   
 6     V   _  -0.458 -2.387 -0.066   
 7     N   H  -0.156 -2.647 -0.185   
 8     E   H   0.565 -2.743 -0.705   
 9     Y   H   1.666 -2.648 -1.554   
 10    E   H   1.686 -2.471 -1.757   
 11    A   H   1.879 -2.324 -2.049   
 12    I   H   1.530 -2.306 -1.791   
 13    A   H   1.205 -2.923 -1.102   
 14    K   H   1.230 -3.308 -0.894   

 IH, IE and IC is a information value of 
 "Helix","Sheet" and "Coil" respectively.
 Predicted stated is the largest information value state.


****/


#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

char Amino_Number();
char Sec_Number();
char Nary_Number();


#define MAX_A_NUM   2048 
#define MAX_LET_MA    21 
#define MAX_SEQ     4000 
#define MAX_ALI      100 

/**** DIRECTORY where Isma.dat and IsA.dat located  **********/
/*
  You must modify following line !!  

static char infodire[] = "/usr/people/takawaba/work/CombiGORJK/Predict/";
*/

static char infodire[] = "./";

/******* Protein Data ************/
char  header[80];     /* Header Comment            */
char  compnd[80];     /* Compound Comment          */
char  source[80];     /* Source Comment            */
char  Seq[MAX_SEQ];      /* Amino Acid Seq            */
char  Sec[MAX_SEQ];      /* Secondary Structure Seq   */
char  NarySeq[MAX_SEQ];  /* N-ary Amino Acid          */
int   pdbrnum[MAX_SEQ];  /* Residue Number of PDB     */
int   all_seq_num;    /* Number of Sequence        */
float SecRatio[3];    /* Secondary Structure Ratio */  
char seqfile[128]; 
char proname[128];
char plistfile[32]; 


/******* Aligned sequence data ********/
int  NALIGN;                 /* The number of aligned sequence  */
int  SEQLENGTH;              /* length of the sequence  */
char ASEQ[MAX_ALI][MAX_SEQ]; /* Aligned Sequence */
char ANAME[MAX_ALI][16];     /* Name of aligned sequence */ 


int LET_ma;
int LET_A;  /* LET-NUMBER */
int N_ARY;
int NLET; /* NLET = N^LET */


/**** N-Symbol Taxonomy Table *********/
char taxfile[128];
int Nary_Tab[20];


/********** Information *************/
float Isma[3][MAX_LET_MA][20];
float ISmA[MAX_LET_MA];
float IsA[3][MAX_A_NUM];
float DC[3];
float Isme[3][MAX_LET_MA];  /* Edge Information */

/********* Information along sequnece *********/
float ISeq[MAX_SEQ][3];


/********** Prediction Data ***********/
char PSec[MAX_SEQ];           /* Predicted Secondary Srtucture */ 
char APSec[MAX_ALI][MAX_SEQ]; /* Aligned Predicted Seconday Structure */

int N;
int Ns[3];   
int Np[3];
int Nsp[3][3];     /* Real ,Predict */

double Lambda;     /* Coupling Prameter */
double DCweight;   /* DC weight (mu)    */
double DCweightM;  /* DC weight (mu) for Multiple */

double HR[4],PR[4],CS[4];  /* 0:h 1:e 2:c 3:total */


char SeqType; /* D :DSSP file  P :Single SEQ  M : Multiple SEQ */ 
char Ftype;   /* T :filter     - :no filter  */ 
char Otype;   /* V :Vertical output */ 
int cgi_len;


main(argc, argv)                                                     /* KOHJI */
  int  argc;                                                         /* KOHJI */
  char  **argv;                                                      /* KOHJI */
{                                                                    /* KOHJI */
    unsigned long  count=0;                                          /* KOHJI */
                                                                     /* KOHJI */
    fprintf(stderr, "\n <BW-MGORf>                              \n");/* KOHJI */
    fprintf(stderr, "\n");                                           /* KOHJI */
    fprintf(stderr, " Program for predicting secondary structure\n");/* KOHJI */
    fprintf(stderr, " using BW-MGOR method\n");                      /* KOHJI */
    fprintf(stderr, "\n");                                           /* KOHJI */
    fprintf(stderr, " Coded by Takeshi Kawabata  \n");               /* KOHJI */
    fprintf(stderr, "   (Email : takawaba@lab.nig.ac.jp)\n\n");      /* KOHJI */
                                                                     /* KOHJI */
    for (;;)                                                         /* KOHJI */
      {                                                              /* KOHJI */
        count++;                                                     /* KOHJI */
        if (count%50 == 0)                                           /* KOHJI */
            fprintf(stderr, "\r%lu", count);                         /* KOHJI */
        former_main(argc, argv);                                     /* KOHJI */
      }                                                              /* KOHJI */
}                                                                    /* KOHJI */


former_main(argc,argv)
 int argc;
 char **argv;
{
 int i,j,k,Nline;
 FILE *fp;
 char buff[128],line[10000],CGItype,type;

 cgi_len = 0;
 
 if (argc<3)
 {
  if (getenv("CONTENT_LENGTH") != NULL) cgi_len = atoi(getenv("CONTENT_LENGTH"));

  if (cgi_len==0)
   { printf("BW-MGORf [SeqType(P,D,M)] [seqfile] (-option value)\n");  
     printf("  P:Plain ascii one-letter file (fasta format),  D:DSSP file,\n");  
     printf("  M:Multiple aligned sequences in ClustalW format\n"); 
     printf(" <Option>\n"); 
     printf("  -O : 'V'ertical output (-)\n");  
     printf("  -mu: Decision Weight Parameter(1.2)\n");  
     printf("  -mm: Decision Weight Parameter for multiple alignment(1.2)\n");  
     printf("  -l : Coupling parameter of BW information(0.45)\n");  
     exit(1); }
  }

 /*** READING INFORMATION ***/
 sprintf(buff,"%sIsma.dat",infodire); Input_Isma(buff);
 sprintf(buff,"%sIsA.dat",infodire);  Input_IsA(buff);


 /**** DEFAULT VALUE ****/
 sprintf(buff,"01111010000010000010"); 
 for (i=0;i<20;++i) if (buff[i]=='0') Nary_Tab[i] = 0; else Nary_Tab[i] = 1; 
 
 DCweightM = 1.2; DCweight = 1.2;
 Lambda = 0.45; Ftype = 'T'; Otype = 'O'; CGItype = '-';

 if (argc>1) SeqType   = argv[1][0];    else CGItype = 'C';
 if (argc>2) sprintf(seqfile,argv[2]); 
        else { sprintf(seqfile,"HTTP"); sprintf(proname,"HTTP"); }

 /*** OPTION VALUE TRANSLATION ***/
 k = 3;
 while (k<argc)
   {
    if (argv[k][0]=='-')
    {
          if (argv[k][1]=='O') { ++k; Otype = argv[k][0]; }
     else if (argv[k][1]=='F') { ++k; Ftype = argv[k][0]; }
     else if ((argv[k][1]=='m')&&(argv[k][2]=='u')) {  ++k; DCweight  = atof(argv[k]); }
     else if ((argv[k][1]=='m')&&(argv[k][2]=='m')) {  ++k; DCweightM = atoi(argv[k]); }
     else if (argv[k][1]=='l') {  ++k; Lambda  = atof(argv[k]); }
     else {printf("#Can't understand option %s\n",argv[k]); exit(1);}
    }
   ++k;
   } /* if '-' */

 /**** READING SEQUENCE ***/
 for (i=0;i<MAX_SEQ;++i) Seq[i] = Sec[i] = ' '; 

#if 0                                                                /* KOHJI */
     if (CGItype=='C') 
  { Read_CGIBIN_Data_Line(line,&Nline,&SeqType); 
    if (SeqType=='P') Input_CGIBIN_Sequence_Data(line,Nline); 
    if (SeqType=='M') Input_CGIBIN_ClustalW_Data(line,Nline); }
 else
#endif                                                               /* KOHJI */
      if (SeqType=='D')                                              /* KOHJI */
        {                                                            /* KOHJI */
          fprintf(stderr, "exit Read_DSSP_File(argv[2]);\n");        /* KOHJI */
          exit(-11);                                                 /* KOHJI */
        }                                                            /* KOHJI */
 else if (SeqType=='P') Input_Sequence_Data(argv[2]); 
 else if (SeqType=='M') Read_ClustalW_File(argv[2]);
 else exit(1);


 /**** PREDICTION ****/ 
 Initialize_N(); 

 if (SeqType=='M') { Group_Prediction(); 
                     Each_Seq_Group_Prediction();  }
              else  Prediction(all_seq_num,Seq,PSec);
 
 if (Ftype == 'T') Filtering3dash(PSec); 

  if (SeqType=='D')
  { Count_RealPred(Nsp,Sec,PSec);
    Cal_HitRatio(Nsp,HR);
    Cal_Predicted_HitRatio(Nsp,PR);
    Cal_Correlation(Nsp,CS); }

 /**** OUTPUT ****/
 Output_Header();

 if (Otype =='V') 
 { 
  if (SeqType=='P') Output_Pre_Data_Vertical();
  if (SeqType=='M') Output_Pre_Data_Vertical_Multiple();
  if (SeqType=='D') Output_Pre_Data_Vertical();
 } 
 else 
 { 
  if (SeqType=='P') Output_Pre_Data_Horizontal();
  if (SeqType=='M') Output_Pre_Data_Horizontal_Multiple();
  if (SeqType=='D') Output_Pre_Data_Horizontal();
 }

 if (SeqType =='D') Output_Prediction_Accuracy(); 

 if (CGItype=='C') printf("</PRE></BODY></HTML>\n");

} /* end of main() */


Input_CGIBIN_Sequence_Data(line,Nline)
 char *line;
 int Nline;
{
 char head[5],sym;
 int i,j,ok;

 i = 0; all_seq_num = 1;
 while (i<Nline)
  { sym = line[i]; ++i; 
   
    if (sym=='>'){  i = line_gets(proname,128,line,i); proname[strlen(proname)-1] = '\0';} 
    if (isalpha(sym)!=0) {Seq[all_seq_num] = sym; ++all_seq_num;}   
    
    if (all_seq_num>MAX_SEQ)
     { printf("Content-type: text/html\n\n");
       printf("<BODY><H1>ERROR : LENGTH OF DATA IS TOO LONG !!</H1>\n"); 
       exit(1); }  
  }
  all_seq_num -= 1;
  Set_NarySeq();

} /* end of Input_CGIBIN_Sequence_Data() */ 


Input_CGIBIN_ClustalW_Data(line,Nline)
 char *line;
 int Nline;
{
 char Fname[128];
 char Line[128];
 int read_on,read_seq,read_group;
 int i,j,k,len,matchline,Llen;
 int m; 

 Llen = 60;
 
 for (i=0;i<MAX_SEQ;++i)
  for (k=0;k<MAX_ALI;++k) ASEQ[k][i] = ' ';

 all_seq_num = 1; NALIGN = 0;
 read_on = 0;  read_seq = read_group = 0;
 m = 0;
 m = line_gets(Line,100,line,m);  /* Read "CLUSTAL W ...." Comment */
 if (Line[0]!='C') m = 0; 
 while (m<Nline)
   {
     m = line_gets(Line,100,line,m);
     len = strlen(Line);
     if (strncmp(Line,"          ",9)==0) matchline = 1; else matchline = 0;
     if ((len>11)&&(matchline==0)&&(read_on==0)) {read_on = 1; read_group = 0; }
     if ((matchline==1)&&(read_on==1)) {read_on = 0; read_group = 0; read_seq += Llen ;}

     if (read_on==1)
      {  
         Get_Part_Of_Line(ANAME[read_group],Line,0,9);
         i = 1;
         while ((i<=Llen)&&(i<(len-16))&&(iscntrl(Line[i+15])==0))
          { k = i + read_seq;
            ASEQ[read_group][k] = Line[i+15];
            if (k>=all_seq_num) all_seq_num = k;
            ++i; }

        ++read_group;
        if (read_group>NALIGN) NALIGN = read_group;
       }
   }


} /* Read_ClustalW_File() */


#if 0                                                                /* KOHJI */
Read_CGIBIN_Data_Line(line,Nline,type)
 char *line;
 int *Nline;
 char *type;
{
 int i,j;
 char q,r,s;

 /* HEADER */ 
 printf("Content-type: text/html\n\n");
 printf("<HTML><BODY BGCOLOR=\"#FFF0E0\">\n");
 printf("<H2>Result of BW-MGORf secondary structure prediction</H2>\n"); 
 
 for (i=0;i<5;++i) getchar();
 *type = getchar();
 for (i=0;i<5;++i) getchar();
  
 i = 11; *Nline = 0;
 while (i<cgi_len)
 {
  q = getchar();
  if (q=='+') q = ' ';
  else if (q=='%')
      { r = getchar(); s = getchar();
        if((48<=r)&&(r<=57))  r -= 48; else if ((65<=r)&&(r<=70)) r -= 55;
        if ((48<=s)&&(s<=57)) s -= 48; else if ((65<=s)&&(s<=70)) s -= 55;
        q = 16*r + s; i +=2;}
  line[*Nline] = q; 
  /* printf("%c",line[*Nline]); */
  ++(*Nline); 
  ++i;
 }
 line[*Nline] = '\0'; 
 printf("<PRE>\n");

} /* end of Read_CGIBIN_DATA_Line() */
#endif                                                               /* KOHJI */


Output_Header()
{
 printf(">%s\n",proname);
 printf("#TYPE %c SEQ %s all_seq_num %d NALIGN %d \n",SeqType,seqfile,all_seq_num,NALIGN);   
 printf("#LET_ma %d LET_A %d Lambda %.2lf mu %.2lf muM %.2lf Ftype %c\n",
   LET_ma,LET_A,Lambda,DCweight,DCweightM,Ftype);   
 printf("#Naa %3d NH %3d (%4.1lf %%) NE %3d (%4.1lf %%) NC %3d (%4.1lf %%)\n\n",
   N,Np[0],(double)Np[0]/N*100.0,Np[1],(double)Np[1]/N*100.0,Np[2],(double)Np[2]/N*100.0);
}

Output_Prediction_Accuracy()
{
   printf("# Prediction Accuracy\n"); 
   printf("#  N %3d NH %3d (%4.1lf) NE %3d (%4.1lf) NC %3d (%4.1lf)\n",
    N,Ns[0],(double)Ns[0]/N*100.0,Ns[1],(double)Ns[1]/N*100.0,Ns[2],(double)Ns[2]/N*100.0);
   printf("#    H      E      C      Total\n");       
   printf("# HR %6.2lf %6.2lf %6.2lf %6.2lf\n",100*HR[0],100*HR[1],100*HR[2],100*HR[3]); 
   printf("# PR %6.2lf %6.2lf %6.2lf %6.2lf\n",100*PR[0],100*PR[1],100*PR[2],100*PR[3]); 
   printf("# CS %6.3lf %6.3lf %6.3lf %6.3lf\n",CS[0],CS[1],CS[2],CS[3]); 
}


Cal_HitRatio(K,HR)
 int K[3][3];
 double HR[4];
{
 int s,p;
 int sum;
 int N,Ns[3];

 N = Ns[0] = Ns[1] = Ns[2] = sum = 0;
 
 for (s=0;s<3;++s) sum += K[s][s];
 
 for (s=0;s<3;++s)  
  for (p=0;p<3;++p)  
   {Ns[s] += K[s][p]; N += K[s][p]; }

 HR[3] =(double)sum/(double)N;

 for (s=0;s<3;++s) 
  if (Ns[s] > 0) HR[s] = (double)K[s][s]/(double)Ns[s];
            else HR[s] = -0.01;

} /* end of Cal_HitRatio() */ 


Cal_Predicted_HitRatio(K,HR)
 int K[3][3];
 double HR[4];
{
 int s,p;
 int sum;
 int Mp[3];
 
 sum = 0;
 for (s=0;s<3;++s) sum += K[s][s];
 HR[3] =(double)sum/(double)N;

 Mp[0]=Mp[1]=Mp[2] = 0;
 
 for (p=0;p<3;++p) 
  for (s=0;s<3;++s) Mp[p] += K[s][p];


 for (p=0;p<3;++p)
   if (Mp[p]>0) HR[p] = (double)K[p][p]/(double)Mp[p];
          else  HR[p] = -0.01;

} /* end of Cal_Predicted_HitRatio() */ 


Cal_Correlation(K,Cs)
 int K[3][3]; 
 double Cs[4]; 
{ 
  int P[2][2];
  int Fr[3],Fp[3]; 
  int s,r,p,sum;
  int m,n; 
  double bunbo;


  /**** Calculate Cs *********/

  for (s=0;s<3;++s) 
   { m = (s+1)%3;  n = (s+2)%3;  
    
     P[0][0] = K[s][s];
     P[1][1] = K[m][m] + K[n][n] + K[m][n] + K[n][m];
     P[0][1] = K[m][s] + K[n][s]; 
     P[1][0] = K[s][m] + K[s][n];

     bunbo = 1;
     if ((P[1][1]>0)||(P[1][0]>0)) bunbo *= sqrt((double)(P[1][1]+P[1][0]));
     if ((P[1][1]>0)||(P[0][1]>0)) bunbo *= sqrt((double)(P[1][1]+P[0][1]));
     if ((P[0][0]>0)||(P[1][0]>0)) bunbo *= sqrt((double)(P[0][0]+P[1][0]));
     if ((P[0][0]>0)||(P[0][1]>0)) bunbo *= sqrt((double)(P[0][0]+P[0][1]));
     Cs[s] = (double)(P[0][0]*P[1][1]-P[1][0]*P[0][1])/bunbo;
   }

  /****** Calculate Ctotal *******/

  for (r=0;r<3;++r) 
   { Fr[r] = 0; 
      for (p=0;p<3;++p) Fr[r] += K[r][p]; } 

  for (p=0;p<3;++p) 
   { Fp[p] = 0; 
     for (r=0;r<3;++r) Fp[p] += K[r][p]; } 

  Cs[3] = 0.0; 
  for (r=0;r<3;++r) 
   for (p=0;p<3;++p) 
     if ((Fr[r]>0)&&(Fp[p]>0)) 
      Cs[3] += (double)K[r][p]*K[r][p]/(double)Fr[r]/(double)Fp[p];
 
   Cs[3] = sqrt(Cs[3]-1.0);


} /* end of Cal_Correlation() */


Initialize_N()
{
 int s,ss,e,ee;

 N = 0;
 for (s=0;s<3;++s) Ns[s] = 0;
 
 for (s=0;s<3;++s)
  for (ss=0;ss<3;++ss) Nsp[s][ss]  = 0;

} /* end of Initialize_N() */


Prediction(all_seq_num,Seq,Pre)
 int all_seq_num;
 char *Seq;
 char *Pre;
{
 int i,j,m,n;
 int a,s,A,A_ok;
 int R_s,P_s; 
 int h_ma,h_A;
 int B[11]; 
 FILE *fp;

 float ISma[3],ISA[3],ISmaA[3];
 float Imax; 
 char Break[MAX_LET_MA];  /* 'A':Normal Amino 'N':None 'G':Edge */
 char SeqW[MAX_LET_MA];  

 /*
 fp = fopen("sing.dat","w");
 */

 for (i=0;i<MAX_SEQ;++i) Pre[i] = ' '; 

 h_ma = (LET_ma-1)/2; h_A = (LET_A-1)/2;
 for (i=0;i<11;++i) B[i] = 0;  

 for (i=1;i<=all_seq_num;++i)
  if (Seq[i]!='!')  
   {
     R_s = Number_Sec(Sec[i]);
    
     for (s=0;s<3;++s) ISma[s] = ISA[s] = ISmaA[s] = 0.0; 
     A_ok = 1;  
     ++N; ++Ns[R_s];

     /****** Range ,Break Check ******/
    
    for (m=-h_ma;m<=h_ma;++m) 
     {  j = i+m; 
             if ((j==0)||(j==all_seq_num+1)) { Break[m+h_ma] = 'G'; SeqW[m+h_ma] = '*';  }
        else if ((j<0)||(j>all_seq_num))     { Break[m+h_ma] = 'N'; SeqW[m+h_ma] = 'n';  } 
        else if (Seq[j]=='!')                { Break[m+h_ma] = 'G'; SeqW[m+h_ma] = '!';  } 
        else                                 { Break[m+h_ma] = 'A'; SeqW[m+h_ma] = Seq[j]; }     
     }

    for (m=-h_ma;m<0    ;++m)  if (Break[m+h_ma]=='G') for (n=m-1;n>=-h_ma;--n) Break[n+h_ma] = 'N';  


    for (m=-h_A;m<=h_A;++m) 
     {  j = i+m; 
        if ((j<1)||(j>all_seq_num)) A_ok = 0;  
        else if (Seq[j]=='!')       A_ok = 0; 
        else B[m+h_A] = Nary_Tab[Number_Amino(Seq[j])];  }


    /******* Cal Information ******/

     for (m=-h_ma;m<=h_ma;++m) 
       {  a = Number_Amino(SeqW[m+h_ma]);
          if (Break[m+h_ma]=='A') for (s=0;s<3;++s) ISma[s] += Isma[s][m+h_ma][a]; 
          if (Break[m+h_ma]=='G') for (s=0;s<3;++s) ISma[s] += Isme[s][m+h_ma];    }

      if (A_ok == 1) 
        {  A  = Nary_to_Index(B[0],B[1],B[2],B[3],B[4],B[5],B[6],B[7],B[8],B[9],B[10]); 
           for (s=0;s<3;++s) ISA[s] +=IsA[s][A]; }  

      for (s=0;s<3;++s) ISma[s] = ISma[s] + DCweight*DC[s]; 

      for (s=0;s<3;++s) ISeq[i][s] = ISmaA[s] =(1.0-Lambda)*ISma[s]+ Lambda*ISA[s];
       
      Imax =-1000.0; P_s = 2;
      for (s=2;s>=0;--s) if (Imax<ISmaA[s]) { Imax=ISmaA[s]; P_s = s; };   
      Pre[i] = Sec_Number(P_s); 
      ++Np[P_s];
    /*
    fprintf(fp,"%d %f %f %f %f %f %f\n",i,ISma[0],ISma[1],ISma[2],ISA[0],ISA[1],ISA[2]);
    */
   } /* i loop */
  
  /*
  fclose(fp);
  */

} /* end of Prediction() */


Each_Seq_Group_Prediction()
{
 char seq[MAX_SEQ],pre[MAX_SEQ];
 int g,i,j,Naa,n; 

 for (g=0;g<NALIGN;++g)
 {
  Naa = 0;
  for (i=1;i<=all_seq_num;++i) if (ASEQ[g][i]!='-') { ++Naa; seq[Naa] = ASEQ[g][i];}
  /*
  printf("%d %s Naa %d\n",g,ANAME[g],Naa); 
  */
  Prediction(Naa,seq,pre);
  if (Ftype == 'T') Filtering3dash(pre);
  n = 0;
  for (i=1;i<=all_seq_num;++i) 
   if (ASEQ[g][i]!='-') { ++n; APSec[g][i] = pre[n];} else APSec[g][i] = '-';
 }

} /* end of Each_Seq_Group_Prediction() */


Group_Prediction()
{
 int i,j,m,n,g;
 int a,s,A,A_ok;
 int R_s,P_s;
 int h_ma,h_A;
 int B[11];
 float ISma[3],ISA[3],ISmaA[3];
 float Imax;
 char Break[MAX_LET_MA];  /* 'A':Normal Amino 'N':None 'G':Edge */
 FILE *fp;
 
 /*
 fp = fopen("mult.dat","w");
 Lambda = 0.35; DCweight = 0.45; 
 */

 for (i=0;i<MAX_SEQ;++i) PSec[i] = ' ';

 h_ma = (LET_ma-1)/2; h_A = (LET_A-1)/2;
 for (i=0;i<11;++i) B[i] = 0;

 for (i=1;i<=all_seq_num;++i)
  if (Seq[i]!='!')
   {
     R_s = Number_Sec(Sec[i]);

     for (s=0;s<3;++s) ISma[s] = ISA[s] = ISmaA[s] = 0.0;

     /****** Range ,Break Check ******/

    for (m=-h_ma;m<=h_ma;++m)
     {  j = i+m;
          if ((j==0)||(j==all_seq_num+1)) Break[m+h_ma] = 'G';
     else if ((j<0)||(j>all_seq_num))     Break[m+h_ma] = 'N';
     else if (Seq[j]=='!')       Break[m+h_ma] = 'G';
     else                              Break[m+h_ma] = 'A';
     }

    for (m=-h_ma;m<0 ;++m) if (Break[m+h_ma]=='G') for (n=m-1;n>=-h_ma;--n) Break[n+h_ma] = 'N';
    for (m=1;m<=h_ma;++m) if (Break[m+h_ma]=='G') for (n=m+1;n<= h_ma;++n) Break[n+h_ma] = 'N';


    for (g=0;g<NALIGN;++g)
    {
      A_ok = 1;
      for (m=-h_A;m<=h_A;++m)
       {  j = i+m;
          if ((j<1)||(j>all_seq_num)) A_ok = 0;
          else if ((ASEQ[g][j]==' ')||(ASEQ[g][j]=='.'))   A_ok = 0;
          else B[m+h_A] = Nary_Tab[Number_Amino(ASEQ[g][j])];  }

  /******* Cal Information ******/

      for (m=-h_ma;m<=h_ma;++m)
        {  j=i+m;
           if (Break[m+h_ma]=='A')
             { if ((ASEQ[g][j]!=' ')&&(ASEQ[g][j]!='.')&&(ASEQ[g][j]!='-'))
               for (s=0;s<3;++s)
                  ISma[s] += Isma[s][m+h_ma][Number_Amino(ASEQ[g][j])];
              }
           else if (Break[m+h_ma]=='G') for (s=0;s<3;++s) ISma[s] += Isme[s][m+h_ma];
 }

       if (A_ok == 1)
         {  A  = Nary_to_Index(B[0],B[1],B[2],B[3],B[4],B[5],B[6],B[7],B[8],B[9],B[10]);

            for (s=0;s<3;++s) ISA[s] +=IsA[s][A]; }


      for (s=0;s<3;++s) ISma[s] += DCweightM*DC[s];

      } /* g loop */

      for (s=0;s<3;++s) ISeq[i][s] =ISmaA[s] =(1.0-Lambda)*ISma[s]+ Lambda*ISA[s];
      
      ++N; ++Ns[R_s];

      Imax =-10000.0; P_s = 2;
      for (s=2;s>=0;--s) if (Imax<ISmaA[s]) { Imax=ISmaA[s]; P_s = s; } 
      PSec[i] = Sec_Number(P_s); 
      ++Np[P_s];
    /*  
    fprintf(fp,"%d %f %f %f %f %f %f\n",i,ISma[0],ISma[1],ISma[2],ISA[0],ISA[1],ISA[2]);
    */ 
 } /* i loop */

  /*
  fclose(fp);
  */
} /* end of Group_Prediction() */


Output_Pre_Data_Horizontal()
{
 int i,j,k,max;

 max = all_seq_num/50 + 1;
 
 for (i=0;i<max;++i)
 {
   if ((i%2)==0)
     printf("              1         2         3         4      %4d\n",(i+1)*50);
   if ((i%2)==1)
     printf("              6         7         8         9      %4d\n",(i+1)*50);
   
   printf("AmAc:"); 
  for (j=1;j<=50;++j) printf("%c",Seq[50*i+j]); printf("\n");
   
   if (SeqType=='D') 
    { printf("Real:"); for (j=1;j<=50;++j) printf("%c",Sec[50*i+j]);  printf("\n"); }
  
  printf("Pred:");  
  for (j=1;j<=50;++j) printf("%c",PSec[50*i+j]);  printf("\n");
   printf("\n");
  }
 printf("\n\n");

} /* end of Output_Pre_Data_Horizontal() */


Output_Pre_Data_Vertical()
{
 int i;

 printf("#    AA  Pr  IH     IE     IC     Re\n"); 
 for (i=1;i<=all_seq_num;++i)
 {
  printf("%-4d  %c   %c  %6.3f %6.3f %6.3f  %c\n",
   i,Seq[i],PSec[i],ISeq[i][0],ISeq[i][1],ISeq[i][2],Sec[i]); 
  }
 printf("\n\n");

} /* end of Output_Pre_Data_Vertical() */


Output_Pre_Data_Vertical_Multiple()
{
 int i,g;
 char space[50];

 for (i=0;i<(NALIGN-2);++i) space[i] = ' ';
 space[NALIGN-2] = '\0';
 
 printf("#    AA%s Pr%s    IH     IE     IC     \n",space,space); 
 for (i=1;i<=all_seq_num;++i)
 {
  printf("%-4d ",i);
  
  /* 
  printf("%-4d  %c   %c  %6.3f %6.3f %6.3f  %c\n",
   i,Seq[i],PSec[i],ISeq[i][0],ISeq[i][1],ISeq[i][2],Sec[i]); 
  */

 for (g=0;g<NALIGN;++g) printf("%c",ASEQ[g][i]); 
 printf(" ");
 for (g=0;g<NALIGN;++g) printf("%c",APSec[g][i]); 
 printf(" %c ",PSec[i]); 
  printf(" %6.3f %6.3f %6.3f  %c\n",ISeq[i][0],ISeq[i][1],ISeq[i][2],Sec[i]); 
  }
 printf("\n\n");

} /* end of Output_Pre_Data_Vertical_Multile() */


Output_Pre_Data_Horizontal_Multiple()
{
 int i,j,k,max,g,Llen;

 Llen = 60;
 max = all_seq_num/Llen + 1;

 for (i=0;i<max;++i)
 {
   /*
   if ((i%2)==0)
     printf("                    1         2         3         4      %4d\n",(i+1)*Llen);
   if ((i%2)==1)
     printf("                    6         7         8         9      %4d\n",(i+1)*Llen);
   */ 
 
  for (g=0;g<NALIGN;++g) 
   {  printf("%10s:",ANAME[g]); 
      for (j=1;j<=Llen;++j) 
       if ((Llen*i+j)<=all_seq_num)  printf("%c",ASEQ[g][Llen*i+j]); 
      printf("\n");
    }
 
  printf("Pred Multi:");  
  for (j=1;j<=Llen;++j) 
   if ((Llen*i+j)<=all_seq_num)  printf("%c",PSec[Llen*i+j]);  
  printf("\n");
   printf("\n");
  
  }
 printf("\n\n");

 printf("#  PREDICTION FOR EACH SEQUENCE  \n\n\n");
 for (i=0;i<max;++i)
 {
  for (g=0;g<NALIGN;++g) 
   {  printf("%10s:",ANAME[g]); 
      for (j=1;j<=Llen;++j) 
      if ((Llen*i+j)<=all_seq_num) printf("%c",APSec[g][Llen*i+j]); 
      printf("\n");
    }
 
  printf("Pred Multi:");  
  for (j=1;j<=Llen;++j) 
    if ((Llen*i+j)<=all_seq_num) printf("%c",PSec[Llen*i+j]);  
    printf("\n");
   printf("\n");
  
  }
 printf("\n\n");

} /* end of Output_Pre_Data_Horizontal_Align() */


Input_Isma(fname)
 char *fname;
{
 int s,m,a;
 char buff[8]; 
 FILE *fp;

 fp = fopen(fname,"r");
 if (fp==NULL) {printf("Can'f find %s\n",fname); exit(1); } 
 fscanf(fp,"LIST %s LET_ma %d\n",plistfile,&LET_ma); 
 for (s=0;s<3;++s)
  { fscanf(fp,"%s\n",buff); 
    for (a=0;a<20;++a)
    { fscanf(fp,"%c ",buff);  
      for (m=0;m<LET_ma;++m)  fscanf(fp,"%f ",&Isma[s][m][a]); 
     fscanf(fp,"\n"); }
    fscanf(fp,"%s ",buff); for (m=0;m<LET_ma;++m)  fscanf(fp,"%f ",&Isme[s][m]); 
    fscanf(fp,"\n"); 
  }

 fscanf(fp,"DC %f %f %f\n",&DC[0],&DC[1],&DC[2]); 

 fclose(fp);

} /* end of Input_Isma() */


Input_IsA(fname)
 char *fname;
{
 int s,A,j;
 FILE *fp;
 char B[11];
 char buff[8];
 int bnum1,bnum2;

 fp = fopen(fname,"r");
 if (fp==NULL) {printf("Can'f find %s\n",fname); exit(1); } 
 fscanf(fp,"LIST %s LET_A %d taxfile %s\n",plistfile,&LET_A,taxfile); 
 N_ARY = 2; 
 NLET = (int)pow((double)N_ARY,(double)LET_A);
 for (A=0;A<NLET;++A)
   fscanf(fp,"%d %s %f %f %f %d\n",
         &bnum1,buff,&IsA[0][A],&IsA[1][A],&IsA[2][A],&bnum2);
 fclose(fp);

} /* end of Input_IsA() */


#if  0                                                               /* KOHJI */
Read_DSSP_File(fname)
 char *fname; 
{
 FILE *fp;
 char Line[200];
 int end;
 char buff[10],key[20];
 int READDATA;
 int Rnum,num;
 char dfname[128];

 /*
 sprintf(dfname,"/users/kawabata/DATA/DSSP/%s",fname);
 */
 
 fp = fopen(fname,"r");

 if (fp==NULL) {printf(" Can't find %s\n",fname); exit(1); }
 else
 {
  end = 0; READDATA = 0;
  while ((end ==0))
   {
      if (fgets(Line,200,fp)==NULL) end = 1;

     /* if (strncmp(Line,"HEADER",4)==0) Get_Part_Of_Line(header,Line,0,70); */
      if (strncmp(Line,"HEADER",4)==0) Get_Part_Of_Line(proname,Line,0,70); 
      if (strncmp(Line,"COMPND",4)==0) Get_Part_Of_Line(compnd,Line,0,70);
      if (strncmp(Line,"SOURCE",4)==0) Get_Part_Of_Line(source,Line,0,70);


      if (READDATA==1)
        {
         Change_Line_to_Rnum(Line,&Rnum,&num);
         if ((Rnum>=1)&&(Rnum<MAX_SEQ))
          { pdbrnum[Rnum] = num;
            Seq[Rnum] = Line[13];
            if (Sec[Rnum]==' ') Sec[Rnum]='-';
            Sec[Rnum] = Sec_Number(Number_Sec(Line[16]));
            all_seq_num = Rnum;
 
          }
        }

      if ((Line[0]=='#')||(Line[2]=='#')) READDATA = 1;

    } /* while */

    Set_NarySeq();
    Set_SecRatio();
    fclose(fp);
  }

} /* end of Read_DSSP_File() */
#endif                                                               /* KOHJI */


Input_Sequence_Data(fname)
 char *fname;
{
int ok,count;
char sym;

 count = 1;

 while (feof(stdin)==0)
  {
   sym = getc(stdin);
        if (sym=='>') { fgets(proname,127,stdin); proname[strlen(proname)-1] = '\0';} 
   else                                                              /* KOHJI */
     {                                                               /* KOHJI */
       if (sym == '\n')                                              /* KOHJI */
           break;                                                    /* KOHJI */
       if (isalpha(sym)!=0) { Seq[count] = sym; ++count;  }          /* KOHJI */
     }                                                               /* KOHJI */
   Sec[count] = ' ';                                                 /* KOHJI */
  }

 if (feof(stdin))                                                    /* KOHJI */
     exit(-12);                                                      /* KOHJI */

 all_seq_num = count-1;
 Set_NarySeq();

} /* end of Input_Sequence_Data() */


int line_gets(part,Npart,Line,st) 
 char *part,*Line;
 int Npart,st;
{ int i,j,N,end;

 N = strlen(Line); 
 j = st; i = 0; end = 0; 
 while ((i<Npart)&&(j<N)&&(end==0))
 {
   if (Line[j]=='\n') { end = 1; ++j;}
   else 
   { part[i] = Line[j];
     ++i; ++j;  }
 }
 part[i] = '\0';
 return(j);

} /* end of line_gets() */


Read_ClustalW_File(fname)
 char *fname;
{
 FILE *fp;
 char Fname[128];
 char Line[128];
 int read_on,read_seq,read_group;
 int i,j,k,len,matchline,Llen;


 Llen = 60;
 fp = fopen(fname,"r");
 if (fp==NULL) {printf(" I can't find %s.\n",fname); exit(1);}
 
 for (i=0;i<MAX_SEQ;++i)
  for (k=0;k<MAX_ALI;++k) ASEQ[k][i] = ' ';


 all_seq_num = 1; NALIGN = 0;
 read_on = 0;  read_seq = read_group = 0;

 fgets(Line,100,fp); /* Read "CLUSTAL W ...." Comment */

 while (feof(fp)==0)
   {
     fgets(Line,100,fp);
     len = strlen(Line);
     if (strncmp(Line,"          ",9)==0) matchline = 1; else matchline = 0;
     if ((len>11)&&(matchline==0)&&(read_on==0)) {read_on = 1; read_group = 0; }
     if ((matchline==1)&&(read_on==1)) {read_on = 0; read_group = 0; read_seq += Llen ;}

     if (read_on==1)
      {  
         Get_Part_Of_Line(ANAME[read_group],Line,0,9);
         i = 1;
         while ((i<=Llen)&&(i<(len-16))&&(iscntrl(Line[i+15])==0))
          { k = i + read_seq;
            ASEQ[read_group][k] = Line[i+15];
            if (k>=all_seq_num) all_seq_num = k;
            ++i; }

        ++read_group;
        if (read_group>NALIGN) NALIGN = read_group;
       }
   }

  fclose(fp);

} /* Read_ClustalW_File() */


Filtering3dash(pred)
 char *pred;
{
 char predf[MAX_SEQ];
 int i,j;
 char word[3][3][3]; /* 0:OK 1:Prohibited */
 int s,t,u,ns,nt,nu;
 int so,uo;
 int max_s,max_t,max_u;
 float I,Imax;

 word[0][0][0] = 0; word[0][0][1] = 0; word[0][0][2] = 0;
 word[0][1][0] = 1; word[0][1][1] = 0; word[0][1][2] = 1;
 word[0][2][0] = 0; word[0][2][1] = 0; word[0][2][2] = 0;
 word[1][0][0] = 0; word[1][0][1] = 1; word[1][0][2] = 1;
 word[1][1][0] = 0; word[1][1][1] = 0; word[1][1][2] = 0;
 word[1][2][0] = 0; word[1][2][1] = 0; word[1][2][2] = 0;
 word[2][0][0] = 0; word[2][0][1] = 1; word[2][0][2] = 1;
 word[2][1][0] = 1; word[2][1][1] = 0; word[2][1][2] = 1;
 word[2][2][0] = 0; word[2][2][1] = 0; word[2][2][2] = 0;

for (i=1;i<=all_seq_num;++i) predf[i] = pred[i];
 for (i=1;i<=(all_seq_num-2);++i)
  {
 
   if (i>=2) so = Number_Sec(predf[i-1]); else so = 0;
   s = Number_Sec(predf[i]);
   t = Number_Sec(predf[i+1]);
   u = Number_Sec(predf[i+2]);
   if (i<=(all_seq_num-3)) uo = Number_Sec(predf[i+3]); else uo = 0;

 if ((word[s][t][u]==1)&&(Seq[i]!='!')&&(Seq[i+1]!='!')&&(Seq[i+2]!='!'))
   {
     Imax = -100.0;

     for (j=1;j<=2;++j)
      {
        ns = (s+j)%3; nt = t; nu = u;
        I = ISeq[i][ns] + ISeq[i+1][nt] + ISeq[i+2][nu];
        if ((I>Imax)&&(word[ns][nt][nu]==0)
            &&((i==1)||(word[so][ns][nu]==0))
            &&((i>(all_seq_num-3))||(word[nt][nu][uo]==0)))
         { Imax = I; max_s = ns; max_t = nt; max_u = nu; }

        ns = s; nt = (t+j)%3; nu = u;
        I = ISeq[i][ns] + ISeq[i+1][nt] + ISeq[i+2][nu];
        if ((I>Imax)&&(word[ns][nt][nu]==0)
            &&((i==1)||(word[so][ns][nu]==0))
            &&((i>(all_seq_num-3))||(word[nt][nu][uo]==0)))
         { Imax = I; max_s = ns; max_t = nt; max_u = nu; }

        ns = s; nt = t; nu = (u+j)%3;
        I = ISeq[i][ns] + ISeq[i+1][nt] + ISeq[i+2][nu];
        if ((I>Imax)&&(word[ns][nt][nu]==0)
            &&((i==1)||(word[so][ns][nu]==0))
            &&((i>(all_seq_num-3))||(word[nt][nu][uo]==0)))
         { Imax = I; max_s = ns; max_t = nt; max_u = nu; }
      }

    predf[i]   = Sec_Number(max_s);
    predf[i+1] = Sec_Number(max_t);
    predf[i+2] = Sec_Number(max_u);
   }
 }

 for (i=1;i<=all_seq_num;++i) pred[i] = predf[i];

  /*
 for (i=1;i<=all_seq_num;++i)
   ++Ns_F[Number_Sec(Sec[i])][Number_Sec(predf[i])];
 */

} /* end of Filtering3dash() */


Count_RealPred(N,sec,pre)
  int N[3][3];
  char *sec; 
  char *pre; 
{  int i;
   for (i=1;i<=all_seq_num;++i) ++N[Number_Sec(sec[i])][Number_Sec(pre[i])];
}


Set_SecRatio()
{
 int NS[3]; 
 int i,s; 

 for (s=0;s<3;++s) NS[s] = 0;
 for (i=1;i<=all_seq_num;++i) ++NS[Number_Sec(Sec[i])];

 for (s=0;s<3;++s) SecRatio[s] = (float)NS[s]/(float)all_seq_num; 

} /* end of Set_SecRatio() */


Change_Line_to_Rnum(line,Rnu,rnu)
 char *line;
 int *Rnu,*rnu;
{
  char Rchar[10],rchar[10];

  Get_Part_Of_Line(Rchar,line,0,4);
  Get_Part_Of_Line(rchar,line,5,9);

  *Rnu = atoi(Rchar);
  *rnu = atoi(rchar);

} /* end of Change_Line_to_Rnum() */


Get_Part_Of_Line(part,line,s,e)
 char *part,*line;
 int s,e;
{
 int i;

 for (i=s;i<=e;++i)
   part[i-s] = line[i];

 part[e-s+1] = '\0';

} /* end of Get_Part_of_Line() */


/***********************************************
 >> CAUTION <<

 Index_to_Nary()
 int Nary_to_Index()

  These function is described up-side-down order.
  ex) if (Nary = 3)
        3 = 01  9 = 001  30 = 0101

***********************************************/

Index_to_Nary(ind,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11)
 int ind;
 int *q1,*q2,*q3,*q4,*q5,*q6,*q7,*q8,*q9,*q10,*q11;
{
 int v1,v2,v3,v4,v5,v6,v7,v8,v9,v10;
 int base[11],i;
 int val;

 val = ind;
 base[0] = 1;
 for (i=1;i<11;++i) base[i] = base[i-1]*N_ARY;

 *q11 = val /base[10]; val -= base[10]* (*q11);
 *q10 = val/base[9];   val -= base[9] * (*q10);
 *q9  = val/base[8];   val -= base[8] * (*q9);
 *q8  = val/base[7];   val -= base[7] * (*q8);
 *q7  = val/base[6];   val -= base[6] * (*q7);
 *q6  = val/base[5];   val -= base[5] * (*q6);
 *q5  = val/base[4];   val -= base[4] * (*q5);
 *q4  = val/base[3];   val -= base[3] * (*q4);
 *q3  = val/base[2];   val -= base[2] * (*q3);
 *q2  = val/base[1];   val -= base[1] * (*q2);
 *q1  = val/base[0];

} /* end of Index_to_Nary() */


Index_to_Nary_Char(ind,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11)
 int ind;
 char *q1,*q2,*q3,*q4,*q5,*q6,*q7,*q8,*q9,*q10,*q11;
{
 int a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11;

 Index_to_Nary(ind,&a1,&a2,&a3,&a4,&a5,&a6,&a7,&a8,&a9,&a10,&a11);

  *q1 = Nary_Number(a1);
  *q2 = Nary_Number(a2);
  *q3 = Nary_Number(a3);
  *q4 = Nary_Number(a4);
  *q5 = Nary_Number(a5);
  *q6 = Nary_Number(a6);
  *q7 = Nary_Number(a7);
  *q8 = Nary_Number(a8);
  *q9 = Nary_Number(a9);
  *q10 = Nary_Number(a10);
  *q11 = Nary_Number(a11);

} /* end of Index_to_Nary() */


int Nary_to_Index(q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11)
 int q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11;
{
 int base[11],i;
 int ret;

 base[0] = 1;
 for (i=1;i<11;++i) base[i] = base[i-1]*N_ARY;

 ret = 0;
 ret += base[0]*q1;
 ret += base[1]*q2;
 ret += base[2]*q3;
 ret += base[3]*q4;
 ret += base[4]*q5;
 ret += base[5]*q6;
 ret += base[6]*q7;
 ret += base[7]*q8;
 ret += base[8]*q9;
 ret += base[9]*q10;
 ret += base[10]*q11;

 return(ret);

} /* end of Nary_to_Index() */


Set_NarySeq()
{
 int i;

 for (i=1;i<=all_seq_num;++i)
  NarySeq[i] = Nary_Number(Nary_Tab[Number_Amino(Seq[i])]);

} /* end of Set_NarySeq() */


char Nary_Number(num)
 int num;
{
 char ret;

  switch(num)
   { case  0:ret = '0'; break;
     case  1:ret = '1'; break;
     case  2:ret = 'c'; break;
     case  3:ret = 'd'; break;
     case  4:ret = 'e'; break;
     case  5:ret = 'f'; break;
     case  6:ret = 'g'; break;
     case  7:ret = 'h'; break;
     case  8:ret = 'i'; break;
     case  9:ret = 'j'; break;
     case 10:ret = 'k'; break;
     case 11:ret = 'l'; break;
     case 12:ret = 'm'; break;
     case 13:ret = 'n'; break;
     case 14:ret = 'o'; break;
     case 15:ret = 'p'; break;
     case 16:ret = 'q'; break;
     case 17:ret = 'r'; break;
     case 18:ret = 's'; break;
     case 19:ret = 't'; break;
     default:ret = 'a';
   }
 return(ret);
} /* end of Nary_Number() */


int  Number_Nary(sec)
 char sec;
{
 int  ret;

  switch(sec)
   { case '0':ret = 0; break;
     case '1':ret = 1; break;
     case 'c':ret = 2; break;
     case 'd':ret = 3; break;
     case 'e':ret = 4; break;
     case 'f':ret = 5; break;
     case 'g':ret = 6; break;
     case 'h':ret = 7; break;
     case 'i':ret = 8; break;
     case 'j':ret = 9; break;
     case 'k':ret = 10; break;
     case 'l':ret = 11; break;
     case 'm':ret = 12; break;
     case 'n':ret = 13; break;
     case 'o':ret = 14; break;
     case 'p':ret = 15; break;
     case 'q':ret = 16; break;
     case 'r':ret = 17; break;
     case 's':ret = 18; break;
     case 't':ret = 19; break;
     default :ret = 0; break;
   }
 return(ret);

} /* end of Number_Nary() */


Number_Amino(rsym)
 char rsym;
{
 int i;
 switch(rsym)
  {
   case 'A':i = 0; break;
   case 'I':i = 1; break;
   case 'L':i = 2; break;
   case 'M':i = 3; break;
   case 'F':i = 4; break;
   case 'P':i = 5; break;
   case 'V':i = 6; break;
   case 'R':i = 7; break;
   case 'D':i = 8; break;
   case 'E':i = 9; break;
   case 'K':i = 10; break;
   case 'N':i = 11; break;
   case 'C':i = 12; break;
   case 'Q':i = 13; break;
   case 'H':i = 14; break;
   case 'S':i = 15; break;
   case 'T':i = 16; break;
   case 'W':i = 17; break;
   case 'Y':i = 18; break;
   case 'G':i = 19; break;
   case 'B':i =  8; break;
   case 'Z':i =  9; break;
   case 'X':i = 0; break; 
   default:i = 0; break;
  }
 return(i);

} /* end of Number_Amino() */


char Amino_Number(num)
 int num;
{
 char i;
 switch(num)
  {
   case 0:i = 'A'; break;
   case 1:i = 'I'; break;
   case 2:i = 'L'; break;
   case 3:i = 'M'; break;
   case 4:i = 'F'; break;
   case 5:i = 'P'; break;
   case 6:i = 'V'; break;
   case 7:i = 'R'; break;
   case 8:i = 'D'; break;
   case 9:i = 'E'; break;
   case 10:i ='K'; break;
   case 11:i = 'N'; break;
   case 12:i = 'C'; break;
   case 13:i = 'Q'; break;
   case 14:i = 'H'; break;
   case 15:i = 'S'; break;
   case 16:i = 'T'; break;
   case 17:i = 'W'; break;
   case 18:i = 'Y'; break;
   case 19:i = 'G'; break;
  }
 return(i);

} /* end of Amino_Number() */


char Sec_Number(num)
 int num;
{
 char ret;

  switch(num)
   { case 0:ret = 'H'; break;
     case 1:ret = 'E'; break;
     case 2:ret = '_'; break;
     default:ret = '_';  break;  
  }
 return(ret);
}


int  Number_Sec(sec)
 char sec;
{
 int  ret;

  switch(sec)
   { case 'H':ret = 0; break;
     case 'E':ret = 1; break;
     case 'C':ret = 2; break;
     default:ret = 2; break; 
   }

 return(ret);

} /* end of Number_Sec() */
