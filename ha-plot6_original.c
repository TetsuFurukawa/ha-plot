/*************************************************************************************************************************************************************/
/*** Program for ha-plot6                                                                                                                                 ****/
/*** Compile: cc -o HPLT6 ha-plot6.c                                                                                                                      ****/
/*** Execution: ./HPLT6 -infile test3.dat -outfile test3.gp -pltdata test3.plt -ginf1 abs.ginf  -ginf2 ord.inf -window 10 -thr 0.9 -shift 250 -optaa 0    ****/
/*** ginf                                                                                                                                                 ****/
/***      num: number of lines for the information to draw a gene structure                                                                               ****/
/***      ach   ast   aed  : start and end positions of non-coding exon (ach = n), intoron (ach = i) , or coding exon (ach = c)                           ****/
/*** gnuplot:    load "test3.gp"                                                                                                                          ****/
/***             test3.dat.ps is generated                                                                                                                ****/
/***      by Hiroyuki Toh     06/11/2013                                                                                                                  ****/
/*************************************************************************************************************************************************************/
/*************************************************************************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void  Get_2DChar(char ***, int, int);
char  SimAA(char);
char  Cn(char);
void  Mk_1DInt(int **, int);
void  Mk_1DChar(char **, int);
void  Get_Gstr(char *, char **, int *, int **, int **, int *, int *, int *);

int main(int argc, char **argv)  {
  char    **Seq, **Name, dummy[100000], infile[100], outfile[50], ginf1[50], ginf2[50],
    pltdata[50], ginf[50], *ach, *och, *state1, *state2, fname1[50], fname2[50], fname3[50], fname4[50];
  int     i, j, k, l, nc, slen, thr, icnt, optaa, L1, L2, x, y, anum, onum,
    n1, i1, c1, n2, i2, c2, *abst, *aben, *orst, *oren, stn1, stn2, shift = 0;
  double  fthr;
  FILE    *fp, *fp1, *fp2, *fp3;

  Get_2DChar(&Seq, 2, 100000);

  Get_2DChar(&Name, 2, 100);

  slen = 20;
  fthr = 0.8;
  optaa = 0;

  for (i = 0 ; i < argc ; ++i)  {
    if (strcmp(*(argv + i), "-infile") == 0) {
      strcpy(infile, *(argv + i + 1));
      i++;
    }
    else if (strcmp(*(argv + i), "-outfile") == 0) {
      strcpy(outfile, *(argv + i + 1));
      i++;
    }
    else if (strcmp(*(argv + i), "-pltdata") == 0) {
      strcpy(pltdata, *(argv + i + 1));
      i++;
    }
    else if (strcmp(*(argv + i), "-ginf1") == 0) {
      strcpy(ginf1, *(argv + i + 1));
      i++;
    }
    else if (strcmp(*(argv + i), "-ginf2") == 0) {
      strcpy(ginf2, *(argv + i + 1));
      i++;
    }
    else if (strcmp(*(argv + i), "-optaa") == 0) {
      optaa = atoi(*(argv + i + 1));
      i++;
    }
    else if (strcmp(*(argv + i), "-window") == 0) {
      slen = atoi(*(argv + i + 1));
      i++;
    }
    else if (strcmp(*(argv + i), "-shift") == 0) {
      shift = atoi(*(argv + i + 1));
      i++;
    }
    else if (strcmp(*(argv + i), "-thr") == 0) {
      fthr = atof(*(argv + i + 1));
      i++;
    }
  }
  thr = (int)(fthr * slen);

  //  printf("%d %lf %s\n", slen, fthr, infile);


  Get_Gstr(ginf1, &state1, &stn1, &abst, &aben, &n1, &i1, &c1);

  Get_Gstr(ginf2, &state2, &stn2, &orst, &oren, &n2, &i2, &c2);

  /*CHK 
  printf(">> 001\n"); 
  printf("%5d %5d %5d %5d\n", stn1, n1, i1, c1);
  for (i = 0 ; i < stn1 ; ++i)  {
    printf("%5d %5d %c\n",*(abst + i), *(aben + i), *(state1 + i)); 
  }
  printf("\n\n");
  printf("%5d %5d %5d %5d\n", stn2, n2, i2, c2);
  for (i = 0 ; i < stn2 ; ++i)  {
    printf("%5d %5d %c\n",*(orst + i), *(oren + i), *(state2 + i)); 
  }  
   exit(0);
  */

  strcpy(fname1, pltdata);
  strcpy(fname2, pltdata);
  strcpy(fname3, pltdata);
  strcpy(fname4, pltdata);
  strcpy(fname1 + strlen(fname1), "n");
  strcpy(fname2 + strlen(fname2), "i");
  strcpy(fname3 + strlen(fname3), "c"); 
  strcpy(fname4 + strlen(fname4), "d");
  if ((n1 + n2) != 0) fp1 = fopen(fname1, "w");
  if ((i1 + i2) != 0) fp2 = fopen(fname2, "w");
  if ((c1 + c2) != 0) fp3 = fopen(fname3, "w");
 
  for (i = 0 ; i < stn1 ; ++i)  {
    switch(*(state1 + i)) {
    case 'n' :fprintf(fp1, "%5d %5d\n", *(abst + i), shift);
              fprintf(fp1, "%5d %5d\n", *(aben + i), shift);
              fprintf(fp1, "\n");
              break;
    case 'i' :fprintf(fp2, "%5d %5d\n", *(abst + i), shift);
              fprintf(fp2, "%5d %5d\n", *(aben + i), shift);
              fprintf(fp2, "\n");
              break;
    case 'c' :fprintf(fp3, "%5d %5d\n", *(abst + i), shift);
              fprintf(fp3, "%5d %5d\n", *(aben + i), shift);
	      fprintf(fp3, "\n");
              break;
    }
  }
  //  exit(0);
  for (i = 0 ; i < stn2 ; ++i)  {
    switch(*(state2 + i)) {
    case 'n' :fprintf(fp1, "%5d %5d\n", shift, *(orst + i));
              fprintf(fp1, "%5d %5d\n", shift, *(oren + i));
              fprintf(fp1, "\n");
              break;
    case 'i' :fprintf(fp2, "%5d %5d\n", shift, *(orst + i));
              fprintf(fp2, "%5d %5d\n", shift, *(oren + i));
              fprintf(fp2, "\n");
              break;
    case 'c' :fprintf(fp3, "%5d %5d\n", shift, *(orst + i));
              fprintf(fp3, "%5d %5d\n", shift, *(oren + i));
	      fprintf(fp3, "\n");
               break;
    } 
  }

  //  exit(0);

  if ((n1 + n2) != 0) fclose(fp1);
  if ((i1 + i2) != 0) fclose(fp2);
  if ((c1 + c2) != 0) fclose(fp3);

  //  exit(0);

  /*CHK 
  printf(">> 002\n"); 
  printf("%5d %5d %5d %5d\n", stn1, n1, i1, c1);
  for (i = 0 ; i < stn1 ; ++i)  {
    printf("%5d %5d %c\n",*(abst + i), *(aben + i), *(state1 + i)); 
  }
  printf("\n\n");
  printf("%5d %5d %5d %5d\n", stn2, n2, i2, c2);
  for (i = 0 ; i < stn2 ; ++i)  {
    printf("%5d %5d %c\n",*(orst + i), *(oren + i), *(state2 + i)); 
  }  
    exit(0);
  */

  nc = 0;
  fp1 = fopen(infile, "r");
  if (fp1 == NULL)  {
    printf("MAIN1:No Such File %s\n", infile);
    exit(0);
  }

  /*CHK 
  printf(">> 002.5\n"); 
  printf("%5d %5d %5d %5d\n", stn1, n1, i1, c1);
  for (i = 0 ; i < stn1 ; ++i)  {
    printf("%5d %5d %c\n",*(abst + i), *(aben + i), *(state1 + i)); 
  }
  printf("\n\n");
  printf("%5d %5d %5d %5d\n", stn2, n2, i2, c2);
  for (i = 0 ; i < stn2 ; ++i)  {
    printf("%5d %5d %c\n",*(orst + i), *(oren + i), *(state2 + i)); 
  }  
  exit(0);
  */

  k = 0;

  while(fgets(dummy, 99990, fp1) != NULL)  {
    //  puts(dummy);
  /*CHK 
  printf(">> 002.6\n"); 
  printf("%5d %5d %5d %5d\n", stn1, n1, i1, c1);
  for (i = 0 ; i < stn1 ; ++i)  {
    printf("%5d %5d %c\n",*(abst + i), *(aben + i), *(state1 + i)); 
  }
  printf("\n\n");
  printf("%5d %5d %5d %5d\n", stn2, n2, i2, c2);
  for (i = 0 ; i < stn2 ; ++i)  {
    printf("%5d %5d %c\n",*(orst + i), *(oren + i), *(state2 + i)); 
  }  
  exit(0);
  */
  /*
    for (i = strlen(dummy) - 1 ; i > 0 ; --i) {
      if(dummy[i] == '\n') {
	dummy[i] = '\0';
        break;
      }
    }
  */
    if (dummy[0] == '>')  {
      strcpy(*(Name + nc), dummy);
      nc++;
      k = 0;
    }
    else {
      //      strcpy(*(Seq + nc - 1)+strlen(*(Seq + nc -1)),dummy); 
      for (i = 0 ; i < strlen(dummy) ; ++i)  {
	if (dummy[i] == ' ' || dummy[i] == '\n')  continue;
        else {
	  *(*(Seq + nc - 1) + k) = dummy[i];
          k++;
	}
      }
      *(*(Seq + nc - 1) + k) = '\0';
    }
  }

  /* */
  for(i = 0 ;i < nc ; ++i)  {
    puts(*(Name + i));
    puts(*(Seq + i));
  }
  /* */

  /*
  printf("%5d %5d\n", (int)strlen(*(Seq)), (int)strlen(*(Seq + 1)));

  for (i = 0 ; i < strlen(*(Seq)); ++i)  {
    if (*(*(Seq)+ i) < 'A' && *(*(Seq) + i) > 'Z')  {
      printf("%5d: %5d %c\n", (int)strlen(*(Seq)), i, *(*(Seq) + i));
    }  
  }
  for (i = 0 ; i < strlen(*(Seq + 1)); ++i)  {
    if (*(*(Seq + 1)+ i) < 'A' && *(*(Seq + 1) + i) > 'Z')  {
      printf("%5d: %5d %c\n", (int)strlen(*(Seq + 1)), i, *(*(Seq + 1) + i));
    }  
  }
  */

  fclose(fp1);

  /*CHK 
  printf(">> 003\n"); 
  printf("%5d %5d %5d %5d\n", stn1, n1, i1, c1);
  for (i = 0 ; i < stn1 ; ++i)  {
    printf("%5d %5d %c\n",*(abst + i), *(aben + i), *(state1 + i)); 
  }
  printf("\n\n");
  printf("%5d %5d %5d %5d\n", stn2, n2, i2, c2);
  for (i = 0 ; i < stn2 ; ++i)  {
    printf("%5d %5d %c\n",*(orst + i), *(oren + i), *(state2 + i)); 
  }  
  exit(0);
   */

  fp = fopen(outfile, "w");

  /*CHK 
  printf(">> 004\n"); 
  printf("%5d %5d %5d %5d\n", stn1, n1, i1, c1);
  for (i = 0 ; i < stn1 ; ++i)  {
    printf("%5d %5d %c\n",*(abst + i), *(aben + i), *(state1 + i)); 
  }
  printf("\n\n");
  printf("%5d %5d %5d %5d\n", stn2, n2, i2, c2);
  for (i = 0 ; i < stn2 ; ++i)  {
    printf("%5d %5d %c\n",*(orst + i), *(oren + i), *(state2 + i)); 
  }  
    exit(0);
   */

  fprintf(fp, "unset key\n");
  fprintf(fp, "set size ratio %f\n", (double)strlen(*(Seq + 1))/(double)strlen(*(Seq)));
  fprintf(fp, "set xrange [1:%d]\n", (int)strlen(*(Seq)));
  fprintf(fp, "set yrange [1:%d]\n", (int)strlen(*(Seq + 1)));
  fprintf(fp, "plot \"%s\" with lines", pltdata);
  if ((n1 + n2) != 0) fprintf(fp, ", \"%s\" linecolor rgbcolor \"blue\" lt 1 lw 7 with lines", fname1);
  if ((i1 + i2) != 0) fprintf(fp, ", \"%s\" linecolor rgbcolor \"blue\" lt 1 lw 1 with lines", fname2);
  if ((c1 + c2) != 0) fprintf(fp, ", \"%s\" linecolor rgbcolor \"blue\" lt 1 lw 15 with lines", fname3);
  if ((c1 + c2) != 0) fprintf(fp, ", \"%s\" linecolor rgbcolor \"blue\" lt 12 lw 1 with lines", fname4);
  fprintf(fp, "\n");
  fprintf(fp, "set terminal postscript enhanced color\n");
  fprintf(fp, "set output \"%s.ps\"\n", infile);
  fprintf(fp, "replot\n"); 

  fclose(fp); 

  /*CHK 
  printf(">> 005\n"); 
  printf("%5d %5d %5d %5d\n", stn1, n1, i1, c1);
  for (i = 0 ; i < stn1 ; ++i)  {
    printf("%5d %5d %c\n",*(abst + i), *(aben + i), *(state1 + i)); 
  }
  printf("\n\n");
  printf("%5d %5d %5d %5d\n", stn2, n2, i2, c2);
  for (i = 0 ; i < stn2 ; ++i)  {
    printf("%5d %5d %c\n",*(orst + i), *(oren + i), *(state2 + i)); 
  }  
  exit(0);
  */


  L1 = strlen(*(Seq));
  L2 = strlen(*(Seq + 1));

  fp1 = fopen(fname4, "w");

  if ((c1 + c2) != 0) {
   for (i = 0 ; i < stn1 ; ++i)  {
     if (*(state1 + i) == 'c')   {
       fprintf(fp1, "%5d %5d\n", *(abst + i), shift);
       fprintf(fp1, "%5d %5d\n", *(abst + i), L2);
       fprintf(fp1, "\n");
       fprintf(fp1, "%5d %5d\n", *(aben + i), shift);
       fprintf(fp1, "%5d %5d\n", *(aben + i), L2);
       fprintf(fp1, "\n");
     }
   }
   for (i = 0 ; i < stn2 ; ++i)  {
     if (*(state2 + i) == 'c')   {
       fprintf(fp1, "%5d %5d\n", shift, *(orst + i));
       fprintf(fp1, "%5d %5d\n", L1, *(orst + i));
       fprintf(fp1, "\n");
       fprintf(fp1, "%5d %5d\n", shift, *(oren + i));
       fprintf(fp1, "%5d %5d\n", L1, *(oren + i));
       fprintf(fp1, "\n");
     }    
   }
  }

  fclose(fp1);

  fp = fopen(pltdata, "w");

  /* 順方向 1 */
  for (i = 0 ; i < L1 - slen ; ++i)   {
    for (j = 0 ; ; j++) {
      x = i + j;
      y = j;
      if (x >= L1 - slen || y >= L2 - slen)  break;
      if (j == 0)  {
	icnt = 0;
	for (k = 0 ; k < slen ; ++k)  {
	  if (optaa == 0)  {
	    if (*(*(Seq) + x + k) == *(*(Seq + 1) + y + k) && *(*(Seq) + x + k) != 'N')  icnt++;
	  }
	  else {
	    if (SimAA(*(*(Seq) + x + k)) == SimAA(*(*(Seq + 1) + y + k)))  icnt++;
	  }
	}
      }
      else {
	if (optaa == 0)  {
	  if (*(*(Seq) + x - 1) == *(*(Seq + 1) + y - 1) && *(*(Seq) + x - 1) != 'N')  icnt--;
          if (*(*(Seq) + x + slen - 1) == *(*(Seq + 1) + y + slen - 1) && *(*(Seq) + x + slen - 1) != 'N') icnt++;
	}
        else  {
	  if (SimAA(*(*(Seq) + x - 1)) == SimAA(*(*(Seq + 1) + y - 1)))  icnt--;
          if (SimAA(*(*(Seq) + x + slen - 1)) == SimAA(*(*(Seq + 1) + y + slen - 1))) icnt++;
	}
      }
      if (icnt >= thr)  {
        fprintf(fp, "%d %d\n", x+1, y+1);
        fprintf(fp, "%d %d\n", x+slen, y+slen);
        fprintf(fp, "\n");
      }
    }
  }

  /* 順方向 2 */
  for (i = 1 ; i < L2 - slen ; ++i)   {
    for (j = 0 ; ; j++) {
      x = j;
      y = i + j;
      if (x >= L1 - slen || y >= L2 - slen)  break;
      if (j == 0)  {
	icnt = 0;
	for (k = 0 ; k < slen ; ++k)  {
	  if (optaa == 0)  {
	    if (*(*(Seq) + x + k) == *(*(Seq + 1) + y + k) && *(*(Seq) + x + k) != 'N')  icnt++;
	  }
	  else {
	    if (SimAA(*(*(Seq) + x + k)) == SimAA(*(*(Seq + 1) + y + k)))  icnt++;
	  }
	}
      }
      else {
	if (optaa == 0)  {
	  if (*(*(Seq) + x - 1) == *(*(Seq + 1) + y - 1) && *(*(Seq) + x - 1) != 'N')  icnt--;
          if (*(*(Seq) + x + slen - 1) == *(*(Seq + 1) + y + slen - 1) && *(*(Seq) + x + slen - 1) != 'N') icnt++;
	}
        else  {
	  if (SimAA(*(*(Seq) + x - 1)) == SimAA(*(*(Seq + 1) + y - 1)))  icnt--;
          if (SimAA(*(*(Seq) + x + slen - 1)) == SimAA(*(*(Seq + 1) + y + slen - 1))) icnt++;
	}
      }
      if (icnt >= thr)  {
        fprintf(fp, "%d %d\n", x+1, y+1);
        fprintf(fp, "%d %d\n", x+slen, y+slen);
        fprintf(fp, "\n");
      }
    }
  }

  /* 逆方向 1 */
  for (i = L1 - 1 ; i >= slen - 1 ; --i)   {                                                                           
    for (j = 0 ; ; j++) {                                                                                              
      x = i - j;                                                                                                       
      y = j;                                                                                                           
      if (x < slen - 1 || y > L2 - slen)  break;                                                                       
      icnt = 0;
      for (k = 0 ; k < slen ; ++k)  {
          if (Cn(*(*(Seq) + x - k)) == *(*(Seq + 1) + y + k) && *(*(Seq) + x - k) != 'N')  icnt++;                     
      }      
      if (icnt >= thr)  {                                                                                              
        fprintf(fp, "%d %d\n", x+1, y+1);                                                                              
        fprintf(fp, "%d %d\n", x- slen + 2, y+slen);                                                                   
        fprintf(fp, "\n");                                                                                             
      }                                                                                                                
    }                                                                                                                  
  }
                                         
  /*
  for (i = L1 - 1 ; i >= slen - 1 ; --i)   {
    for (j = 0 ; ; j++) {
      x = i - j;
      y = j;
      if (x < slen - 1 || y > L2 - slen)  break;
      if (j == 0)  {
	icnt = 0;
 	for (k = 0 ; k < slen ; ++k)  {
	  if (Cn(*(*(Seq) + x - k)) == *(*(Seq + 1) + y + k) && *(*(Seq) + x - k) != 'N')  icnt++;
	}
      }
      else {
	if (Cn(*(*(Seq) + x + 1)) == *(*(Seq + 1) + y - 1) && *(*(Seq) + x + 1) != 'N')  icnt--;
	if (Cn(*(*(Seq) + x - slen + 1)) == *(*(Seq + 1) + y + slen - 1) && *(*(Seq) + x - slen + 1) != 'N') icnt++;
      }
      if (icnt >= thr)  {
        fprintf(fp, "%d %d\n", x+1, y+1);
        fprintf(fp, "%d %d\n", x- slen + 2, y+slen);
        fprintf(fp, "\n");
      }
    }
  }
  */
  
  /* 逆方向 2*/
  for (i = 1 ; i <= L2 - slen ; ++i)  {
    for (j = 0 ; ; j++)  {
      x = L1 - 1 - j;
      y = i + j;
      if ( x < slen - 1 || y > L2 - slen  )  break;
      icnt = 0;
      for (k = 0 ; k < slen ; ++k)  {
        if (Cn(*(*(Seq) + x - k)) == *(*(Seq + 1) + y + k) && *(*(Seq) + x - k) != 'N')  icnt++;
      } 
      if (icnt >= thr)   {
        fprintf(fp, "%d %d\n", x + 1, y + 1);
        fprintf(fp, "%d %d\n", x - slen + 2, y + slen);
        fprintf(fp, "\n");
      }
    }
  }
  /*  */
  /*
  for (i = 1 ; i <= L2 - slen ; ++i)   {
    for (j = 0 ; ; j++) {
      x = L1 - 1 - j;
      y = i + j;
      if (x < slen - 1 || y > L2 - slen)  break;
      if (j == 0)  {
	icnt = 0;
	for (k = 0 ; k < slen ; ++k)  {
	  if (Cn(*(*(Seq) + x - k)) == *(*(Seq + 1) + y + k))  icnt++;
	}
      }
      else {
	if (Cn(*(*(Seq) + x + 1)) == *(*(Seq + 1) + y - 1))  icnt--;
	if (Cn(*(*(Seq) + x - slen + 1)) == *(*(Seq + 1) + y + slen - 1)) icnt++;
      }
      if (icnt >= thr)  {
        fprintf(fp, "%d %d\n", x + 1, y + 1);
        fprintf(fp, "%d %d\n", x - slen + 2, y + slen);
        fprintf(fp, "\n");
      }
    }
  }
  */

  fclose(fp);

}

void  Get_2DChar(char ***Arr, int Row, int Col)  {
  int  i;

  *Arr = (char **)calloc(Row, sizeof(char *));
  if (*Arr == NULL)  {
    printf("GET2DCHAR1:Not Enough Memory\n");
    exit(0);
  }

  for (i = 0 ; i < Row ; ++i)  {
    *(*Arr + i) = (char *)calloc(Col,sizeof(char));
    if (*(*Arr + i) == NULL)  {
      printf("GET2DCHAR2:Not Enough Memory\n");
      exit(0);
    } 
    *(*(*Arr + i)) = '\0';
  }
}

char  SimAA(char AA)  {
  char Amino[] = {'L', 'I', 'M', 'V', 'S', 'P', 'T', 'A', 'G', 'F', 'Y', 'W', 'D', 'E', 'N', 'Q', 'K', 'R', 'H', 'C', 'X'};
  char  SimA[]  = {'L', 'L', 'L', 'L', 'S', 'S', 'S', 'S', 'S', 'F', 'F', 'F', 'D', 'D', 'D', 'D', 'K', 'K', 'K', 'C', 'X'};
    int i;

    for (i = 0 ; i < 21 ; ++i)  {
      if (AA == Amino[i])  {
	return SimA[i];
      }
    }

    return 'X';
}

char  Cn(char  nn)  {
  char an;

  switch(nn)  {
     case 'A' : an = 'T';
       break;
     case 'T' : an = 'A';
       break;
     case 'G' : an = 'C';
       break;
     case 'C' : an = 'G';
       break;
  }

  return  an;
}

void  Mk_1DInt(int **arr,  int size)  {
  *arr = (int *)calloc(size, sizeof(int));
  if (*arr == NULL)  {
    printf("MK1DINT1:Not Enough Memory\n");
    exit(0);
  }
}

void  Mk_1DChar(char **arr, int size) {
  *arr = (char *)calloc(size, sizeof(char));
  if (*arr == NULL)  {
    printf("MK1DCHAR1:Noty Enough Memory\n");
    exit(0);
  }
}

void  Get_Gstr(char *fname, char **state, int *stn, int **st, int **en, int *nc,int *in, int *cd) {
  FILE  *fp;
  char  dummy[110];
  int   i;

  fp = fopen(fname, "r");
  if (fp == NULL)  {
    printf("GETGSTR1:No Such File\n");
    exit(0);
  }

  fgets(dummy, 100, fp);
  sscanf(dummy, "%d", stn);

  *st = (int *)calloc(*stn, sizeof(int));
  if (*st == NULL)  {
    printf("GETGSTR2:Not Enough Memory\n");
    exit(0);
  }
  *en = (int *)calloc(*stn, sizeof(int));
  if (*en == NULL)  {
    printf("GETGSTR3:Not Enough Memory\n");
    exit(0);
  }
  *state = (char *)calloc(*stn, sizeof(char));
  if (*state == NULL)  {
    printf("GETGSTR4:Not Enough Memory\n");
    exit(0);
  }

  *nc = 0;
  *in = 0;
  *cd = 0;

  for (i = 0 ; i < *stn ; ++i)  {
    fgets(dummy, 100, fp);
    sscanf(dummy, "%c %d %d", (*state + i), (*st + i), (*en + i));
    switch(*(*state + i))  {
      case 'n' :  *nc += 1;
	break;
      case 'i' :  *in += 1;
	break;
      case 'c' :  *cd += 1;
	break;
    } 
  }

  fclose(fp);
}
