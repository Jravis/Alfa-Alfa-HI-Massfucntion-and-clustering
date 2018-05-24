#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>

#include <omp.h>
#define deg2rad  M_PI/180.
#define rad2deg  180./M_PI
#define Ho 70.0 //km/Mpc/s
#define c_km_s 3e5

#define omp_get_thread_num() 

//#define MIN(a,b) (((a)<(b))?(a):(b))
//#define MAX(a,b) (((a)>(b))?(a):(b))

FILE *fp =NULL;
FILE *in =NULL;

char ch_pos;
int f_read(FILE *p, long n1);
double Max(double *element, int n);
double Min(double *element, int n);
double Signal(double mass, double dl);

int i, j, k, m, l, iter;
int Mbin = 25, Wbin=20; // Number of bins in log(M_HI) and lo(W50) plane
int N;




main()
{

  double min_M, max_M, min_W, max_W, width_M, width_W;
  double Sint, Slim, Sum; 
  double *logHI_mass, *log_W, *zCMB, *D_Como;

  double *axis_M = (double *)malloc(Mbin*sizeof(double));
  double *axis_W = (double *)malloc(Wbin*sizeof(double));
  double Survey_Volume = 1.562383e+06;

  for (i = 0; i < Mbin; i++) 
    axis_M[i] = 0.0;
    
  double ** phi_jk = (double **)malloc(Mbin*sizeof(double *));

  for (i = 0; i < Mbin; i++) 
    phi_jk[i] = (double *)malloc(Wbin*sizeof(double));
   
  char *name = (char *) malloc(100*sizeof(char));
  sprintf(name, "HI_data_Veff.dat");
  printf("%s\n", name);
  in = fopen(name,"r");

  N = f_read(in, N);
  fseek(in, SEEK_SET,0);

  printf("number of data points = %d\n", N);
  printf("number of bins in log(M_HI) points = %d\n", Mbin);
  printf("number of bins in log(W50)  points = %d\n", Wbin);
  printf("number of data  points = %d\n", N);

//========================================================

  double ***H_ijk; //See Zwaan et. al. 2003 for defination
    
  H_ijk = (double***)malloc(N*sizeof(double**));

  for (i = 0; i< N; i++) 
  {
     H_ijk[i] = (double**) malloc(Mbin*sizeof(double*));
     for (j = 0; j < Mbin; j++) 
     {
       H_ijk[i][j] = (double*)malloc(Wbin*sizeof(double));
     }
  }

  for (i = 0; i < Wbin; i++) 
    axis_W[i] = 0.0;

  logHI_mass = (double*)malloc(N*sizeof(double));
  log_W      = (double*)malloc(N*sizeof(double));
  D_Como     = (double*)malloc(N*sizeof(double));
  double *Slim_gal   = (double*)malloc(N*sizeof(double));
  double *Sint_gal   = (double*)malloc(N*sizeof(double));
  double *Ra         = (double*)malloc(N*sizeof(double));
  double *Dec        = (double*)malloc(N*sizeof(double));
  zCMB        = (double*)malloc(N*sizeof(double));
  double *zHelio        = (double*)malloc(N*sizeof(double));

  for(j = 0; j < Mbin; j++)
    for(k = 0; k < Wbin; k++)
      phi_jk[j][k]=0.0;

  for(i = 0; i < N; i++)
    for(j = 0; j < Mbin; j++)
      for(k = 0; k < Wbin; k++)
        H_ijk[i][j][k] = 0;

//========================================================
//Read data file entries//

  for(i= 0; i< N; i++)
  {
    fscanf(in,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &Ra[i], &Dec[i], &logHI_mass[i], &log_W[i], &zCMB[i], 
            &D_Como[i], &zHelio[i],&Slim_gal[i], &Sint_gal[i] );
  }

  fclose(in);
  sprintf(name, "H_ijk.txt");
  printf("%s\n", name);
  in = fopen(name,"r");


  for(i= 0; i < N; i++)
    for(j= 0; j< Mbin; j++)
      for(k= 0; k< Wbin; k++)
        fscanf(in, "%*d\t%*d\t%*d\t%lf\n", &H_ijk[i][j][k]);

  fclose(in);

  sprintf(name, "Norm_phi_final_iter-25.txt");
  printf("%s\n", name);
  in = fopen(name,"r");

  for(j= 0; j< Mbin; j++)
    for(k= 0; k< Wbin; k++)
      fscanf(in, "%*d\t%*d\t%*d\t%lf\n", &phi_jk[j][k]);

  fclose(in);

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  double *Veff=(double *)malloc(N*sizeof(double));
  for(i=0; i < N; i++)
    Veff[i]= 0.; 
    
  
  #ifdef _OPENMP
  double tmp =0.;
  int omp_i=0,omp_l=0, ll, mm;
  double summ;
    
  #pragma omp parallel for schedule(dynamic) num_threads(50) default(shared) private(omp_l, tmp,summ, ll, mm)
  for(omp_i= 0; omp_i < N; omp_i++)
  {
    tmp=0.0;
    for(omp_l= 0; omp_l < N; omp_l++)
    {
      if (Sint_gal[omp_l] > Slim_gal[omp_i])
      {
        summ=0.0;
        for(ll=0; ll < Mbin; ll++)
        {
          for(mm=0; mm < Wbin; mm++)
          {
            summ =summ+ H_ijk[omp_l][ll][mm]*phi_jk[ll][mm];
          }
        }

        tmp = tmp+(1./summ);
      }
    }
    Veff[omp_i] = tmp;
  }
    #endif
   printf("Done\n") ;

   char *name1 = (char *) malloc(100*sizeof(char));
   sprintf(name1, "HI_Veff.dat");
   printf("%s\n", name1);
   in = fopen(name1,"w");

    fprintf(in, "Ra\tDec\tlogHI_mass\tlog_W\tzCMB\tD_Como\tVeff\n");

    for(i= 0; i < N; i++)
        fprintf(in, "%0.6e\t%0.6e\t%0.6e\t%0.6e\t%0.6e\t%0.6e\t%0.6e\n", Ra[i], Dec[i], logHI_mass[i], log_W[i], zCMB[i], D_Como[i], Veff[i]);

   printf("Done\n") ;
    fclose(in);
//    free(H_ijk);
    free(logHI_mass);
    free(log_W);
    free(zCMB);
    free(Sint_gal);
    free(Slim_gal);
    free(D_Como);
    free(axis_M);
    free(axis_W);
    free(Ra);
    free(Dec);
    free(Veff);
    free(name);
    free(name1);
   
}

double Min(double *element, int n)
{
    double small;
    small = element[0];
    int nn;
    for(nn=1; nn<n;nn++)
    {
        if(element[nn] < small)
            small = element[nn];
    }
    return small;
}

double Max(double *element, int n)
{
    double maximum;
    maximum = element[0];
    int nn;
    for(nn=1; nn<n;nn++)
    {
        if(element[nn] > maximum)
            maximum = element[nn];
    }
    return maximum;
}

int f_read(FILE *p, long n1)      
{
      while(1)
      {
            ch_pos=fgetc(p);
            if(ch_pos==EOF)
                  break;
            if(ch_pos=='\n')
                  n1++;
      }
      return n1;
}


double Signal(double mass, double dl)
{
    double S;
//    S = pow(10.0,mass) * (pow((1+z),2) /( pow(dl, 2)* (2*pow(10.0,5)) ) );
    S = pow(10.0,mass) /( pow(dl, 2)* (2*pow(10.0,5)) ) ;
    return S;
}

