/******************************************
Authur@ Sandeep Rana                      *
This particular code for computing mass   *
function from HI survey using 2DSWL Method* 
for detail see Zwaan et. al. (2003, 2005) *
and Marting et. al. 2010.                 *
******************************************/

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram2d.h>

#define pi 3.14159265
//#include "Param.h"
#include <omp.h>

FILE *in = NULL;
FILE *out = NULL;
FILE *out1 = NULL;
FILE *out2 = NULL;
FILE *out3 = NULL;
//FILE *out4 = NULL;

char ch_pos;
int f_read(FILE *p, long n1);
double Max(double *element, int n);
double Min(double *element, int n);
double Signal(double mass, double dl);
//double Signal(double mass, double z, double dl);
double denom(double **H, double **phi);

int i, j, k, m, l, iter;
int Mbin = 25, Wbin=20; // Number of bins in log(M_HI) and lo(W50) plane
int N;

int main(int argc, char **argv){
  if(argc != 3){
    printf("\n./a.out filename MaxIter\n");
    printf("\n./a.out HI_data 25 \n");
    return(-1);
  }//if

    double ** hist = (double **)malloc(Mbin*sizeof(double *));

    double min_M, max_M, min_W, max_W, width_M, width_W;
    double Sint, Slim, Sum; 
    double *logHI_mass, *log_W, *zCMB, *D_Como;

    double *axis_M = (double *)malloc(Mbin*sizeof(double));
    double *axis_W = (double *)malloc(Wbin*sizeof(double));
    double Survey_Volume = 1.562383e+06;

    for (i = 0; i < Mbin; i++) 
    {
        axis_M[i] = 0.0;
        hist[i] = (double *)malloc(Wbin*sizeof(double));
    }
    
    double ** phi_jk = (double **)malloc(Mbin*sizeof(double *));

    for (i = 0; i < Mbin; i++) 
    {
        phi_jk[i] = (double *)malloc(Wbin*sizeof(double));
    }
   
    char *name = (char *) malloc(100*sizeof(char));
    sprintf(name, "/jbodstorage/data_sandeep/sandeep/Guo_analysis/%s.dat", argv[1]);
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
    {
        for(k = 0; k < Wbin; k++)
        {
            phi_jk[j][k]=0.0;
            hist[j][k]=0.0;
        }
    }


    for(i = 0; i < N; i++)
    {
        for(j = 0; j < Mbin; j++)
        {
            for(k = 0; k < Wbin; k++)
            {   
                H_ijk[i][j][k] = 0;
            }
        }
    }

//========================================================
//Read data file entries//

    for(i= 0; i< N; i++)
    {
        fscanf(in,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &Ra[i], &Dec[i],&logHI_mass[i], &log_W[i], &zCMB[i], &D_Como[i], &Sint_gal[i], &zHelio[i]);
        //if(i<20)
       //     printf("%lf\t%lf\n", D_Como[i], Sint_gal[i]);
        
    }
    fclose(in);
    free(name);
//Maximum and minimum range in log(M_HI) and log(W50) plane 

    min_M = 6.0; //Min(logHI_mass, N);
    max_M = 11.0; //Max(logHI_mass, N);
    min_W = Min(log_W, N);
    max_W = Max(log_W, N);

    printf("Minmum of log(M_HI) bins=%f\n", min_M);
    printf("Maximum of log(M_HI) bins=%f\n", max_M);
    printf("Minmum of log(W50) bins=%f\n", min_W);
    printf("Maximum of log(W50) bins=%f\n", max_W);
    
    width_M = (max_M-min_M)/25;
    width_W = (max_W-min_W)/20;

    printf("Width of log(W50) bins=%f\n", width_M);
    printf("Width of log(MHI) bins=%f\n", width_W);

// 2D histogram in log(M_HI) and log(W50) plane which gives
// number density and statrting bivariate distribution phi_jk
// Which we iterate in next part

    for(j= 0; j< Mbin; j++)
    {
        for(k= 0; k< Wbin; k++)
        {
            for(i= 0; i < N; i++)
            {
                if( ( logHI_mass[i] >= min_M +(j)*width_M)&&(logHI_mass[i]< min_M+(j+1)*width_M) )
                {
                    if( ( log_W[i] >= min_W +(k)*width_W)&&(log_W[i]< min_W+(k+1)*width_W) )
                    {
                        hist[j][k] += 1.0;
                        phi_jk[j][k] += 1.0;
                    } 
                }
            }
        }
    }
 
    out = fopen("/jbodstorage/data_sandeep/sandeep/Guo_analysis/phi_initial.txt","w");

    char *fname1 = (char*)malloc(100*sizeof(char));                                 
    char *fname2 = (char*)malloc(100*sizeof(char));                                 
    char *fname3 = (char*)malloc(100*sizeof(char));                                 
                                                                               
    for(j= 0; j< Mbin; j++)
    {
        for(k= 0; k< Wbin; k++)
        {
            fprintf(out,"%d\t%d\t%0.8e\n",j,k, phi_jk[j][k]);
        }
    }
    fclose(out);
    

    for(i= 0; i< Mbin; i++)
        axis_M[i]= min_M+i*width_M;
   
    for(i= 0; i< Wbin; i++)
        axis_W[i]= min_W+i*width_W;
    
    printf("==========================================\n");

    int Max_iter=atoi(argv[2]);
    printf("Max_iter=%d\n", Max_iter);

//    sprintf(fname1,"H_ijk.txt");
//    printf("%s\n",fname1);

//    out1 = fopen(fname1,"w");

//=============================================================
    
/* 
Evaluating H_ijk for each galaxy point given DL_i 
that is the luminosity distance of ith galaxy using that 
one compute intgrated flux for all HI mass in Mass bin 
and for each mass in Mbin we check for which W50 
Intergrated fulx is greater than equal 50% completness limit
So for each galaxy i we have j,k plane and H_ijk matrix 
*/

    for(i= 0; i < N; i++)
    {
        for(j= 0; j< Mbin; j++)
        {
            Sint = Signal(axis_M[j], D_Como[i]);

            for(k= 0; k< Wbin; k++)
            {
             
                if (axis_W[k] < 2.5) 
                    Slim = pow( 10.0,  (0.5*axis_W[k]-1.14-0.067) );
                if (axis_W[k]>=2.5) 
                    Slim = pow( 10.0, (axis_W[k]-2.39-0.067) );

                if (Sint > Slim)
                {
                    H_ijk[i][j][k]=1.0;
                }
                else
                {
                    H_ijk[i][j][k]=0.0;
                }
            }
        }
    }

//    for(i= 0; i < N; i++)
//        for(j= 0; j< Mbin; j++)
//            for(k= 0; k< Wbin; k++)
//                fprintf(out1,"%d\t%d\t%d\t%0.3f\n",i, j, k, H_ijk[i][j][k]);
     
//    fclose(out1);

    sprintf(fname1,"Norm_phi_final_iter-%d.txt",atoi(argv[2]));
    printf("%s\n",fname1);

    out1 = fopen(fname1,"w");


//========== Computing Bivariate density distribution iteratively ===========================//


     printf("==== Computing Bivariate density distribution iteratively ====\n");

     for(iter= 0; iter< Max_iter; iter++)
     {
        printf("iter=%d\n",iter);

        for(j= 0; j< Mbin; j++)
        {    
            for(k= 0; k< Wbin; k++)
            {
                Sum = 0.0;
                for(i= 0; i< N; i++)
                {

                    Sum=Sum+H_ijk[i][j][k]/denom(H_ijk[i], phi_jk);
                }
                if(Sum!=0)
                    phi_jk[j][k] = hist[j][k]/Sum;
            }
        }
    }
    printf("==== Computing Non normalized phi_j ====\n");
    sprintf(fname2,"/jbodstorage/data_sandeep/sandeep/Guo_analysis/phi_jk_Massfunction_bin_%d-%d_iter-%d.dat",Mbin, Wbin, atoi(argv[2]));

    printf("%s\n",fname2);
    out2 = fopen(fname2,"w");

    for(j= 0; j< Mbin; j++)
    {
        Sum = 0.0;
        for(k= 0; k< Wbin; k++)
        { 
            Sum=Sum+phi_jk[j][k];

            //fprintf(out1,"%d\t%d\t%0.8e\n",j,k, phi_jk[j][k]);
        }

        fprintf(out2,"%0.8e\t%0.8e\n",axis_M[j], Sum);
    }
    
    fclose(out2);


//====== Noramlising Mass function =======================//

    printf("==== Noramlising Mass function  ====\n");

    double avg_number_density, Normalising_factor, Sum1=0.0; 

//first normalise phi_jk to unity
 
    Sum = 0.0;

    for(j= 0; j< Mbin; j++)
    {
        for(k= 0; k< Wbin; k++)
        { 
            Sum=Sum+(phi_jk[j][k]);
        }
        Sum1=Sum1+Sum;
    }
    
    Normalising_factor = 1./Sum1;
    printf("Normalising factor = %0.4e\n", Normalising_factor);

    for(j= 0; j< Mbin; j++)
    {
        for(k= 0; k< Wbin; k++)
        { 
            phi_jk[j][k] = phi_jk[j][k]*Normalising_factor;
        }
    }

//Computing n_bar
    Sum = 0.0;
    double temp = 0.;
    for(i= 0; i< N; i++)
    {   temp=  (1./(denom(H_ijk[i], phi_jk)));
        //if (temp!=0)
        Sum = Sum + temp ;
    }  
    printf("Effective Volume = %0.4e\n", Sum);

    avg_number_density = (Sum/Survey_Volume)*0.7*0.7*0.7 ;
    printf("avg_number_density = %0.6f\n", avg_number_density);


    sprintf(fname3,"/jbodstorage/data_sandeep/sandeep/Guo_analysis/\
Norm_Massfunc_bin_%d-%d_iter-%d_1.dat",Mbin, Wbin, Max_iter);
    printf("%s\n",fname3);
    out3 = fopen(fname3,"w");

    for(j= 0; j< Mbin; j++)
    {
        for(k= 0; k< Wbin; k++)
        { 
            phi_jk[j][k] = phi_jk[j][k]*avg_number_density;

            fprintf(out1,"%d\t%d\t%0.8e\n",j,k, phi_jk[j][k]);
        }
    }
    fclose(out1);
    for(j= 0; j< Mbin; j++)
    {
        Sum = 0.0;
        for(k= 0; k< Wbin; k++)
        { 
            Sum = Sum + phi_jk[j][k]; 
        }

        fprintf(out3,"%0.8e\t%0.8e\t%.8e\n",axis_M[j], log10(Sum/0.2), Sum);
    }
    
    //=========================Veffective==========================
    
    sprintf(fname2,"HI_data_Veff.dat");
    printf("%s\n",fname2);
    out2 = fopen(fname2,"w");


    for(i= 0; i < N; i++)
    {
        if (log_W[i] < 2.5) 
            Slim_gal[i] = pow( 10.0,  (0.5*log_W[i]-1.14-0.067) );
        if (axis_W[i]>=2.5) 
            Slim_gal[i] = pow( 10.0, (log_W[i]-2.39-0.067) );
        
        fprintf(out2,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", Ra[i], Dec[i], logHI_mass[i], log_W[i], zCMB[i], D_Como[i], zHelio[i],
                 Slim_gal[i], Sint_gal[i] );
    }
    
//======= Free Mem and file pointers ==================================//

    printf("==== Free Mem and file pointers ====\n");

    for (i = 0; i < Mbin; i++){  
       free(phi_jk[i]);  
       free(hist[i]);  
    }     
    free(phi_jk);
    free(hist);

    for (i = 0; i < N; i++){  
        for (j = 0; j < Mbin; j++){  
            free(H_ijk[i][j]);  
        }
    }
    free(H_ijk);
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
//    free(Veff);

    fclose(out3);
    fclose(out2);
    printf("Done\n");
    free(fname1);
    printf("Done\n");
    free(fname2);
    printf("Done\n");
    free(fname3);
    printf("Done\n");

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

double denom(double **H, double **phi)
{
    double summ;
    summ=0.0;
    for(l=0; l < Mbin; l++)
    {
        for(m=0; m < Wbin; m++)
        {
            summ =summ+ H[l][ m]*phi[l][ m];
        }
    }
    return summ;
}
