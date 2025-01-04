#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<complex.h>
#include<fcntl.h>
#include<limits.h>


//#include<gsl/gsl_statistics_double.h>

#define SIZE_R 6000




void main()
{

double x,y;
int i,j;
double n, fx;
double arr[SIZE_R];
char *file_name;
FILE *ifp=NULL;

 FILE *fp1,*fp2;
    fp1 =fopen("mov_avg.dat","w");
    fp2 =fopen("st_d.dat","w");
    
    file_name = "t_series.dat";
    ifp =fopen(file_name,"r");
    
     if(ifp!=NULL)
        { 
         for(i=0; i <= SIZE_R; i++)
           { 
             
               fscanf(ifp,"%lf",&n);
               arr[i] = n; 
               printf("%20.12lf\n", arr[i]); 
               }
          }
          
          
          
          for(i = 0; i <= SIZE_R; i++)
           {
             fx = sqrt(i);
             fx =fx/(float) 2500;
             fprintf(fp1,"%20.12lf\n", fx );
          
           double diff = ((arr[i]-fx)/(float)2.0) + fx ;
                   fprintf(fp2,"%20.12lf\n", diff );
     
           
          }
          
          
          
          
          
          
           }
           
           
   
   
   
   
   
   
   
   
      
      






