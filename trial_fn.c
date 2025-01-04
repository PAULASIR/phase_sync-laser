#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<complex.h>
#include<fcntl.h>
#include<limits.h>


//#include<gsl/gsl_statistics_double.h>

#define SIZE_R 500
#define SIZE_C 20



void main()
{

double x,y;
int i,j;
double n;
double arr[SIZE_C][SIZE_R];
char *file_name;
FILE *ifp=NULL;

 FILE *fp1;
    fp1 =fopen("mov_avg.dat","w");
    
    
    file_name = "t_series.dat";
    ifp =fopen(file_name,"r");
    
     if(ifp!=NULL)
        { 
         for(i=0; i <= SIZE_C; i++)
           { 
             for(j=0; j<= SIZE_R; j++){
               fscanf(ifp,"%lf",&n);
               arr[i][j] = n; 
               printf("%lf\n", arr[i][j]); 
               }
          }
           }
           
           
   
   
   
   
   
   
   
   
      
      





}
