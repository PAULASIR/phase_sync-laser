#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<complex.h>


void main()
  {
     double r, theta;
     int i;
     FILE *fp;
     fp = fopen("rt.dat","w");
     
     for(i = -1.0;  i <=1.0; i+=0.1){
        for(theta = 0.0; theta <= 2.0*M_PI; theta+=0.5*M_PI){
            fprintf(fp, "%f   %20.12lf", i, theta);
        
        }

        }
      printf("process over !!");    
  
  }
