
/****************laser time-delay model***************/

#include<stdio.h>
#include<math.h>
#include<time.h>
#include<stdlib.h>
#include<float.h>
#include<complex.h>



#define N 500
double lambda, beta,tau, omega, a,b,c, om, gama, delta, alpha,mu;
double y[20],x,xx;


void RK4(int,double,double,double[],void (*DGL)(double,double[],double[]));
void DGL(double, double[],double[]);


double phi(double x)
{
  double x1;

  x1 = y[3];

return x;
}
/****************/

void main()
 {
   long int i,j;
   int nn=3; 
   double t,h,pi;
   double a,b,c;
//clrscr();
FILE *fp;
fp=fopen("trial1.dat","w");

srand(time(0));

beta = 1.0;
omega =0.1;
gama = 2.0;
lambda = 0.8; 
alpha = 0.5;
mu = 0.01;


tau = -2.5;


 y[1] = (float) rand()/(double)RAND_MAX*-2.0+4.0;  
 y[2] = (float) rand()/(double)RAND_MAX*-2.0+4.0; 
 y[3] = tau;

   
  h= 0.01;


for(i=1;i<=30000;i++)
 {
   t=h*(double)(i);
   RK4(nn,h,t,y,DGL);
  // if(i>=5000)
   fprintf(fp," %f  %f %f \n", y[1],y[2],t);
}

 printf("process over!!\n");
 //getch();
}
//***********************EULER SUBROUTINE*********************************//
void RK4(int nn,double h,double t,double y[20],
	   void (*DGL)(double,double[],double[]))
{
	   int i;
	   double k1[10];
	   double yaux[20];
 
	   DGL(t+h,y,k1);
	   for(i=1;i<=nn;i++)
	   {
	      y[i]=y[i]+h*k1[i];
	   }
}

//*********************FUNCTION SUBROUTINE********************************//
void DGL(double t,double y[20],double F[10])
{


  F[1] = y[2];
  F[2] = -omega*y[2]  -(y[1] -phi(x)) - alpha*pow(y[1],3) + 0.2*sin(omega*t);
  F[3] =  tau;
}

//-0.5*y-x*x*x+x+0.42*sin(t);

//$$$$$$$$$$$$$$$$$$   END  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$//

//-2.0*gama*(y[1]-phi(x))- 2.0*lambda*(y[1]-phi(x))
//-2.0*mu*y[1]
