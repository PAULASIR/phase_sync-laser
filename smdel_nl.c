/* BVP network nonlinear coupling (memductance) */
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<complex.h>


int i,j,k,nn=20,n=5;
double y[101][501], x_tau, r_tau, J_0,g , R_sp, s, G, alpha;
double complex phi, N_0;
double sm1, sm, eps, sum, sum1;
double gama;
double acf;
double omega, gama_n;
double syn_er[30000], z1[30000];

double tau_r, tau_phi;

int p;  //# of nearest neighbours
double z[35000],zz[35000];


void RK4(int,int,double,double,double[101][501],
                   void (*DGL)(double,double[101][501],double[101][501]));
void DGL(double, double[101][501],double[101][501]);



/******************************/
double getsum(double z[], int nn){
    int i; 
    double cumul =0.0;
    for(i = 1; i<=nn; i++) 
        cumul += z[i];
       cumul = cumul/(float) nn;
      //  printf("%f \n",cumul);
       // fprintf(fp,"%d  %f \n", i, z[i]);}
 return 0; 
}
/*********************/

/*******************************************/

double sn_sum(double syn_er[], int nn){
    int i;     double avg_er=0.0; 
    for(i = 1; i<=nn; i++) 
       avg_er = syn_er[i]/(float)nn;
      // printf("%20.16lf  \n", avg_er);
 return avg_er; 
}
/********************************/

double getsum1(double zz[], int nn){
    int i; 
    double cm1 =0.0;
    for(i = 1; i<=nn; i++) 
       cm1 += zz[i];
        cm1 = cm1/(float)nn;
       // fprintf(fp,"%d  %f \n", i, z[i]);}
 return 0; 
}


/**********history function (phi)*******/
double delay(double x_tau )
{
   //double x_tau;

   x_tau = y[j][4];

 return x_tau;
}
/****************************************/


/**********history function (r_spatial coordinate)*******/
double delay1(double r_tau)
{
  // double r_tau;

   r_tau = y[j][5];

 return r_tau;
}
/****************************************/


/************ function for gaussian noise******************/
double rand_normal(double mean, double stddev)
 {//Bo262 muller method
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
     {
        double x, y, r;
         do
          {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;

            r = x*x + y*y;
          }

        while (r == 0.0 || r > 1.0);
         {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
         }
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}

/*--------------------------------------------------------------------*/
void main()
{
//nn= no 0f oscillators, n=dimension of the model//
double t,h,pi;
double amplitude,vmin,vmax,c_max,c_min;
double x_max,x_min;
float cnt;
srand(time(0));

FILE *fp1,*fp2,*fp3, *fp4, *fp5;
fp1=fopen("ts.dat","w");
fp2=fopen("chm.dat","w");
fp3=fopen("ldl_cc.dat","w");
 fp4 = fopen("syn_avg.dat","w");
 fp5 =fopen("chm1.dat","w");
 
 for(j=1;j<=nn;j++)
 {
  y[j][1]=(float) rand()/(double)RAND_MAX*-4.0+2.0 ;
  y[j][2]= (float) rand()/(double)RAND_MAX*-4.0-2.0;
  y[j][3]= (float) rand()/(double)RAND_MAX*-4.0-2.0;
  y[j][4] =  tau_r;
  y[j][5] =  tau_phi;
 }


//parameters
scanf("%lf",&eps);

//eps = 0.5;
 omega= 12.0;   //5.0
 G = 1.5e-8;
 N_0 = 1.5e8;
 s = 2.0e-7;
 gama = 0.5;
 R_sp = 10.0;
 
 J_0 = 4.0*gama*(N_0 - (gama/G));
 tau_r   =  -1.0;  //-0.5 
 tau_phi = -1.5; // -0.8


 
//***time step***//
h=0.01; t=0.0;
x_max=x_min=0.0;


sm = 0.0;  sm1 = 0.0; syn_er[k] = 0.0;  
for(k=1;k<=6000;k++)
	{               
        t=h*(double)(k);
    	RK4(nn,n,h,t,y,DGL);   
    	 cnt =0.0;          
   for(j=1;j<=nn;j++)
	{  
/*    leaving transients*/               
	//if(k>=4900)
	 {          
/*----------------max & min of the series--------
       if(y[j][1]>x_max){x_max=y[j][1];} 
       if(y[j][1]<x_min){x_min=y[j][1];}   */
                                                     
        z[k] = y[j][2];    zz[k] = delay(x_tau);
        
        sm += z[k];   sm1 += zz[k];
   
     acf = (z[k] - getsum(z,nn)) * (zz[k] - getsum1(zz,nn))/  ( (z[k] - getsum(z,nn)) *  (z[k] - getsum(z,nn))) ;
        // getsum(z,nn);
        
  
   //syn_er[k]  += fabs(z[k] - zz[k])/100; 
   
    z1[k]  =  cexp(I*y[j][2]);
    z1[k] =   carg(z1[k]);//z2[k] = y[j-1][2];
   
   
     syn_er[k] += fabs(y[j][1] -y[j+1][1]);
// printf("%f \n", acf);
   
/* -----------------Printing time series-----------------*/        
       if(k == 3000)
          fprintf(fp3,"%d %f   %f %f\n",j, z1[k], y[j][2]/(float) 1e6, y[j][1]/(float) 1e6);
   
         
         for(i =1;i<=nn;i++)
          fprintf(fp2,"  %f  %f  \t", y[i][1], y[i][2]); 
            fprintf(fp2," \n"); 
            
   
   
   fprintf(fp5,"%f   %f\n", t ,y[12][1]);
   
   
         fprintf(fp1,"%f \t",t);  
           for(i=1; i<=nn; i++)
             fprintf(fp1,"%f \t",y[i][1]/(float)1e6);    
             fprintf(fp1," \n"); 
             
           
              fprintf(fp4,"%f  %f\n", sum+=((y[j][1]/(float) 10000)/nn), sum1+=((y[j][2]/(float) 10000) /nn));
    
/**************************************************/    
         // }

	}    
}

     //fprintf(fp4,"%f  %f\n",t, sum);
}

   printf("%f \n",sn_sum(syn_er,nn));

    // getsum(z,nn);
   //  getsum1(zz,nn);

  printf("process over!!\n");
}
//************************RK4 SUBROUTINE*********************************//
void RK4(int nn,int n,double h,double t,double y[101][501],
	   void (*DGL)(double,double[101][501],double[101][501]))
{
	   
	   double k1[101][501],k2[101][501],k3[101][501],k4[101][501];
	   double yaux[101][501];

	   DGL(t,y,k1);
	   for(j=1;j<=nn;j++)
	   {
            for(i=1;i<=n;i++)
	    yaux[j][i]=y[j][i]+h*k1[j][i]/2.0;
	   }
	   
	   DGL(t+h/2.0,yaux,k2);
	   for(j=1;j<=nn;j++)
	   {
            for(i=1;i<=n;i++)
	    yaux[j][i]=y[j][i]+h*k2[j][i]/2.0;
	   }
	   
	   DGL(t+h/2.0,yaux,k3);
	   for(j=1;j<=nn;j++)
	   {
            for(i=1;i<=n;i++)
	    yaux[j][i]=y[j][i]+h*k3[j][i];
	   }
	   
	   DGL(t+h,yaux,k4);
	   for(j=1;j<=nn;j++)
	   {
             for(i=1;i<=n;i++)
	      y[j][i]=y[j][i]+h*((k1[j][i]+2*k2[j][i]+2*k3[j][i]+k4[j][i])/6.0);
	   }
}
//*********************FUNCTION SUBROUTINE********************************//
void DGL(double t,double y[101][501],double F[101][501])
{

    double meanx,A;
    double aa[600][600];  
     double zeta;  
      p=8;
      double eps1;

    meanx=0.0; 

    for(j=1;j<=nn;j++)
     {
       meanx=meanx+y[j][1];
     }

    meanx=meanx/nn; 
    
    
    phi = cexp(y[j][2]*I);
    phi = carg(phi);
    
    
   for(j=1;j<=nn;j++)
     { 
     for(i=1;i<=nn;i++)
      {
         if(i != j && abs(i-j)%(nn-p)<=p) 
           aa[j][i]=1;
          else
            aa[j][i]=0; 
          
            
                
         
      // printf("%f %f \n", fabs(acf), eps);   
         
          zeta = rand_normal(0.0,1.0);
 
                if ( fabs(acf) <= 0.2)
                 {
                      eps1 = eps;}
                 //  printf("%12.20f   %12.20f\n", fabs(acf), eps1);}
                else 
                   { eps1 = 0.0; }  
 
   
F[j][1]= 0.5*(G*(y[j][3]- N_0)/(1+s*y[j][1]*y[j][1]))*y[j][1] + sqrt(R_sp)*zeta*creal(phi) + eps1/(2.0*p)*(y[j][1] - delay1(r_tau))*aa[j][i]*(cos(y[j][2]) - cos(delay(x_tau)));
F[j][2]= alpha/2.0*((g*y[j][3] - N_0)/(1+s*y[j][1]*y[j][1]) - gama)+omega -(sqrt(R_sp)/ y[j][1])*zeta*cimag(phi)+eps1/(2.0*p)*aa[j][i]*y[j][1]-(delay1(r_tau)/ y[j][1])*aa[j][i]*(sin(delay(x_tau)- y[j][2])); 
F[j][3]= J_0 - gama*y[j][3] - (G*(y[j][3]- N_0)/(1+s*y[j][1]*y[j][1]))*y[j][1]*y[j][1] + sqrt(gama_n*y[j][3]*zeta);
    F[j][4]= tau_phi;
    F[j][5] =  tau_r;
     }
}
}

