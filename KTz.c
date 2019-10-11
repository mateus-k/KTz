#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
 
double Ilin(long long int t); int delta(long long int t,long long int t0); int theta(double x);
  
main(void){

  double x0, y0, z0, *x, *y, *z, **par, xAUX, **X, **Xmax, J, Jx, Iel, I0, I, *f, *g, tf, tg;
  long long int t, t0, tMAX, tTRANS, dt, n, i, j, k, l, **tk, m, s, dphi, T;
  FILE *datai, *dataj, *J_dphi, *phiIJ, *picosI, *picosJ, *J_dphiINT;
  
  //Tempo maximo
   tMAX =  500000;
  //Tempo transiente
   tTRANS= 400000;

  //dimensionando as matrizes
   n=2;//==> Quantidade de neuronios
   x=malloc(n*sizeof(double)); 
   y=malloc(n*sizeof(double)); 
   z=malloc(n*sizeof(double)); 
   g=malloc(n*sizeof(double)); f=malloc(n*sizeof(double));
   if (x==NULL||y==NULL||z==NULL||f==NULL||g==NULL){
     printf("Erro ao dimensionar matrizes\n");
     exit(EXIT_FAILURE);
   }
   X=malloc(tMAX*sizeof(double *)); par=malloc(5*sizeof(double *)); Xmax=malloc(tMAX*sizeof(double *)); 
   tk=malloc(tMAX*sizeof(long long int *));
   if (X==NULL||par==NULL||Xmax==NULL||tk ==NULL){
     printf("ERRO NA 1ª ALOCACAO DE MEMÓRIA DE ''array's'' !!\n");
     exit(1);
   }
   for (t=0;  t<5; t++){
    par[t]=malloc(n*sizeof(double));
    if (par[t]==NULL){
      printf("ERRO NA 2ª ALOCACAO DE MEMÓRIA NA MATRIZ dos PARAMETROS!!\n");
    }
  }
  for (t=0 ; t<=tMAX; t++){
   X[t]=malloc((n+1)*sizeof(double)); Xmax[t]=malloc((2*n)*sizeof(double));
   tk[t]=malloc(n*sizeof(long long int));
   if (X[t]==NULL||Xmax[t]==NULL||tk[t]==NULL){
     printf("ERRO NA 2ª ALOCACAO DE MEMÓRIA DE ''array's'' !!\n");
     exit(1);
   }
  }
  //---------------------------------------------------------------------
  
  //Condições Iniciais
   x0 = -0.6971564118917724;           
   y0 = -0.6971564118917724;           
   z0 = -0.022748704865822533;        
   
  //MATRIZ dos PARAMETROS    
  /* K */      par[0][0] =  0.6;      par[0][1] =  0.6;       
  /* T */      par[1][0] =  0.35;     par[1][1] =  0.35; 
  /* delta */  par[2][0] =  0.004;    par[2][1] =  0.004; 
  /* lambda */ par[3][0] =  0.003;    par[3][1] =  0.003;
  /* xR */     par[4][0] = -0.632;    par[4][1] = -0.632; 
 
  //Intensidade do estimulo inicial 
   I0=0.0;
  
  //Inicio do estimulo inicial
   t0=10;
   dt=10;
  
  //synaptic delay s & memory m
   m=6;
   s=4;
  
  //Paramentro de acoplamento para salvar arquivos 
   Jx = 0.00273;
   
  //Imprimindo parametros e condições inciais
   printf("Parametros do pré\nki = %lf\nTi = %lf\nδi = %lf\nλi = %lf\nxRi = %lf\n\n",par[0][0],par[1][0],par[2][0],par[3][0],par[4][0]);
   printf("Parametros do pós\nkj = %lf\nTj = %lf\nδj = %lf\nλj = %lf\nxRj = %lf\n\n",par[0][1],par[1][1],par[2][1],par[3][1],par[4][1]);
   printf("Condições Iniciais\nx0 = %.15lf\ny0 = %.15lf\nz0 = %.15lf\n\n",x0,y0,z0);
  //--------------------------------------------------------- 
 
  //Abrindo arquivos
  J_dphi = fopen("Jxdphi.dat","w");
  phiIJ = fopen("phiIJ.dat","w"); 
  J_dphiINT = fopen("JxphiINT.dat","w");
  fprintf(phiIJ,"J        phiI    phiJ    dphi\n"); 
  if (J_dphi==NULL||phiIJ==NULL||J_dphiINT==NULL){
   printf("ERRO!\nNão foi possível gerar dados nos arquivos.Certifique se há espaço suficiente.\n\n");
    exit(EXIT_FAILURE);
  }
  
  
  
  for (J=0.0; J<=0.015; J = J + 0.00001){
    
    x[0]=x[1]=x0; y[0]=y[1]=y0; z[0]=z[1]=z0;  //Inicializa sempre com as mesmas condições iniciais
   
   //Zerando variaveis
    for (l=0; l<n; l++){
      f[l] = g[l] = 0.0;
      for (t=0; t<tMAX; t++){
	X[t][l] = X[t][n+l] = 0.0;
	Xmax[t][l] = Xmax[t][n+l] = 0.0;
	tk[t][l] = 00; 
      }
    }
    Iel = 0.0;
    //---------------------

    
    
    for (t=0; t<=tMAX; t++){
     
      X[t][0]=t;     //armazena na coluna 0 os tempos
      X[t][1]=x[0];  //armazena na coluna 1 o potencial do pre
      X[t][2]=x[1];  //armazena na coluna 2 o potencial do pos
      
      //Sinapse Elétrica
      if ((t-s)<0 || (t-m)<0){
	Iel=0.0;
      } else {
	Iel=J*(X[t-s][1]-X[t-m][2]);
      }
      
      I=I0*(delta(t,t0) + delta(t,t0+dt));  
     
      /*//Sinpases Químicas
      f[0] = (1.0 - 1.0/tf)*f[0] + g[0];
        g[0] = (1.0 - 1.0/tg)*g[0] + J*theta(x[1]);
     
      f[1] = (1.0 - 1.0/tf)*f[1] + g[1];
      g[1] = (1.0 - 1.0/tg)*g[1] + J*theta(x[0]);
      //-------------------------------------*/
      
      //Neuronios KTz
      xAUX=x[0];
      x[0] = tanh((x[0] - par[0][0]*y[0] + z[0] + I)/par[1][0]);
      y[0] = xAUX;
      z[0] = (1.0-par[2][0])*z[0] - par[3][0]*(xAUX-par[4][0]);
     
      xAUX=x[1];
      x[1]=tanh((x[1] - par[0][1]*y[1] + z[1] + Iel)/par[1][1]);
      y[1]=xAUX;
      z[1] = (1.0-par[2][1])*z[1] - par[3][1]*(xAUX-par[4][1]);
    }
   
    //Acha TODOS os pontos de máximo locais 
    for (l=0; l<n; l++){
     i=0;
      for (t=tTRANS+1; t<tMAX; t++){
        if (X[t-1][l+1] < X[t][l+1] && X[t][l+1] > X[t+1][l+1]){
	  Xmax[i][l]=X[t][0];        //armazena o t dos máximos na coluna 0 (Xi) ou 1 (Xj)
	  Xmax[i][n+l]=X[t][l+1];    //armazena o X dos máximos na coluna 2 (Xi) ou 3 (Xj)
	  i++;
        }
      }
    }
    //-------------------------------------------
   
    //Acha os pontos de spike para calculo da fase
    for (l=0; l<n; l++){
      k=0;
      for (j=1; j<tMAX; j++){
	if (Xmax[j][n+l] > 0 && Xmax[j-1][n+l] < 0){
	  tk[k][l]=Xmax[j][l];            //armazena os t dos spikes de fase na coluna 0 (Xi) ou 1 (Xj) 
	  Xmax[k][n+l]=Xmax[j][n+l];      //Recupera o valor do Xmax para a linha k
	  //printf("tk[%Ld][%Ld]=%Ld\n",k,l,tk[k][l]);
	  ++k;
	}
      }
    }
    
    for (j=1; j<tMAX; j++){
      if (tk[j][0] !=0 && tk[j][1] !=0){
	T = tk[j+1][0]-tk[j][0];
	dphi=tk[j][1]-tk[j][0];
	//printf("Ti=%Ld Tj=%Ld\n",Ti,Tj);
	//fprintf(phiIJ,"#\n# %lf %Ld %Ld %Ld %Ld %Ld\n#\n",J,tk[j][0],tk[j][1],tk[j][1]-tk[j][0]);
	fprintf(J_dphiINT,"%lf %Ld\n",J,tk[j][1]-tk[j][0]);
	if (llabs(dphi)<=T){
	fprintf(phiIJ,"%lf %Ld %Ld %Ld %Ld\n",J,tk[j][0],tk[j][1],dphi,T);
	fprintf(J_dphi,"%lf %Ld\n",J,tk[j][1]-tk[j][0]);
	}
      }
    }
  
    printf("J=%lf\n",J);
  }
  fclose(J_dphi); fclose(phiIJ); fclose(J_dphiINT);
  
  printf("dados gerados no arquivo Jxdphi.dat\n");   
}
//============================================================================================

//FUNÇÕES UTILIZADAS

//Linear varying current
double Ilin(long long int t){
  
  double I, a, I0;
  long long int t0, t1; 
  
  //Parametros
   I0 = 0.3;
   a  = 0.1;
   t0 = 10;
   t1 = 20;
  
  if (t >= t0 && t <= t1) {   
   I= a*(t-t0) + I0;
  } else {
   I=0; 
  }  
  
  return I;
}


//Delta de Kronecker
 int delta(long long int t,long long int t0){
  
  int d;
   
  if (t == t0){
    d=1;
  } else {
    d=0;
  } 
  
  return d;
}

//Theta de Heaviside
int theta(double x){
  
  int theta;
  
  if (x < 0){
    theta = 0;
  } else {
    theta = 1;
  }
  
  return theta;
}
