#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "diff2d.h"
#include <time.h> // inclusão da biblioteca time para calculos temporais de funçãos
#include <x86intrin.h> // inclusão da biblioteca intrin para calculos de ciclos de funçãos
#define FREQ 3500000000

#pragma intrinsic(__rdtsc)

  // variaveis
  double tempo_quinta_fun;
  double tempo_sexta_fun;
  double tempo_setima_fun;
  double tempo_oitava_fun;
  double tempo_anterior;
  double tempo_posterior;
  unsigned __int64 __rdtsc();
  float Lut[256];


/*--------------------------------------------------------------------------*/


float dco (float v,         /* value at one point */
           float w,         /* value at the other point */
           float lambda)    /* contrast parameter */

/* diffusivity */

{
    float result = 0.0, temp_result = 0.0;

    temp_result = (float)fabs(v-w);
    temp_result = pow(temp_result,0.2);
    temp_result = temp_result/lambda;
    if(temp_result != 0.0){
        temp_result = -(temp_result/5.0);
    }
    result = exp(temp_result);
    //result = exp( - (pow(fabs(v-w),0.2)/lambda) / 5.0);

    return (result);
}


/*--------------------------------------------------------------------------*/


void diff2d

     (float    ht,        /* time step size, >0, e.g. 0.5 */
      float    lambda,    /* contrast parameter */
      long     nx,        /* image dimension in x direction */
      long     ny,        /* image dimension in y direction */
      float    **f)       /* input: original image ;  output: smoothed */


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*             NONLINEAR TWO DIMENSIONAL DIFFUSION FILTERING                */
/*                                                                          */
/*                       (Joachim Weickert, 7/1994)                         */
/*                                                                          */
/*--------------------------------------------------------------------------*/


/* Explicit scheme with 9-point stencil and exponential stabilization.      */
/* Conservative, conditionally consistent to the discrete integration       */
/* model, unconditionally stable, preserves maximum-minimum principle.      */


{

long    i, j;                                     /* loop variables */
float   qC, qN, qNE, qE, qSE, qS, qSW, qW, qNW;   /* weights */
float   **g;                                      /* work copy of f */


/* ---- allocate storage for g ---- */

g = (float **) malloc ((nx+2) * sizeof(float *));
if (g == NULL)
   {
     printf("not enough storage available\n");
     exit(1);
   }

  tempo_anterior = (double) clock();

for (i=0; i<=nx+1; i++)
    {
      g[i] = (float *) malloc ((ny+2) * sizeof(float));
      if (g[i] == NULL)
         {
           printf("not enough storage available\n");
           exit(1);
         }
    }

    tempo_posterior = (double) clock();
    tempo_quinta_fun = (tempo_posterior - tempo_anterior);
    i = __rdtsc();
    printf("\n%I64d ticks\n", i);
    printf("\no tempo levado na quinta fun%c%co foi de: %.30lf segundos", 135, 198, tempo_quinta_fun/FREQ);

/* ---- copy f into g ---- */

    tempo_anterior = (double) clock();

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     g[i][j] = f[i-1][j-1];

    tempo_posterior = (double) clock();
    tempo_sexta_fun = (tempo_posterior - tempo_anterior);
    i = __rdtsc();
    printf("\n%I64d ticks\n", i);
    printf("\no tempo levado na sexta fun%c%co foi de: %.30lf segundos", 135, 198, tempo_sexta_fun/FREQ);


/* ---- create dummy boundaries ---- */

    tempo_anterior = (double) clock();

for (i=1; i<=nx; i++)
    {
     g[i][0]    = g[i][1];
     g[i][ny+1] = g[i][ny];
    }

for (j=0; j<=ny+1; j++)
    {
     g[0][j]    = g[1][j];
     g[nx+1][j] = g[nx][j];
    }

    tempo_posterior = (double) clock();
    tempo_setima_fun = (tempo_posterior - tempo_anterior);
    i = __rdtsc();
    printf("\n%I64d ticks\n", i);
    printf("\no tempo levado na setima fun%c%co foi de: %.30lf segundos", 135, 198, tempo_setima_fun/FREQ);

/*------------------------------*/

for(int i=0; i<256;i++){
    Lut[i] = (1.0 - exp(-8.0 * ht * dco(i, 0, lambda))) / 8.0;
}

/* ---- diffusive averaging ---- */

    tempo_anterior = (double) clock();

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)

     {
       /* calculate weights */

       qN  = Lut[(int)fabs(g[i][j]- g[i][j+1])]; /*(1.0 - exp(-8.0 * ht * dco(g[i][j], g[i  ][j+1], lambda))) / 8.0;*/
       qNE = Lut[(int)fabs(g[i][j]- g[i+1][j+1])];/*(1.0 - exp(-8.0 * ht * dco(g[i][j], g[i+1][j+1], lambda))) / 8.0;*/
       qE  = Lut[(int)fabs(g[i][j]- g[i+1][j])];/*(1.0 - exp(-8.0 * ht * dco(g[i][j], g[i+1][j  ], lambda))) / 8.0;*/
       qSE = Lut[(int)fabs(g[i][j]- g[i+1][j-1])];/*(1.0 - exp(-8.0 * ht * dco(g[i][j], g[i+1][j-1], lambda))) / 8.0;*/
       qS  = Lut[(int)fabs(g[i][j]- g[i][j-1])];/*(1.0 - exp(-8.0 * ht * dco(g[i][j], g[i  ][j-1], lambda))) / 8.0;*/
       qSW = Lut[(int)fabs(g[i][j]- g[i-1][j-1])];/*(1.0 - exp(-8.0 * ht * dco(g[i][j], g[i-1][j-1], lambda))) / 8.0;*/
       qW  = Lut[(int)fabs(g[i][j]- g[i-1][j])];/*(1.0 - exp(-8.0 * ht * dco(g[i][j], g[i-1][j  ], lambda))) / 8.0;*/
       qNW = Lut[(int)fabs(g[i][j]- g[i-1][j+1])];/*(1.0 - exp(-8.0 * ht * dco(g[i][j], g[i-1][j+1], lambda))) / 8.0;*/
       qC  = 1.0 - qN - qNE - qE - qSE - qS - qSW - qW - qNW;


       /* weighted averaging */

       f[i-1][j-1] = qNW * g[i-1][j+1] + qN * g[i][j+1] + qNE * g[i+1][j+1] +
                     qW  * g[i-1][j  ] + qC * g[i][j  ] + qE  * g[i+1][j  ] +
                     qSW * g[i-1][j-1] + qS * g[i][j-1] + qSE * g[i+1][j-1];

     }  /* for */

    tempo_posterior = (double) clock();
    tempo_oitava_fun = (tempo_posterior - tempo_anterior);
    i = __rdtsc();
    printf("\n%I64d ticks\n", i);
    printf("\no tempo levado na oitava fun%c%co foi de: %.30lf segundos", 135, 198, tempo_oitava_fun/FREQ);


/* ---- disallocate storage for g ---- */

for (i=0; i<=nx+1; i++)
    free(g[i]);
free(g);

return;

} /* diff */


/*--------------------------------------------------------------------------*/

