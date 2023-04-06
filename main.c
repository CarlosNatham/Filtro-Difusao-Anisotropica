#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "pgmfiles.h"
#include "diff2d.h"
#include <time.h> // inclusão da biblioteca time para calculos temporais de funçãos
#include <x86intrin.h> // inclusão da biblioteca intrin para calculos de ciclos de funçãos
#define FREQ 3500000000

#pragma intrinsic(__rdtsc)

//gcc -o fda pgmtolist.c pgmfiles.c diff2d.c main.c -lm

void main (int argc, char **argv) {
  char   row[80];
  float  **matrix;
  int i, j;
  FILE   *inimage, *outimage;
  long   imax;
  float  lambda;
  int result;
  eightBitPGMImage *PGMImage;


  // variaveis
  double tempo_primeira_fun;
  double tempo_segunda_fun;
  double tempo_terceira_fun;
  double tempo_quarta_fun;
  double tempo_anterior;
  double tempo_posterior;
  unsigned __int64 __rdtsc();

  /* ---- read image name  ---- */

  PGMImage = (eightBitPGMImage *) malloc(sizeof(eightBitPGMImage));

  if (!argv[1])
  {
    printf("name of input PGM image file (with extender): ");
    scanf("%s", PGMImage->fileName);
  }
  else
  {
    strcpy(PGMImage->fileName, argv[1]);
  }

  result = read8bitPGM(PGMImage);

  if(result < 0)
    {
      printPGMFileError(result);
      exit(result);
    }

  /* ---- allocate storage for matrix ---- */

  tempo_anterior = (double)clock(); // variavel para calcular tempo inicial da função

  matrix = (float **) malloc (PGMImage->x * sizeof(float *));
  if (matrix == NULL)
    {
      printf("not enough storage available\n");
      exit(1);
    }
  for (i=0; i<PGMImage->x; i++)
    {
      matrix[i] = (float *) malloc (PGMImage->y * sizeof(float));
      if (matrix[i] == NULL)
        {
	  printf("not enough storage available\n");
	  exit(1);
        }
    }

  tempo_posterior = (double) clock(); // variavel para calcular o tempo final da função
  tempo_primeira_fun = (tempo_posterior - tempo_anterior); // varivael para salvar o tempo final de execução
  i = __rdtsc();
  printf("\n%I64d ticks\n", i);
  printf("\no tempo levado na primeira fun%c%co foi de: %.30lf segundos", 135, 198, tempo_primeira_fun/FREQ);

  /* ---- read image data into matrix ---- */

  tempo_anterior = (double) clock(); // variavel para calcular tempo inicial da função

 for (i=0; i<PGMImage->x; i++)
    for (j=0; j<PGMImage->y; j++)
      matrix[i][j] = (float) *(PGMImage->imageData + (i*PGMImage->y) + j);

  tempo_posterior = (double) clock(); // variavel para calcular o tempo final da função
  tempo_segunda_fun = (tempo_posterior - tempo_anterior); // varivael para salvar o tempo final de execução
  i = __rdtsc();
  printf("\n%I64d ticks\n", i);
  printf("\no tempo levado na segunda fun%c%co foi de: %.30lf segundos", 135, 198, tempo_segunda_fun/FREQ);

  /* ---- process image ---- */

  printf("\ncontrast paramter lambda (>0) : ");
  //~ gets(row);  sscanf(row, "%f", &lambda);
  scanf("%f", &lambda);
  printf("number of iterations: ");
  //~ gets(row);  sscanf(row, "%ld", &imax);
  scanf("%ld", &imax);

  tempo_anterior = (double) clock(); // variavel para calcular tempo inicial da função

  for (i=1; i<=imax; i++)
    {
      printf("\niteration number: %3ld \n", i);
      diff2d (0.5, lambda, PGMImage->x, PGMImage->y, matrix);
    }

  tempo_posterior = (double) clock(); // variavel para calcular o tempo final da função
  tempo_terceira_fun = (tempo_posterior - tempo_anterior); // varivael para salvar o tempo final de execução
  i = __rdtsc();
  printf("\n%I64d ticks\n", i);
  printf("\no tempo levado na terceira fun%c%co foi de: %.30lf segundos", 135, 198, tempo_terceira_fun/FREQ);

  /* copy the Result Image to PGM Image/File structure */

  tempo_anterior = (double) clock(); // variavel para calcular tempo inicial da função

  for (i=0; i<PGMImage->x; i++)
    for (j=0; j<PGMImage->y; j++)
      *(PGMImage->imageData + i*PGMImage->y + j) = (char) matrix[i][j];

  tempo_posterior = (double) clock(); // variavel para calcular o tempo final da função
  tempo_quarta_fun = (tempo_posterior - tempo_anterior); // varivael para salvar o tempo final de execução
  i = __rdtsc();
  printf("\n%I64d ticks\n", i);
  printf("\no tempo levado na quarta fun%c%co foi de: %.30lf segundos", 135, 198, tempo_quarta_fun/FREQ);

  /* ---- write image ---- */

  printf("\nname of output PGM image file (with extender): ");
  scanf("%s", PGMImage->fileName);

  write8bitPGM(PGMImage);

  /* ---- disallocate storage ---- */

  for (i=0; i<PGMImage->x; i++)
    free(matrix[i]);
  free(matrix);

  free(PGMImage->imageData);
  free(PGMImage);
}


