#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>

#define EXECS 10


FILE *arqIn;

typedef struct {
	int index;
	double psUt;
}psUtOrd;

void AlocaMatriz();
void inserir_Sep(FILE *);
void ini_Colonia();
void LiberaMatriz();
void checa_Validade();
void calcula_Fitness(int);
void prep();
void tournament_Select();
void movement(int);
void discretize(int);
void size(int);
void psUtCalculate();
int cmpFunc(const void*, const void*);
bool isFeasible(int);
void dropAdd(int);
bool canAdd(int, int);
void writeBestResult(FILE*, char*, double, int);

int *tp_Recurso, **p_Recurso, *lim_Recurso;
double **colonia;
float *fric_surf;
float *atual_Size;
float *update_cosc;
int num_Itens, num_Rec, best_Sol, total_Val;
float *fitness;
int parent;
int num_Colonias;
float shear_Force = 2;
double *col_Aux;
psUtOrd *psUtOrder;
float energyLoss = 0.3;
double *energy;
double energyInt;
int **rc;
bool isStarve;
int *starvation;


int main(int argc, char *argv[]){

  int i;

  int cont = 0;

  int desl = 0;

  char entrada[40];

  strcpy(entrada, argv[1]);
  num_Colonias = atoi(argv[2]);
  energyInt = atoi(argv[3]);

  arqIn = fopen(entrada, "r");

  if(arqIn == NULL){
    printf("Arquivo de entrada não encontrado");
    exit(1);
  }
  
  fscanf(arqIn, "%d %d %d", &num_Itens, &num_Rec, &best_Sol);

  inserir_Sep(arqIn);

  AlocaMatriz();

  psUtCalculate();

  while(cont < EXECS){

    ini_Colonia();

    prep();


    for(int k = 0; k < num_Colonias; k++){

      isStarve = true;

      while(energy[k] > 0.0){

        calcula_Fitness(k);

        size(k);

        if(energy[k] < 0){
          continue;
        }

        tournament_Select();

        movement(k);

        if(isStarve){
          starvation[k]++;
        }
      }
    }
    cont++;
  }

  LiberaMatriz();

  fclose(arqIn);

  return 0;

}

void inserir_Sep(FILE *arquivo){ // função está dando falha na segmentação

  tp_Recurso = (int *) calloc(num_Itens, sizeof(int));
  p_Recurso = (int **) calloc((num_Itens), sizeof(int*));
  for(int i = 0; i < num_Itens; i++){
    p_Recurso[i] = (int*) calloc(num_Rec, sizeof(int));
  }

  lim_Recurso = (int *) calloc(num_Rec, sizeof(int));

  for (int j = 0; j < num_Itens; j++)
	  fscanf(arquivo, "%d", &(tp_Recurso[j]));

  //Leitura da matriz de pesos
  for (int j = 0; j < num_Rec; j++){
	  for (int k = 0; k < num_Itens; k++){
	    fscanf(arquivo, "%d", &(p_Recurso[k][j]));
    }
  }

  //Leitura do vetor de restrições
  for (int j = 0; j < num_Rec; j++){
	  fscanf(arquivo, "%d", &(lim_Recurso[j]));
  }
}

void AlocaMatriz(){

  int i, j;

  fitness = (float*) calloc(num_Colonias, sizeof(float));
  for(i = 0; i < num_Colonias; i++){
    fitness[i] = 0;
  }
  
  starvation = (int*) calloc(num_Colonias, sizeof(int));

  for(i = 0; i < num_Colonias; i++){
    starvation[i] = 0;
  }

  atual_Size = (float *) calloc(num_Colonias, sizeof(float));
  update_cosc = (float *) calloc(num_Colonias, sizeof(float));
  fric_surf = (float *) calloc(num_Colonias, sizeof(float));

  col_Aux = (double *) calloc(num_Itens, sizeof(double));

  colonia = (double**) calloc(num_Colonias, sizeof(double*));

  for(i = 0; i < num_Colonias; i++){
    colonia[i] = (double*) calloc(num_Itens, sizeof(double));
  }

  rc = (int **) calloc(num_Colonias, sizeof(int*));

  for(i = 0; i < num_Colonias; i++){
    rc[i] = (int*) calloc(num_Rec, sizeof(int));
  }

  energy = (double *) calloc(num_Colonias, sizeof(double));
  psUtOrder = (psUtOrd *) calloc(num_Itens, sizeof(psUtOrd));

  for(int k = 0; k < num_Colonias; k++){
    energy[k] = energyInt;
  }

}

void LiberaMatriz(){


  //Libera cada linha da matriz
  for (int i = 0; i < num_Colonias; i++)
	  free(colonia[i]);

  //Libera a matriz
  free(colonia);
  free(starvation);

  free(tp_Recurso);
  for (int i = 0; i < num_Rec; i++)
    free(p_Recurso[i]);
  free(p_Recurso);
  free(lim_Recurso);
  free(col_Aux);
  free(energy);
  free(psUtOrder);

  free(fric_surf);
  free(atual_Size);
  free(update_cosc);
  free(fitness);

  for (int i = 0; i < num_Colonias; i++)
    free(rc[i]); 
  free(rc);
  
}

void ini_Colonia(){
  int i, j, l;

  srand(time(NULL));

  for(i = 0; i < num_Colonias; i++){
    for(j = 0; j < num_Rec; j++){
      rc[i][j] = lim_Recurso[j];
    }
  }

  for (i = 0; i < num_Colonias; i++){
	  for (j = 0; j < num_Itens; j++){
      colonia[i][j] = rand() % 2;
      if(colonia[i][j] == 1.0){
        for(l = 0; l < num_Rec; l++){
          rc[i][l] = rc[i][l] - p_Recurso[l][j];
        }
      }
      if(!isFeasible(i)){
        dropAdd(i);
      }

    }
  calcula_Fitness(i);
  }
}

void calcula_Fitness(int index){

  int i;

  fitness[index] = 0;

  for(i = 0; i < num_Itens; i++){
    fitness[index] = fitness[index] + colonia[index][i] * tp_Recurso[i];
  }
}

void prep(){

  int i;

  for(i = 0; i < num_Colonias; i++){
    atual_Size[i] = 1.0;
  }
}

void tournament_Select(){
  int i, poolSize;
  int *pool;

  poolSize = 2;

  pool = (int*) calloc(poolSize, sizeof(int));

  for (i = 0; i < poolSize; i++){
    pool[i] = rand() % num_Colonias;
  }

  if(fitness[pool[0]]> fitness[pool[1]]){
    parent = pool[0];
  }
  else {
    parent = pool[1];
  }

  free(pool);

}

void movement(int index){
  double p;
  double alpha, beta;
  int m, k, l;
  double fit;

  // float rand1 = (rand() / (float) RAND_MAX);

  p = ((double)rand() / (RAND_MAX / 2) - 1);
  alpha = ((double)rand() / (RAND_MAX / (2 * M_PI)));
  beta = ((double)rand() / (RAND_MAX / (2 * M_PI)));

  m = rand() % num_Itens;
  k = rand() % num_Itens;
  l = rand() % num_Itens;

  fit = fitness[index];

  for(int i = 0; i < num_Itens; i++){
    col_Aux[i] = colonia[index][i];
  }

  colonia[index][m] = colonia[index][m] + (colonia[parent][m] - colonia[index][m]) * (shear_Force - fric_surf[index] * p);
	colonia[index][k] = colonia[index][k] + (colonia[parent][k] - colonia[index][k]) * (shear_Force - fric_surf[index] * cos(alpha));
  colonia[index][l] = colonia[index][l] + (colonia[parent][l] - colonia[index][l]) * (shear_Force - fric_surf[index] * sin(beta));

  discretize(index);

  dropAdd(index);

  energy[index] = energy[index] - energyLoss;

  calcula_Fitness(index);

  if(fit > fitness[index]){
    for(int i = 0; i < num_Itens; i++){
      colonia[index][i] = col_Aux[i];
    }
    isStarve = false;
    fitness[index] = fit;
  }
  else{
    energy[index] = energy[index] - energyLoss;
  }
}

void discretize(int index){

  double random;
  double gx;

  random = ((int)rand() / (RAND_MAX / 2)) - 1;

  for(int i = 0; i < num_Itens; i++){
    if(fmod(colonia[index][i], 1) != 0.0){
      gx = tanh(colonia[index][i]);
       if(gx < random){
        colonia[index][i] = 0;
       }
       else {
        colonia[index][i] = 1;
       }
    }
  }
}

void psUtCalculate(){
  int i, j;
  double delta;

  for(i = 0; i < num_Itens; i++){
    delta = 0.0;
    
    for(j = 0; j < num_Rec; j++){
      delta += p_Recurso[j][i];
    }
    psUtOrder[i].psUt = tp_Recurso[i] / delta;		
		psUtOrder[i].index = i;
  }
  qsort(psUtOrder, (size_t) num_Itens, sizeof(psUtOrd), cmpFunc);
}

int cmpFunc(const void *a, const void *b){
   //return ( *(int*)a - *(int*)b );
   //return (((AuxPeso*)a)->peso - ((AuxPeso*)b)->peso );
   
  if (((psUtOrd*)a)->psUt > ((psUtOrd*)b)->psUt)
	  return 1;

	if (((psUtOrd*)a)->psUt < ((psUtOrd*)b)->psUt)
		return -1;
		
	return 0;
}

void dropAdd(int j){

  int i, minIndex, maxIndex;
  int dispRc;
  bool found;

  dispRc = j * num_Itens;

  minIndex = 0;

  do{
    found = false;
    while(minIndex < num_Itens && !found){
      if(colonia[j][psUtOrder[minIndex].index] == 1){
        found = true;
      }
      else{
        minIndex++;
      }

    }
    colonia[j][psUtOrder[minIndex].index] = 0;
    for(i = 0; i < num_Rec; i++){

      rc[j][i] = rc[j][i] + p_Recurso[i][psUtOrder[minIndex].index];

    }
  } while(!isFeasible(j));

  for(maxIndex = num_Itens; maxIndex >= 0; maxIndex--){
    if(canAdd(j , psUtOrder[maxIndex].index)){

      colonia[j][psUtOrder[maxIndex].index] = 1;
      for(i = 0; i < num_Rec; i++){

        rc[j][i] = rc[j][i] - p_Recurso[i][psUtOrder[maxIndex].index];

      }
    }
    else{
      colonia[j][psUtOrder[maxIndex].index] = 0;
    }
  }
}

bool canAdd(int j,int index){
  int i;

  for(i = 0; i < num_Rec; i++){
    if(rc[j][i] - p_Recurso[i][index] < 0){
      return false;
    }
    return true;
  }

}

bool isFeasible(int index){
  int i, displRc;

  displRc = index * num_Rec;

  for(i = 0; i < num_Rec; i++){
    if(rc[index][i] < 0){
      return false;
    }
    return true;
  }

}

void size(int index){

  double raiz, aux;

  update_cosc[index] =  (atual_Size[index] + 4*(fitness[index]))/
                      (atual_Size[index] + 2*(fitness[index]));
  aux = (3 * atual_Size[index])/(4 * M_PI);
  raiz = pow(aux, 1.0/3.0);
  fric_surf[index] = 2 * M_1_PI * (raiz);

  atual_Size[index] = update_cosc[index] * atual_Size[index];
}