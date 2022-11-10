#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>


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
void discretize();
void size();
void psUtCalculate();
int cmpFunc(const void*, const void*);
bool isFeasible(int);
void dropAdd(int);
bool canAdd(int, int);

int *tp_Recurso, *p_Recurso, *lim_Recurso;
double **colonia;
float *fric_surf;
float *atual_Size;
float *update_cosc;
int num_Itens, num_Rec, best_Sol, total_Val;
float *fitness;
int parent;
int num_Colonias;
float *m, *k, *l;
float shear_Force = 2;
double *col_Aux;
psUtOrd *psUtOrder;
float energyLoss = 0.3;
double *energy;
double energyInt;
int **rc;


int main(int argc, char *argv[]){

  int i;

  int cont = 0;

  int desl = 0;

  char entrada[40];

  strcpy(entrada, argv[1]);
  int num_Testes = atoi(argv[2]);
  num_Colonias = atoi(argv[3]);
  energyInt = atoi(argv[4]);

  arqIn = fopen(entrada, "r");

  if(arqIn == NULL){
    printf("Arquivo de entrada não encontrado");
    exit(1);
  }
  
  fscanf(arqIn, "%d %d %d", &num_Itens, &num_Rec, &best_Sol);

  total_Val = num_Itens * num_Rec + num_Itens + num_Rec;

  m = (float *) malloc(num_Colonias * sizeof(float));
  k = (float *) malloc(num_Colonias * sizeof(float));
  l = (float *) malloc(num_Colonias * sizeof(float));
  energy = (double *) malloc(num_Colonias * sizeof(double));
  psUtOrder = (psUtOrd *) malloc(num_Itens * sizeof(psUtOrd));

  for(int k = 0; k < num_Colonias; k++){
    energy[k] = energyInt;
  }

  inserir_Sep(arqIn);

  AlocaMatriz();

  psUtCalculate();

  while(cont < num_Testes){

    ini_Colonia();

    prep();

    //printf("teste\n");

    do{

      if(energy[desl] < 0){
        desl++;
        continue;
      }

      //printf("%d\n", desl);
      //printf("%lf\n", energy[desl]);

      tournament_Select();

      movement(desl);
      
      energy[desl] = energy[desl] - energyLoss;
      
    }while(energy[desl] > 0.0);

    cont++;
  }

  LiberaMatriz();

  fclose(arqIn);

}

void inserir_Sep(FILE *arquivo){ // função está dando falha na segmentação

  tp_Recurso = (int *) malloc(num_Itens * sizeof(int));
  p_Recurso = (int *) malloc((num_Itens * num_Rec) * sizeof(int));
  lim_Recurso = (int *) malloc(num_Rec * sizeof(int));

  for (int j = 0; j < num_Itens; j++)
	  fscanf(arquivo, "%d", &(tp_Recurso[j]));

  //Leitura da matriz de pesos
  for (int j = 0; j < num_Rec; j++){
	  for (int k = 0; k < num_Itens; k++){
	    fscanf(arquivo, "%d", &(p_Recurso[j * num_Itens + k]));
    }
  }

  //Leitura do vetor de restrições
  for (int j = 0; j < num_Rec; j++)
	  fscanf(arquivo, "%d", &(lim_Recurso[j]));
}

void AlocaMatriz(){

  int i, j;

  fitness = (float*) malloc(num_Colonias * sizeof(float));
  for(i = 0; i < num_Colonias; i++){
    fitness[i] = 0;
  }
  
  atual_Size = (float *) malloc(num_Colonias * sizeof(float));
  update_cosc = (float *) malloc(num_Colonias * sizeof(float));
  fric_surf = (float *) malloc(num_Colonias * sizeof(float));

  col_Aux = (double *) malloc(num_Itens * sizeof(double));

  colonia = (double**) malloc(num_Colonias * sizeof(double*));

  for(i = 0; i < num_Colonias; i++){
    colonia[i] = (double*) malloc(num_Itens * sizeof(double));
  }

  rc = (int **) malloc(num_Colonias * sizeof(int*));

  for(i = 0; i < num_Colonias; i++){
    rc[i] = (int*) malloc(num_Rec * sizeof(int));
  }

}

void LiberaMatriz(){


  //Libera cada linha da matriz
  for (int i = 0; i < num_Itens; i++)
	  free(colonia[i]);

  //Libera a matriz
  free(colonia);

  free(l);
  free(k);
  free(m);

  free(tp_Recurso);
  free(p_Recurso);
  free(lim_Recurso);
  free(col_Aux);
  free(energy);
  free(psUtOrder);

  free(fric_surf);
  free(atual_Size);
  free(update_cosc);
  free(fitness);

  for (int i = 0; i < num_Rec; i++)
    free(rc[i]);
    
  free(rc);
  
}

void ini_Colonia(){
  int i, j, l;

  srand(time(NULL));

  for(i = 0; i < num_Colonias; i++){
    for(j = 0; j < num_Rec; j++){
      rc[i][j] = lim_Recurso[j];
      //printf("%d ", rc[i][l]);
    }
  }

  for (i = 0; i < num_Colonias; i++){
	  for (j = 0; j < num_Itens; j++){
      colonia[i][j] = rand() % 2;
      //printf("%lf ", colonia[i][j]);
      if(colonia[i][j] == 1){
        for(l = 0; l < num_Rec; l++){
          rc[i][l] = rc[i][l] - p_Recurso[l * num_Itens + j];
          //printf("%d ", rc[i][l]);
        }
      }
    }
    //printf("\n");
    if(!isFeasible(i)){
      dropAdd(i);
    }
    calcula_Fitness(i);
  }
}

/*void checa_Validade(){
  int aux[num_Colonias][num_Itens];
  bool valido = true;

  for(int i = 0; i < num_Colonias; i++){
    for(int k = 0; k < num_Rec; k++){
      aux[i][k] = 0;
      }
    }

  for(int i = 0; i < num_Colonias; i++){
    valido = true;

    for(int j = 0; j < num_Itens; j++){
      if(colonia[i][j] == 1){
        for(int k = 0; k < num_Rec; k++){
          aux[i][k] = aux[i][k] + p_Recurso[j * num_Itens + k];
        }
      }
    }
    for(int k = 0; k < num_Rec; k++){
      if(aux[i][k] > lim_Recurso[k]){
        valido = false;
      }
    }
  }

}*/

/*void calcula_Fitness(int index, int num_Itens){

  fit_Colonia = (int *) malloc(num_Colonias * sizeof(int));

  int i, j;

  for(i = 0; i < num_Colonias; i++){
    for(j = 0; j < num_Itens; j++){
      if(colonia[i][j] == 1){
        fit_Colonia[i] = fit_Colonia[i] + tp_Recurso[j];
      }
      //printf("%d ", fit_Colonia[i]);  
    }
    //printf("%d ", fit_Colonia[i]);  
  }
}*/

void calcula_Fitness(int index){

  int i;

  for(i = 0; i < num_Itens; i++){
    fitness[index] = fitness[index] + colonia[index][i] * p_Recurso[i];
    //printf("%lf ", fitness[index]);
  }

}

void prep(){

  int i, j;
  float aux;
  float raiz;

  for(i = 0; i < num_Colonias; i++){
    atual_Size[i] = 1.0;
  }

  for(i = 0; i < num_Colonias; i++){
    update_cosc[i] =  (atual_Size[i] + 4*(fitness[i]))/
                      (atual_Size[i] + 2*(fitness[i]));
    aux = (3 * atual_Size[i])/(4 * M_PI);
    raiz = pow(aux, 1.0/3.0);
    fric_surf[i] = 2 * M_1_PI * (raiz);
    //printf("%lf", fric_surf[i]);
  }
}

void tournament_Select(){
  int i, poolSize;
  int *pool;

  poolSize = 2;

  //printf("teste\n");

  pool = (int*) malloc(poolSize * sizeof(int));

  for (i = 0; i < poolSize; i++){
    pool[i] = rand() % num_Colonias;
  }

  if(fitness[pool[0]]> fitness[pool[1]]){
    parent = pool[0];
  }
  else {
    parent = pool[1];
  }

  //printf("teste2\n");

  free(pool);

}

void movement(int index){
  double p;
  double alpha, beta;
  int m, k, l;
  float rand1 = (rand() / (float) RAND_MAX);
  double fit;

  // (rand() % (upper - lower + 1)) + lower;
  // float rand1 = (rand() / (float) RAND_MAX);

  p = (rand() % (1 - (-1) + 1)) + (-1);
  alpha = (rand() / (double) (2 * M_PI));
  beta = (rand() / (double) (2 * M_PI));

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
  
  discretize();

  //printf("teste\n");

  dropAdd(index);

  //printf("teste\n");

  calcula_Fitness(index);

  if(fit > fitness[index]){
    for(int i = 0; i < num_Itens; i++){
      col_Aux[i] = colonia[index][i];
      //printf("%lf ", colonia[index][i]);
    }
    fitness[index] = fit;
  }
}

void discretize(){

  double random;
  double gx;

  random = (rand() % (1 - (-1) + 1)) + (-1);

  for(int i = 0; i < num_Itens; i++){
    if(fmod(col_Aux[i], 1) != 0.0){
      gx = tanh(col_Aux[i]);
      //printf("%lf\n", gx);
       if(gx < random){
        col_Aux[i] = 0;
       }
       else {
        col_Aux[i] = 1;
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
      delta += p_Recurso[i + j * num_Itens];
    }
    psUtOrder[i].psUt = tp_Recurso[i] / delta;		
		psUtOrder[i].index = i;
  }
  //printf("%lf ", delta);
  //printf("\nAntes qsort");
  qsort(psUtOrder, (size_t) num_Itens, sizeof(psUtOrd), cmpFunc);
  //printf("\nDepois qsort");
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

  dispRc = j * num_Rec;

  minIndex = 0;

  //printf("teste\n");

  do{
    found = false;
    while(minIndex < num_Itens && !found){
      //printf("%d\n", minIndex);
      //printf("%d\n ",psUtOrder[minIndex].index);
      //printf("depois da colonia");
      if(colonia[j][psUtOrder[minIndex].index] == 1){
        found = true;
        //printf("teste\n");
      }
      else{
        minIndex++;
      }
      //printf("teste\n");
    }
    //printf("teste\n");
    colonia[j][psUtOrder[minIndex].index] = 0;
    //printf("teste\n");
    for(i = 0; i < num_Rec; i++){
      //printf("teste\n");
      rc[dispRc + i] = rc[dispRc + i] + p_Recurso[i * num_Itens + psUtOrder[minIndex].index];
      //printf("teste\n");
    }
    //printf("teste\n");
  } while(!isFeasible(j));

  for(maxIndex = num_Itens; maxIndex >= 0; maxIndex--){
    if(canAdd(dispRc , psUtOrder[maxIndex].index)){
      colonia[j][psUtOrder[maxIndex].index] = 1;
      for(i = 0; i < num_Rec; i++){
        rc[j][i] = rc[j][i] - p_Recurso[i * num_Itens + psUtOrder[maxIndex].index];
      }
    }
  }
  //printf("teste\n");
  
}

bool canAdd(int dispRc,int index){
  int i;

  for(i = 0; i < num_Rec; i++){
    if(rc[dispRc + i] - p_Recurso[index * num_Itens + i] < 0){
      //printf("alou\n");
      return false;
    }
    return true;
  }

}

bool isFeasible(int index){
  int i, displRc;

  displRc = index * num_Rec;

  for(i = 0; i < num_Rec; i++){
    if(rc[displRc + i] < 0){
      //printf("teste\n");
      return false;
    }
    //printf("index\n ");
    return true;
  }

}

void size(){

  for(int i = 0; i < num_Colonias; i++){
    atual_Size[i] = update_cosc[i] * atual_Size[i];
    //printf("%f ", atual_Size[i]);
  }

}