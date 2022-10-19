#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>


FILE *arqIn;

void AlocaMatriz();
void inserir_Sep(FILE *);
void ini_Colonia();
void LiberaMatriz();
void checa_Validade();
void calcula_Fitness(int);
void prep();
void size();

int *tp_Recurso, *p_Recurso, *lim_Recurso;
int **colonia;
int *fit_Colonia;
float *fric_surf;
float *atual_Size;
float *update_cosc;
int num_Itens, num_Rec, best_Sol, total_Val;
float *fitness;


int main(int argc, char *argv[]){

  int **colonia;

  int cont = 0;

  char entrada[40];

  strcpy(entrada, argv[1]);
  int num_Testes = atoi(argv[2]);
  int num_Colonias = atoi(argv[3]);

  arqIn = fopen(entrada, "r");

  if(arqIn == NULL){
    printf("Arquivo de entrada não encontrado");
    exit(1);
  }

  fscanf(arqIn, "%d %d %d", &num_Itens, &num_Rec, &best_Sol);


  total_Val = num_Itens * num_Rec + num_Itens + num_Rec;

  //printf("teste");
  inserir_Sep(arqIn, num_Itens, num_Rec);

  AlocaMatriz();
  
  ini_Colonia();

  while(cont < num_Testes){

    checa_Validade();

      prep();

    cont++;
  }

  LiberaMatriz();

  fclose(arqIn);

}

void inserir_Sep(FILE *arquivo){ // função está dando falha na segmentação

  tp_Recurso = (int *) malloc(num_Itens * sizeof(int));
  p_Recurso = (int *) malloc((num_Itens * n_Rec) * sizeof(int));
  lim_Recurso = (int *) malloc(n_Rec * sizeof(int));

  for (int j = 0; j < num_Itens; j++)
	  fscanf(arquivo, "%d", &(tp_Recurso[j]));

  //Leitura da matriz de pesos
  for (int j = 0; j < n_Rec; j++){
	  for (int k = 0; k < num_Itens; k++){
	    fscanf(arquivo, "%d", &(p_Recurso[j * num_Itens + k]));
    }
  }

  //Leitura do vetor de restrições
  for (int j = 0; j < n_Rec; j++)
	  fscanf(arquivo, "%d", &(lim_Recurso[j]));
}

void AlocaMatriz(){

  int i, j;

  colonia = (int**) malloc(num_Colonias * sizeof(int*));

  for(i = 0; i < num_Colonias; i++){
    colonia[i] = (int*) malloc(num_Itens * sizeof(int));
  }

}

void LiberaMatriz(){


  //Libera cada linha da matriz
  for (int i = 0; i < num_Colonias; i++)
	free(colonia[i]);

  //Libera a matriz
  free(colonia);

  free(tp_Recurso);
  free(p_Recurso);
  free(lim_Recurso);

}

void ini_Colonia(){
  int i, j;

  fitness = (float*) malloc(num_Colonias * sizeof(float));

  srand(time(NULL));

  for (i = 0; i < tamanho; i++){
	  for (j = 0; j < num_Itens; j++){
      colonia[i][j] = rand() % 2;
      //printf("%d ", cabeca[i][j]);
    }
    calcula_Fitness(i);
    //printf("\n");
  }
}

void checa_Validade(){
  int aux[tamanho][num_Itens];
  bool valido = true;

  for(int i = 0; i < tamanho; i++){
    for(int k = 0; k < n_Rec; k++){
      aux[i][k] = 0;
      }
    }

  for(int i = 0; i < tamanho; i++){
    valido = true;

    for(int j = 0; j < num_Itens; j++){
      if(colonia[i][j] == 1){
        for(int k = 0; k < n_Rec; k++){
          aux[i][k] = aux[i][k] + p_Recurso[j + k * num_Itens];
        }
      }
    }
    for(int k = 0; k < n_Rec; k++){
      if(aux[i][k] > lim_Recurso[k]){
        valido = false;
      }
    }
  }

}

/*void calcula_Fitness(int index, int num_Itens){

  fit_Colonia = (int *) malloc(tamanho * sizeof(int));

  int i, j;

  for(i = 0; i < tamanho; i++){
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

  int i, desl;

  desl = index + num_Itens;

  fitness[index] = 0;

  for(i = 0; i < num_Itens; i++){
    fitness[index] += colonia[desl + i] * p_Recurso[i];
  }

}

void prep(){

  atual_Size = (float *) malloc(num_Colonias * sizeof(float));
  update_cosc = (float *) malloc(num_Colonias * sizeof(float));
  fric_surf = (float *) malloc(num_Colonias * sizeof(float));
  int i;
  float aux;

  for(i = 0; i < num_Colonias; i++){
    atual_Size[i] = 1;
  }

  for(i = 0; i < num_Colonias; i++){
    update_cosc[i] =  (atual_Size[i] + 4*(fit_Colonia[i]))/
                      (atual_Size[i] + 2*(fit_Colonia[i]));
    aux = (3 * atual_Size[i])/(4 * M_1_PI);
    fric_surf[i] = 2 * M_1_PI * (pow(aux, 1.0/3.0));
    //printf("%lf", fric_surf[i]);
  }

}

void tournament_Select(){
  int i, poolSize;
  int *pool;

  pool = (int*) malloc(poolSize * sizeof(int));

  for (i = 0; i < poolSize; i++){
    pool[i] = rand() % 
  }

}

void movement(){
	int m[num_Colonias], k[num_Colonias], l[num_Colonias];
	
	for(int i = 1; i < num_Colonias; i++){
		m[i] = 
	}
	
}

void size(){

  for(int i = 0; i < num_Colonias; i++){
    atual_Size[i] = update_cosc[i] * atual_Size[i];
    //printf("%f ", atual_Size[i]);
  }

}