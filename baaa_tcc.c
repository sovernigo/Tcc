#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>


FILE *arqIn;

void AlocaMatriz(int, int);
void inserir_Sep(FILE *, int, int);
void ini_Colonia(int, int);
void LiberaMatriz(int);
void checa_Validade(int, int, int);
void fit(int, int);
void prep(int);
void size(int);

int *tp_Recurso, *p_Recurso, *lim_Recurso;
int **colonia;
int *fit_Colonia;
float *fric_surf;
float *atual_Size;
float *update_cosc;


int main(int argc, char *argv[]){

  int **colonia;

  int num_Itens, num_Rec, best_Sol, total_Val;
  int cont = 0;

  char entrada[40];

  strcpy(entrada, argv[1]);
  int n_Testes = atoi(argv[2]);
  int n_Colonias = atoi(argv[3]);

  arqIn = fopen(entrada, "r");

  if(arqIn == NULL){
    printf("Arquivo de entrada não encontrado");
    exit(1);
  }

  fscanf(arqIn, "%d %d %d", &num_Itens, &num_Rec, &best_Sol);


  total_Val = num_Itens * num_Rec + num_Itens + num_Rec;

  //printf("teste");
  inserir_Sep(arqIn, num_Itens, num_Rec);

  AlocaMatriz(num_Itens, n_Colonias);
  
  ini_Colonia(n_Colonias, num_Itens);

  while(cont < n_Testes){

    checa_Validade(n_Colonias, num_Itens, num_Rec);

      fit(n_Colonias, num_Itens);

      prep(n_Colonias);

    cont++;
  }

  LiberaMatriz(n_Colonias);

  fclose(arqIn);

}

void inserir_Sep(FILE *arquivo, int n_Itens, int n_Rec){ // função está dando falha na segmentação

  tp_Recurso = (int *) malloc(n_Itens * sizeof(int));
  p_Recurso = (int *) malloc((n_Itens * n_Rec) * sizeof(int));
  lim_Recurso = (int *) malloc(n_Rec * sizeof(int));

  for (int j = 0; j < n_Itens; j++)
	  fscanf(arquivo, "%d", &(tp_Recurso[j]));

  //Leitura da matriz de pesos
  for (int j = 0; j < n_Rec; j++){
	  for (int k = 0; k < n_Itens; k++){
	    fscanf(arquivo, "%d", &(p_Recurso[j * n_Itens + k]));
    }
  }

  //Leitura do vetor de restrições
  for (int j = 0; j < n_Rec; j++)
	  fscanf(arquivo, "%d", &(lim_Recurso[j]));
}

void AlocaMatriz(int num_Itens, int n_Colonias){

  int i, j;

  colonia = (int**) malloc(n_Colonias * sizeof(int*));

  for(i = 0; i < n_Colonias; i++){
    colonia[i] = (int*) malloc(num_Itens * sizeof(int));
  }

}

void LiberaMatriz(int n_Colonias){


  //Libera cada linha da matriz
  for (int i = 0; i < n_Colonias; i++)
	free(colonia[i]);

  //Libera a matriz
  free(colonia);

  free(tp_Recurso);
  free(p_Recurso);
  free(lim_Recurso);

}

void ini_Colonia(int tamanho, int n_Itens){
  int i, j;

  srand(time(NULL));

  for (i = 0; i < tamanho; i++){
	  for (j = 0; j < n_Itens; j++){
      colonia[i][j] = rand() % 2;
      //printf("%d ", cabeca[i][j]);
    }
    //printf("\n");
  }
}

void checa_Validade(int tamanho, int n_Itens, int n_Rec){
  int aux[tamanho][n_Itens];
  bool valido = true;

  for(int i = 0; i < tamanho; i++){
    for(int k = 0; k < n_Rec; k++){
      aux[i][k] = 0;
      }
    }

  for(int i = 0; i < tamanho; i++){
    valido = true;

    for(int j = 0; j < n_Itens; j++){
      if(colonia[i][j] == 1){
        for(int k = 0; k < n_Rec; k++){
          aux[i][k] = aux[i][k] + p_Recurso[j + k * n_Itens];
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

void fit(int tamanho, int n_Itens){

  fit_Colonia = (int *) malloc(tamanho * sizeof(int));

  int i, j;

  for(i = 0; i < tamanho; i++){
    for(j = 0; j < n_Itens; j++){
      if(colonia[i][j] == 1){
        fit_Colonia[i] = fit_Colonia[i] + tp_Recurso[j];
      }
      //printf("%d ", fit_Colonia[i]);  
    }
    //printf("%d ", fit_Colonia[i]);  
  }
}

void prep(int n_Colonias){

  atual_Size = (float *) malloc(n_Colonias * sizeof(float));
  update_cosc = (float *) malloc(n_Colonias * sizeof(float));
  fric_surf = (float *) malloc(n_Colonias * sizeof(float));
  int i;
  float aux;

  for(i = 0; i < n_Colonias; i++){
    atual_Size[i] = 1;
  }

  for(i = 0; i < n_Colonias; i++){
    update_cosc[i] =  (atual_Size[i] + 4*(fit_Colonia[i]))/
                      (atual_Size[i] + 2*(fit_Colonia[i]));
    aux = (3 * atual_Size[i])/(4 * M_1_PI);
    fric_surf[i] = 2 * M_1_PI * (pow(aux, 1.0/3.0));
    //printf("%lf", fric_surf[i]);
  }

}

void size(int n_Colonias){

  for(int i = 0; i < n_Colonias; i++){
    atual_Size[i] = update_cosc[i] * atual_Size[i];
    //printf("%f ", atual_Size[i]);
  }

}