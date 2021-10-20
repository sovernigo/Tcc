#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>


FILE *arqIn;

void AlocaMatriz(int, int);
void inserir_Sep(FILE *, int, int);
void ini_Colonia(int **, int, int);
void LiberaMatriz(int **, int);
void checa_Validade(int **, int, int);


int *tp_Recurso, *p_Recurso, *lim_Recurso;
int **colonia;

int main(int argc, char *argv[]){

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


  inserir_Sep(arqIn, num_Itens, num_Rec);

  AlocaMatriz(num_Itens, n_Colonias);

  while(cont < n_Testes){

    ini_Colonia(colonia, n_Colonias, num_Itens);

    cont++;
  }

  LiberaMatriz(colonia, n_Colonias);

  fclose(arqIn);

}

void inserir_Sep(FILE *arquivo, int n_Itens, int n_Rec){ // função está dando falha na segmentação

  tp_Recurso = (int *) malloc(n_Itens * sizeof(int));
  p_Recurso = (int *) malloc((n_Itens * n_Rec) * sizeof(int));
  lim_Recurso = (int *) malloc(n_Rec * sizeof(int));

  for (int j = 0; j < n_Itens; j++)
	fscanf(arquivo, "%d", &(tp_Recurso[j]));

  //Leitura da matriz de pesos
  for (int j = 0; j < n_Rec; j++)
	 for (int k = 0; k < n_Itens; k++)
	  fscanf(arquivo, "%d", &(p_Recurso[j * n_Itens + k]));

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

void LiberaMatriz(int **colonia, int n_Colonias){


  //Libera cada linha da matriz
  for (int i = 0; i < n_Colonias; i++)
	free(colonia[i]);

  //Libera a matriz
  free(colonia);

  free(tp_Recurso);
  free(p_Recurso);
  free(lim_Recurso);

}

void ini_Colonia(int **cabeca, int tamanho, int n_Itens){
  int i, j;

  for (i = 0; i < tamanho; i++){
	  for (j = 0; j < n_Itens; j++){
      cabeca[i][j] = rand() % 2;
      printf("%d ", cabeca[i][j]);
    }
    printf("\n");
  }
}

void checa_Validade(int **cabeca, int tamanho, int n_Itens){

}