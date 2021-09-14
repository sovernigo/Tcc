#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>
#include<time.h>

typedef struct mochila {
  double atributo; // matriz dimenciona a mochila do problema de forma dinamica
  struct mochila *prox;
}knapsack;

FILE *arqIn;

void read(char *entrada); // funcao que realiza a leitura da entrada

void inserir(knapsack **cabeca, FILE *arquivo, int n_Itens, int n_Rec);
void listar (int *vetor, int tamanho);
int **AlocaMatriz(int num_Itens, int n_Colonias);
void ini_Colonia(int *cabeca, int tamanho, int n_Itens);


double *tp_Recurso, *p_Recurso, *lim_Recurso;

int main(int argc, char *argv[]){

  int num_Itens, num_Rec, best_Sol, total_Val;
  int j = 0, cont = 0;

  char *entrada = argv[1];
  int n_Testes = argv[2];
  int n_Colonias = argv[3];

  int *colonia;

  arqIn = fopen(entrada, "r");

  if(arqIn == NULL){
    printf("Arquivo de entrada não encontrado");
    exit(1);
  }

  fscanf(arqIn, "%d %d %d", &num_Itens, &num_Rec, &best_Sol);

  if(arqIn == NULL){
    return 0;
  }

  total_Val = num_Itens * num_Rec + num_Itens + num_Rec;
  colonia = AlocaMatriz(num_Itens, n_Colonias);

  //inserir(&inicio, arqIn, num_Itens, num_Rec);

  inserir_Sep(arqIn, num_Itens, num_Rec);

  while(cont < n_Testes){

    ini_Colonia(colonia, total_Val, num_Itens);

    listar(colonia, total_Val);

    cont++;
  }
  fclose(arqIn);

  LiberaMatriz(colonia, num_Itens);

}

void inserir(knapsack **cabeca, FILE *arquivo, int n_Itens, int n_Rec){
  knapsack *noatual, *novo;
  double valor;

  while(fscanf(arqIn,"%lf", &valor) != EOF){

    if(*cabeca == NULL){
      *cabeca = (knapsack *) malloc(sizeof(knapsack));
      (*cabeca)->atributo = valor;
      (*cabeca)->prox = NULL;
    }
    else {
      noatual = *cabeca;
      while(noatual->prox != NULL){
          noatual = noatual->prox;
      }
      novo = (knapsack *)malloc(sizeof(knapsack));
      novo->atributo = valor;
      novo->prox = NULL;
      noatual->prox = novo;
    }
  }
}

void inserir_Sep(FILE *arquivo, int n_Itens, int n_Rec){
  double valor;

  int cont = 0;
  int i, j = 0;

  tp_Recurso = (double *) malloc(n_Itens * sizeof(double*));
  p_Recurso = (double *) malloc((n_Itens * n_Rec) * sizeof(double*));
  lim_Recurso = (double *) malloc(n_Rec * sizeof(double*));

  while(fscanf(arqIn,"%lf", &valor) != EOF){

    if(cont < n_Itens){
      tp_Recurso[cont] = valor;
    }
    else if(n_Itens <= cont < (n_Itens + 1)* n_Rec) {
      p_Recurso[i] = valor;
      i++;
    }
    else {
      lim_Recurso[j] = valor;
      j++;
    }
    cont++;
  }
}

int **AlocaMatriz(int num_Itens, int n_Colonias){
  int **colonia;
  int i;

  colonia = (int **) malloc(num_Itens * sizeof(int*));
  if(colonia == NULL){
    printf("Memoria insuficiente.\n");
    exit(1);
  }
  for(i = 0; i < num_Itens; i++){
    colonia[i] = (int *)malloc(n_Colonias * sizeof(int));
    if(colonia[i] == NULL){
      printf("Memoria insuficiente.\n");
      exit(1);
    }
  }
  return colonia;
}

void LiberaMatriz(int **colonia, int num_Itens){
  int i;
  for(i = 0; i < num_Itens; i++)
    free(colonia[i]);
  free(colonia);
}

void ini_Colonia(int *cabeca, int tamanho, int n_Itens){
  int i, random;

  while(i < tamanho){
    random = rand();

    if(random % 2 == 0){
      cabeca[i] = 0;
    }
    else {
      cabeca[i] = 1;
    }
    i++;
  }
  return;
}

void listar (int *vetor, int tamanho){
  int i = 0;
  while(i < tamanho) {    /* Enquanto nao chega no fim da lista */

    printf("%d ", vetor[i]);    /* Faz noatual apontar para o proximo no */
    i++;
  }
}
