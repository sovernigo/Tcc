#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>

typedef struct mochila {
  double atributo; // matriz dimenciona a mochila do problema de forma dinamica
  struct mochila *prox;
}knapsack;

FILE *arqIn;

void read(char *entrada); // funcao que realiza a leitura da entrada

void inserir(knapsack **cabeca, FILE *arquivo, int n_Itens, int n_Rec);
void listar (knapsack *noatual);
int **AlocaMatriz(int num_Itens, int n_Colonias);
void ini_Colonia(bool **cabeca, knapsack **inicio);

int main(int argc, char *argv[]){
  knapsack *inicio = NULL;
  knapsack *atual;

  int num_Itens, num_Rec, best_Sol, total_Val;
  int j = 0, cont = 0;

  char *entrada = argv[1];
  int n_Testes = argv[2];
  int n_Colonias = argv[3];

  bool **colonia;

  arqIn = fopen(entrada, "r");

  if(arqIn == NULL){
    printf("Arquivo de entrada não encontrado");
    exit(1);
  }

  fscanf(arqIn, "%d %d %d", &num_Itens, &num_Rec, &best_Sol);

  if(arqIn == NULL){
    return 0;
  }

  total_Val = num_Itens * num_Rec;
  colonia = AlocaMatriz(num_Itens, n_Colonias);

  inserir(&inicio, arqIn, num_Itens, num_Rec);

  listar(inicio);

  while(cont < n_Testes){

    ini_Colonia(colonia, &inicio);

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

int **AlocaMatriz(int num_Itens, int n_Colonias){
  int **colonia;
  int i;

  colonia = (bool **) malloc(num_Itens * sizeof(bool*));
  if(colonia == NULL){
    printf("Memoria insuficiente.\n");
    exit(1);
  }
  for(i = 0; i < num_Itens; i++){
    colonia[i] = (bool *)malloc(n_Colonias * sizeof(bool));
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

void ini_Colonia(bool **cabeca, knapsack **inicio){
  return;
}

void listar (knapsack *noatual){
  int i = 0;
  while(noatual != NULL) {    /* Enquanto nao chega no fim da lista */
    i++;
    printf("%lf ", noatual->atributo);
    noatual = noatual->prox;     /* Faz noatual apontar para o proximo no */
  }
}
