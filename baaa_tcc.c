#include<stdio.h>
#include<stdlib.h>

typedef struct mochila {
  double atributo; // matriz dimenciona a mochila do problema de forma dinamica
  struct mochila *prox;
}knapsack;

FILE *arqIn;

void read(char *entrada); // funcao que realiza a leitura da entrada

void inserir(knapsack **cabeca, FILE *arquivo);
void listar (knapsack *noatual);

int main(int argc, char *argv[]){
  knapsack *inicio = NULL;
  knapsack *atual;

  int num_Itens, num_Rec, best_Sol, total_Val;
  int j = 0;

  char *entrada = argv[1];

  //read(entrada);

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

  printf("%d\n", num_Itens);
  printf("%d\n", num_Rec);
  printf("%d\n", best_Sol);
  printf("%d\n", total_Val);

  inserir(&inicio, arqIn);

  listar(inicio);

  fclose(arqIn);
}

/*void read(char *entrada){

  arqIn = fopen(entrada, "r");

  if(arqIn == NULL){
    printf("Arquivo de entrada não encontrado");
    exit(1);
  }

  fclose(arqIn);
  return;
}*/

void inserir(knapsack **cabeca, FILE *arquivo){
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

void listar (knapsack *noatual){
  int i=0;
  while( noatual != NULL) {    /* Enquanto nao chega no fim da lista */
    i++;
    printf("%lf ", noatual->atributo);
    noatual = noatual->prox;     /* Faz noatual apontar para o proximo no */
  }
}
