#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>
#include<time.h>

/*typedef struct mochila {
  double atributo; // matriz dimenciona a mochila do problema de forma dinamica
  struct mochila *prox;
}knapsack;*/

FILE *arqIn;

// void read(char *entrada); // funcao que realiza a leitura da entrada

void AlocaMatriz(int num_Itens, int n_Colonias);
void inserir_Sep(FILE *arquivo, int n_Itens, int n_Rec);
void ini_Colonia(int **cabeca, int tamanho, int n_Itens);
void LiberaMatriz(int **colonia, int num_Itens, int num_Rec);


int *tp_Recurso, *p_Recurso, *lim_Recurso;
int **colonia;

int main(int argc, char *argv[]){

  int num_Itens, num_Rec, best_Sol, total_Val;
  int cont = 0;

  char *entrada = argv[1];
  int n_Testes = atoi(argv[2]) - 48;
  int n_Colonias = atoi(argv[3]) - 48;

  printf("%d %d\n", n_Testes, n_Colonias);

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
  printf("%d", total_Val);

  inserir_Sep(arqIn, num_Itens, num_Rec);


  AlocaMatriz(num_Itens, n_Colonias);

  while(cont < n_Testes){

    ini_Colonia(colonia, total_Val, num_Itens);

    cont++;
  }

  LiberaMatriz(colonia, num_Itens, num_Rec);

  fclose(arqIn);

}

/*void inserir(knapsack **cabeca, FILE *arquivo, int n_Itens, int n_Rec){
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
}*/

void inserir_Sep(FILE *arquivo, int n_Itens, int n_Rec){
  int valor;

  int cont = 0;
  int i, j, k = 0;

  tp_Recurso = (int *) malloc(n_Itens * sizeof(int));
  p_Recurso = (int *) malloc((n_Itens * n_Rec) * sizeof(int));
  lim_Recurso = (int *) malloc(n_Rec * sizeof(int));

  while(fscanf(arqIn,"%d", &valor) != EOF){

    if(cont < n_Rec){
      continue;
    }
    else if(cont >= n_Rec && cont < n_Itens){
      tp_Recurso[i] = valor;
      i++;
    }
    else if(cont >= n_Itens && cont < (n_Itens * (n_Rec + 1))) {
      p_Recurso[j] = valor;
      j++;
    }
    else if(cont >= (n_Itens * (n_Rec + 1))) {
      lim_Recurso[k] = valor;
      k++;
    }
    cont++;
  }
}

void AlocaMatriz(int num_Itens, int n_Colonias){

  int i, j;

  colonia = malloc(num_Itens * sizeof(int*));

  for(i = 0; i < num_Itens; i++){
    colonia[i] =  malloc(n_Colonias * sizeof(int));
  }

  for(i = 0; i < num_Itens; i++){
    for(j = 0; j < n_Colonias; j++){
      colonia[i][j] = 0;
    }
  }
}

void LiberaMatriz(int **colonia, int num_Itens, int num_Rec){

  free(colonia);
  free(tp_Recurso);
  free(p_Recurso);
  free(lim_Recurso);

}

void ini_Colonia(int **cabeca, int tamanho, int n_Itens){
  int i;

  while(i < tamanho){
    cabeca[i] = rand() % 2;
    i++;
  }
  return;
}

/*checa_Validade(){

}*/
