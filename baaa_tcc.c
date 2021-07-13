#include<stdio.h>

typedef struct knapsack {
  int ;
  float **mochila; // matriz dimenciona a mochila do problema de forma dinamica

};

FILE *arqIn;

void read(char *entrada); // funcao que realiza a leitura da entrada

int main(int argc, char *argv[]){
  int num_Itens, num_Rec, best_Sol;
  int i = 3, j = 0;
  double **matriz;

  char *entrada = argv[1];

  read(entrada);

  num_Itens = arqIn[0];
  num_Rec = arqIn[1];
  best_Sol = arqIn[2];

  while i < num_Itens * num_Rec{
    matriz[j] = arqIn[i]
    j++;
  }

  

}

void read(char *entrada){

  arqIn = fopen(entrada, "rt");

  if(arqIn == NULL){
    printf("Arquivo de entrada não encontrado")
    exit(1);
  }

  fclose(arqIn);
  return;
}
