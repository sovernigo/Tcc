#include<stdio.h>

typedef struct knapsack {
  int ;
  float **mochila; // matriz dimenciona a mochila do problema de forma dinamica

};

void read(); // funcao que realiza a leitura da entrada

int main(){

}

void read(){
  FILE *arqIn;

  arqIn = fopen("argv[1]", "rt");

  if(arqIn == NULL){
    printf("Arquivo de entrada não encontrado")
    exit(1);
  }

  fclose(arqIn);
}
