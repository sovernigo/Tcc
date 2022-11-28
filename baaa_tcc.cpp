
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <limits.h>

#define EXECS 100

struct OrdPsUt
{
	int index;
	float psUt;
};

struct Mochila
{
    int qtdObjs;        // Quantidade de objetos
    int numComparts;   // Numero de compartimentos da mochila
    float* valores;        // Beneficio associado a escolha de um objeto
    int* pesos;             // Onde sera armazenado os pesos referente a cada escolha (matriz em forma de vetor)
    int* capacidades;          // Capacidade limite de cada compartimento da mochila
};

/*
 * Protótipos das funções
 */
void psUtPreProcess();
void generateInitialColonia();
bool isDuplicated(int*, int, int);
bool isFeasible(int*, int);
void dropAdd(int*, int*, int);
bool canAdd(int*, int, int);
void calculaFitness(int);
void tournamentSelection();
int cmpFunc(const void*, const void*);
void movementColonia(int);
void cbInvertedMutation(int*, int*);
int findBestIndividual();
void writeBestResult(FILE*, char*, double, int);
void finalWrite(FILE*, char*, double, double, double, double, int, int);
void size();
void discretize(int);
void sizeMovement();

Mochila m;
int energia;
int tamanhoColonia;
int* dadosColonia; //vetor com tamanho da população, número de genes de cada indivíduo
float* fitness; //vetor contendo o fitness de todos os indivíduos
int* individuoColonia; //vetor contendo os genes de todos os indivíduos
int* rc;//vetor contendo os resíduos dos recursos para cada solução
OrdPsUt* psUtOrder;
int parents; //vetor contendo os índices dos pais
float shear_Force = 2;
double *fric_surf;
double *atual_Size; // vetor contendo o tamanho atual de cada colônia
double *update_cosc; // vetor contendo o coeficiente de atualização de cada colônia
int biggerIndex, smallerIndex = 0;
double *energiaColonia; //vetor contendo energia de cada colônia
float adaptProb = 0.5;
float energyLoss = 0.3;
bool isStarve;
int *starvation; // vetor contendo a fome de cada colônia

int optSol;

int* vetSols;
float* vetGaps;
float* vetTimes;
double avgSol;
double mediumGap;
double mediumTime;
double dpSols;
double dpTimes;
double dpGaps;
int bestOverallSolValue;
int bestSolValue;
double bestOverallTime;
double bestOverallGap;
int numOpt;
int numBetter;
int bestIndex;

/*
 * Função principal
 */
int main(int argc, char *argv[])
{
    if (argc != 4)
    {
		printf("Uso:\n./seq_ag arq_teste tam_pop num_ger\n");

		return 0;
    }

    //Abertura do arquivo de entrada
    FILE* fileIn = NULL;
    FILE* fileOut = NULL;
  
    timeval start, finish;
    double totalTime;
    char nomeArq[120];
  
    strcpy(nomeArq, argv[1]);    

	tamanhoColonia = atoi(argv[2]);
	energia = atoi(argv[3]);
  
    char fileOutName[160];
    sprintf(fileOutName, "%s_cb%s", nomeArq, ".out");

	fileIn = fopen(nomeArq, "r");

	//printf("Arquivo: %s", argv[1]);

	if (fileIn == NULL)
	{
	    printf("Erro na abertura do arquivo de entrada.\n");
	    return 1;
	}

	fscanf(fileIn, "%d", &m.qtdObjs);
	fscanf(fileIn, "%d", &m.numComparts);
	fscanf(fileIn, "%d", &optSol);
	
	//Alocação do vetor de valores
	m.valores = (float*) malloc(m.qtdObjs * sizeof(float));
	
	//Leitura do vetor de valores
	for (int j = 0; j < m.qtdObjs; j++)
		fscanf(fileIn, "%f", &(m.valores[j]));
	
	//Alocação da matriz de pesos como um vetor
	m.pesos = (int*) malloc(m.numComparts * m.qtdObjs * sizeof(int));
	
	//Leitura da matriz de pesos
	for (int j = 0; j < m.numComparts; j++)
		for (int k = 0; k < m.qtdObjs; k++)
			fscanf(fileIn, "%d", &(m.pesos[j * m.qtdObjs + k]));
		
	//Alocação do vetor de restrições
	m.capacidades = (int*) malloc(m.numComparts * sizeof(int));
	
	//Leitura do vetor de restrições
	for (int j = 0; j < m.numComparts; j++)
		fscanf(fileIn, "%d", &(m.capacidades[j]));

    //Criação da população
	fitness = (float*) malloc(tamanhoColonia * sizeof(float));
	individuoColonia = (int*) malloc(tamanhoColonia * m.qtdObjs * sizeof(int));
	rc = (int*) malloc(tamanhoColonia * m.numComparts * sizeof(int));
	psUtOrder = (OrdPsUt*) malloc(m.qtdObjs * sizeof(OrdPsUt));

	fric_surf = (double*) malloc(tamanhoColonia * sizeof(double));
	atual_Size = (double*) malloc(tamanhoColonia * sizeof(double));
	update_cosc = (double*) malloc(tamanhoColonia * sizeof(double));
	energiaColonia = (double*) malloc(tamanhoColonia * sizeof(double));
	starvation = (int*) malloc(tamanhoColonia * sizeof(int));
	
	/*******************************************************************
	 * Laço de execuções                                               *
	 ******************************************************************/
	vetSols = (int*) malloc(EXECS * sizeof(int));
    vetGaps = (float*) malloc(EXECS * sizeof(float));
	vetTimes = (float*) malloc(EXECS * sizeof(float));
    
    avgSol = 0.0;
	mediumTime = 0.0;
	mediumGap = 0.0;
	bestOverallSolValue = 0;
	numOpt = 0;

	srand(time(NULL));
	
	//Invoca a função de pré-processamento das pseudoutilidades uma única vez,
	//pois a entrada é a mesma.
    psUtPreProcess();
	
    for (int exec = 0; exec < EXECS; exec++)
    {
        //Inicio da tomada de tempo
		gettimeofday(&start, NULL);
		
		generateInitialColonia();

		for(int l = 0; l < tamanhoColonia; l++){
			starvation[l] = 0;
		}
		//Laço que controla o número de gerações
		for (int l = 0; l < tamanhoColonia; l++){
			while (energiaColonia[l] > 0.0){
				movementColonia(l);
				if(isStarve == true){
					starvation[l] += 1;
				}
			}
		}//Fim do for de gerações

		sizeMovement();
		
		//Encontra o melhor indivíduo
		bestIndex = findBestIndividual();
				
		//Fim da tomada de tempo
		gettimeofday(&finish, NULL);

		long start_usecs, finish_usecs, diff_usecs;

		finish_usecs = finish.tv_usec + (1000000 * finish.tv_sec);
		start_usecs = start.tv_usec + (1000000 * start.tv_sec);
		diff_usecs = finish_usecs - start_usecs;

		totalTime = (double)(diff_usecs) / 1000000;
		
		//Escreve o resultado final da execução
		bestSolValue = fitness[bestIndex];

		vetGaps[exec] = (float)(optSol - bestSolValue) / optSol * 100;
		vetTimes[exec] = totalTime;
        vetSols[exec] = bestSolValue;
		
        writeBestResult(fileOut, fileOutName, totalTime, exec);
        
        avgSol += bestSolValue;
		mediumTime += totalTime;
		mediumGap += vetGaps[exec];
		
		if (bestSolValue >= bestOverallSolValue)
		{
			bestOverallSolValue = bestSolValue;
			bestOverallTime = totalTime;
		}
		
		if (vetGaps[exec] == 0.0)
			numOpt++;
		else if (vetGaps[exec] < 0.0)
			numBetter++;

    }//Fim do laço de execuções

	//Calcular média e desvio padrão das execuções
	avgSol /= EXECS;
    mediumTime /= EXECS;
	mediumGap /= EXECS;
	
    dpSols = 0.0;
	dpTimes = 0.0;
	dpGaps = 0.0;
	
	for (int i = 0; i < EXECS; i++)
	{
        dpSols += pow((vetSols[i] - avgSol), 2);
		dpTimes += pow((vetTimes[i] - mediumTime), 2);
		dpGaps += pow((vetGaps[i] - mediumGap), 2);
	}
	
    dpSols /= EXECS;
	dpTimes /= EXECS;
	dpGaps /= EXECS;
	
    dpSols = sqrt(dpSols);
	dpTimes = sqrt(dpTimes);
	dpGaps = sqrt(dpGaps);
	
	finalWrite(fileOut, fileOutName, mediumTime, mediumGap, dpTimes, 
			   dpGaps, bestOverallSolValue, numOpt);

	//Libera a memória utilizada
	free(m.valores);
	free(m.pesos);
	free(m.capacidades);

	//free(dadosColonia);
	free(fitness);
	free(individuoColonia);
	free(rc);
	free(psUtOrder);

	free(fric_surf);
	free(atual_Size);
	free(update_cosc);
	free(energiaColonia);
	free(starvation);
	
    free(vetSols);
	free(vetGaps);
	free(vetTimes);
    
    //Fecha o arquivo de entrada
	fclose(fileIn);

    return 0;
} // fim funcao main

/*
 * Função que gera os indivíduos da população inicial. Nesta implementação
 * são utilizadas as ideias de Chu e Beasley.
 */
/*void genFeasibleInitialPopulation()
{    
    //Para cada indivíduo
    for (int j = 0; j < tamanhoColonia; j++)
    {
		fitness[j] = 0;

		int desl = j * m.qtdObjs;
		//printf("Ind: %d Desl: %d\n", init, desl);

		//Gera os genes
		for (int i = 0; i < m.qtdObjs; i++)
			individuoColonia[desl + i] = rand() % 2;

		while (!isFeasible(individuoColonia, desl))
			dropItem(individuoColonia, desl);

		calculaFitness(j);
    }

}*/

/*
 * Função para pré-processamento das pseudoutilidades para auxiliar
 * no processo de adição e remoção de itens, baseada nas ideias de
 * Chu e Beasley.
 */
void psUtPreProcess()
{
	int i, j;
	float sum;
	
	//printf("\n");
	for (i = 0; i < m.qtdObjs; i++)
	{
		//printf("%d ", i);
		sum = 0.0;
		for (j = 0; j < m.numComparts; j++)
			sum += (float) m.pesos[j * m.qtdObjs + i];//tirei a divisão por rc[j]
		
		psUtOrder[i].psUt = (float) m.valores[i] / sum;		
		psUtOrder[i].index = i;		
	}	
	
	//printf("\nAntes qsort");
	qsort(psUtOrder, (size_t) m.qtdObjs, sizeof(OrdPsUt), cmpFunc);
	//printf("\nDepois qsort");

}

/*
 * Função de comparação para qsort ordenar em ordem não decrescente
 */
int cmpFunc(const void* a, const void* b)
{
   //return ( *(int*)a - *(int*)b );
   //return (((AuxPeso*)a)->peso - ((AuxPeso*)b)->peso );
   
    if (((OrdPsUt*)a)->psUt > ((OrdPsUt*)b)->psUt)
	    return 1;

	if (((OrdPsUt*)a)->psUt < ((OrdPsUt*)b)->psUt)
		return -1;
		
	return 0;
}

/*
 * Função que gera os indivíduos da população inicial. Nesta implementação
 * são considerados apenas indivíduos que representam soluções válidas.
 */
void generateInitialColonia()
{    
	int i, j, k, desl, deslRc;
		
    //Para cada indivíduo
    for (j = 0; j < tamanhoColonia; j++)
    {
		fitness[j] = 0;

		deslRc = j * m.numComparts;
		
		desl = j * m.qtdObjs;
		//printf("Ind: %d Desl: %d\n", init, desl);

		//printf("\nAntes de gerar indivíduo %d", j);
		//printf("\nBLAH %d", j);
		do
		{
			//Inicializa o vetor residual
			for (i = 0; i < m.numComparts; i++)
				rc[deslRc + i] = m.capacidades[i];
			
			//Gera os genes
			for (i = 0; i < m.qtdObjs; i++)
			{
				individuoColonia[desl + i] = rand() % 2;
				
				//Precisa diminuir o vetor residual para refletir a adição do item
				if (individuoColonia[desl + i] == 1)
					for (k = 0; k < m.numComparts; k++)
						rc[deslRc + k] -= m.pesos[k * m.qtdObjs + i];
			}
						
			if (!isFeasible(rc, j))
				dropAdd(individuoColonia, rc, j);			
			
		} while(isDuplicated(individuoColonia, j, j - 1));
		
		calculaFitness(j);
    }
}

bool isDuplicated(int* genes, int iDesl, int last)
{
	int i, j, myDesl, curDesl;
	
	myDesl = iDesl * m.qtdObjs;
	
	//Para cada indivíduo presente na população
	for (i = 0; i < last; i++)
	{
		curDesl = i * m.qtdObjs;
		
		j = 0;
		
		//Verifica se todos os itens são iguais
		while (j < m.qtdObjs && (genes[myDesl + j] == individuoColonia[curDesl + j]))
			j++;
		
		if (j == m.qtdObjs)//Todos os itens são iguais
			return true;
			
	}
	
	return false;
}

void dropAdd(int* genes, int* rc, int j)
{
	int i, desl, displRc, minIndex, maxIndex;
	bool found;
	
	desl = j * m.qtdObjs;
	displRc = j * m.numComparts;
	
	//DROP PHASE
	minIndex = 0;
	do
	{
		//Encontra o item com a menor pseudoutilidade presente na solução
		found = false;
		while (minIndex < m.qtdObjs && !found)
		{
			if (genes[desl + psUtOrder[minIndex].index] == 1)//Se o item está na solução
				found = true;
			else	
				minIndex++;
		}
		//Retira da solução
		genes[desl + psUtOrder[minIndex].index] = 0;

		//Atualiza o vetor de capacidades residuais
		for (i = 0; i < m.numComparts; i++)
			rc[displRc + i] += m.pesos[i * m.qtdObjs + psUtOrder[minIndex].index];
	} while (!isFeasible(rc, j));
	
	//ADD PHASE	
	for (maxIndex = m.qtdObjs; maxIndex >= 0; maxIndex--)
	{
		if (canAdd(rc, displRc, psUtOrder[maxIndex].index))
		{
			genes[desl + psUtOrder[maxIndex].index] = 1;
			//Atualiza o vetor de capacidades residuais			
			for (i = 0; i < m.numComparts; i++)
				rc[displRc + i] -= m.pesos[i * m.qtdObjs + psUtOrder[maxIndex].index];
		}
	}
}

bool canAdd(int* rc, int displRc, int index)
{
	int i;	
	
	for (i = 0; i < m.numComparts; i++)
		if (rc[displRc + i] - m.pesos[i * m.qtdObjs + index] < 0)
			return false;
	
	return true;
}

/*
 * Definição da função de cálculo do fitness de um cromossomo em uma populacao.
 * O valor calculado não leva em consideração se o cromossomo representa uma
 * solução válida.
 */
void calculaFitness(int index)
{
	int i, desl;
	
    desl = index * m.qtdObjs;

    fitness[index] = 0;

    for (i = 0; i < m.qtdObjs; i++)
		fitness[index] += individuoColonia[desl + i] * m.valores[i];
}

/*
 * Função que retorna se uma solução é factível
 */
bool isFeasible(int* rc, int j)
{
    int i, displRc;
    
    displRc = j * m.numComparts;
	
    //Para cada restrição
    for (i = 0; i < m.numComparts; i++)
    	if (rc[displRc + i] < 0)
    		return false;
    	
    return true;
}

/*
 * Função para seleção utilizando o método do torneio binário.
 * No método do torneio binário dois pools com 2 elementos cada (por isso, binário) são construídos
 * com elementos escolhidos aleatoriamente da população atual e, então, o melhor elemento de cada
 * pool é selecionado. Esses dois elementos serão os pais dos indivíduos gerados para a próxima
 * geração.
 */
void tournamentSelection()
{
	int i, poolSize;
	int* pool;
	
	poolSize = 2;
	
	pool = (int*) malloc(poolSize * sizeof(int));
	
	//Monta as pools em uma só (0 e 1) e (2 e 3)
	for (i = 0; i < poolSize; i++)
		pool[i] = rand() % tamanhoColonia;

	//Escolhe o primeiro pai vindo da primeira pool	
	if (fitness[pool[0]] >= fitness[pool[1]])
		parents = pool[0];
	else
		parents = pool[1];
	
	free(pool);
}

/*
 * Função que faz o crossover dos pais da nova geração.
 * Nessa implementação, utilizamos o crossover uniforme,
 * de acordo com as ideias de Chu e Beasley.
 */
void movementColonia(int index)
{
	int i, j, desl, deslRc, displP, minIndex;
	int* child;
	int* childRc;

	double p;
  	double alpha, beta;
  	int n, k, l;

    child = (int*) malloc(m.qtdObjs * sizeof(int));
    childRc = (int*) malloc(m.numComparts * sizeof(int));
	desl = m.qtdObjs * index;
    
    do
    {
    	//Escolha dos pais
    	//printf("\nAntes da seleção");
		tournamentSelection();
		//printf("\nTeste");
		//printf("\nDepois da seleção");
		
		//Obtém os deslocamentos dos pais
	    displP = parents * m.qtdObjs;
	    
	    //Inicializa o vetor residual de capacidades para o novo indivíduo
	    for (i = 0; i < m.numComparts; i++)
	    	childRc[i] = m.capacidades[i];
    	
		for (i = 0; i < m.qtdObjs; i++)
		{
			child[i] = individuoColonia[displP + i];
			
			//Atualiza vetor residual de capacidades
			if (child[i] == 1)
				for (j = 0; j < m.numComparts; j++)
					childRc[j] -= m.pesos[j * m.qtdObjs + i];
		}

		p = ((double)rand() / (RAND_MAX / 2) - 1);
  		alpha = ((double)rand() / (RAND_MAX / (2 * M_PI)));
  		beta = ((double)rand() / (RAND_MAX / (2 * M_PI)));

		n = rand() % m.qtdObjs;
		k = rand() % m.qtdObjs;
		l = rand() % m.qtdObjs;

		while(n == k || k == l || n == l){
			if (n == k){
			k = rand() % m.qtdObjs;
			}
			if (n == l){
			l = rand() % m.qtdObjs;
			}
			if (k == l){
			l = rand() % m.qtdObjs;
			}
		}

		individuoColonia[desl + n] = individuoColonia[desl + n] + (individuoColonia[parents + n] - individuoColonia[desl + n]) * (shear_Force - fric_surf[index] * p);
		individuoColonia[desl + k] = individuoColonia[desl + k] + (individuoColonia[parents + k] - individuoColonia[desl + k]) * (shear_Force - fric_surf[index] * cos(alpha));
		individuoColonia[desl + l] = individuoColonia[desl + l] + (individuoColonia[parents + l] - individuoColonia[desl + l]) * (shear_Force - fric_surf[index] * sin(beta));

		discretize(index);
	    
	    //Verifica se é viável o corrige, caso necessário
		if (!isFeasible(childRc, 0)){
			dropAdd(child, childRc, 0);
			energiaColonia[index] -= energyLoss;
			isStarve = true;
		}

    } while(isDuplicated(child, 0, tamanhoColonia));    
    
    //Encontra o menos apto na população
    minIndex = 0;
    for (i = 1; i < tamanhoColonia; i++)
    	if (fitness[i] < fitness[minIndex])
    		minIndex = i;
    
    desl = minIndex * m.qtdObjs;
    
    //Substitui o menos apto: steady-state
    for (i = 0; i < m.qtdObjs; i++)
    	individuoColonia[desl + i] = child[i];
    
    deslRc = minIndex * m.numComparts;
    
    //Substitui o menos apto: steady-state
    for (i = 0; i < m.numComparts; i++)
    	rc[deslRc + i] = childRc[i];

	//Calcula o fitness do novo indivíduo
	calculaFitness(minIndex);

	size();
	
	free(child);
	free(childRc);
}


void discretize(int index){

  double random;
  double gx;
  double diff;

  int desl = index * m.qtdObjs;

  random = ((double)rand() / (RAND_MAX / 2)) - 1;

  for(int i = 0; i < m.qtdObjs; i++){
    if(fmod(individuoColonia[desl + i], 1) != 0.0){
      gx = tanh(individuoColonia[desl + i]);
      //printf("%lf ", random);
      //printf("%lf\n", gx);
      diff = gx - random;
       if(diff > 0){
        individuoColonia[desl + i] = 0.0;
       }
       else {
        individuoColonia[desl + i] = 1.0;
       }
    }
  }
}

void size(){

  double raiz, aux;

  for(int index = 0; index < tamanhoColonia; index++){

	double diff = atual_Size[index] - atual_Size[smallerIndex];

	update_cosc[index] =  (atual_Size[index] + 4*(fitness[index]))/
						(atual_Size[index] + 2*(fitness[index]));
	aux = (3 * atual_Size[index])/(4 * M_PI);
	raiz = pow(aux, 1.0/3.0);
	fric_surf[index] = 2 * M_1_PI * (raiz);

	atual_Size[index] = update_cosc[index] * atual_Size[index];

	if(atual_Size[index] > atual_Size[biggerIndex]){
		biggerIndex = index;
	}
	if(diff < 0.0){
		smallerIndex = index;
	}
  }
}

void sizeMovement(){
  int d;
  double rand2;
  int starveIndex = 0;

  d = rand() % m.qtdObjs;

  rand2 = ((double)rand() / (RAND_MAX / 1));

  if(adaptProb > rand2){
    for(int i = 0; i < tamanhoColonia; i++){
      if(starvation[i] < starvation[starveIndex]){
        starveIndex = i;
      }
    }
    individuoColonia[starveIndex * m.qtdObjs + d] = (individuoColonia[biggerIndex * m.qtdObjs + d] - individuoColonia[starveIndex * m.qtdObjs + d]) * rand2;

    discretize(starveIndex);

    if(!isFeasible(individuoColonia, starveIndex)){
      dropAdd(individuoColonia, rc, starveIndex);
    }
    calculaFitness(starveIndex);

    size();
  }
}

int findBestIndividual()
{
	int i, iMax;
	
	iMax = 0;
	
	for (i = 1; i < tamanhoColonia; i++)
		if (fitness[i] > fitness[iMax])
			iMax = i;
	
	return iMax;
}

/*
 * Função que escreve o melhor resultado de todas as iterações de uma
 * execução.
 */
void writeBestResult(FILE* f, char* nome, double time, int exec)
{
	FILE* fGaps;
	FILE* fSols;
	FILE* fTimes;

    f = fopen(nome, "a+");
	fGaps = fopen("seq_allgaps_nopr.out", "a+");
	fSols = fopen("seq_allsols.out", "a+");
	fTimes = fopen("seq_alltimes_nopr.out", "a+");

	float bestGap = (float)(optSol - bestSolValue) / optSol * 100;

    fprintf(f, "Execução: %d\nTempo Total: %.4f\n Valor da Melhor Solução: %d", 
            exec, time, bestSolValue);
    fprintf(f, "\nNúmero de gerações = %d\nSolução ótima = %d\nMelhor solução = %d\nMelhor índice = %d\nGap: %f%%\n",
		energia, optSol, bestSolValue, bestIndex, bestGap);

	fprintf(fGaps,"%.4f\n", bestGap);
	fprintf(fSols,"%d\n", bestSolValue);
	fprintf(fTimes,"%.4f\n", time);

    fprintf(f, "Solução:\n");

    int displ = bestIndex * m.qtdObjs;
    for (int i = 0; i < m.qtdObjs; i++)
		fprintf(f, "%d ", individuoColonia[displ + i]);

    fprintf(f, "\n\n");

    fclose(f);
    fclose(fGaps);
    fclose(fSols);
	fclose(fTimes);
}

/*
 * Função para escrever no arquivo de saída
 */
void finalWrite(FILE* f, char* nome, double mediumTime, double mediumGap, 
                double dpTimes, double dpGaps, 
                int bestOverallSolValue, int numOpt)
{
    FILE* f2;
	double bestGap = (float)(optSol - bestOverallSolValue) / optSol * 100;
	
    f = fopen(nome, "a+");
    f2 = fopen("seqdados_ag.out", "a+");

    fprintf(f, "Execuções: %d \nTempo Medio: %.4f \nSolução Média: %.4lf \nGap Médio: %.4lf \nDesvio-Padrão (Soluções): %.4lf\nDesvio-Padrão (Gaps): %.4lf\nDesvio-Padrão (Tempos): %.4lf\nMelhor Solução: %d, Melhor Gap: %.4lf, Melhor Tempo: %.4lf\nÓtimas: %d\nMelhores: %d\n",
			EXECS, mediumTime, avgSol, mediumGap, dpSols, dpGaps, dpTimes, 
			bestOverallSolValue, bestGap, bestOverallTime, numOpt, numBetter);

    fprintf(f2, "%d %d %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n", numOpt, 
            bestOverallSolValue, avgSol, dpSols, bestGap, mediumGap, 
            dpGaps, mediumTime, dpTimes);

    fclose(f);
    fclose(f2);
}
