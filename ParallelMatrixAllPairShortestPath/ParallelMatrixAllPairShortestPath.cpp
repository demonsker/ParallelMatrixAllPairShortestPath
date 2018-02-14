#include "stdafx.h"
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#define INF 999999
#define SIZE 2048

void generate(int **);
void initialize(int **, int **);
void findAllPairShortestPath(int **, int **, int **, int);
int getSizePerProcess(int);
int getBeginIndexFromInput();
void printProcess(int **, int n);
void print(int **);
void useExampleData(int **);

int world_size, world_rank;

int main(int argc, char** argv) {

	double start1, end1, start2, end2;

	// Initialize the MPI environment
	MPI_Init(NULL, NULL);

	//Before setup
	start1 = MPI_Wtime();

	// Get the number of processes
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	// Get the rank of the process
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	//Number of row per processor
	int n = getSizePerProcess(world_rank);

	//create memory for copy part of full distance
	int i, j;
	int **partOfDistance;
	partOfDistance = (int**)malloc(n * sizeof(int*));
	for (i = 0; i < n; i++)
	{
		partOfDistance[i] = (int*)malloc(SIZE * sizeof(int));
	}

	//ceate memory for receive path
	int **partOfPath;
	partOfPath = (int**)malloc(n * sizeof(int*));
	for (i = 0; i < n; i++)
	{
		partOfPath[i] = (int*)malloc(SIZE * sizeof(int));
		for (j = 0; j < SIZE; j++)
			partOfPath[i][j] = j;
	}

	//ceate memory for receive full distance
	int **distance;
	distance = (int**)malloc(SIZE * sizeof(int*));
	for (i = 0; i < SIZE; i++)
	{
		distance[i] = (int*)malloc(SIZE * sizeof(int));
	}

	//Master
	if (world_rank == 0)
	{
		//After setup
		end1 = MPI_Wtime();
		
		//declare data for generate
		int **dataGen;
		dataGen = (int**)malloc(SIZE * sizeof(int*));
		int i, j;
		for (i = 0; i < SIZE; i++)
		{
			dataGen[i] = (int*)malloc(SIZE * sizeof(int));
		}

		//Generate data
		generate(dataGen);
		//useExampleData(dataGen);

		//Start calculate
		start2 = MPI_Wtime();

		initialize(dataGen, distance);

	}

	//BroadCast all data to another processor
	for (i = 0; i < SIZE; i++)
			MPI_Bcast(distance[i], SIZE, MPI_INT, 0, MPI_COMM_WORLD);

	//partition data
	int beginIndex = getBeginIndexFromInput();
	int endIndex = beginIndex + n;
	int currentIndex = 0;
	for (i = beginIndex; i < endIndex; i++)
	{
		for (j = 0; j < SIZE; j++)
			partOfDistance[currentIndex][j] = distance[i][j];
		currentIndex++;
	}

	//Find shrotest path
	findAllPairShortestPath(partOfDistance, distance, partOfPath, n);

	//End calculate
	end2 = MPI_Wtime();

	//printf("Process %d : Distance\n", world_rank);
	//printProcess(partOfDistance, n);

	//printf("Process %d : Path\n", world_rank);
	//printProcess(partOfPath, n);

	float diff = (float)(end1 - start1 + end2 - start2);
	printf("Time : %.4f\n", diff);

	if (world_rank > 0)
		MPI_Send(&diff, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
	else
	{
		float temp;
		for (i = 1; i < world_size; i++)
		{
			MPI_Recv(&temp, 1, MPI_FLOAT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (temp > diff)
				diff = temp;
		}
		printf("Finish Time : %.4f\n", diff);

		FILE * fp;
		char fileName[10];
		char filePath[70] = "C:\\Users\\EucliwoodX\\Desktop\\STATS\\Matrix\\stat(SaveMode)\\Parallel\\";

		sprintf(fileName, "%d.txt", SIZE);
		strcat(filePath, fileName);
		fp = fopen(filePath, "a");
		fprintf(fp, "%d Process\n %.4f\n\n", world_size, diff);
		fclose(fp);
	}

	MPI_Finalize();
}

void findAllPairShortestPath(int **partOfDistance, int **distance, int **path, int n)
{
	int **tempGraph;
	int i, j, k, r;
	tempGraph = (int**)malloc(SIZE * sizeof(int*));
	for (i = 0; i < SIZE; i++)
	{
		tempGraph[i] = (int*)malloc(SIZE * sizeof(int));
	}
	
	initialize(distance, tempGraph);
	
	int round = (int)(log10(SIZE) / log10(2));
	for (r = 0; r < round; r++)
	{		
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < SIZE; j++)
			{
				for (k = 0; k < SIZE; k++)
				{
					if (partOfDistance[i][k] + tempGraph[k][j] < partOfDistance[i][j])
					{
						partOfDistance[i][j] = partOfDistance[i][k] + tempGraph[k][j];
						path[i][j] = path[i][k];
					}
				}
			}
		}

		//Integrate partOfDistance to tempGraph
		int row,p;
		if (world_rank != 0)
		{
			for (row = 0; row < n; row++)
				MPI_Send(partOfDistance[row], SIZE, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
		else
		{
			int beginIndex = n;
			for (p = 1; p < world_size; p++)
				for (row = 0; row < getSizePerProcess(p); row++)
					MPI_Recv(tempGraph[beginIndex++], SIZE, MPI_INT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		
		//BroadCast all data to another processor
		for (row = 0; row < SIZE; row++)
			MPI_Bcast(tempGraph[row], SIZE, MPI_INT, 0, MPI_COMM_WORLD);
	}
		
}

void generate(int **data)
{
	int i, j, r;

	for (i = 0; i < SIZE; i++)
	{
		data[i][i] = 0;
		for (j = i + 1; j < SIZE; j++)
		{
			r = (rand() % 20) + 1;
			if (r == 19)
				data[i][j] = INF;
			else
				data[i][j] = r;
			data[j][i] = data[i][j];
		}
	}
}

void initialize(int **sour, int **dest)
{
	int i, j;

	for (i = 0; i < SIZE; i++)
	{
		for (j = 0; j < SIZE; j++)
		{
			dest[i][j] = sour[i][j];
		}
	}
}

void printProcess(int **distance, int n)
{
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < SIZE; ++j)
		{
			if (distance[i][j] == INF)
				printf("%7s", "INF");
			else
				printf("%7d", distance[i][j]);
		}

		printf("\n");
	}
}

void print(int **m)
{
	printf("Shortest distances between every pair of vertices: \n");

	for (int i = 0; i < SIZE; ++i)
	{
		for (int j = 0; j < SIZE; ++j)
		{
			if (m[i][j] == INF)
				printf("%7s", "INF");
			else
				printf("%7d", m[i][j]);
		}

		printf("\n");
	}
}

int getSizePerProcess(int rank)
{
	int n = SIZE / world_size;
	int m = SIZE % world_size;
	if (m != 0)
		if (rank < m)
			n++;
	return n;
}
int getBeginIndexFromInput()
{
	if (world_rank == 0) return 0;
	int begin, end, np, i;
	end = getSizePerProcess(0);
	for (i = 1; i < world_size; i++)
	{
		np = getSizePerProcess(i);
		begin = end;
		end = begin + np;
		if (world_rank == i)
			return begin;
	}
	return -1;
}

void useExampleData(int **data)
{
	int example[8][8] = {
		{ 0,1,9,3,INF,INF,INF,INF },
		{ 1,0,INF,1,INF,3,INF,INF },
		{ 9,INF,0,INF,INF,3,10,INF },
		{ 3,1,INF,0,5,INF,INF,8 },
		{ INF,INF,INF,5,0,2,2,1 },
		{ INF,3,3,INF,2,0,INF,INF },
		{ INF,INF,10,INF,2,INF,0,4 },
		{ INF,INF,INF,8,1,INF,4,0 }
	};

	int i, j;

	for (i = 0; i < SIZE; i++)
	{
		for (j = 0; j < SIZE; j++)
		{
			data[i][j] = example[i][j];
		}
	}
}