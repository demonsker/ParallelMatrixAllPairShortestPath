#include "stdafx.h"
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#define PATH "C:\\Users\\Eucliwood\\Desktop\\stat(SaveMode)\\Parallel\\"
#define INF 999999
#define SIZE 2048

void array_transpose(int[][SIZE]);
void distance_generate(int[][SIZE]);
void distance_useexample(int[][SIZE]);
void find_AllPairShortestPath(int[][SIZE],int[][SIZE], int[][SIZE], int);
int get_datasize_per_process(int);
int get_beginindex_frominput(int);
void process_print(int[][SIZE], int);
void array_print(int[][SIZE]);
void log_save(float);

int world_size, world_rank;

int main(int argc, char** argv) {

	double start_1, end_1, start_2, end_2;

	// Initialize the MPI environment
	MPI_Init(NULL, NULL);

	//Before setup
	start_1 = MPI_Wtime();

	// Get the number of processes
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	// Get the rank of the process
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	//Number of row per processor
	int row_per_process = get_datasize_per_process(world_rank);

	//create memory for receive distance , path
	int i, j;
	int(*part_of_distance)[SIZE], (*part_of_path)[SIZE], (*distance)[SIZE];
	part_of_distance = (int(*)[SIZE]) malloc(row_per_process * sizeof(int[SIZE]));
	part_of_path = (int(*)[SIZE]) malloc(row_per_process * sizeof(int[SIZE]));
	distance = (int(*)[SIZE]) malloc(SIZE * sizeof(int[SIZE]));
	for (i = 0; i < row_per_process; i++)
		for (j = 0; j < SIZE; j++)
		{
			part_of_path[i][j] = j;
		}

	//Master
	if (world_rank == 0)
	{		
		//After setup
		end_1 = MPI_Wtime();

		//Generate data
		distance_generate(distance);
		//distance_useexample(distance);

		//Start calculate
		start_2 = MPI_Wtime();
	}

	//BroadCast all data to another processor
	MPI_Bcast(distance, SIZE*SIZE, MPI_INT, 0, MPI_COMM_WORLD);

	//partition data
	int begin_index = get_beginindex_frominput(world_rank);
	int end_index = begin_index + row_per_process;
	int current_index = 0;
	for (i = begin_index; i < end_index; i++)
	{
		for (j = 0; j < SIZE; j++)
			part_of_distance[current_index][j] = distance[i][j];
		current_index++;
	}

	//Find shrotest path
	find_AllPairShortestPath(part_of_distance, distance, part_of_path, row_per_process);

	//End calculate
	end_2 = MPI_Wtime();

    //printf("Process %d : Distance\n", world_rank);
	//process_print(part_of_distance, row_per_process);

	//printf("Process %d : Path\n", world_rank);
	//process_print(part_of_path, row_per_process);

	float diff = (float)(end_1 - start_1 + end_2 - start_2);
	printf("Time : %.4f\n", diff);

	log_save(diff);

	MPI_Finalize();
}

void find_AllPairShortestPath(int part_of_distance[][SIZE], int distance[][SIZE], int part_of_path[][SIZE], int row_per_process)
{	
	int i, j, k, r;
	int round = (int)(log10(SIZE) / log10(2));
	for (r = 0; r < round; r++)
	{		
		for (i = 0; i < row_per_process; i++)
		{
			for (j = 0; j < SIZE; j++)
			{
				for (k = 0; k < SIZE; k++)
				{
					if (part_of_distance[i][k] + distance[j][k] < part_of_distance[i][j])
					{
						part_of_distance[i][j] = part_of_distance[i][k] + distance[j][k];
						part_of_path[i][j] = part_of_path[i][k];
					}
				}
			}
		}

		//Integrate partOfDistance to tempGraph
		if (world_rank != 0)
		{
			MPI_Send(part_of_distance, row_per_process*SIZE, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
		else
		{
			int p;
			for (p = 1; p < world_size; p++)
			{
				int row_begin = get_beginindex_frominput(p);
				int number_of_row = get_datasize_per_process(p);
				MPI_Recv(distance[row_begin], number_of_row*SIZE, MPI_INT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				array_transpose(distance);
			}
		}
		
		//BroadCast all data to another processor
		MPI_Bcast(distance, SIZE*SIZE, MPI_INT, 0, MPI_COMM_WORLD);
	}
		
}

void distance_generate(int data[][SIZE])
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

void array_transpose(int data[][SIZE])
{
	int i, j;

	for (i = 0; i < SIZE; i++)
	{
		for (j = 0; j < SIZE; j++)
		{
			data[i][j] ^= data[j][i];
			data[j][i] ^= data[i][j];
			data[i][j] ^= data[j][i];
		}
	}
}

void process_print(int distance[][SIZE], int n)
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

void array_print(int m[][SIZE])
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

int get_datasize_per_process(int rank)
{
	int n = SIZE / world_size;
	int m = SIZE % world_size;
	if (m != 0)
		if (rank < m)
			n++;
	return n;
}
int get_beginindex_frominput(int world_rank)
{
	if (world_rank == 0) return 0;
	int begin, end, np, i;
	end = get_datasize_per_process(0);
	for (i = 1; i < world_size; i++)
	{
		np = get_datasize_per_process(i);
		begin = end;
		end = begin + np;
		if (world_rank == i)
			return begin;
	}
	return -1;
}

void distance_useexample(int data[][SIZE])
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

void log_save(float diff)
{
	int i;
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
		char filePath[70] = PATH;

		sprintf(fileName, "%d.txt", SIZE);
		strcat(filePath, fileName);
		fp = fopen(filePath, "a");
		fprintf(fp, "%d Process\n %.4f\n\n", world_size, diff);
		fclose(fp);
	}
}
