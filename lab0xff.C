#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <string.h>
#include <time.h>

//defining values for case loop and pi value

#define PI 3.141
#define CONTINUOUS  true
#define DISPLAYTABLE true
#define EUCLIDEANDISTANCE true

typedef struct Node		 //coordinate for EuclideanDistance method
{
	double x;
	double y;
} Node;

enum Algs				//enum declaring all algos used in code
{
	ALLALGOS,
	GREEDY,
	BRUTEFORCE,
	GREEDYCOMPARE			//change val of line 34 to use diff algos
};

const int NoMatrix = 50; // controls the amount of matrices 
const int code = GREEDYCOMPARE;			// used for case condition to run particular algos
const float DecayRate = 0.01; //setting the decay rate
const float Pheromone = 1; //setting the pheromone factor of the ants
long long unsigned int busyCount; // setting the busy count 
const int R = 100; //setting the Radius
const int TimeTakenteps = 1000; //setting the time steps 
const int NoAnts = 100; // setting the number of ants
const int N = 5;		//busy work count
const int MAXVAL = 10000;
double MatrixCost[NoMatrix][NoMatrix] = { 0 };    // setting Matrix which will hold the weight of a path
int indexes[NoMatrix] = { 0 };                 // array to track different brute force paths
int SP[NoMatrix] = { 0 };			// matrix to store the shortest path
Node points[NoMatrix] = { 0 };                 // array to hold the coordinates for the EuclideanDistance algo
double SPCost = -1;               // var to store cost of shortest path

using namespace std;

void doBusyWork(void)			//counts number of busy count work
{
	for (int k = 0; k < N; k++)
		for (int j = 0; j < N; j++)
			busyCount++;
}

void DisplayMatrix(double MatrixCost[][NoMatrix], long long int Ind)
{
	for (long long i = 0; i < Ind; i++) {
		for (long long k = 0; k < Ind; k++)
			printf("%8.2f", MatrixCost[i][k]);
		printf("\n");
	}
}

void swap(int* x, int* y)
{
	int temp;
	temp = *x;
	*x = *y;
	*y = temp;
}

void TspBruteForce(int indexes[NoMatrix], double MatrixCost[][NoMatrix], int lx, int mx)
{
	if (lx == mx) {
		double cost = 0;
		int Iterate = 0;
		int index = 0;


		while (Iterate <= mx) {
			cost += MatrixCost[index][indexes[Iterate]];
			index = indexes[Iterate];
			Iterate++;
		}


		if (SPCost < 0) {
			SPCost = cost;
			memcpy(SP, indexes, sizeof(indexes) * mx);
		}
		else if (SPCost > cost) {
			SPCost = cost;
			memcpy(SP, indexes, sizeof(indexes) * mx);
		}
	}
	else {
		for (int i = lx; i < mx; i++) {
			swap((indexes[lx]), (indexes[i]));
			TspBruteForce(indexes, MatrixCost, lx + 1, mx);
			swap((indexes[lx]), (indexes[i]));
		}
	}
}

void TspBruteForceWorstCase(int indexes[NoMatrix], double MatrixCost[][NoMatrix], int lx, int mx)
{
	if (lx == mx) {
		double cost = 0;
		int Iterate = 0;
		int index = 0;


		while (Iterate <= mx) {
			cost += MatrixCost[index][indexes[Iterate]];
			index = indexes[Iterate];
			Iterate++;
		}
		if (SPCost < 0) {
			SPCost = cost;
			memcpy(SP, indexes, sizeof(indexes) * mx);
		}
		else if (SPCost < cost) {
			SPCost = cost;
			memcpy(SP, indexes, sizeof(indexes) * mx);
		}
	}
	else {
		for (int i = lx; i < mx; i++) {
			swap((indexes[lx]), (indexes[i]));
			TspBruteForceWorstCase(indexes, MatrixCost, lx + 1, mx);
			swap((indexes[lx]), (indexes[i]));
		}
	}
}

void TspGreedy(double MatrixCost[][NoMatrix], long long int Ind)
{


	double smallest = MAXVAL + 1;
	int index = 0;
	int tempIndex = 0;
	int* visited = new int[Ind];
	for (int i = 0; i < Ind; i++)
		visited[i] = 0;


	//visited[Ind + 1] = { 0 };
	int* path = new int[Ind];
	for (int i = 0; i < Ind; i++)
		path[i] = 0;

	visited[0]++;


	for (long long int j = 0; j < Ind - 1; j++) {
		for (long long int k = 0; k < Ind; k++) {
			if (k != index) {
				if (MatrixCost[index][k] < smallest && visited[k] == 0) {
					smallest = MatrixCost[index][k];
					tempIndex = k;
				}
			}
		}
		index = tempIndex;
		visited[index]++;
		path[j + 1] = index;
		smallest = MAXVAL + 1;
	}
	double cost = 0;
	int pathIndex = 0;
	int ACC = 1;
	while (ACC < Ind) {
		cost += MatrixCost[pathIndex][path[ACC]];
		pathIndex = path[ACC];
		ACC++;
	}
	cost += MatrixCost[pathIndex][0];




	if (CONTINUOUS) {
		printf("Greedy: Best path cost: %f\n", cost);
		printf("Greedy: Best path: ");
		for (int l = 0; l < Ind; l++)
			if (visited[l] > 0)
				printf("%d -> ", path[l]);
		printf("0\n");
	}
}


void RandomMatrixCost(double MatrixCost[][NoMatrix], int Ind, int maxValue)
{
	int RandomNum;
	for (long long int i = 0; i < Ind; i++) {
		for (long long int j = 0; j < Ind; j++) {
			if (i == j)
				MatrixCost[i][j] = 0;
			else {
				RandomNum = rand() % maxValue;
				MatrixCost[i][j] = RandomNum;
				MatrixCost[j][i] = RandomNum;
			}
		}
	}
}

void RandomEuclideanDistanceMatrixCost(double MatrixCost[][NoMatrix], int Ind, int maxValue) //generate random coordinates and then insert EuclideanDistance distance
{
	for (long long int i = 0; i < Ind; i++) {
		points[i].x = rand() % maxValue;
		points[i].y = rand() % maxValue;
	}
	for (long long int i = 0; i < Ind; i++) {
		for (long long int j = 0; j < Ind; j++) {
			if (i != j)
				MatrixCost[i][j] = sqrt(pow(points[j].x - points[i].x, 2) + pow(points[j].y - points[i].y, 2));
		}
	}
}

void RandomCircularGraphMatrixCost(double MatrixCost[][NoMatrix], long long int Ind, int maxValue, int radius)
{
	double Angle = 2 * PI / Ind;
	int SuccNode;
	Node x;
	int* used = new int[Ind + 1];

	for (int i = 0; i < Ind + 1; i++)
		used[i] = 0;
	//used[Ind] = { 0 };
	int* bestPath = new int[Ind];
	for (int i = 0; i < Ind + 1; i++)
		bestPath[i] = 0;
	//bestPath[Ind] = { 0 };
	for (long long int i = 0; i < Ind; i++) {
		points[i].x = radius * sin(i * Angle);
		points[i].y = radius * cos(i * Angle);
	}
	for (long long int i = 0; i < Ind; i++)
		bestPath[i] = i;

	for (long long int i = 1; i < Ind / 2; i += 2) {
		SuccNode = 1 + (rand() % Ind - 1);
		while (used[SuccNode] >= 1 || SuccNode == 0 || i == SuccNode)
			SuccNode = 1 + (rand() % Ind - 1);
		x = points[SuccNode];
		points[SuccNode] = points[i];
		points[i] = x;
		swap(bestPath[SuccNode], bestPath[i]);
		used[SuccNode]++;
		used[i]++;
	}
	double LeastDistance = -1;
	for (long long int i = 0; i < Ind; i++) {
		for (long long int j = 0; j < Ind; j++) {
			if (i != j) {
				MatrixCost[i][j] = sqrt(pow(points[j].x - points[i].x, 2) + pow(points[j].y - points[i].y, 2));
				if (MatrixCost[i][j] < LeastDistance || LeastDistance == -1)
					LeastDistance = MatrixCost[i][j];
			}
		}
	}
	if (CONTINUOUS) {
		printf("\nBest path: ");
		for (long long int i = 0; i < Ind; i++) {
			printf("%d -> ", bestPath[i]);
		}
		printf("0\n");
		printf("Best path points: ");
		for (long long int k = 0; k < Ind - 1; k++)
			printf("(%.2f,%.2f) -> ", points[k].x, points[k].y);
		printf("(%.2f,%.2f)\n", points[Ind].x, points[Ind].y);
		printf("Expected Cost: (%f * %d) = %f\n\n", LeastDistance, Ind, (int)(LeastDistance * Ind));
	}
}

unsigned long long int fact(unsigned int x)
{
	unsigned long long int result = 1;
	for (unsigned long long int i = 2; i <= x; i++)
		result *= i;
	return result;
}



int main(int argc, char** argv)
{
	double MaxTrialTime = .250; // in seconds
	long long int MaxTrialCount = 1,// 1000000,
		MinN = 1,
		trial;
	clock_t TimeTakentampSplit, trialSetStart;
	const long long  Max_Array = 100000;
	double TimeTakenSplit, SetTrialCount, SetTrialTime, SetDummyCount, SetDummyTime, TrialTimeAverage, DummyTrialTimeAverage, TrialEstimatedTime;
	double TimeTaken[NoMatrix] = { 0 };
	int index = 1;
	srand(time(NULL));
	printf("+-------------------------------------------------------------------------------------------+\n");
	printf("| %20s | %20s | %20s | %20s |\n", "N", "Time", "Exprimental Double Ratio", "Expected Double Ratio");
	printf("+-------------------------------------------------------------------------------------------+\n");

	for (long long int HBZ = 3; HBZ < NoMatrix; HBZ++) {
		TimeTakenSplit = 0.0;
		trialSetStart = clock();
		for (trial = 0; trial < MaxTrialCount && TimeTakenSplit < MaxTrialTime; trial++) {
			busyCount = 0;
			switch (code) {
			case GREEDY:
				if (EUCLIDEANDISTANCE)
					RandomEuclideanDistanceMatrixCost(MatrixCost, HBZ, R);
				else
					RandomCircularGraphMatrixCost(MatrixCost, HBZ, MAXVAL, R);
				TspGreedy(MatrixCost, HBZ);
				if (!DISPLAYTABLE && CONTINUOUS)
					puts("---------------------------------------------------------------------------------------------");
				
	
				break;
			case GREEDYCOMPARE:
				for (int STrial = 0; STrial < 10; STrial++)
				{
					if (EUCLIDEANDISTANCE)
						RandomEuclideanDistanceMatrixCost(MatrixCost, HBZ, R);
					else
						RandomCircularGraphMatrixCost(MatrixCost, HBZ, MAXVAL, R);
					if (!DISPLAYTABLE && CONTINUOUS)
					{
						puts("---------------------------------------------------------------------------------------------");
						printf("N: %d\n", HBZ);
						puts("---------------------------------------------------------------------------------------------");
					}
					for (int i = 0; i < HBZ - 1; i++)
					{
						indexes[i] = i + 1;
					}
					TspBruteForce(indexes, MatrixCost, 0, HBZ - 1);
					printf("BruteForce: Best path cost: %f \n", SPCost);
					printf("BruteForce: Best path: 0 -> ");
					for (int i = 0; i < HBZ - 1; i++)
						printf("%d -> ", SP[i]);
					puts("0");


					memset(SP, 0, sizeof(SP));
					SPCost = -1;
					TspGreedy(MatrixCost, HBZ);
				}
				break;

			case BRUTEFORCE:
				if (EUCLIDEANDISTANCE)
					RandomEuclideanDistanceMatrixCost(MatrixCost, HBZ, R);
				else
					RandomCircularGraphMatrixCost(MatrixCost, HBZ, MAXVAL, R);
				for (int i = 0; i < HBZ - 1; i++)
				{
					indexes[i] = i + 1;
				}
				TspBruteForce(indexes, MatrixCost, 0, HBZ - 1);
				if (CONTINUOUS) {
					printf("BrutForce: shortest Path: 0 -> ");
					for (int i = 0; i < HBZ; i++) {
						if (i == HBZ - 1)
							printf("%d\n", SP[i]);
						else
							printf("%d -> ", SP[i]);
					}
					printf("BruteForce: shortest Path Cost: %f\n\n", SPCost);
					DisplayMatrix(MatrixCost, HBZ);
					if (!DISPLAYTABLE)
						puts("---------------------------------------------------------------------------------------------");
				}
				SPCost = -1;
				memset(SP, 0, sizeof(SP));
				break;

			case ALLALGOS:
				if (EUCLIDEANDISTANCE)
					RandomEuclideanDistanceMatrixCost(MatrixCost, HBZ, R);
				else
					RandomCircularGraphMatrixCost(MatrixCost, HBZ, MAXVAL, R);
				printf("-------------------------\n\n");
				
				for (int i = 0; i < HBZ - 1; i++)
				{
					indexes[i] = i + 1;
				}
				TspBruteForce(indexes, MatrixCost, 0, HBZ - 1);
				printf("BruteForce: Best path cost: %f \n", SPCost);
				printf("BruteForce: Best path: 0 -> ");
				for (int i = 0; i < HBZ - 1; i++)
					printf("%d -> ", SP[i]);
				puts("0");
				memset(SP, 0, sizeof(SP));
				SPCost = -1;
				for (int i = 0; i < HBZ - 1; i++)
				{
					indexes[i] = i + 1;
				}
				TspBruteForceWorstCase(indexes, MatrixCost, 0, HBZ - 1);
				printf("WorstCase:  Worst path cost: %f \n", SPCost);
				printf("WorseCase:  Worst path: 0 -> ");
				for (int i = 0; i < HBZ - 1; i++)
					printf("%d -> ", SP[i]);
				puts("0");
				memset(SP, 0, sizeof(SP));
				SPCost = -1;
				TspGreedy(MatrixCost, HBZ);
				break;
			

			default:
				break;
			}

			TimeTakentampSplit = clock();
			TimeTakenSplit = (TimeTakentampSplit - trialSetStart) / (double)CLOCKS_PER_SEC;
		}
		SetTrialCount = trial;
		SetTrialTime = TimeTakenSplit;
		TrialTimeAverage = SetTrialTime / SetTrialCount;
		TimeTakenSplit = 0.0;

		trialSetStart = clock();
		for (trial = 0; trial < SetTrialCount && TimeTakenSplit < MaxTrialTime; trial++) {
			TimeTakentampSplit = clock();
			TimeTakenSplit = (TimeTakentampSplit - trialSetStart) / (double)CLOCKS_PER_SEC;
		}
		SetDummyCount = trial;
		SetDummyTime = TimeTakenSplit;
		DummyTrialTimeAverage = SetDummyTime / SetDummyCount;
		TrialEstimatedTime = TrialTimeAverage - DummyTrialTimeAverage;
		TimeTaken[index] = TrialEstimatedTime;
		if (DISPLAYTABLE) {
			if (HBZ > 3) {
				switch (code) {
				case GREEDY:
					printf("| %20llu | %20.6f | %20.6f | %20.6f |\n", HBZ, TimeTaken[index], TimeTaken[index] / TimeTaken[index - 1], pow(HBZ, 2) / pow(HBZ - 1, 2));
					break;
				case BRUTEFORCE:
					printf("| %20llu | %20.6f | %20.6f | %20.6f |\n", HBZ, TimeTaken[index], TimeTaken[index] / TimeTaken[index - 1], (double)(fact(HBZ) / fact(HBZ - 1)));
					break;
				default:
					break;
				}
			}
		}
		index++;
		busyCount = 0;
	}
}

