#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <omp.h>

int num_threads;

int main()
{
	srand(time(NULL));
	num_threads = 16;
	omp_set_num_threads(num_threads);
	printf("------------------- use thread-------------------\n");
#pragma omp parallel 
	{
#pragma omp for
		for (int i = 0; i < 30; i++)
		{
			printf("thread id : %d   random number : %d\n", omp_get_thread_num(), rand() % 30);
		}

	}

	printf("\n ---------------------- non use thread ----------------\n");
	for(int i = 0; i < 30; i++)
	{
		printf("random number : %d\n", rand() % 30);
	}


	return 0;
}