#include <cstdio>
#include <stdlib.h>
#include <ctime>
#include <pthread.h>

using namespace std;

unsigned int seed_seed = time(NULL);
long long *result_list;

typedef struct arg_struct {
	long long num_tosses;
	long long *result;
	arg_struct(long long num, long long *res) : num_tosses(num), result(res){};
	arg_struct(){};
} arg_struct;

void* monte_carlo_search(void *args)
{
	arg_struct* arg = (arg_struct*) args;
	long long count_inCircle = 0;
	unsigned int seed = rand_r(&seed_seed);
	
	for(long long toss = 0; toss < (arg->num_tosses); toss++)
	{
		double x = ((double)rand_r(&seed) / (double)RAND_MAX);
		double y = ((double)rand_r(&seed) / (double)RAND_MAX);
		if(x * x + y * y <= 1.0)
		{
			count_inCircle++;
		}
	}
	*(arg->result) = count_inCircle;
	pthread_exit(0);
}

int main(const int argc, const char* argv[])
{
	//argv[0]: num of threads, argv[1]: num of tosses
	long long number_of_tosses = atoll(argv[2]);
	int num_thread = atoi(argv[1]);
	long long total_in_circle = 0;		
	arg_struct* args_list= new arg_struct[num_thread];

	result_list = (long long *) new long long [num_thread];
	pthread_t *threads = new pthread_t[num_thread];

	long long quotient = number_of_tosses / num_thread;
	int remainder = number_of_tosses % num_thread;



	for(int i = 0; i < num_thread; i++)
	{
		long long thread_toss = (remainder != 0) ? quotient + (i < remainder) : quotient;
		args_list[i].num_tosses = thread_toss;
		args_list[i].result = &result_list[i];
	}

	for(int i = 0; i < num_thread; i++)
	{
		pthread_create(&threads[i], NULL, monte_carlo_search, &args_list[i]);
	}

	for(int i = 0; i < num_thread; i++)
	{
		pthread_join(threads[i], NULL);
		total_in_circle += result_list[i];
	}

	double pi_estimate = 4 * (double)total_in_circle / (double)number_of_tosses;
	printf("%.8f\n", pi_estimate);
	return 0;
}
