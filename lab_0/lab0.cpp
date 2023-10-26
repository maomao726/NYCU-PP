#include <cstdio>
#include <random>
#include <ctime>

using namespace std;

int main()
{
	default_random_engine generator(time(NULL));
	uniform_real_distribution<double> distribution_xy(-1, 1);
	uniform_int_distribution<long long int> distribution_num(10e6,  0x7FFFFFFFFFFFFFFF);

	long long int number_of_tosses = distribution_num(generator);
        long long int number_in_circle = 0;		
	printf("num of tosses: %lld\n", number_of_tosses);
	
	for(auto toss = 0; toss < number_of_tosses; toss++)
	{
		double x = distribution_xy(generator);
		double y = distribution_xy(generator);
		double distance_squared = x * x + y * y;
		if(distance_squared <= 1)
		{
			number_in_circle++;
		}
	}

	double pi_estimate = 4 * number_in_circle / ((double) number_of_tosses);
	
	printf("pi_estimate: %.4f\n", pi_estimate);
	return 0;
}
