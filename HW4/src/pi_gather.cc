#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>

int main(int argc, char **argv)
{
    // --- DON'T TOUCH ---
    MPI_Init(&argc, &argv);
    double start_time = MPI_Wtime();
    double pi_result;
    long long int tosses = atoi(argv[1]);
    int world_rank, world_size;
    // ---

    // TODO: MPI init
    MPI_Comm_size( MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank( MPI_COMM_WORLD, &world_rank);
    long long local_count[world_size] = {0};
    unsigned int seed = world_rank * time(NULL);
    long long int total_count = 0;
    for(int i = 0; i < tosses / world_size; i++)
    {
        double x = ((double)rand_r(&seed) / (double)RAND_MAX);
        double y = ((double)rand_r(&seed) / (double)RAND_MAX);
        if(x * x + y * y <= 1.0)
        {
            local_count[world_rank]++;
        }
    }
    MPI_Gather(local_count + world_rank, 1, MPI_LONG_LONG, local_count, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    // TODO: use MPI_Gather
    if (world_rank == 0)
    {
        for(int i = 0; i < world_size; i++)
            total_count += local_count[i];
    }

    if (world_rank == 0)
    {
        // TODO: PI result
        pi_result =  4 * (double)total_count / (double)tosses;

        // --- DON'T TOUCH ---
        double end_time = MPI_Wtime();
        printf("%lf\n", pi_result);
        printf("MPI running time: %lf Seconds\n", end_time - start_time);
        // ---
    }
    
    MPI_Finalize();
    return 0;
}
