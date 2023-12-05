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
    long long count = 0;
    unsigned int seed = world_rank * time(NULL);
    long long int total_count = 0;
    for(int i = 0; i < tosses / world_size; i++)
    {
        double x = ((double)rand_r(&seed) / (double)RAND_MAX);
        double y = ((double)rand_r(&seed) / (double)RAND_MAX);
        if(x * x + y * y <= 1.0)
        {
            count++;
        }
    }
    int i = world_rank, level = 1;

    while( i % 2 == 0 && level * 2 <= world_size)
    {
        int temp = 0;
        MPI_Recv(&temp, 1, MPI_LONG_LONG, (i+1) * level, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        count += temp;
        i /= 2;
        level++;
    }
    if(world_rank != 0)
        MPI_Send(&count, 1, MPI_LONG_LONG, (i-1) * level, 0, MPI_COMM_WORLD);

    // TODO: binary tree redunction

    if (world_rank == 0)
    {
        // TODO: PI result
        pi_result =  4 * (double)count / (double)tosses;

        // --- DON'T TOUCH ---
        double end_time = MPI_Wtime();
        printf("%lf\n", pi_result);
        printf("MPI running time: %lf Seconds\n", end_time - start_time);
        // ---
    }

    MPI_Finalize();
    return 0;
}
