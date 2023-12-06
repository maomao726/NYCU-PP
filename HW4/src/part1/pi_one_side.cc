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

    MPI_Win win;

    // TODO: MPI init
    MPI_Comm_size( MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank( MPI_COMM_WORLD, &world_rank);
    
    unsigned int seed = world_rank * time(NULL);
    long long count = 0;
    long long int total_count = 0;
    long long int local_toss = tosses / world_size;
    for(long long int i = 0; i < local_toss; i++)
    {
        double x = ((double)rand_r(&seed) / (double)RAND_MAX);
        double y = ((double)rand_r(&seed) / (double)RAND_MAX);
        if(x * x + y * y <= 1.0)
        {
            count++;
        }
    }

    if (world_rank == 0)
    {
        // Master
        total_count = count;
        MPI_Win_create(&total_count, sizeof(long long), sizeof(long long), MPI_INFO_NULL,
          MPI_COMM_WORLD, &win);
    }
    else
    {
        // Workers
        MPI_Win_create(NULL, 0, 1, MPI_INFO_NULL, MPI_COMM_WORLD, &win);

        // Register with the master
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win);
        MPI_Accumulate(&count, 1, MPI_LONG_LONG, 0, 0, 1, MPI_LONG_LONG, MPI_SUM, win);
        MPI_Win_unlock(0, win);
    }

    MPI_Win_free(&win);

    if (world_rank == 0)
    {
        // TODO: handle PI result
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