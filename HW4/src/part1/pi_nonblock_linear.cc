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

    if (world_rank > 0)
    {
        // TODO: MPI workers
        MPI_Send(&count, 1, MPI_LONG_LONG, 0, 0, MPI_COMM_WORLD);
    }
    else if (world_rank == 0)
    {
        // TODO: non-blocking MPI communication.
        // Use MPI_Irecv, MPI_Wait or MPI_Waitall.
        MPI_Request requests[world_size - 1];
        MPI_Status status[world_size - 1];
        long long int local_count[world_size - 1];
        for(int i = 0; i < world_size - 1; i++)
            MPI_Irecv(local_count+i, 1, MPI_LONG_LONG, i+1, 0, MPI_COMM_WORLD, requests+i);
        MPI_Waitall(world_size - 1, requests, status);

        total_count = count;
        for(int i = 0; i < world_size - 1; i++)total_count += local_count[i];
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
