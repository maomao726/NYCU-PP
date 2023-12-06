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
    int temp = 0;
    int level = 1;
    int temp_rk = world_rank;
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

    
    
    while( temp_rk % 2 == 0 && level < world_size)
    {
        printf("Receive: %d, %d, %d from %d\n", temp_rk, level, world_rank, (temp_rk+1) * level);
        
        MPI_Recv(&temp, 1, MPI_LONG_LONG, (temp_rk+1) * level, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        count += temp;
        //printf("before: %d / rank: %d\n", i, world_rank);
        temp_rk = temp_rk / 2;
        //printf("after: %d / rank: %d\n", i, world_rank);
        level *= 2;
    }
    if(world_rank != 0)
    {
        printf("Send: %d, %d, %d to %d\n", temp_rk, level, world_rank, (temp_rk-1) * level);
        MPI_Send(&count, 1, MPI_LONG_LONG, (temp_rk-1) * level, 0, MPI_COMM_WORLD);
    }
        

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
