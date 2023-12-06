#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void construct_matrices(int *n_ptr, int *m_ptr, int *l_ptr,
                        int **a_mat_ptr, int **b_mat_ptr)
{
    
    int world_rank, world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    

    if(world_rank == 0)
    {
        scanf("%d %d %d", n_ptr, m_ptr,l_ptr);
    }
    if(world_rank == 0)
    {
        *a_mat_ptr = (int*) malloc(sizeof(int) * (*n_ptr) * (*m_ptr));
        *b_mat_ptr = (int*)malloc(sizeof(int) * (*m_ptr) * (*l_ptr));

        for(int i = 0; i < *n_ptr; i++)
        {
            int row = (*m_ptr) * i;
            for(int j = 0; j < *m_ptr; j++)
            {
                scanf("%d", *a_mat_ptr + row + j);
            }
        }
        for(int i = 0; i < *m_ptr; i++)
        {
            int row = (*l_ptr) * i;
            for(int j = 0; j < *l_ptr; j++)
            {
                scanf("%d", *b_mat_ptr + row + j);
            }
        }
    }
    
    MPI_Bcast(n_ptr, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(m_ptr, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(l_ptr, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(world_rank > 0)
    {
        *b_mat_ptr = (int*)malloc(sizeof(int) * (*m_ptr) * (*l_ptr));
    }
    MPI_Bcast(*b_mat_ptr, (*m_ptr) * (*l_ptr), MPI_INT, 0, MPI_COMM_WORLD);
    
    int avg_row = (*n_ptr) / world_size;
    int extra = (*n_ptr) % world_size;
    int row = 0;

    if(world_rank == 0)
    {
        int offset = 0;
        for(int dest = 0; dest < world_size; dest++)
        {
           int send_row = (dest < extra) ?  avg_row+1 : avg_row;
           
            if(dest != 0)
                MPI_Send(*a_mat_ptr+offset, send_row * (*m_ptr), MPI_INT, dest, 0, MPI_COMM_WORLD);
            offset += send_row * (*m_ptr);
        }
    }
    else if(world_rank > 0)
    {
        row = (world_rank < extra) ?  avg_row+1 : avg_row;
        *a_mat_ptr = (int*)malloc(sizeof(int) * row * (*m_ptr));
        MPI_Recv(*a_mat_ptr, row * (*m_ptr), MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

void matrix_multiply(const int n, const int m, const int l,
                     const int *a_mat, const int *b_mat)
{
    int world_rank, world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int avg_row = n / world_size;
    int extra = n % world_size;
    int row = (world_rank < extra) ?  avg_row+1 : avg_row;

    int *local_result = (int*)calloc(row * l, sizeof(int));
    for(int i = 0; i < row; i++)
    {
        for(int k = 0; k < l; k++)
        {
            for(int j = 0; j < m; j++)
            {
                local_result[i*l+k] += a_mat[i*m+j] * b_mat[j*l+k];
            }
        }
    }

    if(world_rank == 0)
    {
        int *result = (int*)malloc(sizeof(int) * n * l);
        memcpy(result, local_result, row * l * sizeof(int));
        int offset = row*l;
        for(int src = 1; src < world_size; src++)
        {
            int send_row = (src < extra) ?  avg_row+1 : avg_row;
            MPI_Recv(result+offset, send_row*l, MPI_INT, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            offset += send_row * l;
        }
        for(int i = 0; i < n; i++)
        {
            int offset = i*l;
            for(int j = 0; j < l; j++)
            {
                printf("%d ", result[offset+j]);
            }
            printf("\n");
        }
        
    }
    if(world_rank > 0)
    {
        MPI_Send(local_result, row * l, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
}

void destruct_matrices(int *a_mat, int *b_mat)
{
    free(a_mat);
    free(b_mat);
}