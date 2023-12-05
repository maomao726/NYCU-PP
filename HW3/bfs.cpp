#include "bfs.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstddef>
#include <omp.h>

#include "../common/CycleTimer.h"
#include "../common/graph.h"

#define ROOT_NODE_ID 0
#define NOT_VISITED_MARKER -1

void vertex_set_clear(vertex_set *list)
{
    list->count = 0;
}

void vertex_set_init(vertex_set *list, int count)
{
    list->max_vertices = count;
    list->vertices = (int *)malloc(sizeof(int) * list->max_vertices);
    vertex_set_clear(list);
}

// Take one step of "top-down" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighboring vertices to the
// new_frontier.
void top_down_step(
    Graph g,
    vertex_set *frontier,
    vertex_set *new_frontier,
    int *distances)
{
    int num_threads = omp_get_max_threads();
    vertex_set *local_frontier = (vertex_set*) malloc(sizeof(vertex_set) * num_threads);
    #pragma omp parallel for
    for(int i = 0; i < num_threads; i++)
    {
        vertex_set_init(&local_frontier[i], g->num_nodes);
    }

    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        vertex_set_clear(&local_frontier[id]);
        #pragma omp for
        for (int i = 0; i < frontier->count; i++)
        {
            int node = frontier->vertices[i];

            int start_edge = g->outgoing_starts[node];
            int end_edge = (node == g->num_nodes - 1)
                            ? g->num_edges
                            : g->outgoing_starts[node + 1];

            // attempt to add all neighbors to the new frontier
            for (int neighbor = start_edge; neighbor < end_edge; neighbor++)
            {
                int outgoing = g->outgoing_edges[neighbor];

                if(distances[outgoing] == NOT_VISITED_MARKER)
                {
                    if (__sync_bool_compare_and_swap(&distances[outgoing], NOT_VISITED_MARKER, distances[node] + 1))
                    {
                        int index = __sync_fetch_and_add(&local_frontier[id].count, 1);
                        local_frontier[id].vertices[index] = outgoing;
                    }
                }
                
            }
        }
    }
    int* start = (int*)malloc(sizeof(int) * num_threads);
    start[0] = 0;
    for(int i = 1; i < num_threads; i++)
    {
        start[i] = (start[i-1] + local_frontier[i-1].count);
    }
    new_frontier->count = start[num_threads-1] + local_frontier[num_threads-1].count;

    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        memcpy(&new_frontier->vertices[start[id]], local_frontier[id].vertices, sizeof(int)*local_frontier[id].count);
    }

}

// Implements top-down BFS.
//
// Result of execution is that, for each node in the graph, the
// distance to the root is stored in sol.distances.
void bfs_top_down(Graph graph, solution *sol)
{

    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set *frontier = &list1;
    vertex_set *new_frontier = &list2;

    // initialize all nodes to NOT_VISITED
    #pragma omp parallel for
    for (int i = 0; i < graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    // setup frontier with the root node
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;

    while (frontier->count != 0)
    {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(new_frontier);

        top_down_step(graph, frontier, new_frontier, sol->distances);

#ifdef VERBOSE
        double end_time = CycleTimer::currentSeconds();
        printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // swap pointers
        vertex_set *tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
    }
}

void bottom_up_step(
    Graph g,
    vertex_set *frontier,
    vertex_set *new_frontier,
    int *distances, int dist
)
{

    int num_threads = omp_get_max_threads();
    vertex_set *local_frontier = (vertex_set*) malloc(sizeof(vertex_set) * num_threads);
    #pragma omp parallel for
    for(int i = 0; i < num_threads; i++)
    {
        vertex_set_init(&local_frontier[i], g->num_nodes);
    }

    #pragma omp parallel 
    {
        int id = omp_get_thread_num();
        vertex_set_clear(&local_frontier[id]);
        
        #pragma omp for schedule(dynamic, 512)
        for(int v = 0; v < g->num_nodes; v++)
        {
            if(distances[v] == NOT_VISITED_MARKER)
            {
                int start_edge = g->incoming_starts[v];
                int end_edge = (v == g->num_nodes - 1)
                                ? g->num_edges
                                : g->incoming_starts[v + 1];
                for(int u = start_edge; u < end_edge; u++)
                {
                    int in_node = g->incoming_edges[u];
                    if(distances[in_node] != NOT_VISITED_MARKER)
                    {
                        int index = local_frontier[id].count++;
                        local_frontier[id].vertices[index] = v;
                        break;
                    }
                }
            }
        }
        #pragma omp barrier
        for(int i = 0; i < local_frontier[id].count; i++)
        {
            distances[local_frontier[id].vertices[i]] = dist;
        }
    }

    int* start = (int*)malloc(sizeof(int) * num_threads);
    start[0] = 0;
    for(int i = 1; i < num_threads; i++)
    {
        start[i] = (start[i-1] + local_frontier[i-1].count);
    }
    new_frontier->count = start[num_threads-1] + local_frontier[num_threads-1].count;

    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        memcpy(&new_frontier->vertices[start[id]], local_frontier[id].vertices, sizeof(int)*local_frontier[id].count);
    }
    
}


void bfs_bottom_up(Graph graph, solution *sol)
{
    // For PP students:
    //
    // You will need to implement the "bottom up" BFS here as
    // described in the handout.
    //
    // As a result of your code's execution, sol.distances should be
    // correctly populated for all nodes in the graph.
    //
    // As was done in the top-down case, you may wish to organize your
    // code by creating subroutine bottom_up_step() that is called in
    // each step of the BFS process.
    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set *frontier = &list1;
    vertex_set *new_frontier = &list2;

    // initialize all nodes to NOT_VISITED
    #pragma omp parallel for
    for (int i = 0; i < graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    // setup frontier with the root node
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;

    int dist = 1;
    while (frontier->count != 0)
    {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(new_frontier);
        bottom_up_step(graph, frontier, new_frontier, sol->distances, dist);

#ifdef VERBOSE
        double end_time = CycleTimer::currentSeconds();
        printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // swap pointers
        vertex_set *tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
        dist++;
    }
}


void bfs_hybrid(Graph graph, solution *sol)
{
    // For PP students:
    //
    // You will need to implement the "hybrid" BFS here as
    // described in the handout.
    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set *frontier = &list1;
    vertex_set *new_frontier = &list2;

    // initialize all nodes to NOT_VISITED
    #pragma omp parallel for
    for (int i = 0; i < graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    // setup frontier with the root node
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;

    int dist = 1;
    while (frontier->count != 0)
    {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(new_frontier);
        if(frontier->count / float(graph->num_nodes) >= 0.1)
            bottom_up_step(graph, frontier, new_frontier, sol->distances, dist);
        else
            top_down_step(graph, frontier, new_frontier, sol->distances);

#ifdef VERBOSE
        double end_time = CycleTimer::currentSeconds();
        printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // swap pointers
        vertex_set *tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
        dist++;
    }
}