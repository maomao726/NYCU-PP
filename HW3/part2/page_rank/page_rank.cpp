#include "page_rank.h"

#include <stdlib.h>
#include <cmath>
#include <omp.h>
#include <utility>

#include "../common/CycleTimer.h"
#include "../common/graph.h"

// pageRank --
//
// g:           graph to process (see common/graph.h)
// solution:    array of per-vertex vertex scores (length of array is num_nodes(g))
// damping:     page-rank algorithm's damping parameter
// convergence: page-rank algorithm's convergence threshold
//
void pageRank(Graph g, double *solution, double damping, double convergence)
{

  // initialize vertex weights to uniform probability. Double
  // precision scores are used to avoid underflow for large graphs

  int numNodes = num_nodes(g);
  double equal_prob = 1.0 / numNodes;
  #pragma omp parallel for
  for (int i = 0; i < numNodes; ++i)
  {
    solution[i] = equal_prob;
  }

  /*
     For PP students: Implement the page rank algorithm here.  You
     are expected to parallelize the algorithm using openMP.  Your
     solution may need to allocate (and free) temporary arrays.

     Basic page rank pseudocode is provided below to get you started:

     // initialization: see example code above
     score_old[vi] = 1/numNodes;

     while (!converged) {

       // compute score_new[vi] for all nodes vi:
       score_new[vi] = sum over all nodes vj reachable from incoming edges
                          { score_old[vj] / number of edges leaving vj  }
       score_new[vi] = (damping * score_new[vi]) + (1.0-damping) / numNodes;

       score_new[vi] += sum over all nodes v in graph with no outgoing edges
                          { damping * score_old[v] / numNodes }

       // compute how much per-node scores have changed
       // quit once algorithm has converged

       global_diff = sum over all nodes vi { abs(score_new[vi] - score_old[vi]) };
       converged = (global_diff < convergence)
     }

   */
  bool converge = false;
  double* solution_new = (double*)malloc(sizeof(double) * numNodes);
  int* outgoing_countList = (int*)malloc(sizeof(int) * numNodes);
  const Vertex** in_start = (const Vertex**)malloc(sizeof(Vertex*) * numNodes);
  const Vertex** in_end = (const Vertex**)malloc(sizeof(Vertex*) * numNodes);
  // double solution_new[numNodes];
  // int outgoing_countList[numNodes];

  #pragma omp parallel for
  for(int i = 0; i < numNodes; i++)
  {
    outgoing_countList[i] = outgoing_size(g, i);
    in_start[i] = incoming_begin(g, i);
    in_end[i] = incoming_end(g, i);
  }

  while(!converge)
  {
    double global_diff = 0.0;
    double no_out = 0.0;
    
    #pragma omp parallel for
    for(int i = 0; i < numNodes; i++)
    {
      solution_new[i] = 0.0;
      if(outgoing_countList[i] == 0)
        #pragma omp atomic
        no_out += damping * solution[i] / numNodes;
    }
   
    #pragma omp parallel for reduction(+:global_diff)
    for(int i = 0; i < numNodes; i++)
    {

      for(const Vertex* vj = in_start[i]; vj != in_end[i]; vj++)
      {  
        double temp = (solution[*vj] / outgoing_countList[*vj]);
        solution_new[i] += temp;
      }

      solution_new[i] = (damping * solution_new[i]) + (1.0 - damping) / numNodes;
      solution_new[i] += no_out;
      
      global_diff += fabs(solution_new[i] - solution[i]);
    }

    #pragma omp parallel for
    for(int i = 0; i < numNodes; i++)
    {
      solution[i] = solution_new[i];
    }
    converge = (global_diff < convergence);
  }
  
}
