//============================================================================
// Name        : Smooth.cpp
// Author      : George Rokos
// Description : 2D Vertex-Smoothing kernel - Smart Laplacian variant
//============================================================================

#include <algorithm>
#include <cmath>
#include <list>
#include <stdio.h>
#include <string.h>
#include <future>

#include "SVD2x2.hpp"
#include "Smooth.hpp"

static bool * vertices_in_neighberhood; //vertices in the neighberhood
static bool * to_examine_all; //vertices to be examined
static int vertices; //number of vertices left to examine

static std::vector<std::vector<int> > slices;

void populate_vertices(Mesh *mesh){
  //put all vertices from the mesh into `to_examine_all`
  to_examine_all = new bool[mesh->NNodes];
  memset(to_examine_all, true, mesh->NNodes);

  vertices_in_neighberhood = new bool[mesh->NNodes];
}

//greedy colouring
void select_vertices(Mesh *mesh, int colour){
  //reset neighberhoods
  memset(vertices_in_neighberhood, false, mesh->NNodes);

  //A bit of optimization
  int pref[4] = {6115, 13024, 14178, 14700};
  if(colour == 0)
    for(size_t i = 0; i < 4; i++){
      size_t it = pref[i];
      if( !vertices_in_neighberhood[it] && to_examine_all[it] ){
        to_examine_all[it] = false;
        vertices--;

        slices[colour].push_back(it);
        vertices_in_neighberhood[it] = true;
        for(std::vector<size_t>::const_iterator nit=mesh->NNList[it].begin(); nit!=mesh->NNList[it].end(); ++nit){
          vertices_in_neighberhood[*nit] = true;
        }
      }
    }

  for(size_t it = 0; it < mesh->NNodes; it++){
    //if a vertex is not in the neighberhood and has to be examined
    if( !vertices_in_neighberhood[it] && to_examine_all[it] ){
      //mark vertex as examined
      to_examine_all[it] = false;
      vertices--;

      //add to a slice
      slices[colour].push_back(it);
      //put in the neighberhood
      vertices_in_neighberhood[it] = true;
      //iterate over adjacent vertices and put them in the neighberhood
      for(std::vector<size_t>::const_iterator nit=mesh->NNList[it].begin(); nit!=mesh->NNList[it].end(); ++nit){
        vertices_in_neighberhood[*nit] = true;
      }
    }
  }
}


//memory write access to - mesh->coords...
void smooth_job(Mesh * mesh, size_t colour, double * quality_cache, bool * vertice_in_cache, int iter, int start_range, int end_range){
  for(size_t vi = start_range; vi < end_range; vi++){
    int vertex = slices[colour][vi];
    // If this is a corner node, it cannot be moved.
    if(mesh->isCornerNode(vertex))
      continue;

    // Find the quality of the worst element adjacent to vid
    double worst_q=1.0;

    //parallelize
    for(std::set<size_t>::const_iterator it=mesh->NEList[vertex].begin();
        it!=mesh->NEList[vertex].end(); ++it){
      double v_quality;

      //compute quality only not in the cache
      if( !vertice_in_cache[*it] ){
        v_quality = quality_cache[*it] = mesh->element_quality(*it);
        vertice_in_cache[*it] = true;
      } else{
        v_quality = quality_cache[*it];
      }

      worst_q = std::min(worst_q, v_quality);
    }

    const double * m0 = &mesh->metric[3*vertex]; //const

    const double x0 = mesh->coords[2*vertex];
    const double y0 = mesh->coords[2*vertex+1];

    double A[4] = {0.0, 0.0, 0.0, 0.0}; //const
    double q[2] = {0.0, 0.0};

    // Iterate over all edges and assemble matrices A and q.
    for(std::vector<size_t>::const_iterator it=mesh->NNList[vertex].begin();
        it!=mesh->NNList[vertex].end(); ++it){
      size_t il = *it;

      const double *m1 = &mesh->metric[3*il]; //const

      // Find the metric in the middle of the edge.
      // Vectorize
      double ml00 = 0.5*(m0[0] + m1[0]); //const
      double ml01 = 0.5*(m0[1] + m1[1]); //const
      double ml11 = 0.5*(m0[2] + m1[2]); //const

      double x = mesh->coords[2*il] - x0;
      double y = mesh->coords[2*il+1] - y0;

      // Calculate and accumulate the contribution of
      // this vertex to the barycentre of the cavity.
      //Vectorize
      q[0] += (ml00*x + ml01*y);
      q[1] += (ml01*x + ml11*y);

      //Vectorize
      if(iter == 0){
        A[0] += ml00; //const
        A[1] += ml01; //const
        A[3] += ml11; //const
      }
    }

    // The metric tensor is symmetric, i.e. ml01=ml10, so A[2]=A[1].
    A[2]=A[1];

    // Displacement vector for vid
    double p[2];

    svd_solve_2x2(vertex, A, p, q);

    /* If this is a surface vertex, restrict the displacement
     * to the surface. The new displacement is the projection
     * of the old displacement on the surface.
     */
    if(mesh->isSurfaceNode(vertex)){
      p[0] -= p[0]*fabs(mesh->normals[2*vertex]);
      p[1] -= p[1]*fabs(mesh->normals[2*vertex+1]);
    }

    // Actually change something
    // Update the coordinates
    mesh->coords[2*vertex] += p[0];
    mesh->coords[2*vertex+1] += p[1];

    double new_worst_q=1.0;

    //parallelize
    for(std::set<size_t>::const_iterator it=mesh->NEList[vertex].begin();
        it!=mesh->NEList[vertex].end(); ++it){

      //store in cache new quality measure
      double v_quality = quality_cache[*it] = mesh->element_quality(*it);
      vertice_in_cache[*it] = true;

      new_worst_q = std::min(new_worst_q, v_quality);
    }

    //Undo the change
    if(new_worst_q < worst_q){
      mesh->coords[2*vertex] -= p[0];
      mesh->coords[2*vertex+1] -= p[1];

      for(std::set<size_t>::const_iterator it=mesh->NEList[vertex].begin();
        it!=mesh->NEList[vertex].end(); ++it){
        vertice_in_cache[*it] = false;
      }
    }
  }
}

//possibly divide into 8-16 threads max
void spawn_threads(Mesh * mesh, size_t colour, double * quality_cache, bool * vertice_in_cache, int iter){
  std::vector<std::future<void> > futures;
  const int MAX_THREADS = 8;
  const int size = slices[colour].size();
  const double segment = size/MAX_THREADS;

  for(int j = 0; j < MAX_THREADS; j++){
    int start_range = (int)j*segment;
    int end_range = std::min( (int)((j+1)*segment), size );

    auto f = std::async( std::launch::async, smooth_job, mesh, colour, quality_cache, vertice_in_cache, iter, start_range, end_range );
    futures.push_back( std::move( f ) );
  }

  //barrier
  std::for_each(futures.begin(), futures.end(), [](std::future<void> & f)
  {
      f.wait();
  });
}

void smooth_parallel(Mesh* mesh, int niter){
  //Colouring phase
  int colour = 0;
  vertices = mesh->NNodes;

  populate_vertices( mesh );

  while( vertices > 0 ){
    std::vector<int> v;
    slices.push_back(v);

    select_vertices( mesh, colour );

    colour++;
  }

  printf("Colours: %d\n", colour);

  delete[] vertices_in_neighberhood;
  delete[] to_examine_all;

  //Execution phase
  svd_init( mesh->NNodes );

  double * quality_cache = new double[mesh->NElements];
  bool * vertice_in_cache = new bool[mesh->NElements];

  memset(quality_cache, 0, mesh->NElements);
  memset(vertice_in_cache, false, mesh->NElements);

  for(int iter = 0; iter < niter; iter++){
    for(size_t c = 0; c < colour; c++){
      spawn_threads( mesh, c, quality_cache, vertice_in_cache, iter );
    }
  }

  svd_teardown( mesh->NNodes );

  delete[] quality_cache;
  delete[] vertice_in_cache;
}

void smooth(Mesh* mesh, size_t niter){
  smooth_parallel(mesh, niter);
}
