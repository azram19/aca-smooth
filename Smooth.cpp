//============================================================================
// Name        : Smooth.cpp
// Author      : George Rokos
// Description : 2D Vertex-Smoothing kernel - Smart Laplacian variant
//============================================================================

#include <algorithm>
#include <cmath>
#include <list>

#include "SVD2x2.hpp"
#include "Smooth.hpp"

bool * vertices_in_neighberhood; //vertices in the neighberhood

std::list<int> to_examine_all; //all vertices to be examined
std::list<int> to_examine_now; //vertices to be examined by a population of threads

void populate_vertices(Mesh *mesh){
  //put all vertices from the mesh into `to_examine_all`
}

void select_vertices(Mesh *mesh){
  //iterate through vertices in 'to_examine_all'
  for (int &v : to_examine_all) {
    //if a vertex is not in the neighberhood
    if( !vertices_in_neighberhood[v] ){
      //add to the 'to_examine_now'
      to_examine_now.push_back(v);
      //put in the neighberhood
      vertices_in_neighberhood[v] = true;
      //iterate over adjacent vertices and put them in the neighberhood
      for( auto &va : mesh->NNList[v] ){
        vertices_in_neighberhood[va] = true;        
      }
    }
  }
}

void remove_vertices(){
  //remove items in `to_examine_now` from `to_examine_all`
}

void spawn_threads(){
  //for each element in to_examine_now spawn new smooth_job
  //barrier
}

void smooth_job(){

}

void smooth_parallel(Mesh* mesh, int iter){
  populate_vertices( mesh );
  
  while( to_examine_all.size() > 0){
    select_vertices();
    remove_vertices();
    spawn_threads();
  }
}

void smooth(Mesh* mesh, size_t niter){
  svd_init( mesh->NNodes );

  double * quality_cache = new double[mesh->NElements];
  bool * vertice_in_cache = new bool[mesh->NElements];

  // For the specified number of iterations, loop over all mesh vertices.
  for(size_t iter=0; iter<niter; ++iter){
    for(size_t vid=0; vid<mesh->NNodes; ++vid){
      // If this is a corner node, it cannot be moved.
      if(mesh->isCornerNode(vid))
        continue;

      // Find the quality of the worst element adjacent to vid
      double worst_q=1.0;

      //parallelize
      for(std::set<size_t>::const_iterator it=mesh->NEList[vid].begin();
          it!=mesh->NEList[vid].end(); ++it){
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

      /* Find the barycentre (centre of mass) of the cavity. A cavity is
       * defined as the set containing vid and all its adjacent vertices and
       * triangles. Since we work on metric space, all lengths have to measured
       * using the metric. The metric tensor is a 2x2 symmetric matrix which
       * encodes the ideal length and orientation of an edge containing vid. As
       * an edge is defined by two vertices, we calculate the edge length using
       * the value of the metric in the middle of the edge, i.e. the average of
       * the two metric tensors of the vertices defining the edge.
       */

      const double * m0 = &mesh->metric[3*vid]; //const

      const double x0 = mesh->coords[2*vid];
      const double y0 = mesh->coords[2*vid+1];

      double A[4] = {0.0, 0.0, 0.0, 0.0}; //const
      double q[2] = {0.0, 0.0};

      // Iterate over all edges and assemble matrices A and q.
      for(std::vector<size_t>::const_iterator it=mesh->NNList[vid].begin();
          it!=mesh->NNList[vid].end(); ++it){
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

      /* The displacement p for vid is found by solving the linear system:
       * ┌─       ─┐   ┌    ┐   ┌    ┐
       * │A[0] A[1]│   │p[0]│   │q[0]│
       * │         │ x │    │ = │    │
       * │A[2] A[3]│   │p[1]│   │q[0]│
       * └─       ─┘   └    ┘   └    ┘
       */
      svd_solve_2x2(vid, A, p, q);

      /* If this is a surface vertex, restrict the displacement
       * to the surface. The new displacement is the projection
       * of the old displacement on the surface.
       */
      if(mesh->isSurfaceNode(vid)){
        p[0] -= p[0]*fabs(mesh->normals[2*vid]);
        p[1] -= p[1]*fabs(mesh->normals[2*vid+1]);
      }

      // Actually change something
      // Update the coordinates
      mesh->coords[2*vid] += p[0];
      mesh->coords[2*vid+1] += p[1];

      /************************************************************************
       * At this point we must also interpolate the metric tensors from all   *
       * neighbouring vertices in order to calculate the new value of vid's   *
       * metric tensor at the new location. This is a quite complex procedure *
       * and has been omitted for simplicity of the exercise. A vertex will   *
       * always use its original metric tensor, no matter whether it has been *
       * relocated or not.                                                    *
       ************************************************************************/

      /* Find the quality of the worst element after smoothing. If an element
       * of the cavity was inverted, i.e. if vid was relocated outside the
       * interior convex hull of the cavity, then the calculated area of that
       * element will be negative and mesh->element_quality() will return a
       * negative number. In such a case, the smoothing operation has to be
       * rejected.
       */
      double new_worst_q=1.0;

      //parallelize
      for(std::set<size_t>::const_iterator it=mesh->NEList[vid].begin();
          it!=mesh->NEList[vid].end(); ++it){

        //store in cache new quality measure
        double v_quality = quality_cache[*it] = mesh->element_quality(*it);
        vertice_in_cache[*it] = true;

        new_worst_q = std::min(new_worst_q, v_quality);
      }

      /* If quality is worse than before, either because of element inversion
       * or just because relocating vid to the barycentre of the cavity does
       * not improve quality, revert the changes.
       */

      //Undo the change
      if(new_worst_q < worst_q){
        mesh->coords[2*vid] -= p[0];
        mesh->coords[2*vid+1] -= p[1];

        for(std::set<size_t>::const_iterator it=mesh->NEList[vid].begin();
          it!=mesh->NEList[vid].end(); ++it){
          vertice_in_cache[*it] = false;
        }
      }
    }
  }

  svd_teardown(mesh->NNodes);
  
  delete quality_cache;
  delete vertice_in_cache;
}
