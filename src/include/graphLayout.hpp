/**
 * @file    graphLayout.hpp
 * @brief   routines to compute ordering of vertices/edges in graph for better performance
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef GRAPH_LAYOUT_HPP
#define GRAPH_LAYOUT_HPP

#include <cassert>
#include <iostream>
#include <queue>

//Own includes
#include "csr.hpp"

namespace psgl
{
  /**
   * @brief                   compute topological sort order using Kahn's algorithm
   * @param[in]   graph       input graph in CSR format
   * @param[out]  finalOrder  vertex ordering (vertex [0 - n-1] to position [0 - n-1] mapping, i.e.
   *                          finalOrder[0] denotes where v_0 should go)
   */
  template <typename VertexIdType, typename EdgeIdType>  
    void topologicalSort(const CSR_container<VertexIdType, EdgeIdType> &graph, std::vector<VertexIdType> &finalOrder)
    {
      assert(finalOrder.size() == graph.numVertices);

      //intermediate data structure 
      std::queue<VertexIdType> Q;  

      std::vector<VertexIdType> in_degree(graph.numVertices, 0);

      //compute in-degree of all vertices in graph
      for (VertexIdType i = 0; i < graph.numVertices; i++)
        for (auto j = graph.offsets[i]; j < graph.offsets[i+1]; j++)
          in_degree[ graph.adjcny[j] ] ++;

      //push 0 in-degree vertices to Q
      for (VertexIdType i = 0; i < graph.numVertices; i++)
        if (in_degree[i] == 0)
          Q.push(i);

      VertexIdType currentOrder = 0;

      while(!Q.empty())
      {
        //pick new vertex from Q
        VertexIdType v = Q.front();  Q.pop();

        //add to vertex order
        finalOrder[v] = currentOrder++;
        
        //remove out-edges of vertex 'v'
        for(auto j = graph.offsets[v]; j < graph.offsets[v+1]; j++)
          if (--in_degree[ graph.adjcny[j] ] == 0)
            Q.push (graph.adjcny[j]);     //add to Q
      }
    }

  /**
   * @brief                   compute maximum distance between connected vertices in the order (a.k.a. 
   *                          directed bandwidth), while noting that each node is a chain of characters
   * @param[in]   graph       input graph in CSR format
   * @param[in]   finalOrder  vertex ordering (vertex to position [0 - n-1] mapping, i.e.
   *                          finalOrder[0] denotes the position where v_0 should go)
   * @return                  directed graph bandwidth
   * @details                 output of this fn decides the maximum count of prior columns we need during DP 
   */
  template <typename VertexIdType, typename EdgeIdType>
    VertexIdType directedBandwidth(const CSR_container<VertexIdType, EdgeIdType> &graph, const std::vector<VertexIdType> &finalOrder)
    {
      assert(finalOrder.size() == graph.numVertices);

      VertexIdType bandwidth = 0;   //temporary value 

      //Sorted position to vertex mapping (reverse order)
      std::vector<VertexIdType> reverseOrder(graph.numVertices);
      for(VertexIdType i = 0; i < graph.numVertices; i++)
        reverseOrder[ finalOrder[i] ] = i;

      //iterate over all vertices in graph to compute bandwidth
      for(VertexIdType i = 0; i < graph.numVertices; i++)
      {
        //iterate over neighbors of vertex i
        for(auto j = graph.offsets[i]; j < graph.offsets[i+1]; j++)
        {
          auto from_pos = finalOrder[i];
          auto to_pos = finalOrder[ graph.adjcny[j] ];

          assert(to_pos > from_pos);

          //now the bandwidth between vertex i and its neighbor equals
          //(to_pos - from_pos) plus the width of intermediate vertices

          auto tmp_bandwidth = to_pos - from_pos;

          for(auto k = from_pos + 1; k < to_pos; k++)
            tmp_bandwidth += graph.vertex_metadata[reverseOrder[k]].length() - 1;

          if(tmp_bandwidth > bandwidth)
            bandwidth = tmp_bandwidth;
        }
      }

      return bandwidth;
    }
}

#endif
