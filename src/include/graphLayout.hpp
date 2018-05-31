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
#include "utils.hpp"

namespace psgl
{
  
  /**
   * @brief                   compute topological sort order using Kahn's algorithm
   * @param[in]   graph       input graph in CSR format
   * @param[out]  finalOrder  vertex ordering (vertex [0 - n-1] to position [0 - n-1] mapping, i.e.
   *                          finalOrder[0] denotes where v_0 should go)
   */
  template <typename VertexIdType, typename EdgeIdType>  
    void topologicalSort( const CSR_container<VertexIdType, EdgeIdType> &graph, 
                          std::vector<VertexIdType> &finalOrder )
    {
      assert(finalOrder.size() == graph.numVertices);

      //intermediate data structure 
      std::list<VertexIdType> Q;  

      std::vector<VertexIdType> in_degree(graph.numVertices, 0);

      //compute in-degree of all vertices in graph
      for (VertexIdType i = 0; i < graph.numVertices; i++)
        in_degree[i] = graph.offsets_in[i+1] -  graph.offsets_in[i];

      //push 0 in-degree vertices to Q
      for (VertexIdType i = 0; i < graph.numVertices; i++)
        if (in_degree[i] == 0)
          Q.push_back(i);

      VertexIdType currentOrder = 0;

      while(!Q.empty())
      {
        //pick new vertex from Q (FIFO)
        VertexIdType v = Q.front();  Q.pop_front();

        //add to vertex order
        finalOrder[v] = currentOrder++;
        
        //remove out-edges of vertex 'v'
        for(auto j = graph.offsets_out[v]; j < graph.offsets_out[v+1]; j++)
          if (--in_degree[ graph.adjcny_out[j] ] == 0)
            Q.push_back (graph.adjcny_out[j]);     //add to Q
      }
    }

  /**
   * @brief                   compute topological sort order using several random runs
   *                          of Kahn's algorithm
   * @param[in]   graph       input graph in CSR format
   * @param[in]   runs        count of random runs to execute
   * @param[out]  finalOrder  vertex ordering (vertex [0 - n-1] to position [0 - n-1] mapping, i.e.
   *                          finalOrder[0] denotes where v_0 should go)
   * @details                 ordering with least directed bandwidth is chosen 
   */
  template <typename VertexIdType, typename EdgeIdType>  
    void topologicalSort( const CSR_container<VertexIdType, EdgeIdType> &graph, 
                          int runs,
                          std::vector<VertexIdType> &finalOrder )
    {
      assert(finalOrder.size() == graph.numVertices);
      assert(runs > 0);

      std::vector<VertexIdType> in_degree(graph.numVertices, 0);

      //compute in-degree of all vertices in graph
      for (VertexIdType i = 0; i < graph.numVertices; i++)
        in_degree[i] = graph.offsets_in[i+1] -  graph.offsets_in[i];

      VertexIdType minbandwidth = std::numeric_limits<VertexIdType>::max();

      for (int i = 0; i < runs; i++)
      {
        std::vector<VertexIdType> tmpOrder(graph.numVertices);
        VertexIdType currentOrder = 0;

        //intermediate data structure 
        std::list<VertexIdType> Q;  

        //copy of degree vector
        auto deg = in_degree;

        //push 0 in-degree vertices to Q
        for (VertexIdType i = 0; i < graph.numVertices; i++)
          if (deg[i] == 0)
            Q.emplace_back(i);

        while(!Q.empty())
        {
          //pick new vertex from Q
          auto it = random::select(Q.begin(), Q.end());
          VertexIdType v = *it;
          Q.erase(it);

          //add to vertex order
          tmpOrder[v] = currentOrder++;

          //remove out-edges of vertex 'v'
          for(auto j = graph.offsets_out[v]; j < graph.offsets_out[v+1]; j++)
            if (--deg[ graph.adjcny_out[j] ] == 0)
              Q.emplace_back (graph.adjcny_out[j]);     //add to Q
        }

        auto currentBandwidth = directedBandwidth(graph, tmpOrder);

#ifdef DEBUG 
        std::cout << "DEBUG, psgl::topologicalSort, Random run #" << i+1 << " , bandwidth = " << currentBandwidth << std::endl;
#endif

        if (minbandwidth > currentBandwidth)
        {
          minbandwidth = currentBandwidth;
          finalOrder = tmpOrder;
        }
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
    std::size_t directedBandwidth( const CSR_container<VertexIdType, EdgeIdType> &graph, 
                                    const std::vector<VertexIdType> &finalOrder )
    {
      assert(finalOrder.size() == graph.numVertices);

      std::size_t bandwidth = 0;   //temporary value 

      //Sorted position to vertex mapping (reverse order)
      std::vector<VertexIdType> reverseOrder(graph.numVertices);
      for(VertexIdType i = 0; i < graph.numVertices; i++)
        reverseOrder[ finalOrder[i] ] = i;

      std::pair<VertexIdType, VertexIdType> logFarthestVertices;
      std::pair<VertexIdType, VertexIdType> logFarthestPositions;

      //iterate over all vertices in graph to compute bandwidth
      for(VertexIdType i = 0; i < graph.numVertices; i++)
      {
        //iterate over neighbors of vertex i
        for(auto j = graph.offsets_out[i]; j < graph.offsets_out[i+1]; j++)
        {
          auto from_pos = finalOrder[i];
          auto to_pos = finalOrder[ graph.adjcny_out[j] ];

          //for a valid topological sort order
          assert(to_pos > from_pos);

          //now the bandwidth between vertex i and its neighbor equals
          //(to_pos - from_pos) plus the width of intermediate vertices

          std::size_t tmp_bandwidth = to_pos - from_pos;

          for(auto k = from_pos + 1; k < to_pos; k++)
            tmp_bandwidth += graph.vertex_metadata[reverseOrder[k]].length() - 1;

          if(tmp_bandwidth > bandwidth)
          {
            bandwidth = tmp_bandwidth;
            logFarthestVertices = std::make_pair(i, graph.adjcny_out[j]);
            logFarthestPositions = std::make_pair(from_pos, to_pos);
          }
        }
      }

#ifdef DEBUG
      std::cout << "DEBUG, psgl::directedBandwidth, Bandwidth deciding positions = " << logFarthestPositions.first << ", " << logFarthestPositions.second << std::endl;
#endif

#ifdef DEBUG
      printOrderedSequences(graph, finalOrder, logFarthestVertices.first, logFarthestVertices.second);
#endif

      return bandwidth;
    }

  /**
   * @brief                   compute a lower bound for bandwidth (not necessarily tight)
   * @param[in]   graph       input graph in CSR format
   * @details                 lower bound is computed by checking in-neighbors 
   *                          and out-neighbors of each vertex
   * @return                  lower bound for directed bandwidth
   */
  template <typename VertexIdType, typename EdgeIdType>
    std::size_t lowerBoundBandwidth(const CSR_container<VertexIdType, EdgeIdType> &graph)
    {
      std::size_t lbound = 0;

      //based on out-neighbors
      {
        for (VertexIdType i = 0; i < graph.numVertices; i++)
        {
          std::size_t minimum_dist = 1;   //since minimum width required is 1
          std::size_t max_width = 0;

          for (auto j = graph.offsets_out[i]; j < graph.offsets_out[i+1]; j++)
          {
            max_width = std::max(graph.vertex_metadata[ graph.adjcny_out[j] ].length(), max_width);
            minimum_dist += graph.vertex_metadata[ graph.adjcny_out[j] ].length();    //sum up widths
          }

          minimum_dist -= max_width; //subtract neighbor with max width
          lbound = std::max(minimum_dist, lbound);
        }
      }

      //based on in-neighbors
      {
        for (VertexIdType i = 0; i < graph.numVertices; i++)
        {
          std::size_t minimum_dist = 1;   //since minimum width required is 1
          std::size_t max_width = 0;

          for (auto j = graph.offsets_in[i]; j < graph.offsets_in[i+1]; j++)
          {
            max_width = std::max(graph.vertex_metadata[ graph.adjcny_in[j] ].length(), max_width);
            minimum_dist += graph.vertex_metadata[ graph.adjcny_in[j] ].length();    //sum up widths
          }

          minimum_dist -= max_width; //subtract neighbor with max width
          lbound = std::max(minimum_dist, lbound);
        }
      }

      return lbound;
    }

  /**
   * @brief                   print all vertex sequences in their topological sorted order
   * @param[in]   graph       input graph in CSR format
   * @param[in]   finalOrder  vertex ordering (vertex to position [0 - n-1] mapping, i.e.
   *                          finalOrder[0] denotes the position where v_0 should go)
   */
  template <typename VertexIdType, typename EdgeIdType>
    VertexIdType printOrderedSequences( const CSR_container<VertexIdType, EdgeIdType> &graph, 
                                        const std::vector<VertexIdType> &finalOrder )
    {
      //Sorted position to vertex mapping (reverse order)
      std::vector<VertexIdType> reverseOrder(graph.numVertices);
      for(VertexIdType i = 0; i < graph.numVertices; i++)
        reverseOrder[ finalOrder[i] ] = i;

      std::cout << "Complete order\n";

      for(auto &e: reverseOrder)
        std::cout << "[v_id = " << e << "] " << graph.vertex_metadata[e]  <<  " -> \n";

      std::cout << "End" << std::endl;
    }

  /**
   * @brief                   print vertex sequences in their topological sorted order
   *                          b/w two given vertices (inclusive)
   * @param[in]   graph       input graph in CSR format
   * @param[in]   finalOrder  vertex ordering (vertex to position [0 - n-1] mapping, i.e.
   *                          finalOrder[0] denotes the position where v_0 should go)
   * @param[in]   v_from      first vertex desired in the printed output
   * @param[in]   v_to        last vertex desired in the printed output
   * @note                    vertex 'v_from' should precede vertex 'v_to' 
   */
  template <typename VertexIdType, typename EdgeIdType>
    VertexIdType printOrderedSequences( const CSR_container<VertexIdType, EdgeIdType> &graph, 
                                        const std::vector<VertexIdType> &finalOrder,
                                        VertexIdType v_from, VertexIdType v_to )
    {
      auto from_pos = finalOrder[v_from];
      auto to_pos = finalOrder[v_to];

      //sanity check
      assert(to_pos >= 0 && to_pos < graph.numVertices);
      assert(from_pos >= 0 && from_pos < graph.numVertices);
      assert(to_pos > from_pos);

      //Sorted position to vertex mapping (reverse order)
      std::vector<VertexIdType> reverseOrder(graph.numVertices);
      for(VertexIdType i = 0; i < graph.numVertices; i++)
        reverseOrder[ finalOrder[i] ] = i;

      std::cout << "Partial order between vertices " << v_from << " and " << v_to << "\n";

      for(VertexIdType i = from_pos; i <= to_pos; i++)
        std::cout << "[v_id = " << reverseOrder[i] << "] " << graph.vertex_metadata[ reverseOrder[i] ]  <<  " -> \n";

      std::cout << "End" << std::endl;
    }

}

#endif
