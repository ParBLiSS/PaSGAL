/**
 * @file    csr.hpp
 * @brief   routines to store graph in CSR format 
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef CSR_CONTAINER_HPP
#define CSR_CONTAINER_HPP

#include <cassert>
#include <iostream>

//Own includes
#include "utils.hpp"
#include "stream.hpp"
#include "vg.pb.h"

namespace psgl
{

  /**
   * @brief     class to support storage of directed graph in CSR format
   * @details   CSR is a popular graph storage format-
   *            - vertex numbering starts from 0
   *            - both outgoing and incoming edges are store separately 
   *              - this is redundant but convenient for analysis
   *              - Adjacency (out/in) list of vertex i is stored in array 
   *                'adjcny_' starting at index offsets_[i] and 
   *                ending at (but not including) index offsets_[i+1]
   *
   *            - Procedure to create CSR_container
   *              addVertexCount() : specify # vertices
   *              initVertexSequence() : specify all vertex sequences
   *              initEdges() : specify edges
   *              sort() : topological sort
   *              verify() : verify correctness of graph
   */
  template <typename VertexIdType, typename EdgeIdType>
    class CSR_container
    {
      public:

        //Count of edges and vertices in the graph
        VertexIdType numVertices;
        EdgeIdType numEdges;

        //contiguous adjacency list of all vertices, size = numEdges
        std::vector<VertexIdType> adjcny_in;  
        std::vector<VertexIdType> adjcny_out;  

        //cumulative prefix sequence length till any vertex, size = numVertices
        std::vector<std::size_t> cumulativeSeqLength;

        //offsets in adjacency list for each vertex, size = numVertices + 1
        std::vector<EdgeIdType> offsets_in;
        std::vector<EdgeIdType> offsets_out;

        //Container to hold metadata (e.g., DNA sequence) of all vertices in graph
        std::vector<std::string> vertex_metadata;

        /**
         * @brief     constructor
         */
        CSR_container()
        {
          numVertices = 0;
          numEdges = 0;
        }

        /**
         * @brief     sanity check for correctness of graph storage in CSR format
         */
        void verify() const
        {
          assert(this->numVertices > 0);
          assert(this->numEdges > 0);

          //sequences
          {
            assert(vertex_metadata.size() == this->numVertices);

            for(auto &seq : vertex_metadata)
            {
              //should be non-empty
              assert(seq.length() > 0);

              //all characters should be upper case
              for(auto &c : seq)
              {
                assert(std::isupper(c));
                assert(c == 'A' || c == 'T' || c == 'G' || c == 'C' || c == 'N');
              }
            }
          }

          //adjacency list
          {
            assert(adjcny_in.size() == this->numEdges);
            assert(adjcny_out.size() == this->numEdges);

            for(auto vId : adjcny_in)
              assert(vId >=0 && vId < this->numVertices);

            for(auto vId : adjcny_out)
              assert(vId >=0 && vId < this->numVertices);
          }

          //offset array
          {
            assert(offsets_in.size() == this->numVertices + 1);
            assert(offsets_out.size() == this->numVertices + 1);

            for(auto off : offsets_in)
              assert(off >=0 && off <= this->numEdges); 

            for(auto off : offsets_out)
              assert(off >=0 && off <= this->numEdges); 

            assert(std::is_sorted(offsets_in.begin(), offsets_in.end()));
            assert(std::is_sorted(offsets_out.begin(), offsets_out.end()));

            assert(offsets_in.front() == 0 && offsets_in.back() == this->numEdges);
            assert(offsets_out.front() == 0 && offsets_out.back() == this->numEdges);
          }

          //topologically sorted order
          {
            for(VertexIdType i = 0; i < this->numVertices; i++)
              for(auto j = offsets_out[i]; j < offsets_out[i+1]; j++)
                assert( i < adjcny_out[j] );
          }

          //prefix sequence lengths
          {
            for(VertexIdType i = 1; i < this->numVertices; i++)
              assert( cumulativeSeqLength[i] > cumulativeSeqLength[i-1] );
            assert( cumulativeSeqLength[ this->numVertices - 1 ] == this->totalRefLength() );
          }
        }

        /**
         * @brief           add vertices in graph
         * @param[in]   n   count of vertices to add
         */
        void addVertexCount(VertexIdType n)
        {
          assert(n > 0);  

          this->numVertices += n;
          this->vertex_metadata.resize(this->numVertices);
          this->cumulativeSeqLength.resize(this->numVertices);
        }

        /**
         * @brief             add vertex sequence in the graph
         * @param[in]   id    vertex id
         * @param[in]   seq   sequence
         */
        void initVertexSequence(VertexIdType id, const std::string &seq)
        {
          assert(id >= 0 && id < vertex_metadata.size());  
          assert(vertex_metadata[id].length() == 0);

          vertex_metadata[id] = seq;
        }

        /**
         * @brief                   add edges from a vector
         * @param[in] diEdgeVector  vector containing directed edges as pairs <from, to> 
         */
        void initEdges(std::vector<std::pair<VertexIdType, VertexIdType>> &diEdgeVector)
        {
          for(auto &edge: diEdgeVector)
          {
            assert(edge.first >= 0 && edge.first < this->numVertices);
            assert(edge.second >=0 && edge.second < this->numVertices);
          }

          this->numEdges = diEdgeVector.size();

          //out-edges
          {
            adjcny_out.reserve(this->numEdges);
            offsets_out.resize(this->numVertices + 1);

            //Sort the edge vector before adding it into adjcny
            std::sort(diEdgeVector.begin(), diEdgeVector.end());

            //Compute offsets and adjcny vectors
            offsets_out[0] = 0;
            auto it_b = diEdgeVector.begin();

            for(VertexIdType i = 0; i < this->numVertices; i++)
            {
              //Range for adjacency list of vertex i
              auto it_e = std::find_if(it_b, diEdgeVector.end(), [i](std::pair<VertexIdType, VertexIdType> &e) { return e.first > i; });

              offsets_out[i+1] = std::distance(diEdgeVector.begin(), it_e);

              for(auto it = it_b; it != it_e; it++)
                adjcny_out.push_back(it->second);

              it_b = it_e;   
            }
          }

          //for in-edges  
          {
            adjcny_in.reserve(this->numEdges);
            offsets_in.resize(this->numVertices + 1);

            //reverse the edge vector: <from, to> -> <to, from>
            for(auto &e: diEdgeVector)
              e = std::make_pair(e.second, e.first);

            //Sort the edge vector before adding it into adjcny
            std::sort(diEdgeVector.begin(), diEdgeVector.end());

            //Compute offsets and adjcny vectors
            offsets_in[0] = 0;
            auto it_b = diEdgeVector.begin();

            for(VertexIdType i = 0; i < this->numVertices; i++)
            {
              //Range for adjacency list of vertex i
              auto it_e = std::find_if(it_b, diEdgeVector.end(), [i](std::pair<VertexIdType, VertexIdType> &e) { return e.first > i; });

              offsets_in[i+1] = std::distance(diEdgeVector.begin(), it_e);

              for(auto it = it_b; it != it_e; it++)
                adjcny_in.push_back(it->second);

              it_b = it_e;   
            }
          }

        }

        /**
         * @brief             check existence of edge from vertex u to v
         * @param[in]   u     vertex id
         * @param[in]   v     vertex id
         * @return            true if edge exists, false otherwise
         */
        bool edgeExists(VertexIdType u, VertexIdType v) const
        {
          assert(u >= 0 && u < this->numVertices);
          assert(v >= 0 && v < this->numVertices);

          bool returnVal = false;

          for (auto i = offsets_out[u]; i < offsets_out[u+1]; i++)
          {
            if ( v == adjcny_out[i] )
              returnVal = true;
          }

          return returnVal;
        }

        /**
         * @brief     print the loaded directed graph to stderr 
         * @details   Format details (assuming n = no. of vertices):
         *              Print n+1 rows in total
         *              First row: value of n 
         *              Subsequent rows contain OUT-neighbors of vertices and their DNA sequences,
         *              one row per vertex
         *            This function is implemented for debugging purpose
         */
        void printGraph() const
        {
          std::cerr << "DEBUG, psgl::CSR_container::printGraph, Printing complete graph" << std::endl;
          std::cerr << this->numVertices << std::endl;

          for (VertexIdType i = 0; i < this->numVertices; i++)
          {
            for (EdgeIdType j = offsets_out[i]; j < offsets_out[i+1]; j++)
              std::cerr << adjcny_out[j] << " ";

            std::cerr << vertex_metadata[i] << "\n";
          }

          std::cerr << "DEBUG, psgl::CSR_container::printGraph, Printing done" << std::endl;
        }

        /**
         * @brief     total reference sequence length represented as graph
         * @return    the length
         */
        std::size_t totalRefLength() const
        {
          std::size_t totalLen = 0;

          for (auto &i : vertex_metadata)
          {
            assert(i.length() > 0);

            totalLen += i.length();
          }

          return totalLen;
        }

        /**
         * @brief     total reference sequence length between 
         *            vertices v1 and V2 (both inclusive) in the sorted order
         * @return    the length
         */
        std::size_t totalRefLength (VertexIdType v1, VertexIdType v2) const
        {
          assert(v1 >= 0 && v1 < this->numVertices);
          assert(v2 >= 0 && v2 < this->numVertices);
          assert(v1 <= v2);

          std::size_t totalLen = 0;

          for (auto i = v1; i <= v2; i++)
          {
            assert(vertex_metadata[i].length() > 0);

            totalLen += vertex_metadata[i].length();
          }

          return totalLen;
        }

        /**
         * @brief     relabel graph vertices in the topologically sorted order
         */
        void sort()
        {
          std::vector<VertexIdType> order(this->numVertices);
          const int runs = 5;

          topologicalSort(runs, order); 

          std::cout << "INFO, psgl::CSR_container::sort, topological sort [rand" << runs << "] computed, bandwidth = " << directedBandwidth(order) << std::endl;
          std::cout << "INFO, psgl::CSR_container::sort, a loose lower bound on bandwidth = " << lowerBoundBandwidth() << std::endl;
          std::cout << "INFO, psgl::CSR_container::sort, relabeling graph based on the computed order" << std::endl;

          //Sorted position to vertex mapping (reverse order)
          std::vector<VertexIdType> rOrder(this->numVertices);
          for(VertexIdType i = 0; i < this->numVertices; i++)
            rOrder[ order[i] ] = i;

          //Relabel the graph completely in this order
          {
            //meta-data
            {
              std::vector<std::string> vertex_metadata_new(this->numVertices);

              for (VertexIdType i = 0; i < this->numVertices; i++)
                vertex_metadata_new[i] = vertex_metadata[ rOrder[i] ];

              vertex_metadata = vertex_metadata_new;
            }

            //adjacency lists
            {
              std::vector<VertexIdType> adjcny_in_new;  
              std::vector<VertexIdType> adjcny_out_new;  

              
              for(VertexIdType i = 0; i < this->numVertices; i++)
              {
                std::vector<VertexIdType> tmp;

                for(auto j = offsets_in[ rOrder[i] ]; j < offsets_in[ rOrder[i] + 1 ]; j++)
                  tmp.push_back( order[adjcny_in[j]] );

                //insert adjacency elements in sorted order
                std::sort(tmp.begin(), tmp.end());

                adjcny_in_new.insert(adjcny_in_new.end(), tmp.begin(), tmp.end());
              }

              for(VertexIdType i = 0; i < this->numVertices; i++)
              {
                std::vector<VertexIdType> tmp;

                for(auto j = offsets_out[ rOrder[i] ]; j < offsets_out[ rOrder[i] + 1 ]; j++)
                  tmp.push_back( order[adjcny_out[j]] );

                //insert adjacency elements in sorted order
                std::sort(tmp.begin(), tmp.end());

                adjcny_out_new.insert(adjcny_out_new.end(), tmp.begin(), tmp.end());
              }

              adjcny_in = adjcny_in_new;
              adjcny_out = adjcny_out_new;
            }

            //offsets
            {
              std::vector<EdgeIdType> offsets_in_new (this->numVertices + 1, 0);
              std::vector<EdgeIdType> offsets_out_new (this->numVertices + 1, 0);

              for(VertexIdType i = 0; i < this->numVertices; i++)
                offsets_in_new[i + 1] = offsets_in_new[i] + ( offsets_in[ rOrder[i] + 1] - offsets_in[ rOrder[i] ] ); 

              for(VertexIdType i = 0; i < this->numVertices; i++)
                offsets_out_new[i + 1] = offsets_out_new[i] + ( offsets_out[ rOrder[i] + 1] - offsets_out[ rOrder[i] ] ); 

              offsets_in  = offsets_in_new;
              offsets_out = offsets_out_new;
            }

          }

          //compute prefix sequence length
          this->cumulativeSeqLength[0] = vertex_metadata[0].length();

          for(VertexIdType i = 1; i < this->numVertices; i++)
          {
            cumulativeSeqLength[i] = cumulativeSeqLength[i-1] + vertex_metadata[i].length();
          }
        }

        /**
         * @brief                   compute maximum distance between connected vertices in the graph (a.k.a. 
         *                          directed bandwidth), while noting that each node is a chain of characters
         * @return                  directed graph bandwidth
         * @details                 output of this fn decides the maximum count of prior columns we need during DP 
         */
        std::size_t directedBandwidth() const
        {
          std::size_t bandwidth = 0;   //temporary value 

          std::pair<VertexIdType, VertexIdType> logFarthestVertices;

          //iterate over all vertices in graph to compute bandwidth
          for(VertexIdType i = 0; i < this->numVertices; i++)
          {
            //iterate over neighbors of vertex i
            for(auto j = offsets_out[i]; j < offsets_out[i+1]; j++)
            {
              auto from_pos = i;
              auto to_pos = adjcny_out[j];

              //for a valid topological sort order
              assert(to_pos > from_pos);

              //now the bandwidth between vertex i and its neighbor equals
              //(to_pos - from_pos) plus the width of intermediate vertices

              std::size_t tmp_bandwidth = to_pos - from_pos;

              for(auto k = from_pos + 1; k < to_pos; k++)
                tmp_bandwidth += vertex_metadata[k].length() - 1;

              if(tmp_bandwidth > bandwidth)
              {
                bandwidth = tmp_bandwidth;
                logFarthestVertices = std::make_pair(i, adjcny_out[j]);
              }
            }
          }

#ifdef DEBUG
          std::cerr << "DEBUG, psgl::CSR_container::directedBandwidth, Bandwidth deciding vertices = " << logFarthestVertices.first << ", " << logFarthestVertices.second << std::endl;
#endif

          return bandwidth;
        }

        /**
         * @brief                   get all in-neighbor vertices of a vertex
         * @param[in]   v
         * @param[out]  vec
         */
        void getInNeighbors(VertexIdType v, std::vector<VertexIdType> &vec) const
        {
          assert(vec.size() == 0);
          assert(v >= 0 && v < this->numVertices);

          for(auto i = offsets_in[v]; i < offsets_in[v+1]; i++)
            vec.push_back( adjcny_in[i] );
        }

        /**
         * @brief                   get all out-neighbor vertices of a vertex
         * @param[in]   v
         * @param[out]  vec
         */
        void getOutNeighbors(VertexIdType v, std::vector<VertexIdType> &vec) const
        {
          assert(vec.size() == 0);
          assert(v >= 0 && v < this->numVertices);

          for(auto i = offsets_out[v]; i < offsets_out[v+1]; i++)
            vec.push_back( adjcny_out[i] );
        }

        /**
         * @brief                   get prefix sequence sum offsets for all in-neighbor vertices
         * @param[in]   v
         * @param[out]  vec
         * @details                 useful during DP execution- to access left neighboring reference 
         *                          cells
         */
        void getInSeqOffsets(VertexIdType v, std::vector<std::size_t> &vec) const
        {
          assert(vec.size() == 0);
          assert(v >= 0 && v < this->numVertices);

          for(auto i = offsets_in[v]; i < offsets_in[v+1]; i++)
            vec.push_back( cumulativeSeqLength[adjcny_in[i]] - 1 );
        }

        /**
         * @brief             compute farthest reachable vertex on left side
         *                    while accounting for length of vertex sequences
         * @param[in]   v     id of starting vertex (where alignment ends)
         * @param[in]   dist  alignment length bound
         * @details           useful to determining traceback range in DP
         * TODO               make the bound more precise
         */
        VertexIdType computeLeftMostReachableVertex(VertexIdType v, std::size_t dist) const
        {
          assert(v >= 0 && v < this->numVertices);

          //trivial case
          if (v == 0)
            return 0;

          //Initialize a distance vector for vertices {0, 1,... , v-1, v}
          std::vector<std::size_t> distvec(v+1, dist);

          //Initialize distance of starting vertex
          distvec[v] = 0;

          //Traverse vertices in reverse of topological sorted order
          for (VertexIdType i = v; i != (VertexIdType) -1; i--)
          {
            //Update in-neigbors of vertex i
            for(auto j = offsets_in[i]; j < offsets_in[i+1]; j++)
            {
              assert( adjcny_in[j] < i );

              if (i == v)
                distvec[ adjcny_in[j] ] = 1;      //immediate neighbors of v are just one unit away
              else
                distvec[ adjcny_in[j] ] = std::min( distvec[ adjcny_in[j] ] , distvec[i] + vertex_metadata[i].length() );
            }
          }

          //Check the leftmost vertex inside dist
          for (VertexIdType i = 0; i <= v; i++)
          {
            if (distvec[i] < dist)
              return i;
          }

          return 0;
        }

      private:

        /**
         * @brief                   compute topological sort order using several runs
         *                          of Kahn's algorithm
         *                          ties are decided randomly
         * @param[in]   runs        count of random runs to execute
         * @param[out]  finalOrder  vertex ordering (vertex [0 - n-1] to position [0 - n-1] 
         *                          mapping, i.e. finalOrder[0] denotes where v_0 should go)
         * @details                 ordering with least directed bandwidth is chosen 
         */
        void topologicalSort(int runs, std::vector<VertexIdType> &finalOrder) const
        {
          assert(finalOrder.size() == this->numVertices);
          assert(runs > 0);

          std::vector<VertexIdType> in_degree(this->numVertices, 0);

          //compute in-degree of all vertices in graph
          for (VertexIdType i = 0; i < this->numVertices; i++)
            in_degree[i] = offsets_in[i+1] -  offsets_in[i];

          VertexIdType minbandwidth = std::numeric_limits<VertexIdType>::max();

          for (int i = 0; i < runs; i++)
          {
            std::vector<VertexIdType> tmpOrder(this->numVertices);
            VertexIdType currentOrder = 0;

            //intermediate data structure 
            std::list<VertexIdType> Q;  

            //copy of degree vector
            auto deg = in_degree;

            //push 0 in-degree vertices to Q
            for (VertexIdType i = 0; i < this->numVertices; i++)
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
              for(auto j = offsets_out[v]; j < offsets_out[v+1]; j++)
                if (--deg[ adjcny_out[j] ] == 0)
                  Q.emplace_back (adjcny_out[j]);     //add to Q
            }

            auto currentBandwidth = directedBandwidth(tmpOrder);

#ifdef DEBUG 
            std::cerr << "DEBUG, psgl::CSR_container::topologicalSort, Random run #" << i+1 << " , bandwidth = " << currentBandwidth << std::endl;
#endif

            if (minbandwidth > currentBandwidth)
            {
              minbandwidth = currentBandwidth;
              finalOrder = tmpOrder;
            }
          }
        }

        /**
         * @brief                   compute maximum distance between connected vertices given the new ordering (a.k.a. 
         *                          directed bandwidth), while noting that each node is a chain of characters
         * @param[in]   finalOrder  vertex ordering (vertex [0 - n-1] to position [0 - n-1] 
         *                          mapping, i.e. finalOrder[0] denotes where v_0 should go)
         * @return                  directed graph bandwidth
         * @details                 output of this fn decides the maximum count of prior columns we need during DP 
         */
        std::size_t directedBandwidth(const std::vector<VertexIdType> &finalOrder) const
        {
          assert(finalOrder.size() == this->numVertices);

          std::size_t bandwidth = 0;   //temporary value 

          //Sorted position to vertex mapping (reverse order)
          std::vector<VertexIdType> reverseOrder(this->numVertices);
          for(VertexIdType i = 0; i < this->numVertices; i++)
            reverseOrder[ finalOrder[i] ] = i;

          std::pair<VertexIdType, VertexIdType> logFarthestVertices;
          std::pair<VertexIdType, VertexIdType> logFarthestPositions;

          //iterate over all vertices in graph to compute bandwidth
          for(VertexIdType i = 0; i < this->numVertices; i++)
          {
            //iterate over neighbors of vertex i
            for(auto j = offsets_out[i]; j < offsets_out[i+1]; j++)
            {
              auto from_pos = finalOrder[i];
              auto to_pos = finalOrder[ adjcny_out[j] ];

              //for a valid topological sort order
              assert(to_pos > from_pos);

              //now the bandwidth between vertex i and its neighbor equals
              //(to_pos - from_pos) plus the width of intermediate vertices

              std::size_t tmp_bandwidth = to_pos - from_pos;

              for(auto k = from_pos + 1; k < to_pos; k++)
                tmp_bandwidth += vertex_metadata[reverseOrder[k]].length() - 1;

              if(tmp_bandwidth > bandwidth)
              {
                bandwidth = tmp_bandwidth;
                logFarthestVertices = std::make_pair(i, adjcny_out[j]);
                logFarthestPositions = std::make_pair(from_pos, to_pos);
              }
            }
          }

#ifdef DEBUG
          std::cerr << "DEBUG, psgl::CSR_container::directedBandwidth, Bandwidth deciding positions = " << logFarthestPositions.first << ", " << logFarthestPositions.second << std::endl;
#endif

          return bandwidth;
        }

        /**
         * @brief                   compute a lower bound for bandwidth (loose bound)
         * @details                 lower bound is computed by checking in-neighbors 
         *                          and out-neighbors of each vertex
         * @return                  lower bound for directed bandwidth
         */
        std::size_t lowerBoundBandwidth() const
        {
          //lower bound
          std::size_t lbound = 0;

          //vertex where bound is obtained
          VertexIdType lbound_v;

          //based on out-neighbors
          {
            for (VertexIdType i = 0; i < this->numVertices; i++)
            {
              std::size_t minimum_dist = 1;   //since minimum width required is 1
              std::size_t max_width = 0;

              for (auto j = offsets_out[i]; j < offsets_out[i+1]; j++)
              {
                max_width = std::max(vertex_metadata[ adjcny_out[j] ].length(), max_width);
                minimum_dist += vertex_metadata[ adjcny_out[j] ].length();    //sum up widths
              }

              minimum_dist -= max_width; //subtract neighbor with max width

              if(minimum_dist > lbound)
                lbound_v = i;

              lbound = std::max(minimum_dist, lbound);

            }
          }

          //based on in-neighbors
          {
            for (VertexIdType i = 0; i < this->numVertices; i++)
            {
              std::size_t minimum_dist = 1;   //since minimum width required is 1
              std::size_t max_width = 0;

              for (auto j = offsets_in[i]; j < offsets_in[i+1]; j++)
              {
                max_width = std::max(vertex_metadata[ adjcny_in[j] ].length(), max_width);
                minimum_dist += vertex_metadata[ adjcny_in[j] ].length();    //sum up widths
              }

              minimum_dist -= max_width; //subtract neighbor with max width

              if(minimum_dist > lbound)
                lbound_v = i;

              lbound = std::max(minimum_dist, lbound);
            }
          }

          //based on single insertion variation
          {
            for (VertexIdType i = 0; i < this->numVertices; i++)
            {
              std::size_t minimum_dist = 1;   //since minimum width required is 1

              if (offsets_out[i+1] - offsets_out[i] == 2)
              {
                auto j = offsets_out[i];

                //two neigbors of vertex i
                auto u = adjcny_out[j];
                auto v = adjcny_out[j + 1];

                //edge from u -> v
                if (edgeExists(u,v))
                {
                  minimum_dist += vertex_metadata[u].length();
                }
                else if (edgeExists(v,u))  //edge from v -> u
                {
                  minimum_dist += vertex_metadata[v].length();
                }
              }

              if(minimum_dist > lbound)
                lbound_v = i;

              lbound = std::max(minimum_dist, lbound);
            }
          }

#ifdef DEBUG
          std::cerr << "DEBUG, psgl::CSR_container::lowerBoundBandwidth, lower bound obtained at vertex id = " << lbound_v << std::endl; 
#endif
          return lbound;
        }

    };
}

#endif
