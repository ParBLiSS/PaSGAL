/**
 * @file    graphLoad.hpp
 * @brief   routines to load DNA sequence graph from various file formats 
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef GRAPH_LOADER_HPP
#define GRAPH_LOADER_HPP

#include <cassert>
#include <iostream>

//Own includes
#include "csr.hpp"
#include "stream.hpp"
#include "vg.pb.h"

namespace psgl
{

  /**
   * @brief                       supports loading of sequence graphs from VG file format
   * @tparam[in]    VertexIdType  data type to store vertex ids (tbd by total #vertices)
   * @tparam[in]    EdgeIdType    data type to store offsets in CSR graph (tbd by total #edges)
   */
  template <typename VertexIdType, typename EdgeIdType>
    class graphLoader
    {
      public:

        //initialize an empty di-graph 
        CSR_container<VertexIdType, EdgeIdType> diGraph;

        /**
         * @brief     load graph from VG graph format
         * @details   VG tool (https://github.com/vgteam/vg) uses .vg format to save graphs
         */
        void loadFromVG(const std::string &filename)
        {
          //Read vertices in the graph 
          {
            std::ifstream graphFile {filename, std::ios::in | std::ios::binary};

            std::function<void(vg::Graph&)> lambda_v = [this](vg::Graph& g) {

              //Get total count of vertices
              diGraph.addVertexCount(g.node_size());

              for (int i = 0; i < g.node_size(); i++)
              {
                auto vg_vertex = g.node(i);

                //add vertex to diGraph
                //vertex numbering in vg starts from 1, so adjust accordingly
                diGraph.initVertexSequence(vg_vertex.id() - 1, vg_vertex.sequence());
              }
            };

            stream::for_each(graphFile, lambda_v);
          }

          //Read edges in the graph 
          {
            std::ifstream graphFile {filename, std::ios::in | std::ios::binary};

            std::vector <std::pair <VertexIdType, VertexIdType> > edgeVector;

            std::function<void(vg::Graph&)> lambda_e = [&edgeVector](vg::Graph& g) {
              for (int i = 0; i < g.edge_size(); i++)
              {
                auto vg_edge = g.edge(i);

                //todo: add support for bi-directed graphs
                assert(("Bi-directed graph not supported yet", vg_edge.from_start() == false));
                assert(("Bi-directed graph not supported yet", vg_edge.to_end() == false));
                assert(("Graph overlaps not supported yet", vg_edge.overlap() == 0));

                //add edge to diGraph
                //vertex numbering in vg starts from 1, so adjust accordingly
                edgeVector.emplace_back(vg_edge.from() - 1, vg_edge.to() - 1);
              }
            };

            stream::for_each(graphFile, lambda_e);

            diGraph.initEdges(edgeVector);
          }

#ifndef NDEBUG
          //verify correctness of CSR container
          diGraph.verify();
#endif
        }

        /**
         * @brief     print the loaded directed graph to stderr 
         */
        void printGraph()
        {
          diGraph.printGraph();
        }
    };

}

#endif
