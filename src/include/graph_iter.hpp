/**
 * @file    graph_iter.hpp
 * @brief   allows convenient iteration over ref. graph during alignment
 * @author  Chirag Jain <cjain7@gatech.edu>
 */


#ifndef GRAPH_ITERATOR_HPP
#define GRAPH_ITERATOR_HPP

namespace psgl
{
  template <typename VertexIdType, typename EdgeIdType>
    class graphIter
    {
      private:

        //vertex id for current position
        VertexIdType currentVid;

        //sequence offset within vertex id for current position
        std::size_t seqOffset;

        //global sequence offset
        std::size_t globalOffset;

        //underlying graph container
        const CSR_container<VertexIdType, EdgeIdType> &graph;

      public:

        /**
         * @brief     public constructor
         */
        graphIter(const CSR_container<VertexIdType, EdgeIdType> &g) :
          graph(g)
        {
          assert(graph.numVertices > 0);

          this->currentVid = 0;
          this->seqOffset = 0;
          this->globalOffset = 0;
        }

        /**
         * @brief     get current character
         */
        char curChar() const
        {
          assert(this->seqOffset >= 0);
          assert(this->seqOffset < graph.vertex_metadata[currentVid].length());

          return graph.vertex_metadata[currentVid].at(seqOffset);
        }

        /**
         * @brief     shift ahead by one character
         */
        void next()
        {
          if (this->seqOffset < graph.vertex_metadata[currentVid].length() - 1)
          {
            this->seqOffset++;
          }
          else
          {
            this->currentVid ++;  this->seqOffset = 0;
          }

          this->globalOffset ++;
        }

        /**
         * @brief                   get offsets in DP matrix associated with in-neighbors of current character
         * @param[out]    offsets   
         * @details                 offsets start from 0 to total reference sequence length
         */
        void getNeighborOffsets(std::vector<std::size_t> &offsets) const
        {
          assert(offsets.size() == 0);

          if(this->seqOffset > 0)
          {
            offsets.push_back( this->globalOffset - 1 );
          }
          else
          {
            std::vector<VertexIdType> inVertices;
            graph.getInNeighbors(currentVid, inVertices);

            for(auto i = inVertices.begin(); i != inVertices.end(); i++)
            {
              std::size_t off = 1;
              
              for(auto j = *i + 1; j < currentVid; j++)
                off += graph.vertex_metadata[j].length();

              offsets.push_back( this->globalOffset - off);
            }
          }
        }

        /**
         * @brief     check if iterator reached end
         */
        bool end() const
        {
          return (currentVid == graph.numVertices);
        }

        /**
         * @brief     return global character offset in ordered graph
         */
        std::size_t getGlobalOffset()
        {
          return this->globalOffset;
        }
    };
}

#endif
