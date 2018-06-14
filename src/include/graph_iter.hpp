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
         * @brief             public constructor
         * @param[in]   g
         */
        graphIter(const CSR_container<VertexIdType, EdgeIdType> &g) :
          graph(g)
        {
          this->currentVid = 0;
          this->seqOffset = 0;
          this->globalOffset = 0;
        }

       /**
         * @brief             public constructor
         *                    begin iterator from a given vertex id
         * @param[in]   g
         * @param[in]   v
         */
        graphIter(const CSR_container<VertexIdType, EdgeIdType> &g, VertexIdType v) :
          graph(g)
        {
          assert(v >= 0 && v < graph.numVertices);

          this->currentVid = v;
          this->seqOffset = 0;

          if(v == 0)
            this->globalOffset = 0;
          else
            this->globalOffset = graph.cumulativeSeqLength[v-1];
        }

        /**
         * @brief     get current character
         */
        char curChar() const
        {
          assert(this->currentVid >= 0 && this->currentVid < graph.numVertices);
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
            graph.getInSeqOffsets(currentVid, offsets);
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
        std::size_t getGlobalOffset() const
        {
          return this->globalOffset;
        }

        /**
         * @brief     return current vertex id
         */
        VertexIdType getCurrentVertexId() const
        {
          return this->currentVid;
        }

        /**
         * @brief     return current sequence offset
         */
        std::size_t getCurrentSeqOffset()
        {
          return this->seqOffset;
        }
    };
}

#endif
