/**
 * @file    graph_iter.hpp
 * @brief   allows convenient iteration over ref. graph during alignment
 * @author  Chirag Jain <cjain7@gatech.edu>
 */


#ifndef GRAPH_ITERATOR_HPP
#define GRAPH_ITERATOR_HPP

#include "base_types.hpp"

namespace psgl
{
  /**
   * @brief     class to allow convenient iteration over reference characters 
   *            in the graph while executing DP
   */
  template <typename VertexIdType, typename EdgeIdType>
    class graphIterFwd
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

        graphIterFwd() = delete;

        /**
         * @brief             public constructor
         * @param[in]   g
         */
        graphIterFwd(const CSR_container<VertexIdType, EdgeIdType> &g) :
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
        graphIterFwd(const CSR_container<VertexIdType, EdgeIdType> &g, VertexIdType v) :
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

  /**
   * @brief     class to allow convenient iteration in reverse direction 
   *            over reference graph while executing traceback in DP
   */
  template <typename VertexIdType, typename EdgeIdType>
    class graphIterRev
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

        graphIterRev() = delete;

        /**
         * @brief               public constructor
         * @param[in]   g
         * @param[in]   best    best score coordinates, i.e. starting location of iterator
         */
        template <typename ScoreType>
          graphIterRev(const CSR_container<VertexIdType, EdgeIdType> &g, const BestScoreInfo<ScoreType, VertexIdType> &best) :
            graph(g)
        {
          assert(best.vid >= 0 && best.vid < graph.numVertices);
          this->currentVid = best.vid;

          assert(best.vertexSeqOffset >= 0 && best.vertexSeqOffset < graph.vertex_metadata[currentVid].length());
          this->seqOffset = best.vertexSeqOffset;

          if (currentVid == 0)
            assert(best.refColumn == seqOffset);
          else
            assert(best.refColumn == seqOffset + graph.cumulativeSeqLength[currentVid-1]);
          this->globalOffset = best.refColumn;
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
         * @brief     return global character offset in ordered graph
         */
        std::size_t getGlobalOffset() const
        {
          return this->globalOffset;
        }

        /**
         * @brief     jump to a prior offset 
         */
        void jump (std::size_t globalPos)
        {
          assert(globalPos >= 0 && globalPos < this->globalOffset);

          while(this->globalOffset != globalPos)
            this->prev();
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

      private:

        /**
         * @brief     shift back by one character
         */
        void prev()
        {
          if (this->seqOffset > 0)
          {
            this->seqOffset--;
          }
          else
          {
            this->currentVid --;  this->seqOffset = graph.vertex_metadata[currentVid].length() - 1;
          }

          this->globalOffset --;
        }
    };
}

#endif
