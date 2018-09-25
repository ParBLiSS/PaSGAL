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
  class graphIterFwd
  {
    private:

      //vertex id for current position
      int32_t currentVid;

      //sequence offset within vertex id for current position
      std::size_t seqOffset;

      //global sequence offset
      std::size_t globalOffset;

      //underlying graph container
      const CSR_container &graph;

    public:

      graphIterFwd() = delete;

      /**
       * @brief             public constructor
       * @param[in]   g
       */
      graphIterFwd(const CSR_container &g) :
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
      graphIterFwd(const CSR_container &g, int32_t v) :
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
      void getInNeighborOffsets(std::vector<int32_t> &offsets) const
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
       * @brief                   get offsets in DP matrix associated with out-neighbors of current character
       * @param[out]    offsets   
       * @details                 offsets start from 0 to total reference sequence length
       */
      void getOutNeighborOffsets(std::vector<int32_t> &offsets) const
      {
        assert(offsets.size() == 0);

        if (this->seqOffset == graph.vertex_metadata[currentVid].length() - 1 )
        {
          graph.getOutSeqOffsets(currentVid, offsets);
        }
        else
        {
          offsets.push_back( this->globalOffset + 1 );
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
      int32_t getCurrentVertexId() const
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
