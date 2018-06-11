/**
 * @file    main.cpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */
#include <iostream>

#include "clipp.h"
#include "graphLoad.hpp"
#include "align.hpp"
#include "utils.hpp"
#include "base_types.hpp"

int main(int argc, char **argv)
{
  std::string infile = "";

  auto cli = (
      clipp::required("-i") & clipp::value("input file", infile)
      );

  if(!clipp::parse(argc, argv, cli)) 
  {
    std::cout << clipp::make_man_page(cli, argv[0]);
    exit(1);
  }

  std::cout << "INFO, psgl::main, input file = " << infile << std::endl;

  psgl::graphLoader<uint32_t, uint32_t> g;
  g.loadFromVG(infile);

  std::cout << "INFO, psgl::main, graph ready in CSR format, n = " << g.diGraph.numVertices << ", m = " << g.diGraph.numEdges << std::endl;

  std::vector<std::string> reads = {"GATCACAGGTCTATCACCCTATTAACCACTCACCGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTATCTGTCTTTGATTCCTGCCTCATCCTATTATTTATCGCACCTACGTTCAATATTACAGGCGAACATACTTACTAAAGTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCACTTTCCACACAGACATCATAACAAAAAATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACACATCTCTGCCAAACCCCAAAAACAAAGAACCCTAACACCAGCCTAACCAGATTTCAAATTTTATCTTTTGGCGGTATGCACTTTTAACAGTCACCCCCCAACAACACATTATTTTCCCCTCCCA"};
  
  psgl::alignToDAG<int>(reads, g.diGraph, psgl::MODE::LOCAL);  

#ifdef DEBUG
  g.printGraph();
#endif

}
