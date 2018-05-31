/**
 * @file    main.cpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */
#include <iostream>

#include "clipp.h"
#include "graphLoad.hpp"
#include "graphLayout.hpp"

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

  psgl::graphLoader<uint32_t, uint32_t> g;
  g.loadFromVG(infile);

#ifdef DEBUG
  g.printGraph();
#endif

  std::vector<uint32_t> order(g.diGraph.numVertices);
  psgl::topologicalSort(g.diGraph, 1000, order); 

  std::cout << "directed bandwidth distance = " << psgl::directedBandwidth(g.diGraph, order) << std::endl;
}