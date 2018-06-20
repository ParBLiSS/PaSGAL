/**
 * @file    main.cpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */
#include <iostream>

#include "graphLoad.hpp"
#include "align.hpp"
#include "utils.hpp"
#include "base_types.hpp"
#include "clipp.h"

int main(int argc, char **argv)
{
  std::string rfile = "", qfile = "";

  auto cli = (
      clipp::required("-r") & clipp::value(".vg reference graph file", rfile),
      clipp::required("-q") & clipp::value("an input query file (fasta/fastq)[.gz]", qfile)
      );

  if(!clipp::parse(argc, argv, cli)) 
  {
    clipp::print ( clipp::make_man_page(cli, argv[0]) );
    exit(1);
  }

  std::cout << "INFO, psgl::main, reference file = " << rfile << std::endl;
  std::cout << "INFO, psgl::main, query file = " << qfile << std::endl;

  psgl::graphLoader<> g;
  g.loadFromVG(rfile);

  std::cout << "INFO, psgl::main, graph ready in CSR format, n = " << g.diGraph.numVertices << ", m = " << g.diGraph.numEdges << ", len = " << g.diGraph.totalRefLength() << std::endl;
  
  psgl::alignToDAG<int>(qfile, g.diGraph, psgl::MODE::LOCAL);  

#ifdef DEBUG
  g.printGraph();
#endif
}
