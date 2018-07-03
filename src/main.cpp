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
  std::string rfile = "", qfile = "", mode = "";

  auto cli = (
      clipp::required("-m") & clipp::value("mode", mode).doc("reference graph format [vg or txt]"),
      clipp::required("-r") & clipp::value("ref", rfile).doc("reference graph file"),
      clipp::required("-q") & clipp::value("query", qfile).doc("query file (fasta/fastq)[.gz]")
      );

  if(!clipp::parse(argc, argv, cli)) 
  {
    clipp::print ( clipp::make_man_page(cli, argv[0]) );
    exit(1);
  }

  std::cout << "INFO, psgl::main, reference file = " << rfile << " (in " << mode  << " format) " << std::endl;
  std::cout << "INFO, psgl::main, query file = " << qfile << std::endl;

  psgl::graphLoader<> g;

  if (mode.compare("vg") == 0)
    g.loadFromVG(rfile);
  else if(mode.compare("txt") == 0)
    g.loadFromTxt(rfile);
  else
  {
    std::cerr << "Invalid format " << mode << std::endl;
    exit(1);
  }

  std::cout << "INFO, psgl::main, graph ready in CSR format, n = " << g.diGraph.numVertices << ", m = " << g.diGraph.numEdges << ", len = " << g.diGraph.totalRefLength() << std::endl;
  
  psgl::alignToDAG<int>(qfile, g.diGraph, psgl::MODE::LOCAL);  

#ifdef DEBUG
  g.printGraph();
#endif
}
