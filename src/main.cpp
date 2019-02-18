/**
 * @file    main.cpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */
#include <iostream>


#ifdef VTUNE_SUPPORT
#include <ittnotify.h>
#endif

#include "parseCmdArgs.hpp"
#include "graphLoad.hpp"
#include "align.hpp"
#include "utils.hpp"
#include "base_types.hpp"

int main(int argc, char **argv)
{
#ifdef VTUNE_SUPPORT
  __itt_pause();
#endif

  //parse command line arguments   
  psgl::Parameters parameters;        
  psgl::parseandSave(argc, argv, parameters);   

  //buffer for results
  std::vector< psgl::BestScoreInfo > bestScoreVector;

  //execute alignment
  if (psgl::alignToDAG (parameters, psgl::MODE::LOCAL, bestScoreVector) == PSGL_STATUS_OK)
    std::cout << "INFO, psgl::main, run finished" << std::endl;

  return 0;
}
