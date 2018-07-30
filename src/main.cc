//=============================================================================
// DAFMPB: DASHMM Accelerated Adaptive Fast Multipole Poisson-Boltzmann Solver 
//
// Portions Copyright (c) 2014, Institute of Computational Mathematics, CAS
// Portions Copyright (c) 2014, Oak Ridge National Laboratory
// Portions Copyright (c) 2017, Trustees of Indiana University,
//
// All rights reserved.
//
// This program is a free software; you can redistribute it and/or modify it
// uner the terms of the GNU General Public License version 3 as published by 
// the Free Software Foundation. 
//=============================================================================

#include "dafmpb.h"

int main(int argc, char **argv) {
  auto parameter = dafmpb::init(argc, argv); 

  if (parameter) {
    dafmpb::DAFMPB system(std::move(parameter)); 

    auto status = system.computePotential(); 

    system.computeEnergy(status); 

    system.finalize(status); 
  } 

  auto err = dafmpb::finalize(); 
  assert(err == 0); 

  return 0;
}

