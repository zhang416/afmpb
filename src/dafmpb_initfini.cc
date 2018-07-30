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

#include <cstdio>
#include <getopt.h>
#include "dafmpb.h"

namespace dafmpb {

void usage(char *program) {
  fprintf(stdout, "usage: %s --pqr-file=FILE [options] \n"
          "  Mesh options:\n"
          "  --mesh-format=[0|1|2|3]\n"
          "    0: built-in mesh, 1: MSMS, 2: TMSMesh\n"
          "  --mesh-file=FILE\n"
          "    required if --mesh-format=[1|2|3]\n"
          "  --mesh-density=num\n"
          "    optional if built-in mesh is used, default value is 40.0\n"
          "  --probe-radius=num\n"
          "    optional if built-in mesh is used, default value is 0.0\n\n"
          "  Molecular system options:\n"
          "  --dilectric-interior=num\n"
          "    interior dilectric constant, default value is 2.0\n"
          "  --dilectric-exterior=num\n"
          "    exterior dilectric constant, default value is 80.0\n"
          "  --ion-concentration=num\n"
          "    ionic concentration, default value is 150.0\n"
          "  --temperature=num\n"
          "    temperature, default value is 300.0\n"
          "  --surface-tension=num\n"
          "    surface tension coefficient, default value is 0.005\n"
          "  --pressure=num\n"
          "    pressure, default value is 0.035\n\n"
          "  Solver options:\n"
          "  --accuracy=[3|6]\n"
          "    accuracy of the multipole method, default value is 3-digits\n"
          "  --rel-tolerance=num\n"
          "    GMRES relative tolerance, default value is 1e-2 \n"
          "  --abs-tolerance=num\n"
          "    GMRES absolute tolerance, default value is 0\n" 
          "  --restart=num\n"
          "    maxium dimension of the Krylov space, default value is 50\n"
          "  --max-restart=num\n"
          "    maximum number of times GMRES can restart, default value is 1\n\n"
          "  Output options:\n"
          "  --log-file=FILE\n"
          "    name of the log file, default value is output.txt\n"
          "  --potential-file=FILE\n"
          "    name of the potential file, no output written if unspecified\n"
          , program); 
}

std::unique_ptr<Configuration> init(int argc, char **argv) {

  if (HPX_SUCCESS != hpx_init(&argc, &argv)) 
    return std::unique_ptr<Configuration>{nullptr}; 

  // Parse command line arguments 
  static struct option long_options[] = {
    {"pqr-file", required_argument, 0, 'q'}, 
    {"mesh-file", required_argument, 0, 'm'}, 
    {"log-file", required_argument, 0, 'l'}, 
    {"potential-file", required_argument, 0, 's'}, 
    {"mesh-format", required_argument, 0, 'f'}, 
    {"mesh-density", required_argument, 0, 'd'}, 
    {"probe-radius", required_argument, 0, 'r'}, 
    {"dielectric-interior", required_argument, 0, 'i'}, 
    {"dielectric-exterior", required_argument, 0, 'e'}, 
    {"ion-concentration", required_argument, 0, 'c'}, 
    {"temperature", required_argument, 0, 't'}, 
    {"surface-tension", required_argument, 0, 'g'}, 
    {"pressure", required_argument, 0, 'p'}, 
    {"accuracy", required_argument, 0, 'a'}, 
    {"rel-tolerance", required_argument, 0, 'y'}, 
    {"abs-tolerance", required_argument, 0, 'z'}, 
    {"restart", required_argument, 0, 'k'}, 
    {"max-restart", required_argument, 0, 'n'}, 
    {"help", no_argument, 0, 'h'}, 
    {0, 0, 0, 0}
  };

  Configuration *p = new Configuration; 

  int opt = 0; 
  int long_index = 0; 
  while ((opt = getopt_long(argc, argv, 
                            "q:m:l:s:f:d:r:i:e:c:t:g:p:a:k:n:y:z:h", 
                            long_options, &long_index)) != -1) {
    switch (opt) {
    case 'q':
      p->pqr_file = optarg; 
      break;
    case 'm':
      p->mesh_file = optarg;
      break;
    case 'l':
      p->log_file = optarg;
      break;
    case 's':
      p->potential_file = optarg;
      break;
    case 'f':
      p->mesh_format = atoi(optarg); 
      break;
    case 'd':
      p->mesh_density = atof(optarg); 
      break;
    case 'r':
      p->probe_radius = atof(optarg); 
      break;
    case 'i':
      p->dielectric_interior = atof(optarg); 
      break;
    case 'e':
      p->dielectric_exterior = atof(optarg); 
      break;
    case 'c':
      p->ion_concentration = atof(optarg); 
      break;
    case 't':
      p->temperature = atof(optarg); 
      break;
    case 'g':
      p->surface_tension = atof(optarg); 
      break;
    case 'p':
      p->pressure = atof(optarg); 
      break;
    case 'a':
      p->accuracy = atoi(optarg); 
      break; 
    case 'k':
      p->restart = atoi(optarg); 
      break;
    case 'n':
      p->max_restart = atoi(optarg); 
      break;
    case 'y':
      p->rel_tolerance = atof(optarg);
      break;
    case 'z':
      p->abs_tolerance = atof(optarg);
      break; 
    case 'h':
    case '?':
      usage(argv[0]); 
      delete p;
      return std::unique_ptr<Configuration>{nullptr}; 
    }
  }

  // Check if the command line arguments are valid 
  if (p->pqr_file.empty() || 
      (p->mesh_format && p->mesh_file.empty()) ||
      !(p->accuracy == 3 || p->accuracy == 6)) {
    usage(argv[0]); 
    delete p; 
    p = nullptr; 
  } 

  return std::unique_ptr<Configuration>{p}; 
}
 
int finalize() {
  hpx_finalize(); 
  return 0;
}

DAFMPB::DAFMPB(std::unique_ptr<Configuration> p) {
  // Process the input stream on rank 0 only 
  if (hpx_get_my_rank() == 0) {
    pqr_.open(p->pqr_file); 

    if (p->mesh_format) 
      mesh_.open(p->mesh_file); 

    log_.open(p->log_file); 

    if (p->potential_file.empty() == false) 
      potential_.open(p->potential_file); 
  }

  mesh_format_ = p->mesh_format; 
  mesh_density_ = p->mesh_density; 
  probe_radius_ = p->probe_radius; 
  dielectric_exterior_ = p->dielectric_exterior; 
  dielectric_interior_ = p->dielectric_interior; 
  dielectric_ = p->dielectric_exterior / p->dielectric_interior; 
  kap_ = sqrt(2.528639884 * std::max(p->ion_concentration, 1e-10) / 
              p->dielectric_exterior / p->temperature); 
  surface_tension_ = p->surface_tension; 
  pressure_ = p->pressure; 
  accuracy_ = p->accuracy; 
  restart_ = p->restart; 
  max_restart_ = p->max_restart; 
  rel_tolerance_ = p->rel_tolerance; 
  abs_tolerance_ = p->abs_tolerance; 
  refine_limit_ = 50; 

  if (mesh_format_ == 0) {
    cut1_ = 0.1; 
    cut2_ = 0.0001; 
  } else if (mesh_format_ == 1) {
    cut1_ = 0.4; 
    cut2_ = 0.4;
  } else {
    cut1_ = 0.0; 
    cut2_ = 0.3;
  } 
  sigma_ = 0.001; 

  if (mesh_format_) {
    xi_ = new double[7]{0.101286507323456, 0.797426958353087, 
                        0.101286507323456, 0.470142064105115, 
                        0.059715871789770, 0.470142064105115, 1.0/3.0};
    eta_ = new double[7]{0.101286507323456, 0.101286507323456, 
                         0.797426958353087, 0.470142064105115, 
                         0.470142064105115, 0.059715871789770, 1.0/3.0};
    weight_ = new double[7]{0.125939180544827, 0.125939180544827, 
                            0.125939180544827, 0.132394152788506, 
                            0.132394152788506, 0.132394152788506, 0.225};
  }

  hess_.resize(restart_ * (restart_ + 1) / 2 + 1); 
  cosine_.resize(restart_); 
  sine_.resize(restart_); 
  residual_.resize(restart_ + 2); 

  hpx_run(&allocate_reducer_, &reducer_); 
  t_dag_ = 0.0; 
  t_exec_ = 0.0; 
  t_gmres_ = 0.0; 
  t_inner_ = 0.0; 
  n_inner_ = 0; 

  natoms_ = 0; 
  nnodes_ = 0; 
  ngauss_ = 0; 
  area_ = 0.0; 
  volume_ = 0.0; 

  setup(); 
}

} // namespace dafmpb
