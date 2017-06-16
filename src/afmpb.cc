#include <algorithm>
#include <cstdio>
#include <stdexcept>
#include <getopt.h>
#include <hpx/hpx.h>
#include "afmpb.h"

namespace afmpb {

void usage(char *program) {
  fprintf(stdout, "usage: %s\n"
          "  (--pqr-file=FILE)\n"
          "  (--mesh-format=0 --mesh-density=num --probe-radius=num | \n"
          "   --mesh-format=[1|2] --mesh-file=FILE)\n"
          "  [--dielectric-interior=num     [default:        2.0]]\n"
          "  [--dielectric-exterior=num     [default:       80.0]]\n"
          "  [--ion-concentration=num       [default:      150.0]]\n"
          "  [--temperature=num             [default:      300.0]]\n"
          "  [--surface-tension-coeff=num   [default:      0.005]]\n"
          "  [--pressure=num                [default:      0.035]]\n"
          "  [--log-file=FILE               [default: output.txt]]\n"
          "  [--surface-potential-file=FILE [default:  surfp.txt]]\n"
          "  [--accuracy=num                [default:          3]]\n"
          , program);
}

std::unique_ptr<Configuration> init(int argc, char **argv) {
  if (HPX_SUCCESS != hpx_init(&argc, &argv)) 
    return std::unique_ptr<Configuration>{nullptr}; 

  // Parse command line arguments 
  static struct option long_options[] = {
    {"dielectric-interior", required_argument, 0, 'i'}, 
    {"dielectric-exterior", required_argument, 0, 'e'}, 
    {"ion-concentration", required_argument, 0, 'c'}, 
    {"temperature", required_argument, 0, 't'}, 
    {"surface-tension-coeff", required_argument, 0, 'g'}, 
    {"pressure", required_argument, 0, 'p'}, 
    {"mesh-format", required_argument, 0, 'f'}, 
    {"mesh-density", required_argument, 0, 'd'}, 
    {"probe-radius", required_argument, 0, 'r'}, 
    {"pqr-file", required_argument, 0, 'q'}, 
    {"mesh-file", required_argument, 0, 'm'}, 
    {"log-file", required_argument, 0, 'l'}, 
    {"surface-potential-file", required_argument, 0, 's'}, 
    {"accuracy", required_argument, 0, 'a'}, 
    {"help", no_argument, 0, 'h'}, 
    {0, 0, 0, 0}
  };

  Configuration *p = new Configuration; 

  int opt = 0; 
  int long_index = 0; 
  while ((opt = getopt_long(argc, argv, 
                            "i:e:c:t:g:p:f:d:r:q:m:l:s:a:h", 
                            long_options, &long_index)) != -1) {
    switch (opt) {
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
      p->surface_tension_coeff = atof(optarg); 
      break;
    case 'p':
      p->pressure = atof(optarg); 
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
      p->surface_potential_file = optarg;
      break;
    case 'a':
      p->accuracy = atoi(optarg); 
      break; 
    case 'h':
    case '?':
      usage(argv[0]); 
      delete p;
      return std::unique_ptr<Configuration>{nullptr}; 
    }
  }

  // Check if the command line arguments are valid 
  if (// missing pqr file 
      p->pqr_file.empty() || 
      // missing mesh file when not using built-in meshing
      (p->mesh_format &&     
       p->mesh_file.empty()) ||
      // missing mesh density or probe radius when using built-in meshing
      (!p->mesh_format && 
       (p->mesh_density < 0 || p->probe_radius < 0))
      ) {
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

AFMPB::AFMPB(std::unique_ptr<Configuration> param) {
  dielectric_ = param->dielectric_exterior / param->dielectric_interior; 
  kap_ = sqrt(2.528639884 * std::max(param->ion_concentration, 1e-10) / 
              dielectric_ / param->temperature); 

  // Parameters for node-patch
  if (param->mesh_format == 0) {
    cut1_ = 0.1; 
    cut2_ = 0.0001;
  } else if (param->mesh_format == 1) {
    cut1_ = 0.4; 
    cut2_ = 0.4;
  } else {
    cut1_ = 0.0; 
    cut2_ = 0.3;
  }
  sigma_ = 0.001; 

  gamma_ = param->surface_tension_coeff; 
  p_ = param->pressure; 

  // Open input and output streams
  if (!param->mesh_format) 
    mesh_.open(param->mesh_file); 
  log_.open(param->log_file); 
  surface_potential_.open(param->surface_potential_file); 

  // Process the input PQR file 
  processPQRFile(param->pqr_file); 
  std::cout << "exit afmpb constructor\n";
}

void AFMPB::processPQRFile(std::string &pqr_file) {
  char buffer[200]; 
  sprintf(buffer, "grep \'ATOM\\|HETATM\' %s | cut -c 30- > %s; wc -l < %s", 
          pqr_file.c_str(), "processed.txt", "processed.txt"); 
  std::shared_ptr<FILE> pipe(popen(buffer, "r"), pclose); 
  if (!pipe) throw std::runtime_error("popen() failed!"); 
  
  if (fgets(buffer, sizeof(buffer), pipe.get()) != nullptr) 
    sscanf(buffer, "%d", &natoms_); 

  pqr_.open("processed.txt"); 
  std::cout << natoms_ << "\n";
}

void AFMPB::readPQRFile() {
  Atom *molecule = new Atom[natoms_]; 
  
  int i = 0; 
  double x, y, z, q, r; 
  while (pqr_ >> x >> y >> z >> q >> r) {
    molecule[i].position = dashmm::Point{x, y, z}; 
    molecule[i].charge = q; 
    molecule[i].radius = r; 
    i++; 
  }

  int err = atoms_.allocate(natoms_); 
  assert(err == dashmm::kSuccess); 
  err = atoms_.put(0, natoms_, molecule); 
  assert(err == dashmm::kSuccess); 
  delete [] molecule; 
}

void AFMPB::readMeshFile() {

} 

void AFMPB::generateMesh() {

}


} // namespace afmpb
