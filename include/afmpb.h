#ifndef __AFMPB_H__
#define __AFMPB_H__

#include <iostream> //debug
#include <fstream>
#include <memory>
#include <string>
#include "dashmm/dashmm.h"

namespace afmpb {

struct Configuration; 
struct Atom;  
struct Node;          
struct Element; 

// Usage of AFMPB
void usage(char *program); 

// Initialize AFMPB
std::unique_ptr<Configuration> init(int argc, char **argv); 

// Finalize AFMPB
int finalize(); 

struct Configuration {
  ~Configuration() {
    std::cout << "destruct configuration\n";
  }
  std::string pqr_file;
  std::string mesh_file; 
  std::string log_file{"output.txt"};
  std::string surface_potential_file{"surfp.txt"};
  int mesh_format = 0; 
  double mesh_density = -1.0;
  double probe_radius = -1.0; 
  double dielectric_interior = 2.0; 
  double dielectric_exterior = 80.0;
  double ion_concentration = 150.0; 
  double temperature = 300.0; 
  double surface_tension_coeff = 0.005; 
  double pressure = 0.035; 
  int accuracy = 3; 
}; 

struct Atom {
  dashmm::Point position; 
  double charge;  
  double radius;
}; 

struct Node {
  dashmm::Point position; // Position of the node 
  dashmm::Point normal_i; // Inner normal derivative of the node 
  dashmm::Point normal_o; // Outer normal derivative of the node 
  double area;            // Area of the patch for the node 
  double projected;       // Projected area of the patch for the node 
  double solution[2];     // Solution value 
  double rhs[2];          // Right-hand side value 
  int index;              // Index of the node
}; 

struct Element {
  dashmm::Point normal;   // Normal derivative of the element 
  double area;            // Area of the element 
  int index[3];           // Index of the nodes of the element
}; 

class AFMPB {
public:
  AFMPB(std::unique_ptr<Configuration> param); 
  ~AFMPB() {
    pqr_.close(); 
    log_.close(); 
    surface_potential_.close(); 
    if (mesh_.is_open()) 
      mesh_.close(); 
  }

  void setup() {
    readPQRFile(); 
    if (mesh_.is_open()) {
      readMeshFile(); 
    } else {
      generateMesh();
    }
  }

  double solvation_energy() const {
    return gamma_ * area_ + p_ * volume_ + b_;
  }

private: 
  void processPQRFile(std::string &pqr_file); 
  void readPQRFile(); 
  void readMeshFile(); 
  void generateMesh(); 

private: 
  double dielectric_; 
  double kap_; 
  double cut1_; 
  double cut2_; 
  double sigma_; 

  int natoms_; 
  dashmm::Array<Atom> atoms_; 

  double gamma_;  // nonpolar solvation energy 
  double p_;      // gamma_ * area_ + p_ * volume_ + b_ 
  double b_; 
  double area_;   // surface area
  double volume_; // and volume of the cavity created by the molecule
 
  std::ifstream pqr_; 
  std::ifstream mesh_; 
  std::ofstream log_; 
  std::ofstream surface_potential_; 
}; 




} // namespace afmpb 


/*
struct Atom {
  dashmm::Point position;
  double charge; 
  double radius;
}; 

*/

#endif
