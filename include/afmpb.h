#ifndef __AFMPB_H__
#define __AFMPB_H__

#include <iostream> //debug
#include <fstream>
#include <memory>
#include <string>
#include "dashmm/dashmm.h"

namespace afmpb {

struct Configuration {
  std::string pqr_file; 
  std::string mesh_file; 
  std::string log_file{"output.txt"}; 
  std::string potential_file{"potential.txt"}; 
  int mesh_format = 0; 
  double mesh_density = 40.0; 
  double probe_radius = 0.0; 
  double dielectric_interior = 2.0; 
  double dielectric_exterior = 80.0; 
  double ion_concentration = 150.0; 
  double temperature = 300.0; 
  double surface_tension = 0.005; 
  double pressure = 0.035; 
  int accuracy = 3;
}; 

// Initialize AFMPB 
std::unique_ptr<Configuration> init(int argc, char **argv); 

// Finalize AFMPB 
int finalize(); 

struct Atom {
  dashmm::Point position; 
  double charge;  
  double radius;
}; 

struct Patch {
  dashmm::Point position; 
  dashmm::Point normal;   // Normal direction of the patch 
  double weight;          // Quadrature weight 
}; 

struct Node {
  int index;                // Index of the node
  dashmm::Point position;   // Position of the node 
  dashmm::Point normal_i;   // Inner normal derivative of the node
  dashmm::Point normal_o;   // Outer normal derivative of the node 
  std::vector<Patch> patch; // Node-patch 
  double area;              // Area of the patch for the node 
  double projected;         // Projected area 
  double solution[2];       // Solution value 
  double rhs[2];            // Right-hand side value 
  std::map<int, std::vector<double>> cached; // Cached values for S_to_T
}; 

struct Element {
  int index;                // Index of the element
  dashmm::Point normal;     // Normal direction of the element
  std::vector<int> nodes;   // Index of the nodes of the element
}; 

// Gaussian quadrature points inside each element
struct GNode {
  int index; 
  dashmm::Point position;   
  dashmm::Point normal; 
  double value[2];          
}; 

class AFMPB {
public:
  AFMPB(std::unique_ptr<Configuration> p); 
  ~AFMPB() {
    pqr_.close(); 
    log_.close(); 
    potential_.close(); 
    if (mesh_.is_open())
      mesh_.close();
  }
  
  void setup(); 

  double totalFreeEnergy() const {
    return surface_tension_ * area_ + pressure_ * volume_ + b_; 
  }

private: 
  void processPQRFile(std::string &pqr_file); 
  std::vector<Atom> readAtoms(); 
  void generateMesh(int s, std::vector<Atom> &molecule, 
                    std::vector<Node> &nodes, 
                    std::vector<GNode> &gauss); 
  void readMesh(std::vector<Node> &nodes); 
  double tetrahedronVolume(dashmm::Point &A, dashmm::Point &B, 
                           dashmm::Point &C); 

private: 
  std::ifstream pqr_; 
  std::ifstream mesh_; 
  std::ofstream log_; 
  std::ofstream potential_; 
  int mesh_format_; 
  double mesh_density_; 
  double probe_radius_; 
  double dielectric_; 
  double kap_; 
  double surface_tension_; 
  double pressure_; 
  int accuracy_; 

  // Parameters for node-patch 
  double cut1_; 
  double cut2_; 
  double sigma_; 

  int natoms_; 
  dashmm::Array<Atom> atoms_; 
  std::vector<Element> elements_; 
  dashmm::Array<Node>  nodes_; 
  dashmm::Array<GNode> gauss_; 
  int ngauss_per_element_; 

  double area_; 
  double volume_; 
  double b_ = 0.0; 
}; 

  /*    
class AFMPB {
public:

  void setup() {
    readPQRFile(); 
    if (mesh_.is_open()) {
      readMeshFile(); 
    } else {
      generateMesh();
    }
  }

private: 
  void readPQRFile(); 
  void readMeshFile(); 
  void generateMesh(); 

}; 

  */


} // namespace afmpb 



#endif
