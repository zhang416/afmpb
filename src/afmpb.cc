#include <iomanip>
#include <hpx/hpx.h>
#include "afmpb.h"
#include "afmpb_rhs.h"

namespace afmpb {

dashmm::Evaluator<Atom, GNode, dashmm::AFMPBRHS, dashmm::FMM97> interp{}; 

void AFMPB::setup() {
  auto molecule = readAtoms(); 
  int err = atoms_.allocate(natoms_, molecule); 
  assert(err == dashmm::kSuccess); 

  std::vector<Node> nodes; 
  std::vector<GNode> gauss; 

  if (!mesh_format_) {
    for (int i = 0; i < natoms_; ++i) 
      generateMesh(i, molecule, nodes, gauss);
  } else {
    readMesh(nodes);
    removeIsolatedNodes(nodes); 
    processElementGeometry(nodes);  
    generateGaussianPoint(nodes, gauss); 
  }

  nnodes_ = nodes.size(); 
  err = nodes_.allocate(nnodes_); 
  assert(err == dashmm::kSuccess); 
  err = nodes_.put(0, nnodes_, nodes.data()); 
  assert(err == dashmm::kSuccess); 

  ngauss_ = gauss.size(); 
  err = gauss_.allocate(ngauss_); 
  assert(err == dashmm::kSuccess); 
  err = gauss_.put(0, ngauss_, gauss.data()); 
  assert(err == dashmm::kSuccess); 

  log_ << "\n----------------------------------------------------------------\n"
       << "*      Adaptive Fast Multipole Poisson Boltzmann Solver        *\n"
       << "----------------------------------------------------------------\n\n"
       << "Problem parameters:\n"
       << std::setw(50) << std::left << "... n_atoms:" 
       << std::setw(14) << std::right << natoms_ << "\n"
       << std::setw(50) << std::left << "... n_elements:" 
       << std::setw(14) << std::right << elements_.size() << "\n"
       << std::setw(50) << std::left << "... n_nodes:" 
       << std::setw(14) << std::right << nnodes_ << "\n"
       << std::setw(50) << std::left << "... Area: " 
       << std::setw(14) << std::right << std::setprecision(5) 
       << std::scientific << area_ << "\n"
       << std::setw(50) << std::left << "... Volume: "
       << std::setw(14) << std::right << std::setprecision(5) 
       << std::scientific << volume_ << "\n"
       << std::setw(50) << std::left << "... Kap: "
       << std::setw(14) << std::right << std::setprecision(5) 
       << std::scientific << kap_ << "\n";

  dashmm::FMM97<Atom, GNode, dashmm::AFMPBRHS> method{};
  std::vector<double> kparam{}; 
  err = interp.evaluate(atoms_, gauss_, refine_limit_, &method, 
                        accuracy_, &kparam);
  assert(err == dashmm::kSuccess);  
}

void AFMPB::collect() {
  auto gauss = gauss_.collect(); 
  auto nodes = nodes_.collect(); 
  
  if (gauss) {
    std::sort(&gauss[0], &gauss[ngauss_], 
              [] (const GNode &a, const GNode &b) -> bool {
                return (a.index < b.index);
              });
  }

  if (nodes) {
    std::sort(&nodes[0], &nodes[nnodes_], 
              [] (const Node &a, const Node &b) -> bool {
                return (a.index < b.index);
              });
  }

  // Compute total free energy 
  double energy = totalFreeEnergy(gauss.get(), ngauss_, 
                                  nodes.get(), nnodes_); 

  /*
  // Write potentials 
  potential_.precision(5); 
  potential_ << std::scientific; 
  for (int i = 0; i < nnodes_; ++i) {
    const Node &n = nodes[i]; 
    potential_ << n.position.x() << " " 
               << n.position.y() << " "
               << n.position.z() << " " 
               << n.normal_o.x() << " "
               << n.normal_o.y() << " " 
               << n.normal_o.z() << " " 
               << n.value[0]  << " " 
               << n.value[1]  << "\n";
  }

  if (!mesh_format_) {
    for (auto && e : elements_) {
      potential_ << e.nodes[0] << " " << e.nodes[1] << " " << e.nodes[2] << "\n";
    }
  } 
  */
}


void AFMPB::evaluateGaussianPoint() {

}

double AFMPB::totalFreeEnergy(const GNode *gauss, int ngauss, 
                              const Node *nodes, int nnodes) const {
  double b = 0; 

  /*
  if (mesh_format_) {
    for (auto && e : elements_) {
      int i1 = e.nodes[0]; 
      int i2 = e.nodes[1]; 
      int i3 = e.nodes[2]; 
      int index = e.index; 
      double temp = 0; 
      for (int j = 0; j < 7; ++j) {
        double zeta = 1.0 - xi_[j] - eta_[j]; 
        double f = nodes[i1].value[0] * zeta + 
          nodes[i2].value[0] * xi_[j] + nodes[i3].value[0] * eta_[j];
        double h = nodes[i1].value[1] * zeta + 
          nodes[i2].value[1] * xi_[j] + nodes[i3].value[1] * eta_[j]; 
        temp += (gauss[index + j].rhs[0] * h * dielectric_ - 
                 gauss[index + j].rhs[1] * f) * weight_[j]; 
      }
      b += temp * e.area; 
    }
      
    b /= 8 * M_PI;     
  } else {
    // When using built-in mesh, the number of Gaussian quadrature points is the
    // same as the number of nodes of the surface mesh 
    for (int i = 0; i < ngauss; ++i) {
      b += (gauss[i].rhs[0] * nodes[i].value[1] * dielectric_ -
            gauss[i].rhs[1] * nodes[i].value[0] * nodes[i].area) / 2; 
    }
    
    b /= 4 * M_PI / 0.985; 
  }
  */
  return surface_tension_ * area_ + pressure_ * volume_ + b;
}

} // namespace afmpb
