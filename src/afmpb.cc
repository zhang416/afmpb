#include <iomanip>
#include <hpx/hpx.h>
#include "afmpb.h"
#include "afmpb_rhs.h"

namespace afmpb {

dashmm::Evaluator<Atom, GNode, dashmm::AFMPBRHS, dashmm::FMM97> interp{};

void AFMPB::setup() {
  //std::vector<Atom> molecule; 
  Atom *molecule{nullptr}; 
  std::vector<Node> nodes;
  //std::vector<GNode> gauss;

  if (hpx_get_my_rank() == 0) {
    // Input file is read from rank 0 only
    //readAtoms(molecule); 
    molecule = readAtoms(); 

    if (!mesh_format_) {
      for (int i = 0; i < natoms_; ++i)
        //generateMesh(i, molecule, nodes, gauss);
        generateMesh(i, molecule, nodes); 
    } else {
      readMesh(nodes);
      removeIsolatedNodes(nodes);
      processElementGeometry(nodes);
      //generateGaussianPoint(nodes, gauss);
    }
  
    log_
      << "\n----------------------------------------------------------------\n"
      << "*      Adaptive Fast Multipole Poisson Boltzmann Solver        *\n"
      << "----------------------------------------------------------------\n\n"
      << "Problem parameters:\n"
      << std::setw(50) << std::left << "... n_atoms:"
      << std::setw(14) << std::right << natoms_ << "\n" << std::flush;

    if (mesh_format_) {
      log_ << std::setw(50) << std::left << "... n_elements:"
           << std::setw(14) << std::right << elements_.size() 
           << "\n" << std::flush;
    }

    log_ << std::setw(50) << std::left << "... n_nodes:"
         << std::setw(14) << std::right << nodes.size() << "\n"
         << std::setw(50) << std::left << "... Area: "
         << std::setw(14) << std::right << std::setprecision(5)
         << std::scientific << area_ << "\n"
         << std::setw(50) << std::left << "... Volume: "
         << std::setw(14) << std::right << std::setprecision(5)
         << std::scientific << volume_ << "\n"
         << std::setw(50) << std::left << "... Kap: "
         << std::setw(14) << std::right << std::setprecision(5)
         << std::scientific << kap_ << "\n" << std::flush;
  } 

  //natoms_ = molecule.size(); 
  nnodes_ = nodes.size();
  //ngauss_ = gauss.size();

  auto err = atoms_.allocate(natoms_, molecule); 
  //  auto err = atoms_.allocate(natoms_); 
  //assert(err == dashmm::kSuccess); 
  //err = atoms_.put(0, natoms_, molecule.data()); 
  //assert(err == dashmm::kSuccess); 

  err = nodes_.allocate(nnodes_);
  assert(err == dashmm::kSuccess);
  err = nodes_.put(0, nnodes_, nodes.data());
  assert(err == dashmm::kSuccess);

  /*
  err = gauss_.allocate(ngauss_);
  assert(err == dashmm::kSuccess);
  err = gauss_.put(0, ngauss_, gauss.data());
  assert(err == dashmm::kSuccess);

  dashmm::FMM97<Atom, GNode, dashmm::AFMPBRHS> method{};
  std::vector<double> kparam{};

  err = interp.evaluate(atoms_, gauss_, refine_limit_, &method,
                        accuracy_, &kparam);
  assert(err == dashmm::kSuccess);
  */
}

void AFMPB::computeEnergy(bool status) {
  // Compute the energy on if the GMRES has converged
  if (!status) 
    return; 

  // Collect the all the nodes to create the Gaussian 

  // Setup the Gaussian interpolation points on each element of the mesh
  


  const double unitfactor = 4171.8; 
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

  if (hpx_get_my_rank())
    return; 

  for (int i = 0; i < nnodes_; ++i) {
    nodes[i].gmres[0] *= unitfactor;
    nodes[i].gmres[1] *= unitfactor;
  }

  // Compute total free energy
  double nonpolar = surface_tension_ * area_ + pressure_ * volume_;
  double polar = polarEnergy(gauss.get(), ngauss_, nodes.get(), nnodes_);

  log_ << "\nResults:\n"
       << std::setw(50) << std::left << "Total solvation energy:"
       << std::setw(14) << std::right << std::setprecision(5)
       << std::scientific << nonpolar + polar << "\n"
       << std::setw(50) << std::left << "... Polar part:"
       << std::setw(14) << std::right << std::setprecision(5)
       << std::scientific << polar << "\n"
       << std::setw(50) << std::left << "... Nonpolar part:"
       << std::setw(14) << std::right << std::setprecision(5)
       << std::scientific << nonpolar << "\n" << std::flush;   
}


void AFMPB::collect() {
  int myrank = hpx_get_my_rank(); 
  const double unitfactor = 4171.8;
  auto nodes = nodes_.collect();
  if (nodes) {
    std::sort(&nodes[0], &nodes[nnodes_],
              [] (const Node &a, const Node &b) -> bool {
                return (a.index < b.index);
              });
  }


  if (mesh_format_)  {
    // Now generate the Gaussian points 
    std::vector<GNode> gauss; 
    if (!mesh_format_) {
    } else {
      generateGaussianPoint(nodes.get(), gauss);      
    }
    
    ngauss_ = gauss.size(); 
    auto err = gauss_.allocate(ngauss_); 
    err = gauss_.put(0, ngauss_, gauss.data()); 
    assert(err == dashmm::kSuccess); 
    
    // Compute values on the Gaussian points
    dashmm::FMM97<Atom, GNode, dashmm::AFMPBRHS> method{}; 
    std::vector<double> kparam{};
    
    err = interp.evaluate(atoms_, gauss_, refine_limit_, &method,
                          accuracy_, &kparam);
    assert(err == dashmm::kSuccess);
  
    // Collect results on the Gaussian points 
    auto gauss1 = gauss_.collect(); 
    
    if (gauss1) {
      std::sort(&gauss1[0], &gauss1[ngauss_],
                [] (const GNode &a, const GNode &b) -> bool {
                  return (a.index < b.index);
                });
    }
    

    if (myrank) 
      return; 

    for (int i = 0; i < nnodes_; ++i) {
      nodes[i].gmres[0] *= unitfactor;
      nodes[i].gmres[1] *= unitfactor;
    }

    // Compute total free energy
    double nonpolar = surface_tension_ * area_ + pressure_ * volume_;
    double polar = polarEnergy(gauss1.get(), ngauss_, nodes.get(), nnodes_);
    log_ << "\nResults:\n"
         << std::setw(50) << std::left << "Total solvation energy:"
         << std::setw(14) << std::right << std::setprecision(5)
         << std::scientific << nonpolar + polar << "\n"
         << std::setw(50) << std::left << "... Polar part:"
         << std::setw(14) << std::right << std::setprecision(5)
         << std::scientific << polar << "\n"
         << std::setw(50) << std::left << "... Nonpolar part:"
         << std::setw(14) << std::right << std::setprecision(5)
         << std::scientific << nonpolar << "\n" << std::flush;
  } else {
    if (myrank) 
      return; 

    for (int i = 0; i < nnodes_; ++i) {
      nodes[i].gmres[0] *= unitfactor;
      nodes[i].gmres[1] *= unitfactor;
    }

    // Compute total free energy
    double nonpolar = surface_tension_ * area_ + pressure_ * volume_;
    double polar = polarEnergy(nullptr, 0, nodes.get(), nnodes_);
    log_ << "\nResults:\n"
         << std::setw(50) << std::left << "Total solvation energy:"
         << std::setw(14) << std::right << std::setprecision(5)
         << std::scientific << nonpolar + polar << "\n"
         << std::setw(50) << std::left << "... Polar part:"
         << std::setw(14) << std::right << std::setprecision(5)
         << std::scientific << polar << "\n"
         << std::setw(50) << std::left << "... Nonpolar part:"
         << std::setw(14) << std::right << std::setprecision(5)
         << std::scientific << nonpolar << "\n" << std::flush;
  }    



  /*
  // Write potentials
  potential_.precision(8);
  potential_ << std::scientific;
  for (int i = 0; i < nnodes_; ++i) {
    const Node &n = nodes[i];
    potential_ << n.position.x() << " " << n.position.y() << " "
               << n.position.z() << " " << n.normal_o.x() << " "
               << n.normal_o.y() << " " << n.normal_o.z() << " "
               << n.gmres[0]  << " " << n.gmres[1]  << "\n";
  }

  if (!mesh_format_) {
    for (auto && e : elements_) {
      potential_ << std::setw(8) << e.nodes[0] << " "
                 << std::setw(8) << e.nodes[1] << " "
                 << std::setw(8) << e.nodes[2] << "\n";
    }
  }
  */
}

double AFMPB::polarEnergy(const GNode *gauss, int ngauss,
                          const Node *nodes, int nnodes) const {
  double b = 0;

  if (mesh_format_) {
    for (auto && e : elements_) {
      int i1 = e.nodes[0];
      int i2 = e.nodes[1];
      int i3 = e.nodes[2];
      int index = e.index * 7;
      double temp = 0;
      for (int j = 0; j < 7; ++j) {
        double zeta = 1.0 - xi_[j] - eta_[j];
        double f = nodes[i1].gmres[0] * zeta +
          nodes[i2].gmres[0] * xi_[j] + nodes[i3].gmres[0] * eta_[j];
        double h = nodes[i1].gmres[1] * zeta +
          nodes[i2].gmres[1] * xi_[j] + nodes[i3].gmres[1] * eta_[j];
        temp += (gauss[index + j].rhs[0] * h * dielectric_ -
                 gauss[index + j].rhs[1] * f) * weight_[j];
      }
      b += temp * e.area;
    }

    b /= 8 * M_PI;
  } else {
    /*
    // When using built-in mesh, the number of Gaussian quadrature
    // points is the same as the number of nodes of the surface mesh
    for (int i = 0; i < ngauss; ++i) {
      b += (gauss[i].rhs[0] * nodes[i].gmres[1] * dielectric_ -
            gauss[i].rhs[1] * nodes[i].gmres[0]) * nodes[i].area / 2.0;
    }
    b /= 4 * M_PI / 0.985;
    */

    // When a built-in mesh has been used in previous computation, the
    // Gaussian points here will be the same, meaning the values are
    // saved in the rhs field of each node. However, these values were
    // scaled by 1.0 / dielectric_exterior_ during the GMRES solve
    // phase and the results need to be adjusted here. 
    for (int i = 0; i < nnodes; ++i) {
      b += (nodes[i].rhs[0] * nodes[i].gmres[1] * dielectric_ - 
            nodes[i].rhs[1] * nodes[i].gmres[0]) * nodes[i].area / 2.0; 
    }

    b = b * dielectric_exterior_ * 0.985 / 4 / M_PI; 
  }

  return b;
}

void AFMPB::finalize() {
  if (hpx_get_my_rank() == 0) {
    pqr_.close(); 
    log_.close(); 
    if (potential_.is_open())
      potential_.close(); 
    if (mesh_.is_open())
      mesh_.close();
  }
  
  assert(atoms_.destroy() == dashmm::kSuccess); 
  //assert(gauss_.destroy() == dashmm::kSuccess); 
  assert(nodes_.destroy() == dashmm::kSuccess); 
}


} // namespace afmpb
