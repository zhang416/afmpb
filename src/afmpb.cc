#include <algorithm>
#include <iomanip>
#include <cmath>
#include <hpx/hpx.h>
#include "afmpb.h"
#include "afmpb_lhs.h"
#include "afmpb_rhs.h"
#include "fmm97NL3_method.h"

namespace afmpb {

dashmm::Evaluator<Atom, GNode, dashmm::AFMPBRHS, dashmm::FMM97> interp{}; 
dashmm::Evaluator<Atom, Node, dashmm::AFMPBRHS, dashmm::FMM97> rhs{};
dashmm::Evaluator<Node, Node, dashmm::AFMPBLHS, dashmm::FMM97NL3> lhs{}; 
dashmm::ArrayMapAction<Node, double> rhs_action{set_rhs}; 
dashmm::ArrayMapAction<Node, double> r0_action{set_r0}; 

void AFMPB::solve() {
  // Solve Ax = b using restarted GMRES 

  // Compute right-hand side b 
  dashmm::FMM97<Atom, Node, dashmm::AFMPBRHS> m_rhs{}; 
  std::vector<double> kparam_rhs{}; 
  auto err = rhs.evaluate(atoms_, nodes_, refine_limit_, &m_rhs, 
                          accuracy_, &kparam_rhs); 
  assert(err == dashmm::kSuccess); 

  // Set initial guess x0 = b 
  nodes_.map(rhs_action, &dielectric_exterior_); 

  // Create DAG for Ax computation
  dashmm::FMM97NL3<Node, Node, dashmm::AFMPBLHS> m_lhs{};
  std::vector<double> kparam_lhs;
  kparam_lhs.push_back(kap_); 
  kparam_lhs.push_back(dielectric_); 
  kparam_lhs.push_back(cut1_); 
  kparam_lhs.push_back(cut2_);
  kparam_lhs.push_back(sigma_); 
  kparam_lhs.push_back(restart_); 
  auto tree = lhs.create_tree(nodes_, nodes_, refine_limit_); 
  auto dag = lhs.create_DAG(tree, accuracy_, &kparam_lhs, &m_lhs); 

  // Compute Ax0 and initial residual r0 = b - Ax0
  err = lhs.execute_DAG(tree, dag.get()); 
  nodes_.map(r0_action, (double *)nullptr); 

  // Compute 2-norm of the rhs, setup tolerance for GMRES 
  double rhs_norm2 = generalizedInnerProduct(-1, -1); 
  double tolerance = rel_tolerance_ * rhs_norm2 + abs_tolerance_; 
  
  // Compute 2-norm of r0, normalize it to q0
  residual_[0] = generalizedInnerProduct(0, 0); 

  bool terminateLoop = false, computeApproxSolution = false; 

  while (true) {
    // Compute A * qk
    err = lhs.execute_DAG(tree, dag.get()); 

    // Orthogonalize the result against q0, ..., qk
    modifiedGramSchmidtReOrth(); 

    // Apply previous Givens rotations on the new column of Hessenberg matrix
    // and generate a new one to eliminate the subdiagonal element
    applyGivensRotation(); 
    
    // Compute the new residual norm 
    double alpha = updateResidualNorm(); 

    if (alpha < tolerance) {
      terminateLoop = true; 
      computeApproxSolution = true;
    } else {
      int nMV = dashmm::builtin_afmpb_table_->increFetchIter(); 
      if (nMV == maxMV_ - 1) {
        // Reach maximum allowed matrix-vector multiply 
        terminateLoop = true;
      } else if (nMV % restart_ == restart_ - 1) {
        computeApproxSolution = true;
      }
    }

    if (computeApproxSolution) {
      if (!terminateLoop) 
        dashmm::builtin_afmpb_table_->increIter(); 
    } 


    if (terminateLoop) 
      break;
  }





  // Cleanup 
  err = lhs.destroy_DAG(tree, std::move(dag)); 
  err = lhs.destroy_tree(tree); 
}


void AFMPB::modifiedGramSchmidtReOrth() {  
  int k = dashmm::builtin_afmpb_table_->s_iter(); 

  // Compute the square of the 2-norm of the input vector A * q_k, without
  // normalizing it
  double Aqk_norm2 = generalizedInnerProduct(k + 1, -(k + 1)); 

  // Threhold for performing re-orthogonalization. If the cosine between the two
  // vectors is greather than 0.99, (0.98 = 0.99^2), re-orthogonalization is
  // performed. The norm of the new input vector is kept in Aqk_norm2 and
  // updated after operating with each vector. 
  double threshold = Aqk_norm2 * 0.98; 

  // Starting index of h_{0, k+1} in the Hessenberg matrix 
  int hidx = (k + 1) * (k + 2) / 2; 

  // Modified Gram-Schmidt loop 
  for (int i = 0; i <= k; ++i) {
    // h_{i, k + 1} = <A * q_k, q_i> 
    // A * q_k = A * q_k - h_{i, k + 1} * q_i
    hess_[hidx] = generalizedInnerProduct(k + 1, i); 

    if (hess_[hidx] * hess_[hidx] > threshold) {
      // Reorthogonalization 
      double temp = generalizedInnerProduct(k + 1, i); 
      hess_[hidx] += temp;
    }

    Aqk_norm2 -= hess_[hidx] * hess_[hidx]; 
    if (Aqk_norm2 < 0.0) 
      Aqk_norm2 = 0.0; 

    threshold = Aqk_norm2 * 0.98; 

    hidx++; 
  }

  // Normalize the resulting vector
  Aqk_norm2 = generalizedInnerProduct(k + 1, k + 1); 
  hess_[hidx] = Aqk_norm2; 
}

void AFMPB::applyGivensRotation() {
  int k = dashmm::builtin_afmpb_table_->s_iter(); 

  // Starting index of h_{*, k + 1} in the Hessenberg matrix 
  int hidx = (k + 1) * (k + 2) / 2; 
  
  // Apply the previous k Givens rotations
  for (int i = 0; i < k; ++i) {
    double c = cosine_[i]; 
    double s = sine_[i]; 
    double alpha = hess_[hidx + i]; 
    double beta = hess_[hidx + i + 1]; 
    hess_[hidx + i] = c * alpha + s * beta; 
    hess_[hidx + i + 1] = -s * alpha + c * beta; 
  }

  // Generate a new Givens rotation 
  generateGivensRotation(hess_[hidx + k], hess_[hidx + k + 1]); 
  
}

void AFMPB::generateGivensRotation(double &x, double &y) {
  int k = dashmm::builtin_afmpb_table_->s_iter(); 

  if (x == 0.0 && y == 0.0) {
    cosine_[k] = 1.0; 
    sine_[k] = 0.0;
  } else if (abs(y) > abs(x)) {
    double t = x / y; 
    x = sqrt(1.0 + t * t); 
    sine_[k] = std::copysign(1.0 / x, y); 
    cosine_[k] = t * sine_[k];
  } else if (abs(y) <= abs(x)) {
    double t = y / x; 
    y = sqrt(1.0 + t * t); 
    cosine_[k] = std::copysign(1.0 / y, x); 
    sine_[k] = t * cosine_[k];
  } else {
    // x or y must be an invalid floating-point number, set both to zero
    x = 0.0; 
    y = 0.0; 
    cosine_[k] = 1.0; 
    sine_[k] = 0.0; 
  }

  x = abs(x * y);  
}

double AFMPB::updateResidualNorm() {
  int k = dashmm::builtin_afmpb_table_->s_iter(); 

  double c = cosine_[k]; 
  double s = sine_[k]; 
  double alpha = residual_[k]; 
  residual_[k] = c * alpha; 
  residual_[k + 1] = -s * alpha; 

  return residual_[k + 1]; 
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



} // namespace afmpb
