#include <iomanip>
#include <cmath>
#include <hpx/hpx.h>
#include "afmpb.h"
#include "afmpb_lhs.h"
#include "afmpb_rhs.h"
#include "fmm97NL3_method.h"

namespace afmpb {

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

  // Compute 2-norm of the rhs, setup tolerance for GMRES 
  double rhs_norm2 = generalizedInnerProduct(-1, -1); 
  double tolerance = rel_tolerance_ * rhs_norm2 + abs_tolerance_; 

  log_ << "Solver status:\n"
       << std::setw(50) << std::left << "... GMRES solver tolerance:"
       << std::setw(14) << std::right << std::setprecision(5) 
       << std::scientific << tolerance << "\n"; 

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
  err = lhs.reset_DAG(dag.get());  
  nodes_.map(r0_action, (double *)nullptr); 
  
  // Compute 2-norm of r0, normalize it to q0
  residual_[0] = generalizedInnerProduct(0, 0); 

  log_ << std::setw(50) << std::left << "... Iteration 0 residual norm:" 
       << std::setw(14) << std::right << std::setprecision(5) 
       << std::scientific << residual_[0] << "\n";

  bool terminateLoop = false, computeSolution = false; 

  while (true) {
    // Compute A * qk
    err = lhs.execute_DAG(tree, dag.get()); 
    err = lhs.reset_DAG(dag.get());  

    // Orthogonalize the result against q0, ..., qk
    modifiedGramSchmidtReOrth(); 

    // Apply previous Givens rotations on the new column of Hessenberg matrix
    // and generate a new one to eliminate the subdiagonal element
    applyGivensRotation(); 
    
    // Compute the new residual norm 
    double alpha = fabs(updateResidualNorm()); 

    int nMV = dashmm::builtin_afmpb_table_->increFetchIter(); 

    log_ << std::setw(50) << std::left << "... Iteration " << nMV
         << " residual norm:" << std::setw(14) << std::right 
         << std::setprecision(5) << std::scientific << alpha << "\n"; 

    if (alpha < tolerance) {
      terminateLoop = true; 
      computeSolution = true;

      log_ << std::setw(50) << std::left << "... GMRES solver has converged\n";
    } else {
      if (nMV == maxMV_ - 1) {
        // Reach maximum allowed matrix-vector multiply 
        terminateLoop = true;

        log_ << std::setw(50) << std::left 
             << "... GMRES solver is terminated without convergence\n"; 

      } else if (nMV % restart_ == restart_ - 1) {
        computeSolution = true;
        
        log_ << std::setw(50) << std::left << "... GMRES solver restarts\n"; 
      }
    }

    if (computeSolution) {
      assert(computeApproxSolution() == 0); 
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

  // [Aq_0, ...., Aq_k] = [q_0, ..., q_k, q_{k + 1}] * hess; 
  
  // Compute the square of the 2-norm of the input vector A * q_k, without
  // normalizing it
  double Aqk_norm2 = generalizedInnerProduct(k + 1, -(k + 1)); 

  // Threhold for performing re-orthogonalization. If the cosine between the two
  // vectors is greather than 0.99, (0.98 = 0.99^2), re-orthogonalization is
  // performed. The norm of the new input vector is kept in Aqk_norm2 and
  // updated after operating with each vector. 
  double threshold = Aqk_norm2 * 0.98; 

  // Starting index of h_{0, k} in the Hessenberg matrix 
  int hidx = k * (k + 1) / 2; 

  // Modified Gram-Schmidt loop 
  for (int i = 0; i <= k; ++i) {
    // h_{i, k} = <Aq_k, q_i> 
    // Aq_k = Aq_k - h_{i, k} * q_i
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

  // [Aq_0, ...., Aq_k] = [q_0, ..., q_k, q_{k + 1}] * hess; 

  // Index of h_{0, k} of the Hessenberg matrix 
  int hidx = k * (k + 1) / 2; 
  
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
  } else if (fabs(y) > fabs(x)) {
    double t = x / y; 
    x = sqrt(1.0 + t * t); 
    sine_[k] = std::copysign(1.0 / x, y); 
    cosine_[k] = t * sine_[k];
  } else if (fabs(y) <= fabs(x)) {
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

  x = fabs(x * y);  
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

int AFMPB::computeApproxSolution() {
  int k = dashmm::builtin_afmpb_table_->s_iter(); 

  // [Aq_0, ..., Aq_{k - 1}] = [q_0, ..., q_{k - 1}, q_k] * hess_; 
  
  // Check the diagonal element of the last column of Hessenberg matrix. If its
  // value is zero, reduce k by 1 so that a smaller triangular system is
  // solved. [It should only happen when the matrix is singular, and at most
  // once.] 
  while (k >= 0) {
    int hidx = k * (k + 1) / 2 - 1; // hess_{k - 1, k - 1} 
    if (hess_[hidx] == 0) {
      k--;
    } else {
      break;
    }
  }

  if (k < 0) // Triangular system has null rank
    return -1;

  // Backward solve of the triangular system, overwrite residual_ 
  for (int i = k - 1; i >= 0; i--) {
    int hidx = i * (i + 1) / 2; // hess_{0, i} 
    residual_[i] /= hess_[hidx + i]; 
    for (int j = 0; j <= i - 1; ++j) {
      residual_[j] -= residual_[i] * hess_[hidx + j]; 
    }
  }

  // Compute linear combination y_0 * q_0 + ... + y_{k - 1} * q_{k - 1}
  linearCombination(k - 1);    

  return 0;
}

} // namespace afmpb
