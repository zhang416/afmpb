#include <iomanip>
#include <cmath>
#include <memory>
#include <chrono>
#include <hpx/hpx.h>
#include "afmpb.h"
#include "afmpb_lhs.h"
#include "afmpb_rhs.h"
#include "afmpb_serializer.h"
#include "fmm97NL3_method.h"

namespace afmpb {

dashmm::Evaluator<Atom, Node, dashmm::AFMPBRHS, dashmm::FMM97> rhs{};
dashmm::Evaluator<Node, Node, dashmm::AFMPBLHS, dashmm::FMM97NL3> lhs{}; 
dashmm::ArrayForEachAction<Node, double> rhs_action{set_rhs}; 
dashmm::ArrayForEachAction<Node, double> r0_action{set_r0}; 

using namespace std::chrono; 
high_resolution_clock::time_point t1, t2; 

bool AFMPB::computePotential() {
  return true;
}

bool AFMPB::solve() {
  int myrank = hpx_get_my_rank(); 

  bool terminate = false, converged = false, compute = false; 

  int n_restarted = 0; 

  // Solve Ax = b using restarted GMRES 
  using Serializer = dashmm::Serializer; 
  using NodeFullSerializer = dashmm::NodeFullSerializer; 
  using NodePartialSerializer = dashmm::NodePartialSerializer; 
  using NodeMinimumSerializer = dashmm::NodeMinimumSerializer; 

  std::unique_ptr<Serializer> m_full{new NodeFullSerializer}; 
  std::unique_ptr<Serializer> m_part{new NodePartialSerializer};
  std::unique_ptr<Serializer> m_min{new NodeMinimumSerializer}; 

  // Compute right-hand side b 
  dashmm::FMM97<Atom, Node, dashmm::AFMPBRHS> m_rhs{}; 
  std::vector<double> kparam_rhs{}; 
 
  auto err = nodes_.set_manager(std::move(m_full)); 
  assert(err == dashmm::kSuccess); 

  err = rhs.evaluate(atoms_, nodes_, refine_limit_, &m_rhs, 
                     accuracy_, &kparam_rhs); 
  assert(err == dashmm::kSuccess); 

  // Scale b, let x0 = b, and copy x0 into q0 slot of Krylov basis
  nodes_.forEach(rhs_action, &dielectric_exterior_); 

  // Compute 2-norm of the rhs, setup tolerance for GMRES 
  double rhs_norm2 = generalizedInnerProduct(-1, -1); 
  double tolerance = rel_tolerance_ * rhs_norm2 + abs_tolerance_; 

  if (!myrank) {
    log_ << "\nSolver status:\n"
         << std::setw(50) << std::left << "... GMRES solver tolerance:"
         << std::setw(14) << std::right << std::setprecision(5) 
         << std::scientific << tolerance << "\n" 
         << std::setw(50) << std::left << "... GMRES restart:" 
         << std::setw(14) << std::right << restart_ + 1<< "\n" 
         << std::setw(50) << std::left << "... GMRES maximum iterations:"
         << std::setw(14) << std::right << (restart_ + 1) * (max_restart_ + 1)
         << "\n" << std::flush; 
  }

  // Create DAG for Ax computation
  dashmm::FMM97NL3<Node, Node, dashmm::AFMPBLHS> m_lhs{};

  std::vector<double> kparam_lhs;
  kparam_lhs.push_back(kap_); 
  kparam_lhs.push_back(dielectric_); 
  kparam_lhs.push_back(cut1_); 
  kparam_lhs.push_back(cut2_);
  kparam_lhs.push_back(sigma_); 
  kparam_lhs.push_back(restart_); 

  t1 = high_resolution_clock::now(); 
  auto tree = lhs.create_tree(nodes_, nodes_, refine_limit_); 
  auto dag = lhs.create_DAG(tree, accuracy_, &kparam_lhs, &m_lhs); 
  t2 = high_resolution_clock::now(); 
  t_dag_ = duration_cast<duration<double>>(t2 - t1).count(); 

  while (terminate == false) {
    dashmm::builtin_afmpb_table_->resetIter(); 

    // Compute Ax0 
    t1 = high_resolution_clock::now(); 
    err = lhs.execute_DAG(tree, dag.get()); 
    assert(err == dashmm::kSuccess); 
    t2 = high_resolution_clock::now(); 
    t_exec_ += duration_cast<duration<double>>(t2 - t1).count(); 

    err = lhs.reset_DAG(dag.get()); 
    assert(err == dashmm::kSuccess); 

    // Compute r0 = b - Ax0
    nodes_.forEach(r0_action, (double *)nullptr); 

    // Compute 2-norm of r0 and normalize r0 to q0 
    residual_[0] = generalizedInnerProduct(0, 0); 
    n_inner_++;

    if (!myrank) {
      if (n_restarted) 
        log_ << "... GMRES solver restarts\n" << std::flush;

      log_ << "... Iteration " << std::setw(3) 
           << n_restarted * (restart_ + 1) 
           << std::setw(33) << std::left << " residual norm:"
           << std::setw(14) << std::right << std::setprecision(5) 
           << std::scientific << residual_[0] << "\n" << std::flush;
    }

    if (n_restarted == 0) {
      // Update the manager once the cached value is computed 
      err = nodes_.set_manager(std::move(m_part)); 
      assert(err == dashmm::kSuccess); 
    }


    double alpha; 

    for (int k = 1; k <= restart_; ++k) {
      // Compute A * q_{k-1}
      t1 = high_resolution_clock::now(); 
      err = lhs.execute_DAG(tree, dag.get()); 
      assert(err == dashmm::kSuccess); 
      t2 = high_resolution_clock::now(); 
      t_exec_ += duration_cast<duration<double>>(t2 - t1).count(); 

      t1 = high_resolution_clock::now(); 
      // Orthogonalize the result against q0, ..., q_{k-1} 
      modifiedGramSchmidtReOrth(); 

      // Apply previous Givens rotations on the new column of the Hessenberg
      // matrix and generate a new one to eliminate the subdiagonal element 
      applyGivensRotation(); 

      // Compute the new residual norm
      alpha = fabs(updateResidualNorm()); 
      t2 = high_resolution_clock::now(); 
      t_gmres_ += duration_cast<duration<double>>(t2 - t1).count();  

      if (!myrank) {
        log_ << "... Iteration " << std::setw(3) 
             << n_restarted * (restart_ + 1) + k 
             << std::setw(33) << std::left << " residual norm:"
             << std::setw(14) << std::right << std::setprecision(5)
             << std::scientific << alpha << "\n" << std::flush;
      }

      if (alpha < tolerance) {
        terminate = true; 
        converged = true; 
        compute = true;
        break;
      } 

      dashmm::builtin_afmpb_table_->increIter(); 

      // Here, alpha is above the tolerance. Reset the DAG if GMRES has not
      // reached maximum allowed matrix-vector multiply
      bool reset_dag = true; 
      if (k == restart_ && n_restarted == max_restart_) 
        reset_dag = false; 

      if (reset_dag) 
        err = lhs.reset_DAG(dag.get()); 
    }

    if (alpha >= tolerance) {      
      if (n_restarted < max_restart_) {
        compute = true;
        n_restarted++;
      } else {
        terminate = true;
      }
    }

    if (compute) {
      t1 = high_resolution_clock::now(); 
      assert(computeApproxSolution(converged) == 0); 
      t2 = high_resolution_clock::now(); 
      t_gmres_ += duration_cast<duration<double>>(t2 - t1).count(); 

      compute = false; 
    }
  }

  if (!myrank) {
    if (converged) {
      log_ << "... GMRES solver has converged\n" << std::flush;
    } else {
      log_ << "... GMRES solver is terminated without convergence\n"
           << std::flush;
    }
  }

  // Cleanup 
  err = lhs.destroy_DAG(tree, std::move(dag)); 
  assert(err == dashmm::kSuccess); 

  err = lhs.destroy_tree(tree); 
  assert(err == dashmm::kSuccess); 

  // Update serialization manager 
  err = nodes_.set_manager(std::move(m_min)); 
  assert(err == dashmm::kSuccess); 

  if (!myrank && converged) {
    log_ << "\nSolver statistics:\n"
         << std::setw(50) << std::left << "... t(dag_construction):"
         << std::setw(14) << std::right << std::setprecision(5)
         << std::scientific << t_dag_ << "\n" 
         << std::setw(50) << std::left << "... t(dag_execution):"
         << std::setw(14) << std::right << std::setprecision(5) 
         << std::scientific << t_exec_ << "\n"
         << std::setw(50) << std::left << "... t(GMRES):" 
         << std::setw(14) << std::right << std::setprecision(5) 
         << std::scientific << t_gmres_ << "\n" 
         << std::setw(50) << std::left << "...... t(inner_product):"
         << std::setw(14) << std::right << std::setprecision(5) 
         << std::scientific << t_inner_ << "\n"
         << std::setw(50) << std::left << "...... n(inner_product):"
         << std::setw(14) << std::right << n_inner_ 
         << "\n" << std::flush;
  }

  return (converged ? true : false); 
}


void AFMPB::modifiedGramSchmidtReOrth() {  
  int k = dashmm::builtin_afmpb_table_->s_iter(); 

  // [Aq_0, ...., Aq_k] = [q_0, ..., q_k, q_{k + 1}] * hess; 
  
  // Compute the square of the 2-norm of the input vector A * q_k, without
  // normalizing it
  double Aqk_norm2 = generalizedInnerProduct(k + 1, -(k + 1)); 
  n_inner_++; 

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
    n_inner_++;

    if (hess_[hidx] * hess_[hidx] > threshold) {
      // Reorthogonalization 
      double temp = generalizedInnerProduct(k + 1, i); 
      hess_[hidx] += temp;
      n_inner_++;
    }

    Aqk_norm2 -= hess_[hidx] * hess_[hidx]; 
    if (Aqk_norm2 < 0.0) 
      Aqk_norm2 = 0.0; 

    threshold = Aqk_norm2 * 0.98; 

    hidx++; 
  }

  // Normalize the resulting vector
  Aqk_norm2 = generalizedInnerProduct(k + 1, k + 1); 
  n_inner_++;
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

int AFMPB::computeApproxSolution(bool converged) {
  int k = (converged ? dashmm::builtin_afmpb_table_->t_iter() : 
           dashmm::builtin_afmpb_table_->s_iter()); 

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
