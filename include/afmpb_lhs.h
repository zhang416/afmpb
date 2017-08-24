//=============================================================================
// AFMPB: Adaptive Fast Multipole Poisson-Boltzmann Solver 
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

#ifndef __AFMPB_LHS_H__
#define __AFMPB_LHS_H__

#include <cmath>
#include <complex>
#include <map>
#include <memory>
#include <vector> 

#include "dashmm/index.h"
#include "builtins/laplace.h"
#include "builtins/yukawa.h"
#include "dashmm/point.h"
#include "dashmm/types.h"
#include "dashmm/viewset.h"

namespace dashmm {

void dlap_s_to_m(Point dist, double q, double scale, Point normal, 
                 dcomplex_t *M); 
void dlap_s_to_l(Point dist, double q, double scale, Point normal, 
                 dcomplex_t *L); 
void dyuk_s_to_m(Point dist, double q, double scale, Point normal, 
                 dcomplex_t *M); 
void dyuk_s_to_l(Point dist, double q, double scale, Point normal, 
                 dcomplex_t *L); 

class AFMPBTable {
public: 
  AFMPBTable(double dielectric, double cut1, double cut2, 
             double sigma, int restart) 
    : dielectric_{dielectric}, cut1_{cut1}, cut2_{cut2}, 
    sigma_{sigma}, restart_{restart} { }
  double dielectric() const {return dielectric_;}
  double cut1() const {return cut1_;}
  double cut2() const {return cut2_;}
  double sigma() const {return sigma_;}
  int s_iter() const {return iter_;} 
  int t_iter() const {return iter_ + 1;} 
  void resetIter() {iter_ = 0;}
  void increIter() {iter_++;} 

private:
  double dielectric_; 
  double cut1_; 
  double cut2_;
  double sigma_;
  int restart_; 
  int iter_ = 0; 
};

extern std::unique_ptr<AFMPBTable> builtin_afmpb_table_; 

void update_afmpb_table(double dielectric, double cut1, double cut2, 
                        double sigma, int restart); 

template <typename Source, typename Target>
class AFMPBLHS {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = AFMPBLHS<Source, Target>;

  AFMPBLHS(ExpansionRole role, double scale = 1.0, Point center = Point{})
    : views_{ViewSet{role, center, scale}} {
    // The AFMPBLHS class handles four kernel expansions
    // Single-layer Laplace potential
    // Double-layer Laplace potential
    // Single-layer Yukawa potential
    // Double-Layer Yukawa potential 
    // The scale saved in the views_ is for Yukawa potential
    
    // Order of the spherical harmonic expansion 
    int p = builtin_laplace_table_->p(); 
    int nsh = (p + 1) * (p + 2) / 2; 

    // Length of the exponential expansion 
    int nexp_l = builtin_laplace_table_->nexp(); 
    int nexp_y = builtin_yukawa_table_->nexp(scale); 

    if (role == kSourcePrimary || role == kTargetPrimary) {
      size_t bytes = sizeof(dcomplex_t) * nsh; 
      for (int i = 0; i < 4; ++i) {
        char *data = new char[bytes](); 
        views_.add_view(i, bytes, data); 
      }
    } else if (role == kSourceIntermediate) {
      size_t bytes_l = sizeof(dcomplex_t) * nexp_l; 
      size_t bytes_y = sizeof(dcomplex_t) * nexp_y;
     
      for (int i = 0; i < 12; ++i) {
        char *data = new char[bytes_l](); 
        views_.add_view(i, bytes_l, data);
      }

      for (int i = 12; i < 24; ++i) {
        char *data = new char[bytes_y](); 
        views_.add_view(i, bytes_y, data);
      }
    } else if (role == kTargetIntermediate) {
      size_t bytes_l = sizeof(dcomplex_t) * nexp_l; 
      size_t bytes_y = sizeof(dcomplex_t) * nexp_y; 

      for (int i = 0; i < 56; ++i) {
        char *data = new char[bytes_l](); 
        views_.add_view(i, bytes_l, data); 
      }

      for (int i = 56; i < 112; ++i) {
        char *data = new char[bytes_y](); 
        views_.add_view(i, bytes_y, data); 
      }
    }
  }

  AFMPBLHS(const ViewSet &views) : views_{views} { }

  ~AFMPBLHS() {
    int count = views_.count();
    if (count) {
      for (int i = 0; i < count; ++i) {
        delete [] views_.view_data(i);
      }
    }
  }

  void release() {views_.clear();}

  bool valid(const ViewSet &view) const {
    // \p view is assumed to be a subset of \p views_ (no range checking
    // performed). The function returns true if and only if each entry in the
    // required subset is associated with some data.
    bool is_valid = true;
    int count = view.count();
    for (int i = 0; i < count; ++i) {
      int idx = view.view_index(i);
      if (views_.view_data(idx) == nullptr) {
        is_valid = false;
        break;
      }
    }
    return is_valid;
  }

  int view_count() const { return views_.count(); }

  ViewSet get_all_views() const {return views_;}

  ExpansionRole role() const {return views_.role();}

  Point center() const {return views_.center();}

  size_t view_size(int view) const {
    return views_.view_bytes(view) / sizeof(dcomplex_t);
  }

  dcomplex_t view_term(int view, size_t i) const {
    dcomplex_t *data = reinterpret_cast<dcomplex_t *>(views_.view_data(view));
    return data[i];
  }

  std::unique_ptr<expansion_t> S_to_M(const Source *first, 
                                      const Source *last) const {
    Point center = views_.center();
    double scale_y = views_.scale(); 
    int level = builtin_yukawa_table_->level(scale_y); 
    double scale_l = builtin_laplace_table_->scale(level);    
    int iter = builtin_afmpb_table_->s_iter(); 

    expansion_t *ret{new expansion_t{kSourcePrimary, scale_y, center}};
    dcomplex_t *M1 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(0));
    dcomplex_t *M2 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(1)); 
    dcomplex_t *M3 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(2)); 
    dcomplex_t *M4 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(3)); 

    for (auto i = first; i != last; ++i) {
      double scale_c = i->area * 0.79577472e-1; 
      Point dist = point_sub(i->position, center); 
      double q0 = i->gmres[2 * iter] * scale_c; 
      double q1 = i->gmres[2 * iter + 1] * scale_c; 
      lap_s_to_m(dist, q1, scale_l, M1); 
      dlap_s_to_m(dist, q0, scale_l, i->normal_i, M2); 
      yuk_s_to_m(dist, q1, scale_y, M3); 
      dyuk_s_to_m(dist, q0, scale_y, i->normal_i, M4); 

    }
    return std::unique_ptr<expansion_t>{ret};
  }

  std::unique_ptr<expansion_t> S_to_L(const Source *first, 
                                      const Source *last) const {
    Point center = views_.center();
    double scale_y = views_.scale();
    int level = builtin_yukawa_table_->level(scale_y); 
    double scale_l = builtin_laplace_table_->scale(level); 
    int iter = builtin_afmpb_table_->s_iter(); 

    expansion_t *ret{new expansion_t{kTargetPrimary}};
    dcomplex_t *L1 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(0));
    dcomplex_t *L2 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(1));
    dcomplex_t *L3 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(2));
    dcomplex_t *L4 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(3));

    for (auto i = first; i != last; ++i) {
      double scale_c = i->area * 0.79577472e-1; 
      Point dist = point_sub(i->position, center); 
      double q0 = i->gmres[2 * iter] * scale_c; 
      double q1 = i->gmres[2 * iter + 1] * scale_c; 
      lap_s_to_l(dist, q1, scale_l, L1);
      dlap_s_to_l(dist, q0, scale_l, i->normal_i, L2); 
      yuk_s_to_l(dist, q1, scale_y, L3); 
      dyuk_s_to_l(dist, q0, scale_y, i->normal_i, L4); 
    }
    return std::unique_ptr<expansion_t>{ret};
  }

  std::unique_ptr<expansion_t> M_to_M(int from_child) const {
    double scale = views_.scale(); 

    dcomplex_t *M1 = reinterpret_cast<dcomplex_t *>(views_.view_data(0));
    dcomplex_t *M2 = reinterpret_cast<dcomplex_t *>(views_.view_data(1));
    dcomplex_t *M3 = reinterpret_cast<dcomplex_t *>(views_.view_data(2));
    dcomplex_t *M4 = reinterpret_cast<dcomplex_t *>(views_.view_data(3));

    expansion_t *ret{new expansion_t{kSourcePrimary}};
    dcomplex_t *W1 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(0));
    dcomplex_t *W2 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(1));
    dcomplex_t *W3 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(2));
    dcomplex_t *W4 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(3));

    lap_m_to_m(from_child, M1, W1); 
    lap_m_to_m(from_child, M2, W2); 
    yuk_m_to_m(from_child, M3, scale, W3); 
    yuk_m_to_m(from_child, M4, scale, W4); 
    return std::unique_ptr<expansion_t>{ret};
  }

  std::unique_ptr<expansion_t> M_to_L(Index s_index, Index t_index) const {
    return std::unique_ptr<expansion_t>{nullptr};
  }

  std::unique_ptr<expansion_t> L_to_L(int to_child) const {
    double scale = views_.scale(); 

    dcomplex_t *L1 = reinterpret_cast<dcomplex_t *>(views_.view_data(0));
    dcomplex_t *L2 = reinterpret_cast<dcomplex_t *>(views_.view_data(1));
    dcomplex_t *L3 = reinterpret_cast<dcomplex_t *>(views_.view_data(2));
    dcomplex_t *L4 = reinterpret_cast<dcomplex_t *>(views_.view_data(3));

    expansion_t *ret{new expansion_t{kTargetPrimary}};
    dcomplex_t *W1 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(0));
    dcomplex_t *W2 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(1));
    dcomplex_t *W3 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(2));
    dcomplex_t *W4 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(3));

    lap_l_to_l(to_child, L1, W1); 
    lap_l_to_l(to_child, L2, W2);
    yuk_l_to_l(to_child, L3, scale, W3);
    yuk_l_to_l(to_child, L4, scale, W4); 
    return std::unique_ptr<expansion_t>{ret};
  }

  void M_to_T(Target *first, Target *last) const {
  }

  void L_to_T(Target *first, Target *last) const {
    double scale_y = views_.scale(); 
    int level = builtin_yukawa_table_->level(scale_y); 
    double scale_l = builtin_laplace_table_->scale(level); 
    double dielectric = builtin_afmpb_table_->dielectric(); 
    double lambda = builtin_yukawa_table_->lambda(); 
    int iter = builtin_afmpb_table_->t_iter(); 

    dcomplex_t *L1 = reinterpret_cast<dcomplex_t *>(views_.view_data(0));
    dcomplex_t *L2 = reinterpret_cast<dcomplex_t *>(views_.view_data(1));
    dcomplex_t *L3 = reinterpret_cast<dcomplex_t *>(views_.view_data(2));
    dcomplex_t *L4 = reinterpret_cast<dcomplex_t *>(views_.view_data(3));
    
    for (auto i = first; i != last; ++i) {
      double f = 0, h = 0; 
      double nx = i->normal_o.x(); 
      double ny = i->normal_o.y(); 
      double nz = i->normal_o.z(); 

      Point dist = point_sub(i->position, views_.center()); 
      // Single layer Laplace 
      auto SL = lap_l_to_t(dist, scale_l, L1, true); 
      f -= 4 * M_PI * SL[0];
      h += 4 * M_PI * (SL[1] * nx + SL[2] * ny +  SL[3] * nz); 

      // Double layer Laplace 
      auto DL = lap_l_to_t(dist, scale_l, L2, true); 
      f -= DL[0] / dielectric * 4 * M_PI;
      h -= 4 * M_PI * (DL[1] * nx + DL[2] * ny + DL[3] * nz) / dielectric;

      // Single layer Yukawa 
      auto SY = yuk_l_to_t(dist, scale_y, L3, true); 
      f += SY[0] * 8 * lambda; 
      h += 8 * lambda * (SY[1] * nx + SY[2] * ny + SY[3] * nz) / dielectric; 

      // Double layer Yukawa
      auto DY = yuk_l_to_t(dist, scale_y, L4, true); 
      f += DY[0] * 8 * lambda; 
      h -= 8 * lambda * (DY[1] * nx + DY[2] * ny + DY[3] * nz) / dielectric;

      i->gmres[2 * iter] += f; 
      i->gmres[2 * iter + 1] += h; 
    }
  }

  void S_to_T(const Source *s_first, const Source *s_last,
              Target *t_first, Target *t_last) const {  
    int key = s_first->index; 
    int s_iter = builtin_afmpb_table_->s_iter(); 
    int t_iter = builtin_afmpb_table_->t_iter(); 

    for (auto i = t_first; i != t_last; ++i) {
      double f = 0, h = 0; 
      auto it = i->cached.find(key); 
      std::vector<double> tbl; 

      if (it != i->cached.end()) {
        tbl = it->second; 
      } else {
        // Generate table 
        generate_direct_table(i, s_first, s_last, tbl); 
        i->cached[key] = tbl;
      }

      int k = 0; 
      for (auto j = s_first; j != s_last; ++j) {
        double q0 = j->gmres[2 * s_iter]; 
        double q1 = j->gmres[2 * s_iter + 1]; 
        f += (tbl[k] * q1 + tbl[k + 1] * q0); 
        h += (tbl[k + 2] * q1 + tbl[k + 3] * q0); 
        k += 4; 
      }

      i->gmres[2 * t_iter] += f; 
      i->gmres[2 * t_iter + 1] += h;
    }  
  }

  std::unique_ptr<expansion_t> M_to_I() const {
    dcomplex_t *M1 = reinterpret_cast<dcomplex_t *>(views_.view_data(0));
    dcomplex_t *M2 = reinterpret_cast<dcomplex_t *>(views_.view_data(1));
    dcomplex_t *M3 = reinterpret_cast<dcomplex_t *>(views_.view_data(2));
    dcomplex_t *M4 = reinterpret_cast<dcomplex_t *>(views_.view_data(3));

    double scale_y = views_.scale(); 
    expansion_t *ret{new expansion_t{kSourceIntermediate, scale_y}}; 
    
    lap_m_to_i(M1, ret->views_, 0); 
    lap_m_to_i(M2, ret->views_, 6); 
    yuk_m_to_i(M3, ret->views_, scale_y, 12); 
    yuk_m_to_i(M4, ret->views_, scale_y, 18); 
    return std::unique_ptr<expansion_t>(ret);
  }

  std::unique_ptr<expansion_t> I_to_I(Index s_index, Index t_index) const {
    ViewSet views{kTargetIntermediate}; 
    double scale = views_.scale(); 

    lap_i_to_i(s_index, t_index, views_, 0, 0, views); 
    lap_i_to_i(s_index, t_index, views_, 6, 28, views); 
    yuk_i_to_i(s_index, t_index, views_, 12, 56, scale, views); 
    yuk_i_to_i(s_index, t_index, views_, 18, 84, scale, views); 

    expansion_t *ret = new expansion_t{views};
    return std::unique_ptr<expansion_t>{ret};
  }

  std::unique_ptr<expansion_t> I_to_L(Index t_index) const {
    // t_index is the index of the child
    expansion_t *ret{new expansion_t{kTargetPrimary}};
    dcomplex_t *L1 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(0));
    dcomplex_t *L2 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(1));
    dcomplex_t *L3 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(2));
    dcomplex_t *L4 = reinterpret_cast<dcomplex_t *>(ret->views_.view_data(3));

    double scale_y = views_.scale(); 
    int level = builtin_yukawa_table_->level(scale_y); 
    double scale_l = builtin_laplace_table_->scale(level); 

    lap_i_to_l(views_, 0, t_index, scale_l * 2.0, L1); 
    lap_i_to_l(views_, 28, t_index, scale_l * 2.0, L2); 
    yuk_i_to_l(views_, 56, t_index, scale_y / 2.0, L3); 
    yuk_i_to_l(views_, 84, t_index, scale_y / 2.0, L4); 

    return std::unique_ptr<expansion_t>(ret);
  }

  void add_expansion(const expansion_t *temp1) {
    // This operation assumes that the views included in \p temp1 is a subset of
    // \p views_. No range checking performed.
    int count = temp1->views_.count();
    for (int i = 0; i < count; ++i) {
      int idx = temp1->views_.view_index(i);
      int size = temp1->views_.view_bytes(i) / sizeof(dcomplex_t);
      dcomplex_t *lhs = reinterpret_cast<dcomplex_t *>(views_.view_data(idx));
      dcomplex_t *rhs =
        reinterpret_cast<dcomplex_t *>(temp1->views_.view_data(i));

      for (int j = 0; j < size; ++j) {
        lhs[j] += rhs[j];
      }
    }
  }

  static void update_table(int n_digits, double domain_size,
                           const std::vector<double> &kernel_params) {
    update_laplace_table(n_digits, domain_size); 
    update_yukawa_table(n_digits, domain_size, kernel_params[0]); 
    update_afmpb_table(kernel_params[1], kernel_params[2], 
                       kernel_params[3], kernel_params[4], 
                       ((int) kernel_params[5])); 
  }

  static void delete_table() { }

  static double compute_scale(Index index) {
    return builtin_yukawa_table_->scale(index.level());
  }

  static int weight_estimate(Operation op,
                             Index s = Index{}, Index t = Index{}) {
    int weight = 0;
    if (op == Operation::MtoI) {
      weight = 6;
    } else if (op == Operation::ItoI) {
      int dx = s.x() - 2 * t.x();
      int dy = s.y() - 2 * t.y();
      int dz = s.z() - 2 * t.z();
      for (int i = 0; i < 3; ++i) {
        int tag = merge_and_shift_table[dx + 2][dy + 2][dz + 2][i];
        if (tag == -1) {
          break;
        }
        weight++;
      }
    } else {
      weight = 1;
    }
    return weight;
  }


 private:
  ViewSet views_;

  void generate_direct_table(Target *t, const Source *s_first, 
                             const Source *s_last, 
                             std::vector<double> &table) const {    
    double dielectric = builtin_afmpb_table_->dielectric();
    double cut1 = builtin_afmpb_table_->cut1(); 
    double cut2 = builtin_afmpb_table_->cut2(); 

    double factor = (1.0 / 2.0 + 1.0 / 2.0 / dielectric) * 4 * M_PI; 

    for (auto i = s_first; i != s_last; ++i) {
      Point dist = point_sub(t->position, i->position); 
      double r = dist.norm(); 
      
      if (r == 0) {
        table.push_back(0.0); 
        table.push_back(factor); 
        table.push_back(factor); 
        table.push_back(0.0);
      } else if (r < cut1) {
        table.push_back(0.0); 
        table.push_back(0.0); 
        table.push_back(0.0); 
        table.push_back(0.0);
      } else if (r < cut2) {
        double A = 0, B = 0, C = 0, D = 0; 
        compute_close_coeff(t, i, A, B, C, D); 
        table.push_back(-A); 
        table.push_back(B); 
        table.push_back(-C);
        table.push_back(D);
      } else {
        double A = 0, B = 0, C = 0, D = 0; 
        compute_nsingular_coeff(t, i, A, B, C, D); 
        table.push_back(-A); 
        table.push_back(B); 
        table.push_back(-C);
        table.push_back(D);
      }        
    }
  }

  void compute_close_coeff(Target *t, const Source *s, double &A, 
                           double &B, double &C, double &D) const {
    double dielectric = builtin_afmpb_table_->dielectric(); 

    auto patch = s->patch; 
    for (auto && p : patch) {
      double a, b, c, d; 
      compute_coeff(t->position, t->normal_o, p.position, p.normal, 
                    a, b, c, d); 
      A += a * p.weight; 
      B += b * p.weight; 
      C += c * p.weight; 
      D += d * p.weight;
    }
    D /= dielectric; 
  }
  
  void compute_nsingular_coeff(Target *t, const Source *s, double &A, 
                               double &B, double &C, double &D) const {
    compute_coeff(t->position, t->normal_o, s->position, s->normal_i, 
                  A, B, C, D); 

    double dielectric = builtin_afmpb_table_->dielectric(); 
    A *= s->area;
    B *= s->projected; 
    C *= s->area; 
    D *= s->projected / dielectric; 
  }

  void compute_coeff(const Point &t, const Point &tn, 
                     const Point &s, const Point &sn, 
                     double &A, double &B, double &C, double &D) const {
    double lambda = builtin_yukawa_table_->lambda(); 
    double sigma = builtin_afmpb_table_->sigma(); 
    double dielectric = builtin_afmpb_table_->dielectric(); 

    Point dx = point_sub(s, t); 
    double r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2] + sigma; 
    double r = sqrt(r2); 
    double r3 = r2 * r; 
    double r4 = r2 * r2; 
    double r5 = r2 * r3; 
    double gr0 = 1.0; 
    double gr1 = gr0 / r; 
    double gr3 = gr0 / r3;
    double gr5 = gr0 / r5; 
    double ur0 = exp(-lambda * r); 
    double ur1 = ur0 / r; 
    double ur2 = ur0 / r2;
    double ur3 = ur0 / r3;
    double ur4 = ur0 / r4;
    double ur5 = ur0 / r5; 
    double ur3ur2 = ur3 + lambda * ur2; 
    double pur4ur3 = lambda * (3 * ur4 + lambda * ur3) + 3 * ur5; 
    double gd[3], ud[3], gdd[3][3], udd[3][3]; 
   
    for (int j = 0; j < 3; ++j) {
      gd[j] = -dx[j] * gr3; 
      ud[j] = -dx[j] * ur3ur2; 
      for (int k = 0; k < 3; ++k) {
        int p0 = (k == j ? 1 : 0); 
        double d0 = dx[j] * dx[k]; 
        gdd[j][k] = 3.0 * d0 * gr5 - p0 * gr3; 
          udd[j][k] = d0 * pur4ur3 - p0 * ur3ur2;
      }
    }
    
    A = gr1 - ur1; 
    B = 0; 
    C = 0; 
    D = 0; 

    for (int j = 0; j < 3; ++j) {
      B += (gd[j] / dielectric - ud[j]) * sn[j]; 
      C -= (gd[j] - ud[j] / dielectric) * tn[j];
      for (int k = 0; k < 3; ++k) {
        D -= tn[k] * (gdd[j][k] - udd[j][k]) * sn[j];
      }
    }  
  }
};

} // namespace dashmm

#endif // __AFMPB_LHS_H__
