#ifndef __AFMPB_RHS_H__
#define __AFMPB_RHS_H__

#include <cassert>
#include <cmath>
#include <complex>
#include <map>
#include <memory>
#include <vector> 

#include "dashmm/index.h"
#include "builtins/laplace.h"
#include "dashmm/point.h"
#include "dashmm/types.h"
#include "dashmm/viewset.h"

namespace dashmm {

template <typename Source, typename Target>
class AFMPBRHS {
 public:
  using source_t = Source;
  using target_t = Target;
  using expansion_t = AFMPBRHS<Source, Target>;

  AFMPBRHS(ExpansionRole role, double scale = 1.0, Point center = Point{})
    : views_{ViewSet{role, center, scale}} {
    // View size for each spherical harmonic expansion
    int p = builtin_laplace_table_->p();
    int nsh = (p + 1) * (p + 2) / 2;

    // View size for each exponential expansion
    int nexp = builtin_laplace_table_->nexp();

    if (role == kSourcePrimary || role == kTargetPrimary) {
      size_t bytes = sizeof(dcomplex_t) * nsh;
      char *data = new char[bytes]();
      views_.add_view(0, bytes, data);
    } else if (role == kSourceIntermediate) {
      size_t bytes = sizeof(dcomplex_t) * nexp;
      for (int i = 0; i < 6; ++i) {
        char *data = new char[bytes]();
        views_.add_view(i, bytes, data);
      }
    } else if (role == kTargetIntermediate) {
      size_t bytes = sizeof(dcomplex_t) * nexp;
      for (int i = 0; i < 28; ++i) {
        char *data = new char[bytes]();
        views_.add_view(i, bytes, data);
      }
    }
  }

  AFMPBRHS(const ViewSet &views) : views_{views} { }

  ~AFMPBRHS() {
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

  std::unique_ptr<expansion_t> S_to_M(Source *first, Source *last) const {
    double scale = views_.scale();
    Point center = views_.center();
    expansion_t *retval{new expansion_t{kSourcePrimary, scale, center}};
    dcomplex_t *M = reinterpret_cast<dcomplex_t *>(retval->views_.view_data(0));
    for (auto i = first; i != last; ++i) {
      Point dist = point_sub(i->position, center);
      lap_s_to_m(dist, i->charge, scale, M); 
    }
   return std::unique_ptr<expansion_t>{retval};
  }

  std::unique_ptr<expansion_t> S_to_L(Source *first, Source *last) const {
    double scale = views_.scale();
    Point center = views_.center();
    expansion_t *retval{new expansion_t{kTargetPrimary}};
    dcomplex_t *L = reinterpret_cast<dcomplex_t *>(retval->views_.view_data(0));
    for (auto i = first; i != last; ++i) {
      Point dist = point_sub(i->position, center); 
      lap_s_to_l(dist, i->charge, scale, L);
    }
    return std::unique_ptr<expansion_t>{retval};
  }

  std::unique_ptr<expansion_t> M_to_M(int from_child) const {
    expansion_t *retval{new expansion_t{kSourcePrimary}};
    dcomplex_t *M = reinterpret_cast<dcomplex_t *>(views_.view_data(0));
    dcomplex_t *W = reinterpret_cast<dcomplex_t *>(retval->views_.view_data(0));
    lap_m_to_m(from_child, M, W); 
    return std::unique_ptr<expansion_t>{retval};
  }

  std::unique_ptr<expansion_t> M_to_L(Index s_index, Index t_index) const {
    return std::unique_ptr<expansion_t>{nullptr};
  }

  std::unique_ptr<expansion_t> L_to_L(int to_child) const {
    expansion_t *retval{new expansion_t{kTargetPrimary}};
    dcomplex_t *L = reinterpret_cast<dcomplex_t *>(views_.view_data(0));
    dcomplex_t *W = reinterpret_cast<dcomplex_t *>(retval->views_.view_data(0));
    lap_l_to_l(to_child, L, W); 
    return std::unique_ptr<expansion_t>{retval};
  }

  void M_to_T(Target *first, Target *last) const {
    double scale = views_.scale();
    dcomplex_t *M = reinterpret_cast<dcomplex_t *>(views_.view_data(0));

    for (auto i = first; i != last; ++i) {
      Point dist = point_sub(i->position, views_.center()); 
      auto result = lap_m_to_t(dist, scale, M, true); 
      i->value[0] += result[0]; 
      i->value[1] += (result[1] * i->normal_o.x() + result[2] * i->normal_o.y() 
                      + result[3] * i->normal_o.z()); 
    }
  }

  void L_to_T(Target *first, Target *last) const {
    double scale = views_.scale();
    dcomplex_t *L = reinterpret_cast<dcomplex_t *>(views_.view_data(0));

    for (auto i = first; i != last; ++i) {
      Point dist = point_sub(i->position, views_.center()); 
      auto result = lap_l_to_t(dist, scale, L, true); 
      i->value[0] += result[0]; 
      i->value[1] += (result[1] * i->normal_o.x() + result[2] * i->normal_o.y() 
                      + result[3] * i->normal_o.z()); 
    }
  }

  void S_to_T(Source *s_first, Source *s_last,
              Target *t_first, Target *t_last) const {
    for (auto i = t_first; i != t_last; ++i) {
      double potential = 0, fx = 0, fy = 0, fz = 0; 
      for (auto j = s_first; j != s_last; ++j) {
        Point s2t = point_sub(i->position, j->position); 
        double r = s2t.norm(); 
        if (r > 0) {
          double temp = j->charge / pow(r, 3); 
          potential += j->charge / r; 
          fx += temp * s2t.x(); 
          fy += temp * s2t.y();
          fz += temp * s2t.z();
        }
      }
      i->value[0] += potential; 
      i->value[1] += (fx * i->normal_o.x() + fy * i->normal_o.y() + 
                      fz * i->normal_o.z()); 
    }
  }

  std::unique_ptr<expansion_t> M_to_I() const {
    expansion_t *retval{new expansion_t{kSourceIntermediate}};
    dcomplex_t *M = reinterpret_cast<dcomplex_t *>(views_.view_data(0));
    lap_m_to_i(M, retval->views_, 0); 
    return std::unique_ptr<expansion_t>(retval);
  }

  std::unique_ptr<expansion_t> I_to_I(Index s_index, Index t_index) const {
    ViewSet views{kTargetIntermediate}; 
    lap_i_to_i(s_index, t_index, views_, 0, views); 
    expansion_t *retval = new expansion_t{views};
    return std::unique_ptr<expansion_t>{retval};
  }

  std::unique_ptr<expansion_t> I_to_L(Index t_index) const {
    // t_index is the index of the child
    expansion_t *retval{new expansion_t{kTargetPrimary}};
    dcomplex_t *L = reinterpret_cast<dcomplex_t *>(retval->views_.view_data(0));
    lap_i_to_l(views_, 0, t_index, L); 
    return std::unique_ptr<expansion_t>(retval);
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
  }

  static void delete_table() { }

  static double compute_scale(Index index) {
    return builtin_laplace_table_->scale(index.level());
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
};

} // namespace dashmm

#endif // __AFMPB_RHS_H__
