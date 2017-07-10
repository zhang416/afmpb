#include "afmpb_lhs.h"

namespace dashmm {

std::unique_ptr<AFMPBTable> builtin_afmpb_table_; 

void dlap_s_to_m(Point dist, double q, double scale, Point normal, 
                 dcomplex_t *M) {
  int p = builtin_laplace_table_->p();
  const double *sqf = builtin_laplace_table_->sqf();

  std::vector<double> legendre((p + 1) * (p + 2) / 2); 
  std::vector<double> powers_r(p + 1); 
  std::vector<dcomplex_t> powers_ephi(p + 1); 

  powers_r[0] = 1.0; 
  powers_ephi[0] = dcomplex_t{1.0, 0.0}; 

  double proj = sqrt(dist.x() * dist.x() + dist.y() * dist.y());
  double r = dist.norm();
  double ctheta = (r <= 1e-14 ? 1.0 : dist.z() / r);

  // Compute exp(-i * phi) for the azimuthal angle phi
  dcomplex_t ephi = (proj / r <= 1e-14 ? dcomplex_t{1.0, 0.0} :
                     dcomplex_t{dist.x() / proj, -dist.y() / proj});

  // Compute powers of r
  r *= scale;
  for (int j = 1; j <= p; ++j) 
    powers_r[j] = powers_r[j - 1] * r;
  
  // Compute powers of exp(-i * phi)
  for (int j = 1; j <= p; ++j) 
    powers_ephi[j] = powers_ephi[j - 1] * ephi;

  // Compute legendre polynomial 
  legendre_Plm(p, ctheta, legendre.data()); 

  dcomplex_t z1{normal[0], -normal[1]}; 
  dcomplex_t z2{normal[0], normal[1]}; 

  for (int n = 1; n <= p - 1; ++n) {
    int m = 1; 
    M[midx(n + 1, m - 1)] -= real(powers_ephi[m] * z2) * 
      q * powers_r[n] * legendre[midx(n, m)] * scale * 
      sqf[n - m + 2] / sqf[n + m]; 
  }

  for (int n = 2; n <= p - 1; ++n) {
    dcomplex_t temp = q * powers_r[n] * scale * z2 / 2.0; 
    for (int m = 2; m <= n; ++m) {
      M[midx(n + 1, m - 1)] -= temp * legendre[midx(n, m)] * powers_ephi[m]
        * sqf[n - m + 2] / sqf[n + m]; 
    }
  }

  for (int n = 0; n <= p - 1; ++n) {
    double temp = q * powers_r[n] * scale * normal[2]; 
    for (int m = 0; m <= n; ++m) {
      M[midx(n + 1, m)] -= temp * legendre[midx(n, m)] * 
        powers_ephi[m] * sqf[n - m + 1] / sqf[n + m] * 
        sqf[n + m + 1] / sqf[n + m]; 
    }
  }

  for (int n = 0; n <= p - 1; ++n) {
    dcomplex_t temp = q * powers_r[n] * scale * z1 / 2.0; 
    for (int m = 0; m <= n; ++m) {
      M[midx(n + 1, m + 1)] += temp * legendre[midx(n, m)] * powers_ephi[m] * 
        sqf[n + m + 2] / sqf[n + m] * sqf[n - m] / sqf[n + m]; 
    }
  }
}


void dlap_s_to_l(Point dist, double q, double scale, Point normal, 
                 dcomplex_t *L) {
  int p = builtin_laplace_table_->p();
  const double *sqf = builtin_laplace_table_->sqf();

  std::vector<double> legendre((p + 1) * (p + 2) / 2); 
  std::vector<double> powers_r(p + 1); 
  std::vector<dcomplex_t> powers_ephi(p + 1); 
  powers_ephi[0] = dcomplex_t{1.0, 0.0};

  double proj = sqrt(dist.x() * dist.x() + dist.y() * dist.y());
  double r = dist.norm();

  // Compute cosine of the polar angle theta
  double ctheta = (r <= 1e-14 ? 1.0 : dist.z() / r);
  
  // Compute exp(-i * phi) for the azimuthal angle phi
  dcomplex_t ephi = (proj / r <= 1e-14 ? dcomplex_t{1.0, 0.0} :
                     dcomplex_t{dist.x() / proj, -dist.y() / proj});

  // Compute powers of 1 / r
  powers_r[0] = 1.0 / r;
  r *= scale;
  for (int j = 1; j <= p; ++j) 
    powers_r[j] = powers_r[j - 1] / r;
  
  // Compute powers of exp(-i * phi)
  for (int j = 1; j <= p; ++j) 
    powers_ephi[j] = powers_ephi[j - 1] * ephi;
  
  // compute local expansion L_n^m
  legendre_Plm(p, ctheta, legendre.data());

  dcomplex_t z1{normal[0], -normal[1]}; 
  dcomplex_t z2{normal[0], normal[1]}; 

  for (int n = 1; n <= p; ++n) {
    int m = 1; 
    L[midx(n - 1, m - 1)] -= real(powers_ephi[m] * z2) * 
      q * powers_r[n] * legendre[midx(n, m)] / scale * 
      sqf[n - m] / sqf[n + m - 2]; 
  }

  for (int n = 2; n <= p; ++n) {
    dcomplex_t temp = q * powers_r[n] / scale * z2 / 2.0; 
    for (int m = 2; m <= n; ++m) {
      L[midx(n - 1, m - 1)] -= temp * legendre[midx(n, m)] * powers_ephi[m] 
        * sqf[n - m] / sqf[n + m - 2]; 
    }
  }

  for (int n = 1; n <= p; ++n) {
    double temp = q * powers_r[n] / scale * normal[2];
    for (int m = 0; m <= n - 1; ++m) {
      L[midx(n - 1, m)] += temp * legendre[midx(n, m)] * powers_ephi[m] * 
        sqf[n - m] / sqf[n - m - 1] * sqf[n - m] / sqf[n + m - 1]; 
    }
  }

  for (int n = 2; n <= p; ++n) {
    dcomplex_t temp = q * powers_r[n] / scale * z1 / 2.0; 
    for (int m = 0; m <= n - 2; ++m) {
      L[midx(n - 1, m + 1)] += temp * legendre[midx(n, m)] * powers_ephi[m] 
        * sqf[n - m] / sqf[n + m] * sqf[n - m] / sqf[n - m - 2]; 
    }
  }
}

void dyuk_s_to_m(Point dist, double q, double scale, Point normal, 
                 dcomplex_t *M) {
  int p = builtin_yukawa_table_->p();
  const double *sqf = builtin_yukawa_table_->sqf();
  double lambda = builtin_yukawa_table_->lambda();

  std::vector<double> legendre((p + 1) * (p + 2) / 2); 
  std::vector<dcomplex_t> powers_ephi(p + 1); 
  std::vector<double> bessel(p + 1); 

  double proj = sqrt(dist.x() * dist.x() + dist.y() * dist.y());
  double r = dist.norm();
  
  // Compute cosine of the polar angle theta
  double ctheta = (r <= 1.0e-14 ? 1.0 : dist.z() / r);

  // Compute exp(-i * phi) for the azimuthal angle phi
  dcomplex_t ephi = (proj / r <= 1.0e-14 ? dcomplex_t{1.0, 0.0} :
                     dcomplex_t{dist.x() / proj, -dist.y() /proj});

  // Compute powers of exp(-i * phi)
  powers_ephi[0] = 1.0;
  for (int j = 1; j <= p; ++j) 
    powers_ephi[j] = powers_ephi[j - 1] * ephi;
  
  // Compute scaled modified spherical bessel function
  bessel_in_scaled(p, lambda * r, scale, bessel.data());

  // Compute legendre polynomial
  legendre_Plm(p, ctheta, legendre.data());
  
  dcomplex_t z1{normal[0], -normal[1]}; 
  dcomplex_t z2{normal[0], normal[1]}; 

  for (int n = 1; n <= p; ++n) {
    int m = 1; 
    dcomplex_t temp = q * bessel[n] * sqf[midx(n, m)] * 
      legendre[midx(n, m)] * lambda / (2.0 * n + 1) * powers_ephi[m]; 
    M[midx(n - 1, m - 1)] += ((double) n * (n + 1)) * scale * 
      real(z2 * temp); 
  }

  for (int n = 1; n <= p - 1; ++n) {
    int m = 1; 
    dcomplex_t temp = q * bessel[n] * sqf[midx(n, m)] * 
      legendre[midx(n, m)] * lambda / (2.0 * n + 1) * powers_ephi[m]; 
    M[midx(n + 1, m - 1)] -= ((double) n * (n + 1)) / scale * 
      real(z2 * temp);
  }

  for (int n = 2; n <= p; ++n) {
    for (int m = 2; m <= n; ++m) {
      dcomplex_t temp = q * bessel[n] * sqf[midx(n, m)] * 
        legendre[midx(n, m)] * lambda / (2.0 * n + 1) * powers_ephi[m]; 
      M[midx(n - 1, m - 1)] += temp * z2 / 2.0 * 
        ((double) (n + m - 1) * (n + m)) * scale; 
    }
  }

  for (int n = 1; n <= p; ++n) {
    for (int m = 0; m <= n - 1; ++m) {
      dcomplex_t temp = q * bessel[n] * sqf[midx(n, m)] * 
        legendre[midx(n, m)] * lambda / (2.0 * n + 1) * powers_ephi[m]; 
      M[midx(n - 1, m)] -= temp * normal[2] * ((double) (n + m)) * scale; 
    }
  }

  for (int n = 2; n <= p; ++n) {
    for (int m = 0; m <= n - 2; ++m) {
      dcomplex_t temp = q * bessel[n] * sqf[midx(n, m)] * 
        legendre[midx(n, m)] * lambda / (2.0 * n + 1) * powers_ephi[m]; 
      M[midx(n - 1, m + 1)] -= temp * z1 / 2.0 * scale; 
    }
  }

  for (int n = 2; n <= p - 1; ++n) {
    for (int m = 2; m <= n; ++m) {
      dcomplex_t temp = q * bessel[n] * sqf[midx(n, m)] * 
        legendre[midx(n, m)] * lambda / (2.0 * n + 1) * powers_ephi[m]; 
      M[midx(n + 1, m - 1)] -= temp * z2 / 2.0 * 
        ((double) (n - m + 1) * (n - m + 2)) / scale; 
    }
  }

  for (int n = 0; n <= p - 1; ++n) {
    for (int m = 0; m <= n; ++m) {
      dcomplex_t temp = q * bessel[n] * sqf[midx(n, m)] * 
        legendre[midx(n, m)] * lambda / (2.0 * n + 1) * powers_ephi[m]; 
      M[midx(n + 1, m)] -= temp * normal[2] * ((double) (n - m + 1)) / scale; 
    }
  }

  for (int n = 0; n <= p - 1; ++n) {
    for (int m = 0; m <= n; ++m) {
      dcomplex_t temp = q * bessel[n] * sqf[midx(n, m)] * 
        legendre[midx(n, m)] * lambda / (2.0 * n + 1) * powers_ephi[m]; 
      M[midx(n + 1, m + 1)] += temp * z1 / 2.0 / scale; 
    }
  }
}

void dyuk_s_to_l(Point dist, double q, double scale, Point normal, 
                 dcomplex_t *L) {
  int p = builtin_yukawa_table_->p();
  const double *sqf = builtin_yukawa_table_->sqf();
  double lambda = builtin_yukawa_table_->lambda();

  std::vector<double> legendre((p + 1) * (p + 2) / 2); 
  std::vector<double> bessel(p + 1); 
  std::vector<dcomplex_t> powers_ephi(p + 1); 

  double proj = sqrt(dist.x() * dist.x() + dist.y() * dist.y());
  double r = dist.norm();

  // Compute cosine of the polar angle theta
  double ctheta = (r <= 1.0e-14 ? 1.0 : dist.z() / r);
  
  // Compute exp(-i * phi) for the azimuthal angle phi
  dcomplex_t ephi = (proj / r <= 1.0e-14 ? dcomplex_t{1.0, 0.0} :
                     dcomplex_t{dist.x() / proj, -dist.y() / proj});
  
  // Compute powers of exp(-i * phi)
  powers_ephi[0] = 1.0;
  for (int j = 1; j <= p; j++) 
    powers_ephi[j] = powers_ephi[j - 1] * ephi;
      
  // Compute scaled modified spherical bessel function
  bessel_kn_scaled(p, lambda * r, scale, bessel.data());

  // Compute legendre polynomial
  legendre_Plm(p, ctheta, legendre.data());

  dcomplex_t z1{normal[0], -normal[1]}; 
  dcomplex_t z2{normal[0], normal[1]}; 

  for (int n = 1; n <= p; ++n) {
    int m = 1; 
    dcomplex_t temp = q * bessel[n] * sqf[midx(n, m)] * legendre[midx(n, m)] * 
      lambda / (2.0 * n + 1) * powers_ephi[m]; 
    L[midx(n - 1, m - 1)] -= ((double) n * (n + 1)) * real(temp * z2) / scale; 
  }

  for (int n = 1; n <= p - 1; ++n) {
    int m = 1; 
    dcomplex_t temp = q * bessel[n] * sqf[midx(n, m)] * legendre[midx(n, m)] * 
      lambda / (2.0 * n + 1) * powers_ephi[m]; 
    L[midx(n + 1, m - 1)] += ((double) n * (n + 1)) * real(temp * z2) * scale; 
  }

  for (int n = 2; n <= p; ++n) {
    for (int m = 2; m <= n; ++m) {
      dcomplex_t temp = q * bessel[n] * sqf[midx(n, m)] * legendre[midx(n, m)] 
        * lambda / (2.0 * n + 1) * powers_ephi[m]; 
      L[midx(n - 1, m - 1)] -= temp * z2 / 2.0 * 
        ((double) (n + m - 1) * (n + m)) / scale; 
    }
  }

  for (int n = 1; n <= p; ++n) {
    for (int m = 0; m <= n - 1; ++m) {
      dcomplex_t temp = q * bessel[n] * sqf[midx(n, m)] * legendre[midx(n, m)] 
        * lambda / (2.0 * n + 1) * powers_ephi[m]; 
      L[midx(n - 1, m)] += temp * normal[2] * ((double) (n + m)) / scale; 
    }
  }

  for (int n = 2; n <= p; ++n) {
    for (int m = 0; m <= n - 2; ++m) {
      dcomplex_t temp = q * bessel[n] * sqf[midx(n, m)] * legendre[midx(n, m)] 
        * lambda / (2.0 * n + 1) * powers_ephi[m]; 
      L[midx(n - 1, m + 1)] += temp * z1 / 2.0 / scale; 
    }
  }

  for (int n = 2; n <= p - 1; ++n) {
    for (int m = 2; m <= n; ++m) {
      dcomplex_t temp = q * bessel[n] * sqf[midx(n, m)] * legendre[midx(n, m)] 
        * lambda / (2.0 * n + 1) * powers_ephi[m]; 
      L[midx(n + 1, m - 1)] += temp * z2 / 2.0 * 
        ((double) (n - m + 1) * (n - m + 2)) * scale; 
    }
  }

  for (int n = 0; n <= p - 1; ++n) {
    for (int m = 0; m <= n; ++m) {
      dcomplex_t temp = q * bessel[n] * sqf[midx(n, m)] * legendre[midx(n, m)] 
        * lambda / (2.0 * n + 1) * powers_ephi[m]; 
      L[midx(n + 1, m)] += temp * normal[2] * ((double) (n - m + 1)) * scale; 
    }
  }

  for (int n = 0; n <= p - 1; ++n) {
    for (int m = 0; m <= n; ++m) {
      dcomplex_t temp = q * bessel[n] * sqf[midx(n, m)] * legendre[midx(n, m)] 
        * lambda / (2.0 * n + 1) * powers_ephi[m]; 
      L[midx(n + 1, m + 1)] -= temp * z1 / 2.0 * scale; 
    }
  }
}

void update_afmpb_table(double dielectric, double cut1, double cut2, 
                        double sigma) {
  if (builtin_afmpb_table_ == nullptr) 
    builtin_afmpb_table_ = std::unique_ptr<AFMPBTable>
      {new AFMPBTable{dielectric, cut1, cut2, sigma}}; 
}

} // namespace dashmm 
