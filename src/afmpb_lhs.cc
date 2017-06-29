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

  dcomplex_t z1{normal[0], -norma[1]}; 
  dcomplex_t z2{normal[0], normal[1]}; 

  // (0, 0) 
  M[midx(1, 1)] += z1 * q * scale / sqf[2]; 
  M[midx(1, 0)] -= q * scale * normal[2]; 

  legendre_Plm(p, ctheta, legendre.data()); 

  // (n, 0), n = 1, ..., p - 1
  for (int n = 1; n <= p - 1; ++n) {
    double cp = q * powers_r[n] * legendre_Plm[midx(n, 0)] * scale; 
    double cpz = cp * sqf[n + 2] / sqf[n] / 2; 
    M[midx(n + 1, 1)] += cpz * z1; 
    M[midx(n + 1, 0)] -= cp * (n + 1) * normal[2];
  }

  // (n, m), n = 1, ..., p - 1, m = 1, ..., n
  for (int n = 1; n <= p - 1; ++n) {
    for (int m = 1; m <= n; ++m) {
      double cp = q * powers_r[n] * sqf[n - m] / sqf[n + m] * scale *
        legendre[midx(n, m)];
      double cpz = cp * powers_ephi[m]; 
      double sr2 = sqf[n - m + 2] / sqf[n - m]; 

      if (m == 1) {
        cp = sr2 * real(cpz * z2);
        M[midx(n + 1, 0)] -= cp;
      } else {
        M[midx(n + 1, m - 1)] -= sr2 / 2.0 * cpz * z2;
      }

      sr2 = sqf[n + m + 2] / sqf[n + m] / 2.0; 
      M[midx(n + 1, m + 1)] += sr2 * cpz * z1; 
      M[midx(n + 1, m)] -= cpz * normal[2] * sqf[n + m + 1] * 
        sqf[n - m + 1] / sqf[n + m] / sqf[n - m];
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

  dcomplex_t z1{normal[0], -norma[1]}; 
  dcomplex_t z2{normal[0], normal[1]}; 

  // (1, 0) 
  L[midx(0, 0)] += q * scale * legendre[midx(1, 0)] * powers_r[1] * normal[2];
  
  // (n, 0), n = 2, ..., p
  for (int n = 2; n <= p; ++n) {
    double cp = q * scale * legendre[midx(n, 0)] * powers_r[n]; 
    doule cp1 = cp * sqf[n] / sqf[n - 2] / 2.0;
    L[midx(n - 1, 1)] += cp1 * z1; 
    L[midx(n - 1, 0)] += cp * n * normal[2];
  }

  // (n, m), n = 1, ..., p, m = 1, ...n
  for (int n = 1; n <= p; ++n) {
    // m = 1, ..., n - 2
    for (int m = 1; m <= n - 2; ++m) {
      double cp = q * scale * legendre[midx(n, m)] * powers_r[n] * 
        sqf[n - m] / sqf[n + m]; 
      dcomplex_t cpz = cp * powers_ephi[m]; 
      double sr2 = sqf[n + m] / sqf[n + m - 2]; 

      if (m == 1) {
        cp = sr2 * real(cpz * z2); 
        L[midx(n - 1, 0)] -= cp;
      } else {
        L[midx(n - 1, m - 1)] -= sr2 / 2.0 * cpz * z2;
      }
      
      sr2 = sqf[n - m] / sqf[n - m -2]; 
      L[midx(n - 1, m + 1)] += sr2 * cpz * z1; 
      L[midx(n - 1, m)] += cpz * normal[2] * sqf[n + m] / sqf[n + m -1] * 
        sqf[n - m] / sqf[n - m - 1];
    }

    // m = n - 1
    {
      m = n - 1; 
      double cp = q * scale * legendre[midx(n, m)] * powers_r[n] * 
        sqf[n - m] / sqf[n + m];
      dcomplex_t cpz = cp * powers_ephi[m]; 
      double sr2 = sqf[n + m] / sqf[n + m - 2]; 
      if (m == 1) {
        cp = sr2 * real(cpz * z2); 
        L[midx(n - 1, 0)] -= cp;
      } else {
        L[midx(n - 1, m - 1)] -= sr2 / 2.0 * cpz * z2;
      }

      L[midx(n - 1, m)] += normal[2] * cpz * sqf[n + m] / sqf[n + m - 1] * 
        sqf[n - m] / sqf[n - m - 1]; 
    }

    // m = n
    {
      m = n; 
      double cp = q * scale * legendre[midx(n, m)] * powers_r[n] * 
        sqf[n - m] / sqf[n + m];
      dcomplex_t cpz = cp * powers_ephi[m]; 
      double sr2 = sqf[n + m] / sqf[n + m - 2]; 

      if (m == 1) {
        cp = sr2 * real(cp * z2); 
        L[midx(n - 1, 0)] -= cp;
      } else {
        L[midx(n - 1, m - 1)] -= sr2 / 2.0 * cpz * z2;
      }
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







  // (1, 0) 
  {
    double cp = q * bessel[0] * lambda / scale; 
    M[midx(1, 1)] += cp / 2.0 * z1;
    M[midx(1, 0)] -= cp * normal[2]; 
  }

  // (n, 0), n = 1, ...., p - 1
  for (int n = 2; n <= p - 1; ++n) {
    double cp = q * lengendre[midx(n, 0)] * bessel[n] * lambda / (2 * n + 1); 
    dcomplex_t cpz = cp / 2.0 * z1; 
    cp *= normal[2]; 

    M[midx(n - 1, 1)] -= cpz * scale; 
    M[midx(n + 1, 1)] += cpz / scale; 
    M[midx(n - 1, 0)] -= cp * scale * n;
    M[midx(n + 1, 0)] -= cp / scale * (n + 1);
  }

  // (p, 0)
  {
    n = p; 
    double cp = q * legendre[midx(n, 0)] * bessel[n] * lambda / 
      (2 * n + 1) * scale;
    M[midx(n - 1, 1)] -= cp / 2.0 * z1;
    M[midx(n - 1, 0)] -= cp * n * normal[2];
  }

  // (1, 1) and (1, -1)
  {
    double cp = q * bessel[1] * legendre[midx(1, 1)] * lambda / 3.0; 
    dcomplex_t cpz = cp * powers_ephi[1]; 
    cp = 2.0 * real(cpz * z2);
    M[midx(0, 0)] += cp * scale;
    M[midx(2, 0)] -= cp / scale;
    cpz /= scale; 
    M[midx(2, 2)] += cpz / 2.0 * z1;
    M[midx(2, 1)] -= cpz * normal[2];
  }

  // (n, m), n = 2, ..., p - 1
  for (int n = 2; n <= p - 1; ++n) {
    // m = 1, ..., n - 2
    for (int m = 1; m <= n - 2; ++m) {
      double cp = q * bessel[n] * sqf[n - m] / sqf[n + m] * 
        legendre[midx(n, m)] * lambda / (2 * n + 1); 
      dcomplex_t cpz = cp * powers_ephi[m]; 
      dcomplex_t cpz2 = cpz * normal[2]; 
      
      if (m == 1) {
        cp = n * (n + 1) * real(cpz2); 
        M[midx(n - 1, 0)] += cp * scale; 
        M[midx(n + 1, 0)] -= cp / scale;
      } else {
        cpz2 /= 2.0; 
        M[midx(n - 1, m - 1)] += cpz2 * (n + m - 1) * (n + m) * scale; 
        M[midx(n + 1, m - 1)] -= cpz2 * (n - m + 1) * (n - m + 2) / scale;
      }

      cpz2 = cpz * z1 / 2.0; 
      M[midx(n - 1, m + 1)] -= cpz2 * scale; 
      M[midx(n + 1, m + 1)] += cpz2 / scale; 
      cpz2 = cpz * normal[2]; 
      M[midx(n - 1, m)] -= cpz2 * (n + m) * scale; 
      M[midx(n + 1, m)] -= cpz2 * (n - m + 1) / scale; 
    }

    // m = n - 1
    {
      int m = n - 1; 
      double cp = q * bessel[n] * sqf[n - m] / sqf[n + m] * lambda * 
        legendre[midx(n, m)] / (2 * n + 1); 
      dcomplex_t cpz = cp * powers_ephi[m]; 
      dcomplex_t cpz2 = cpz * z2; 

      if (m == 1) {
        cp = real(cpz2); 
        M[midx(n - 1, 0)] += cp * n * (n + 1) * scale; 
        M[midx(n + 1, 0)] -= cp * 6.0 / scale; 
      } else {
        cpz2 /= 2.0;
        M[midx(n - 1, m - 1)] += cpz2 * (n + m - 1) * (n + m) * scale; 
        M[midx(n + 1, m - 1)] -= cpz2 * 6.0 / scale;
      }

      cpz2 = cpz * z1 / 2.0; 
      M[midx(n + 1, m + 1)] += cpz2 / scale; 
      cpz2 = cpz * normal[2]; 
      M[midx(n - 1, m)] -= cpz2 * (n + m) * scale; 
      M[midx(n + 1, m)] -= cpz2 * 2.0 / scale;
    }

    // m = n
    {
      int m = n; 
      double cp = q * bessel[n] * sqf[n - m] / sqf[n + m] * lambda * 
        legendre[midx(n, m)] / (2 * n + 1); 
      dcomplex_t cpz = cp * powers_ephi[m]; 
      dcomplex_t cpz2 = cpz * z2 / 2.0;
      M[midx(n - 1, m - 1)] += cpz2 * (n + m - 1) * (n + m) * scale; 
      M[midx(n + 1, m - 1)] -= cpz2 * 2.0 / scale; 
      M[midx(n + 1, m + 1)] += cpz / (2.0 * scale) * z1; 
      M[midx(n + 1, m)] -= cpz * normal[2] / scale;
    }    
  }

  // (n, m), n = p
  {
    int n = p; 

    // (p, m), m = 1, ..., p - 2
    for (int m = 1; m <= n - 2; ++m) {
      double cp = q * bessel[n] * sqf[n - m] / sqf[n + m] * lambda * 
        legendre[midx(n, m)] / (2 * n + 1) * scale; 
      dcomplex_t cpz = cp * powers_ephi[m]; 
      dcomplex_t cpz2 = cpz * normal[2]; 

      if (m == 1) {
        cp = real(cpz2); 
        M[midx(n - 1, 0)] += cp * n * (n + 1); 
      } else {
        M[midx(n - 1, m - 1)] += cpz2 * (n + m - 1) * (n + m) / 2.0;
      } 

      M[midx(n - 1, m + 1)] -= cpz / 2.0 * z1; 
      M[midx(n - 1, m)] -= cpz * (n + m) * normal[2];
    }

    // m = n - 1
    {
      int m = n - 1; 
      double cp = q * bessel[n] * sqf[n - m] / sqf[n + m] * lambda * 
        legendre[midx(n, m)] / (2 * n + 1) * scale; 
      dcomplex_t cpz = cp * powers_ephi[m]; 

      M[midx(n - 1, m - 1)] += cpz * (n + m - 1) * (n + m) / 2.0 * z2; 
      M[midx(n - 1, m)] -= cpz * (n + m) * normal[2]; 
    }

    // m = n
    {
      int m = n; 
      double cp = q * bessel[n] * sqf[n - m] / sqf[n + m] * lambda * 
        legendre[midx(n, m)] / (2 * n + 1) * scale; 
      dcomplex_t cpz = cp * powers_ephi[m]; 

      M[midx(n - 1, m - 1)] += cpz * (n + m - 1) * (n + m) / 2.0 * z2; 
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

  // (0, 0)
  {
    double cp = q * bessel[0] * lambda * scale;
    L[midx(1, 1)] -= cp / 2.0 * z1; 
    L[midx(1, 0)] += cp * normal[2]; 
  }

  // (1, 0) 
  {
    doble cp = q * bessel[1] * lambda * legendre[midx(1, 0)] / 3.0; 
    L[midx(2, 1)] -= cp / 2.0 * scale * z1; 
    L[midx(0, 0)] += cp / scale * normal[2]; 
    L[midx(2, 0)] += cp * scale * 2.0 * normal[2];
  }

  // (n, 0), n = 2, ..., p - 1
  for (int n = 2; n <= p - 1; ++n) {
    double cp = q * bessel[n] * legendre[midx(n, 0)] * lambda / (2 * n + 1); 
    dcomplex_t cpz = cp / 2.0 * z1; 
    cp *= normal[2]; 
    
    L[midx(n - 1, 1)] += cpz / scale; 
    L[midx(n + 1, 1)] -= cpz * scale; 
    L[midx(n - 1, 0)] += cp / scale * n;
    L[midx(n + 1, 0)] += cp * scale * (n + 1); 
  }

  // (n, 0), n = p
  {
    int n = p; 
    double cp = q * legendre[midx(n, 0)] * bessel[n] * lambda / 
      (2 * n + 1) / scale;

    L[midx(n - 1, 1)] += cp / 2.0 * z1; 
    L[midx(n - 1, 0)] += cp * normal[2] * n;    
  }

  // (1, 1) 
  {
    double cp = q * bessel[1] * legendre[midx(1, 1)] * lambda / 3.0 / sqf[2];
    dcomplex_t cpz = cp * powers_ephi[1]; 
    cp = 2.0 * real(cpz * z2);
    L[midx(0, 0)] -= cp / scale; 
    L[midx(2, 0)] += cp * scale; 

    cpz *= scale; 
    L[midx(2, 2)] -= cpz / 2.0 * z1; 
    L[midx(2, 1)] += cpz * normal[2]; 
  }
  
  // (n, m), n = 2, ..., p - 1
  for (int n = 2; n <= p - 1; ++n) {
    // m = 1, ..., n - 2
    for (int m= 1; m <= n - 2; ++m) {
      double cp = q * bessel[n] * sqf[n - m] / sqf[n + m] * 
        legendre[midx(n, m)] * lambda / (2 * n + 1); 
      dcomplex_t cpz = cp * powers_ephi[m]; 
      dcomplex_t cpz2 = cpz * z2; 

      if (m == 1) {
        cp = n * (n + 1) * real(cpz2); 
        L[midx(n - 1, 0)] -= cp / scale; 
        L[midx(n + 1, 0)] += cp * scale; 
      } else {
        cpz2 /= 2.0; 
        L[midx(n - 1, m - 1)] -= cpz2 * (n + m - 1) * (n + m) / scale; 
        L[midx(n + 1, m - 1)] += cpz2 * (n - m + 1) * (n - m + 2) * scale;
      }

      cpz2 = cpz * z1 / 2.0; 
      L[midx(n - 1, m + 1)] += cpz2 / scale; 
      L[midx(n + 1, m + 1)] -= cpz2 * scale; 

      cpz2 *= normal[2]; 
      L[midx(n - 1, m)] += cpz2 * (n + m) / scale; 
      L[midx(n + 1, m)] += cpz2 * (n - m + 1) * scale; 
    }
  
    // m = n - 1
    {
      int m = n - 1; 
      double cp = q * bessel[n] * sqf[n - m] / sqf[n + m] * 
        legendre[midx(n, m)] * lambda / (2 * n + 1); 
      dcomplex_t cpz = cp * powers_ephi[m]; 
      dcomplex_t cpz2 = cpz * z2; 

      if (m == 1) {
        cp = real(cpz2); 
        L[midx(n - 1, 0)] -= cp * n * (n + 1) / scale; 
        L[midx(n + 1, 0)] += cp * 6.0  * scale;
      } else {
        cpz2 /= 2.0; 
        L[midx(n - 1, m - 1)] -= cpz2 * (n + m - 1) * (n + m) / scale; 
        L[midx(n + 1, m - 1)] += cpz2 * 6.0 * scale;
      }

      cpz2 = cpz * z1 / 2.0; 
      L[midx(n + 1, m + 1)] -= cpz2 * scale; 
      
      cpz2 = cpz * normal[2]; 
      L[midx(n - 1, m)] += cpz2 * (n + m) / scale; 
      L[midx(n + 1, m)] += cpz2 * 2.0 * scale;
    }

    // m = n
    {
      int m = n; 
      double cp = q * bessel[n] * sqf[n - m] / sqf[n + m] * 
        legendre[midx(n, m)] * lambda / (2 * n + 1); 
      dcomplex_t cpz = cp * powers_ephi[m]; 
      dcomplex_t cpz2 = cpz / 2.0 * z2; 
      
      L[midx(n - 1, m - 1)] -= cpz2 * (n + m - 1) * (n + m) / scale; 
      L[midx(n + 1, m - 1)] += cpz2 * 2.0 * scale; 
      L[midx(n + 1, m + 1)] -= cpz / 2.0 * scale * z1; 
    }
  }

  // n = p
  {
    int n = p;

    // m = 1, ..., n - 2
    for (int m = 1; m <= n - 2; ++m) {
      double cp = q * bessel[n] * sqf[n - m] / sqf[n + m] * 
        legendre[midx(n, m)] * lambda / (2 * n + 1) / scale; 
      dcomplex_t cpz = cp * powers_ephi[m]; 
      dcomplex_t cp2 = cpz * z2; 
      
      if (m == 1) {
        cp = real(cpz2); 
        L[midx(n - 1, 0)] -= cp * n * (n + 1);
      } else {
        L[midx(n - 1, m - 1)] -= cpz2 / 2.0 * (n + m - 1) * (n + m);
      }

      L[midx(n - 1, m + 1)] += cpz / 2.0 * z1; 
      L[midx(n - 1, m)] += cpz * (n + m) * normal[2];
    }

    // m = n - 1
    {
      int m = n - 1;
      double cp = q * bessel[n] * sqf[n - m] / sqf[n + m] * 
        legendre[midx(n, m)] * lambda / (2 * n + 1) / scale; 
      dcomplex_t cpz = cp * powers_ephi[m]; 

      L[midx(n - 1, m - 1)] -= cpz / 2.0 * (n + m - 1) * (n + m) * z2; 
      L[midx(n - 1, m)] += cpz * (n + m) * normal[2]; 
    }

    // m = n
    {
      int m = n; 
      double cp = q * bessel[n] * sqf[n - m] / sqf[n + m] * 
        legendre[midx(n, m)] * lambda / (2 * n + 1) / scale; 
      dcomplex_t cpz = cp * powers_ephi[m]; 

      L[midx(n - 1, m - 1)] -= cpz / 2.0 * (n + m - 1) * (n + m) * z2;
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
