#include <algorithm>
#include <cstdio>
#include <stdexcept>
#include <getopt.h>
#include <hpx/hpx.h>
#include "afmpb.h"

namespace afmpb {

void usage(char *program) {
  fprintf(stdout, "usage: %s (--pqr-file=FILE)\n"
          "  (--mesh-format=0 --mesh-density=num --probe-radius=num | \n"
          "   --mesh-format=[1|2] --mesh-file=FILE)\n"
          "  [--dielectric-interior=num         [default:            2.0]]\n"
          "  [--dielectric-exterior=num         [default:           80.0]]\n"
          "  [--ion-concentration=num           [default:          150.0]]\n"
          "  [--temperature=num                 [default:          300.0]]\n"
          "  [--surface-tension=num             [default:          0.005]]\n"
          "  [--pressure=num                    [default:          0.035]]\n"
          "  [--log-file=FILE                   [default:     output.txt]]\n"
          "  [--surface-potential-file=FILE     [default:  potential.txt]]\n"
          "  [--accuracy=num                    [default:              3]]\n"
          , program);
}

std::unique_ptr<Configuration> init(int argc, char **argv) {
  if (HPX_SUCCESS != hpx_init(&argc, &argv)) 
    return std::unique_ptr<Configuration>{nullptr}; 

  // Parse command line arguments 
  static struct option long_options[] = {
    {"pqr-file", required_argument, 0, 'q'}, 
    {"mesh-file", required_argument, 0, 'm'}, 
    {"log-file", required_argument, 0, 'l'}, 
    {"potential-file", required_argument, 0, 's'}, 
    {"mesh-format", required_argument, 0, 'f'}, 
    {"mesh-density", required_argument, 0, 'd'}, 
    {"probe-radius", required_argument, 0, 'r'}, 
    {"dielectric-interior", required_argument, 0, 'i'}, 
    {"dielectric-exterior", required_argument, 0, 'e'}, 
    {"ion-concentration", required_argument, 0, 'c'}, 
    {"temperature", required_argument, 0, 't'}, 
    {"surface-tension", required_argument, 0, 'g'}, 
    {"pressure", required_argument, 0, 'p'}, 
    {"accuracy", required_argument, 0, 'a'}, 
    {"help", no_argument, 0, 'h'}, 
    {0, 0, 0, 0}
  };

  Configuration *p = new Configuration; 

  int opt = 0; 
  int long_index = 0; 
  while ((opt = getopt_long(argc, argv, "q:m:l:s:f:d:r:i:e:c:t:g:p:a:h", 
                            long_options, &long_index)) != -1) {
    switch (opt) {
    case 'q':
      p->pqr_file = optarg; 
      break;
    case 'm':
      p->mesh_file = optarg;
      break;
    case 'l':
      p->log_file = optarg;
      break;
    case 's':
      p->potential_file = optarg;
      break;
    case 'f':
      p->mesh_format = atoi(optarg); 
      break;
    case 'd':
      p->mesh_density = atof(optarg); 
      break;
    case 'r':
      p->probe_radius = atof(optarg); 
      break;
    case 'i':
      p->dielectric_interior = atof(optarg); 
      break;
    case 'e':
      p->dielectric_exterior = atof(optarg); 
      break;
    case 'c':
      p->ion_concentration = atof(optarg); 
      break;
    case 't':
      p->temperature = atof(optarg); 
      break;
    case 'g':
      p->surface_tension = atof(optarg); 
      break;
    case 'p':
      p->pressure = atof(optarg); 
      break;
    case 'a':
      p->accuracy = atoi(optarg); 
      break; 
    case 'h':
    case '?':
      usage(argv[0]); 
      delete p;
      return std::unique_ptr<Configuration>{nullptr}; 
    }
  }

  // Check if the command line arguments are valid 
  if (p->pqr_file.empty() || (p->mesh_format && p->mesh_file.empty())) {
    usage(argv[0]); 
    delete p; 
    p = nullptr; 
  } 

  return std::unique_ptr<Configuration>{p}; 
}
 
int finalize() {
  hpx_finalize(); 
  return 0;
}

AFMPB::AFMPB(std::unique_ptr<Configuration> p) {
  processPQRFile(p->pqr_file); 
  log_.open(p->log_file); 
  potential_.open(p->potential_file); 
  if (!p->mesh_format) 
    mesh_.open(p->mesh_file); 

  mesh_format_ = p->mesh_format; 
  mesh_density_ = p->mesh_density; 
  probe_radius_ = p->probe_radius; 
  dielectric_ = p->dielectric_exterior / p->dielectric_interior; 
  kap_ = sqrt(2.528639884 * std::max(p->ion_concentration, 1e-10) / 
              dielectric_ / p->temperature); 
  surface_tension_ = p->surface_tension; 
  pressure_ = p->pressure; 
  accuracy_ = p->accuracy; 

  if (mesh_format_ == 0) {
    cut1_ = 0.1; 
    cut2_ = 0.0001; 
  } else if (mesh_format_ == 1) {
    cut1_ = 0.4; 
    cut2_ = 0.4;
  } else {
    cut1_ = 0.0; 
    cut2_ = 0.3;
  } 
  sigma_ = 0.001; 

  ngauss_per_element_ = (mesh_format_ ? 7 : 1); 
}

void AFMPB::processPQRFile(std::string &pqr_file) {
 char buffer[200]; 
  sprintf(buffer, "grep \'ATOM\\|HETATM\' %s | cut -c 30- > %s; wc -l < %s", 
          pqr_file.c_str(), "processed.txt", "processed.txt"); 
  std::shared_ptr<FILE> pipe(popen(buffer, "r"), pclose); 
  if (!pipe) throw std::runtime_error("popen() failed!"); 
  
  if (fgets(buffer, sizeof(buffer), pipe.get()) != nullptr) 
    sscanf(buffer, "%d", &natoms_); 

  pqr_.open("processed.txt"); 
  std::cout << natoms_ << "\n";
}

std::vector<Atom> AFMPB::readAtoms() {
  std::vector<Atom> molecule(natoms_); 
  int i = 0; 
  double x, y, z, q, r;
  double probe_radius = (!mesh_format_ ? probe_radius_ : 0.0); 
  while (pqr_ >> x >> y >> z >> q >> r) {
    molecule[i].position = dashmm::Point{x, y, z}; 
    molecule[i].charge = q; 
    molecule[i].radius = r + probe_radius; 
    i++; 
  }
  return molecule; 
}

void AFMPB::generateMesh(int s, std::vector<Atom> &molecule, 
                         std::vector<Node> &nodes, 
                         std::vector<GNode> &gauss) {
  using namespace dashmm; 
  const int max_polar_intervals = 600; 
  const int min_polar_intervals = 5; 

  const double scale = 1.0 / 0.985; 
  Atom &S = molecule[s]; 
  Point &center = S.position; 
  double radius = S.radius; 
  std::vector<int> intersected; 

  // Find the list of intersecting atoms
  for (int t = 0; t < natoms_; ++t) {
    if (s == t) 
      continue; 

    Point dist = point_sub(center, molecule[t].position); 
    if (dist.norm() < radius + molecule[t].radius)
      intersected.push_back(t);
  }

  // Number of intervals used to divide polar angle theta
  int n = radius * sqrt(2 * M_PI * mesh_density_); 
  n = std::min(max_polar_intervals, std::max(n, min_polar_intervals)); 

  // Number of intervals used to divide azimuthal angle beta
  int m = 2 * n; 

  double dbeta = 2 * M_PI / m; 
  double area = 4 * radius * M_PI / n / m * scale; 

  std::unique_ptr<double []> theta{new double[n + 1]}; 
  std::unique_ptr<double []> thetaBar{new double[n]}; 
  std::unique_ptr<double []> beta{new double[m]}; 

  theta[0] = -M_PI / 2; 
  for (int i = 1; i <= n; ++i) 
    theta[i] = asin(-1 + 2 * i / n); 

  thetaBar[0] = (theta[0] + 2 * theta[1]) / 3; 
  for (int i = 1; i <= n - 2; ++i) 
    thetaBar[i] = (theta[i] + theta[i + 1]) / 2; 
  thetaBar[n - 1] = (theta[n - 1] + 2 * theta[n]) / 3; 

  for (int i = 0; i < m; ++i) 
    beta[i] = i * dbeta; 

  // Discretize the surface of S 
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      // (ox, oy, oz) is the outer normal direction
      double ox = cos(thetaBar[i]) * cos(beta[j]); 
      double oy = cos(thetaBar[i]) * sin(beta[j]); 
      double oz = sin(thetaBar[i]); 

      Point p{center.x() + radius * ox, center.y() + radius * oy, 
          center.z() + radius * oz}; 
      
      // Check if p is buried inside any of the intersecting atoms 
      bool buried = false; 
      for (auto t : intersected) {
        Point dist = point_sub(p, molecule[t].position); 
        if (dist.norm() < molecule[t].radius) {
          buried = true; 
          break;
        }
      }

      if (buried) 
        continue; 

      area_ += area / scale; // total mesh surface area

      GNode gnode; 
      gnode.position = p; 
      gnode.normal = Point{ox, oy, oz}; 
      gnode.index = gauss.size(); 
      gauss.push_back(gnode); 

      Node node; 
      node.position = p; 
      node.normal_o = Point{ox, oy, oz}; 
      node.area = area; 

      double p1 = area / 4 * (sin(theta[i + 1]) - sin(thetaBar[i])) * n; 
      double p2 = p1; 
      double p3 = area / 4 * (sin(thetaBar[i]) - sin(theta[i])) * n; 
      double p4 = p3; 

      double c1 = cos((thetaBar[i] + theta[i + 1]) / 2); 
      double s1 = sin((thetaBar[i] + theta[i + 1]) / 2); 
      double c2 = cos((thetaBar[i] + theta[i]) / 2); 
      double s2 = sin((thetaBar[i] + theta[i]) / 2); 
      double c3 = cos(beta[j] - M_PI / m); 
      double s3 = sin(beta[j] - M_PI / m); 
      double c4 = cos(beta[j] + M_PI / m); 
      double s4 = sin(beta[j] + M_PI / m); 

      double v11 = c1 * c3; 
      double v21 = c1 * s3; 
      double v31 = s1; 
      
      double v12 = c1 * c4; 
      double v22 = c1 * s4;
      double v32 = s1; 

      double v13 = c2 * c3; 
      double v23 = c2 * s3; 
      double v33 = s2; 
      
      double v14 = c2 * c4; 
      double v24 = c2 * s4; 
      double v34 = s2; 
      
      // (ix, iy, iz) is the inner normal 
      double ix = p1 * v11 + p2 * v12 + p3 * v13 + p4 * v14; 
      double iy = p1 * v21 + p2 * v22 + p3 * v23 + p4 * v24; 
      double iz = p1 * v31 + p2 * v32 + p3 * v33 + p4 * v34; 
      double norm = sqrt(ix * ix + iy * iy + iz * iz); 

      node.projected = norm; 
      node.normal_i = Point{ix / norm, iy / norm, iz / norm}; 

      volume_ += point_dot(node.normal_i, node.position) / 3 * area / scale; 

      node.patch[0].position = Point{center.x() + radius * v11, 
                                     center.y() + radius * v21, 
                                     center.z() + radius * v31}; 
      node.patch[0].normal = Point{v11, v21, v31}; 
      node.patch[0].weight = p1; 

      node.patch[1].position = Point{center.x() + radius * v12, 
                                     center.y() + radius * v22, 
                                     center.z() + radius * v32}; 
      node.patch[1].normal = Point{v12, v22, v32}; 
      node.patch[1].weight = p1; 

      node.patch[2].position = Point{center.x() + radius * v13, 
                                     center.y() + radius * v23, 
                                     center.z() + radius * v33}; 
      node.patch[2].normal = Point{v13, v23, v33}; 
      node.patch[2].weight = p3; 

      node.patch[3].position = Point{center.x() + radius * v14, 
                                     center.y() + radius * v24, 
                                     center.z() + radius * v34}; 
      node.patch[3].normal = Point{v14, v24, v34}; 
      node.patch[3].weight = p3;

      node.index = nodes.size(); 
      nodes.push_back(node); 
    }
  }
}

void AFMPB::readMesh(std::vector<Node> &nodes) {
  using namespace dashmm; 
  int nnodes, nelements; 
  mesh_ >> nnodes >> nelements; 

  nodes.resize(nnodes); 
  elements_.resize(nelements); 
  
  if (mesh_format_ == 1) { 
    // MSMS format, nodes are indexed from 1
    for (int i = 0; i < nnodes; ++i) {
      double x, y, z, nx, ny, nz, u1, u2, u3; 
      mesh_ >> x >> y >> z >> nx >> ny >> nz >> u1 >> u2 >> u3; 
      nodes[i].index = i; 
      nodes[i].position = Point{x, y, z}; 
      nodes[i].normal_o = Point{nx, ny, nz}; 
    }

    for (int i = 0; i < nelements; ++i) {
      int i1, i2, i3, u1, u2; 
      mesh_ >> i1 >> i2 >> i3 >> u1 >> u2; 
      elements_[i].index = i; 
      elements_[i].nodes.push_back(i1 - 1); 
      elements_[i].nodes.push_back(i2 - 1); 
      elements_[i].nodes.push_back(i3 - 1);
    }    
  } else { 
    // OFF format, nodes are index from 0 
    for (int i = 0; i < nnodes; ++i) {
      double x, y, z; 
      mesh_ >> x >> y >> z; 
      nodes[i].index = i; 
      nodes[i].position = dashmm::Point{x, y, z};
    }

    for (int i = 0; i < nelements; ++i) {
      int u, i1, i2, i3; 
      mesh_ >> u >> i1 >> i2 >> i3; 
      elements_[i].index = i; 
      elements_[i].nodes.push_back(i1); 
      elements_[i].nodes.push_back(i2); 
      elements_[i].nodes.push_back(i3); 
    }
  }

  // Make sure the normal direction is the outer direction 
  double volume = 0; 
  for (int i = 0; i < nelements; ++i) {
    int i1 = elements_[i].nodes[0]; 
    int i2 = elements_[i].nodes[1]; 
    int i3 = elements_[i].nodes[2]; 
    Node &n1 = nodes[i1]; 
    Node &n2 = nodes[i2]; 
    Node &n3 = nodes[i3]; 
    volume += tetrahedronVolume(n1.position, n2.position, n3.position); 
  } 

  if (volume < 0) {
    for (int i = 0; i < nelements; ++i) 
      std::swap(elements_[i].nodes[0], elements_[i].nodes[2]); 
  }
}

double AFMPB::tetrahedronVolume(dashmm::Point &A, dashmm::Point &B, 
                                dashmm::Point &C) {
  double x1 = A.x(); 
  double y1 = A.y(); 
  double z1 = A.z(); 
  double x2 = B.x(); 
  double y2 = B.y(); 
  double z2 = B.z(); 
  double x3 = C.x(); 
  double y3 = C.y(); 
  double z3 = C.z(); 
  double x21 = x2 - x1; 
  double y21 = y2 - y1; 
  double z21 = z2 - z1; 
  double x31 = x3 - x1; 
  double y31 = y3 - y1;
  double z31 = z3 - z1; 
  double a = y21 * z31 - y31 * z21; 
  double b = z21 * x31 - z31 * x21; 
  double c = x21 * y31 - x31 * y21; 

  return ((x1 + x2 + x3) * a + (y1 + y2 + y3) * b + (z1 + z2 + z3) * c) / 18.0;
}

void AFMPB::setup() {
  auto molecule = readAtoms(); 

  std::vector<Node> nodes; 
  std::vector<GNode> gauss; 

  if (!mesh_format_) {
    for (int i = 0; i < natoms_; ++i) 
      generateMesh(i, molecule, nodes, gauss);
  } else {
    readMesh(nodes);
  }

  int err = atoms_.allocate(natoms_); 
  assert(err == dashmm::kSuccess); 
  err = atoms_.put(0, natoms_, molecule.data()); 
  assert(err == dashmm::kSuccess); 
}


} // namespace afmpb
