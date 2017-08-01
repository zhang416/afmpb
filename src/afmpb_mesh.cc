#include "afmpb.h"

namespace afmpb {

void AFMPB::readAtoms(std::vector<Atom> &molecule) {
  pqr_ >> natoms_; 

  molecule.resize(natoms_); 
  double x, y, z, q, r; 
  double probe_radius = (!mesh_format_ ? probe_radius_ : 0.0); 
  int i = 0; 
  while (pqr_ >> x >> y >> z >> q >> r) {
    molecule[i].position = dashmm::Point{x, y, z};
    molecule[i].charge = q; 
    molecule[i].radius = r + probe_radius; 
    i++;
  }
}

void AFMPB::generateMesh(int s, const std::vector<Atom> &molecule, //const Atom *molecule, 
                         std::vector<Node> &nodes, 
                         std::vector<GNode> &gauss) {
  using namespace dashmm; 
  const int max_polar_intervals = 600; 
  const int min_polar_intervals = 5; 

  const double scale = 1.0 / 0.985; 
  const Atom &S = molecule[s]; 
  const Point &center = S.position; 
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
  int n = radius * sqrt(2.0 * M_PI * mesh_density_); 
  n = std::min(max_polar_intervals, std::max(n, min_polar_intervals)); 

  // Number of intervals used to divide azimuthal angle beta
  int m = 2 * n; 

  double dbeta = 2.0 * M_PI / m; 
  double area = 4.0 * radius * radius * M_PI / n / m * scale; 

  std::unique_ptr<double []> theta{new double[n + 1]}; 
  std::unique_ptr<double []> thetaBar{new double[n]}; 
  std::unique_ptr<double []> beta{new double[m]}; 

  theta[0] = -M_PI / 2; 
  for (int i = 1; i <= n; ++i) 
    theta[i] = asin(-1 + 2.0 * i / n); 

  thetaBar[0] = (theta[0] + 2 * theta[1]) / 3; 
  for (int i = 1; i <= n - 2; ++i) 
    thetaBar[i] = (theta[i] + theta[i + 1]) / 2; 
  thetaBar[n - 1] = (2 * theta[n - 1] + theta[n]) / 3; 

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

      gauss.emplace_back(gauss.size(), p, Point{ox, oy, oz}); 

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

      // Patch 1
      double p1x = center.x() + radius * v11; 
      double p1y = center.y() + radius * v21;
      double p1z = center.z() + radius * v31; 
      node.patch.emplace_back(Point{p1x, p1y, p1z}, Point{v11, v21, v31}, p1);
      
      // Patch 2
      double p2x = center.x() + radius * v12; 
      double p2y = center.y() + radius * v22;
      double p2z = center.z() + radius * v32; 
      node.patch.emplace_back(Point{p2x, p2y, p2z}, Point{v12, v22, v32}, p1); 

      // Patch 3
      double p3x = center.x() + radius * v13;
      double p3y = center.y() + radius * v23;
      double p3z = center.z() + radius * v33;
      node.patch.emplace_back(Point{p3x, p3y, p3z}, Point{v13, v23, v33}, p3); 

      // Patch 4 
      double p4x = center.x() + radius * v14;
      double p4y = center.y() + radius * v24;
      double p4z = center.z() + radius * v34;
      node.patch.emplace_back(Point{p4x, p4y, p4z}, Point{v14, v24, v34}, p3); 

      node.index = nodes.size(); 
      node.gmres.resize((restart_  + 1) * 2); 
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

    double x1 = n1.position.x(); 
    double y1 = n1.position.y(); 
    double z1 = n1.position.z(); 
    double x2 = n2.position.x(); 
    double y2 = n2.position.y(); 
    double z2 = n2.position.z(); 
    double x3 = n3.position.x(); 
    double y3 = n3.position.y(); 
    double z3 = n3.position.z(); 
    double x21 = x2 - x1; 
    double y21 = y2 - y1; 
    double z21 = z2 - z1; 
    double x31 = x3 - x1; 
    double y31 = y3 - y1;
    double z31 = z3 - z1; 
    double a = y21 * z31 - y31 * z21; 
    double b = z21 * x31 - z31 * x21; 
    double c = x21 * y31 - x31 * y21; 

    volume += ((x1 + x2 + x3) * a + (y1 + y2 + y3) * b + 
               (z1 + z2 + z3) * c) / 18.0;
  } 

  if (volume < 0) {
    for (int i = 0; i < nelements; ++i) 
      std::swap(elements_[i].nodes[0], elements_[i].nodes[2]); 
  }
}

void AFMPB::removeIsolatedNodes(std::vector<Node> &nodes) {
  std::vector<int> nnbr(nodes.size()); 
  std::vector<int> isolated; 

  // Count how many elements each node belongs to 
  for (auto && e : elements_) {
    nnbr[e.nodes[0]]++; 
    nnbr[e.nodes[1]]++; 
    nnbr[e.nodes[2]]++; 
  }

  for (auto && n : nodes) {
    if (!nnbr[n.index])
      isolated.push_back(n.index);
  }

  // For each node, count the number of isolated nodes whose indices are smaller
  for (auto && n : nodes) {
    int count = 0; 
    for (auto && i : isolated) {
      if (n.index > i) {
        count++;
      } else {
        n.index -= count;
        break;
      }
    }
  }

  // Update the node indices for each element 
  for (auto && e : elements_) {
    int i1 = e.nodes[0]; // old indices
    int i2 = e.nodes[1]; 
    int i3 = e.nodes[2]; 

    e.nodes[0] = nodes[i1].index; // new indices
    e.nodes[1] = nodes[i2].index; 
    e.nodes[2] = nodes[i3].index;
  }

  // Remove isolated nodes from the storage 
  for (int i = 0; i < nodes.size(); ++i) {
    int new_index = nodes[i].index; 
    if (new_index == i) 
      continue; 
    nodes[new_index] = nodes[i]; 
  }

  int updated_count = nodes.size() - isolated.size(); 
  nodes.resize(updated_count); 

  // Resize the gmres buffer for each remained node 
  for (int i = 0; i < nodes.size(); ++i) 
    nodes[i].gmres.resize((restart_ + 1) * 2); 
}

void AFMPB::processElementGeometry(std::vector<Node> &nodes) {
  using namespace dashmm; 
  for (auto &&e : elements_) {
    int i1 = e.nodes[0]; 
    int i2 = e.nodes[1]; 
    int i3 = e.nodes[2]; 
    Node &n1 = nodes[i1]; 
    Node &n2 = nodes[i2]; 
    Node &n3 = nodes[i3]; 

    double x1 = n1.position.x(); 
    double y1 = n1.position.y(); 
    double z1 = n1.position.z(); 
    double x2 = n2.position.x(); 
    double y2 = n2.position.y(); 
    double z2 = n2.position.z(); 
    double x3 = n3.position.x(); 
    double y3 = n3.position.y(); 
    double z3 = n3.position.z(); 
    double x21 = x2 - x1; 
    double y21 = y2 - y1; 
    double z21 = z2 - z1; 
    double x31 = x3 - x1; 
    double y31 = y3 - y1;
    double z31 = z3 - z1; 
    double x32 = x3 - x2; 
    double y32 = y3 - y2; 
    double z32 = z3 - z2; 

    double a = y21 * z31 - y31 * z21; 
    double b = z21 * x31 - z31 * x21; 
    double c = x21 * y31 - x31 * y21; 

    // Set the normal direction of the element 
    e.normal = Point{a, b, c}; 

    // Compute the area of the element 
    double area = e.normal.norm() / 2; 
    e.area = area; 
    
    // Update the total surface area of the molecule 
    area_ += area; 

    // Update the volume of the molecule 
    volume_ += ((x1 + x2 + x3) * a + (y1 + y2 + y3) * b + 
                (z1 + z2 + z3) * c) / 18.0;

    if (area > 1e-7) {
      // Normalize the normal 
      e.normal = e.normal.scale(0.5 / area); 

      if (mesh_format_ == 2) { // OFF format 
        double s31 = sqrt(x31 * x31 + y31 * y31 + z31 * z31); 
        double s21 = sqrt(x21 * x21 + y21 * y21 + z21 * z21); 
        double s32 = sqrt(x32 * x32 + y32 * y32 + z32 * z32); 

        n1.normal_o = point_add(n1.normal_o, e.normal.scale(1.0 / s31 / s21)); 
        n2.normal_o = point_add(n2.normal_o, e.normal.scale(1.0 / s32 / s21)); 
        n3.normal_o = point_add(n3.normal_o, e.normal.scale(1.0 / s31 / s32)); 
      }
    }        
      
    // Update the inner normal of the nodes
    n1.normal_i = point_add(n1.normal_i, e.normal.scale(area / 3)); 
    n2.normal_i = point_add(n2.normal_i, e.normal.scale(area / 3)); 
    n3.normal_i = point_add(n3.normal_i, e.normal.scale(area / 3)); 

    n1.area += area / 3; 
    n2.area += area / 3; 
    n3.area += area / 3; 

    // Update the patch of the nodes 
    Point m12 = Point{(x1 + x2) / 2, (y1 + y2) / 2, (z1 + z2) / 2}; 
    Point m13 = Point{(x1 + x3) / 2, (y1 + y3) / 2, (z1 + z3) / 2}; 
    Point m23 = Point{(x2 + x3) / 2, (y2 + y3) / 2, (z2 + z3) / 2}; 
    Point O = Point{(x1 + x2 + x3) / 3, (y1 + y2 + y3) / 3, (z1 + z2 + z3) / 3}; 

    // Patch for n1 
    {
      double p1x = (x1 + m12.x() + m13.x()) / 3;
      double p1y = (y1 + m12.y() + m13.y()) / 3; 
      double p1z = (z1 + m12.z() + m13.z()) / 3; 

      double p2x = (m12.x() + m13.x() + O.x()) / 3; 
      double p2y = (m12.y() + m13.y() + O.y()) / 3; 
      double p2z = (m12.z() + m13.z() + O.z()) / 3; 

      Patch p1{Point{p1x, p1y, p1z}, e.normal, area / 4.0}; 
      Patch p2{Point{p2x, p2y, p2z}, e.normal, area / 12.0}; 

      n1.patch.push_back(p1); 
      n1.patch.push_back(p2);
    }
    
    // Patch for n2 
    {
      double p1x = (x2 + m12.x() + m23.x()) / 3; 
      double p1y = (y2 + m12.y() + m23.y()) / 3; 
      double p1z = (z2 + m12.z() + m23.z()) / 3; 

      double p2x = (m12.x() + m23.x() + O.x()) / 3;
      double p2y = (m12.y() + m23.y() + O.y()) / 3; 
      double p2z = (m12.z() + m23.z() + O.z()) / 3; 

      Patch p1{Point{p1x, p1y, p1z}, e.normal, area / 4.0}; 
      Patch p2{Point{p2x, p2y, p2z}, e.normal, area / 12.0};

      n2.patch.push_back(p1); 
      n2.patch.push_back(p2);
    }
          

    // Patch for n3
    {
      double p1x = (x3 + m23.x() + m13.x()) / 3; 
      double p1y = (y3 + m23.y() + m13.y()) / 3; 
      double p1z = (z3 + m23.z() + m13.z()) / 3; 
      
      double p2x = (m23.x() + m13.x() + O.x()) / 3; 
      double p2y = (m23.y() + m13.y() + O.y()) / 3; 
      double p2z = (m23.z() + m13.z() + O.z()) / 3; 

      Patch p1{Point{p1x, p1y, p1z}, e.normal, area / 4.0}; 
      Patch p2{Point{p2x, p2y, p2z}, e.normal, area / 12.0}; 

      n3.patch.push_back(p1);
      n3.patch.push_back(p2); 
    }
  }

  if (mesh_format_ == 2) {
    // Normalize outer normal of each node 
    for (auto && n : nodes) {
      double norm = n.normal_o.norm(); 
      n.normal_o = n.normal_o.scale(1.0 / norm); 
    }
  }

  // Normalize inner normal of each node 
  for (auto && n : nodes) {
    double norm = n.normal_i.norm(); 
    n.projected = norm; 
    n.normal_i = n.normal_i.scale(1.0 / norm);
  }
}

void AFMPB::generateGaussianPoint(const std::vector<Node> &nodes, 
                                  std::vector<GNode> &gauss) {
  using namespace dashmm; 

  for (auto && e : elements_) {
    int i1 = e.nodes[0]; 
    int i2 = e.nodes[1]; 
    int i3 = e.nodes[2]; 
    const Point &p1 = nodes[i1].position; 
    const Point &p2 = nodes[i2].position; 
    const Point &p3 = nodes[i3].position; 

    for (int j = 0; j < 7; ++j) {
      int index = 7 * e.index + j; 
      double zeta = 1 - xi_[j] - eta_[j]; 
      double x = p1.x() * zeta + p2.x() * xi_[j] + p3.x() * eta_[j]; 
      double y = p1.y() * zeta + p2.y() * xi_[j] + p3.y() * eta_[j]; 
      double z = p1.z() * zeta + p2.z() * xi_[j] + p3.z() * eta_[j]; 
      gauss.emplace_back(index, Point{x, y, z}, e.normal); 
    }
  }
}


} // namespace afmpb 
