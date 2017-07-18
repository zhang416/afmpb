#ifndef __AFMPB_SERIALIZER_H__
#define __AFMPB_SERIALIZER_H__

#include <cstring> 
#include "dashmm/serializer.h"
#include "afmpb.h"
#include "afmpb_lhs.h"

namespace dashmm {

// This serializer only skips the gmres field of the Node class
class NodeFullSerializer : public Serializer {
public: 
  ~NodeFullSerializer() { } 

  size_t size(void *object) const override {
    afmpb::Node *n = reinterpret_cast<afmpb::Node *>(object); 
    size_t retval = 0; 
    retval = sizeof(int) * 2 + // Index and number of elements in the patch
      sizeof(Point) * 2 + // normal_i and normal_o 
      sizeof(afmpb::Patch) * n->patch.size() + // patch 
      sizeof(double) * 4; // area, projected, and rhs
    return retval; 
  }

  void *serialize(void *object, void *buffer) const override {
    afmpb::Node *n = reinterpret_cast<afmpb::Node *>(object); 
    char *dest = reinterpret_cast<char *>(buffer); 
    size_t bytes = 0; 
    int n_patches = n->patch.size(); 

    bytes = sizeof(int); 
    memcpy(dest, &(n->index), bytes); 
    dest += bytes; 
    
    memcpy(dest, &n_patches, bytes); 
    dest += bytes; 

    bytes = sizeof(Point); 
    memcpy(dest, &(n->normal_i), bytes); 
    dest += bytes; 

    memcpy(dest, &(n->normal_o), bytes); 
    dest += bytes; 

    bytes = sizeof(afmpb::Patch) * n_patches; 
    memcpy(dest, n->patch.data(), bytes); 
    dest += bytes; 

    bytes = sizeof(double); 
    memcpy(dest, &(n->area), bytes); 
    dest += bytes; 

    memcpy(dest, &(n->projected), bytes); 
    dest += bytes; 

    memcpy(dest, &(n->rhs[0]), bytes * 2); 
    dest += bytes * 2; 

    return dest;
  }

  void *deserialize(void *buffer, void *object) const override {
    afmpb::Node *n = reinterpret_cast<afmpb::Node *>(object); 
    char *src = reinterpret_cast<char *>(buffer); 
    size_t bytes = 0; 
    int n_patches = 0; 

    bytes = sizeof(int); 
    memcpy(&(n->index), src, bytes); 
    src += bytes; 

    memcpy(&n_patches, src, bytes); 
    src += bytes; 

    bytes = sizeof(Point); 
    memcpy(&(n->normal_i), src, bytes); 
    src += bytes; 

    memcpy(&(n->normal_o), src, bytes); 
    src += bytes; 

    afmpb::Patch *p = reinterpret_cast<afmpb::Patch *>(src); 
    n->patch.assign(p, p + n_patches); 

    src += sizeof(afmpb::Patch) * n_patches; 

    bytes = sizeof(double); 
    memcpy(&(n->area), src, bytes); 
    src += bytes; 

    memcpy(&(n->projected), src, bytes); 
    src += bytes; 

    memcpy(n->rhs, src, bytes * 2); 
    src += bytes * 2; 

    return src; 
  }
}; 

// This serializer skips position, normal, and rhs field 
class NodePartialSerializer : public Serializer {
public: 
  ~NodePartialSerializer() { } 

  size_t size(void *object) const override {
    afmpb::Node *n = reinterpret_cast<afmpb::Node *>(object); 
    size_t retval = 0; 
    retval = sizeof(int) * 2 + // Index and number of elements in the patch 
      sizeof(afmpb::Patch) * n->patch.size() + // patch 
      sizeof(double) * 4; // area, projected, and gmres[0]@2 
    return retval; 
  }

  void *serialize(void *object, void *buffer) const override {
    afmpb::Node *n = reinterpret_cast<afmpb::Node *>(object); 
    char *dest = reinterpret_cast<char *>(buffer); 
    size_t bytes = 0; 
    int n_patches = n->patch.size(); 

    bytes = sizeof(int); 
    memcpy(dest, &(n->index), bytes); 
    dest += bytes; 

    memcpy(dest, &n_patches, bytes); 
    dest += bytes; 

    bytes = sizeof(afmpb::Patch) * n_patches; 
    memcpy(dest, n->patch.data(), bytes); 
    dest += bytes; 

    bytes = sizeof(double); 
    memcpy(dest, &(n->area), bytes); 
    dest += bytes; 

    memcpy(dest, &(n->projected), bytes); 
    dest += bytes; 

    double *gmres = n->gmres.data(); 
    memcpy(dest, gmres, bytes * 2); // gmres[0] and gmres[1]
    dest += bytes * 2; 

    return dest; 
  }

  void *deserialize(void *buffer, void *object) const override {
    afmpb::Node *n = reinterpret_cast<afmpb::Node *>(object); 
    char *src = reinterpret_cast<char *>(buffer); 
    size_t bytes = 0; 
    int n_patches = 0; 

    bytes = sizeof(int); 
    memcpy(&(n->index), src, bytes); 
    src += bytes; 

    memcpy(&n_patches, src, bytes); 
    src += bytes; 

    afmpb::Patch *p = reinterpret_cast<afmpb::Patch *>(src); 
    n->patch.assign(p, p + n_patches); 

    src += sizeof(afmpb::Patch) * n_patches; 

    bytes = sizeof(double); 
    memcpy(&(n->area), src, bytes); 
    src += bytes; 

    memcpy(&(n->projected), src, bytes); 
    src += bytes; 

    double *v = reinterpret_cast<double *>(src); 
    n->gmres[0] = v[0]; 
    n->gmres[1] = v[1]; 

    src += bytes * 2; 
    return src; 
  }
}; 

// This serializer only packs index and gmres
class NodeMinimumSerializer : public Serializer {
public: 
  ~NodeMinimumSerializer() { } 

  size_t size(void *object) const override {
    return sizeof(int) + sizeof(double) * 2; 
  } 

  void *serialize(void *object, void *buffer) const override {
    int iter = builtin_afmpb_table_->s_iter(); 
    afmpb::Node *n = reinterpret_cast<afmpb::Node *>(object); 
    char *dest = reinterpret_cast<char *>(buffer); 
    size_t bytes = 0; 

    bytes = sizeof(int);
    memcpy(dest, &(n->index), bytes); 
    dest += bytes; 

    bytes = sizeof(double) * 2; 
    double *v = &(n->gmres[2 * iter]); 
    memcpy(dest, v, bytes); 
    dest += bytes; 

    return dest; 
  }
   
  void *deserialize(void *buffer, void *object) const override {
    int iter = builtin_afmpb_table_->s_iter(); 
    afmpb::Node *n = reinterpret_cast<afmpb::Node *>(object); 
    char *src = reinterpret_cast<char *>(buffer); 
    size_t bytes = 0; 

    bytes = sizeof(int); 
    memcpy(&(n->index), src, bytes); 
    src += bytes; 

    bytes = sizeof(double) * 2; 
    double *v = reinterpret_cast<double *>(src); 
    n->gmres[2 * iter] = v[0]; 
    n->gmres[2 * iter + 1] = v[1]; 

    src += bytes; 

    return src; 
  }
}; 

} // namespace dashmm


#endif // __AFMPB_SERIALIZER_H__
