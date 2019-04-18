// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "quadrature.h"

template <class T> fe::Quadrature<T>::Quadrature(size_t order) {

  d_element_p = new T(order);
}

template <class T>
std::vector<fe::QuadData>
fe::Quadrature<T>::getQuadPoints(const std::vector<util::Point3> &nodes) {

  return d_element_p->getQuadPoints(nodes);
}

template <class T> size_t fe::Quadrature<T>::getElemType() {
  return d_element_p->getElemType();
}

template <class T> size_t fe::Quadrature<T>::getQuadOrder() {
  return d_element_p->getQuadOrder();
}

template <class T> size_t fe::Quadrature<T>::getNumQuadPoints() {
  return d_element_p->getNumQuadPoints();
}

template <class T>
std::vector<double> fe::Quadrature<T>::getShapes(const util::Point3 &p) {
  return d_element_p->getShapes(p);
}

template <class T>
std::vector<std::vector<double>>
fe::Quadrature<T>::getDerShapes(const util::Point3 &p) {
  return d_element_p->getDerShapes(p);
}

template <class T>
double fe::Quadrature<T>::mapRefElemToElem(
    util::Point3 &p, const std::vector<double> &shapes,
    const std::vector<std::vector<double>> &der_shapes,
    const std::vector<util::Point3> &nodes) {
  return d_element_p->mapRefElemToElem(p, shapes, der_shapes, nodes);
}

/*!
 * Important: Need to explicitly instantiate the template for the Quadrature
 * class to work
 */
template class fe::Quadrature<fe::TriElem>;
template class fe::Quadrature<fe::QuadElem>;