// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "writer.h"
#include "vtkWriter.h"
#include "legacyVtkWriter.h"

//
// VtkWriterInterface
//
rw::writer::VtkWriterInterface::VtkWriterInterface() : d_vtkWriter_p(nullptr) {}

rw::writer::VtkWriterInterface::VtkWriterInterface(
    const std::string &filename, const std::string &compress_type)
    : d_vtkWriter_p(nullptr) {

  d_vtkWriter_p = new rw::writer::VtkWriter(filename, compress_type);
}

void rw::writer::VtkWriterInterface::open(const std::string &filename,
                                          const std::string &compress_type) {
  if (!d_vtkWriter_p)
    d_vtkWriter_p = new rw::writer::VtkWriter(filename, compress_type);
}

rw::writer::VtkWriterInterface::~VtkWriterInterface() {
  delete (d_vtkWriter_p);
}

void rw::writer::VtkWriterInterface::appendNodes(
    const std::vector<util::Point3> *nodes) {
  d_vtkWriter_p->appendNodes(nodes);
}

void rw::writer::VtkWriterInterface::appendNodes(
    const std::vector<util::Point3> *nodes,
    const std::vector<util::Point3> *u) {
  d_vtkWriter_p->appendNodes(nodes, u);
}

void rw::writer::VtkWriterInterface::appendMesh(
    const std::vector<util::Point3> *nodes, const size_t &element_type,
    const std::vector<size_t> *en_con, const std::vector<util::Point3> *u) {

  d_vtkWriter_p->appendMesh(nodes, element_type, en_con, u);
}

void rw::writer::VtkWriterInterface::appendPointData(
    const std::string &name, const std::vector<uint8_t> *data) {
  d_vtkWriter_p->appendPointData(name, data);
}

void rw::writer::VtkWriterInterface::appendPointData(
    const std::string &name, const std::vector<size_t> *data) {
  d_vtkWriter_p->appendPointData(name, data);
}

void rw::writer::VtkWriterInterface::appendPointData(
    const std::string &name, const std::vector<int> *data) {
  d_vtkWriter_p->appendPointData(name, data);
}

void rw::writer::VtkWriterInterface::appendPointData(
    const std::string &name, const std::vector<float> *data) {
  d_vtkWriter_p->appendPointData(name, data);
}

void rw::writer::VtkWriterInterface::appendPointData(
    const std::string &name, const std::vector<double> *data) {
  d_vtkWriter_p->appendPointData(name, data);
}

void rw::writer::VtkWriterInterface::appendPointData(
    const std::string &name, const std::vector<util::Point3> *data) {
  d_vtkWriter_p->appendPointData(name, data);
}

void rw::writer::VtkWriterInterface::appendPointData(
    const std::string &name, const std::vector<util::SymMatrix3> *data) {
  d_vtkWriter_p->appendPointData(name, data);
}

void rw::writer::VtkWriterInterface::appendCellData(
    const std::string &name, const std::vector<float> *data) {
  d_vtkWriter_p->appendCellData(name, data);
}

void rw::writer::VtkWriterInterface::appendCellData(
    const std::string &name, const std::vector<util::SymMatrix3> *data) {
  d_vtkWriter_p->appendCellData(name, data);
}

void rw::writer::VtkWriterInterface::addTimeStep(const double &timestep) {
  d_vtkWriter_p->addTimeStep(timestep);
}

void rw::writer::VtkWriterInterface::close() { d_vtkWriter_p->close(); }

void rw::writer::VtkWriterInterface::appendFieldData(const std::string &name,
                                                     const double &data) {
  d_vtkWriter_p->appendFieldData(name, data);
}

void rw::writer::VtkWriterInterface::appendFieldData(const std::string &name,
                                                     const float &data) {
  d_vtkWriter_p->appendFieldData(name, data);
}

//
// LegacyVtkWriterInterface
//
rw::writer::LegacyVtkWriterInterface::LegacyVtkWriterInterface() : d_vtkWriter_p(nullptr) {}

rw::writer::LegacyVtkWriterInterface::LegacyVtkWriterInterface(
  const std::string &filename, const std::string &compress_type)
  : d_vtkWriter_p(nullptr) {

  d_vtkWriter_p = new rw::writer::LegacyVtkWriter(filename, compress_type);
}

void rw::writer::LegacyVtkWriterInterface::open(const std::string &filename,
                                          const std::string &compress_type) {
  if (!d_vtkWriter_p)
    d_vtkWriter_p = new rw::writer::LegacyVtkWriter(filename, compress_type);
}

rw::writer::LegacyVtkWriterInterface::~LegacyVtkWriterInterface() {
  delete (d_vtkWriter_p);
}

void rw::writer::LegacyVtkWriterInterface::appendNodes(
  const std::vector<util::Point3> *nodes) {
  d_vtkWriter_p->appendNodes(nodes);
}

void rw::writer::LegacyVtkWriterInterface::appendNodes(
  const std::vector<util::Point3> *nodes,
  const std::vector<util::Point3> *u) {
  d_vtkWriter_p->appendNodes(nodes, u);
}

void rw::writer::LegacyVtkWriterInterface::appendMesh(
  const std::vector<util::Point3> *nodes, const size_t &element_type,
  const std::vector<size_t> *en_con, const std::vector<util::Point3> *u) {

  d_vtkWriter_p->appendMesh(nodes, element_type, en_con, u);
}