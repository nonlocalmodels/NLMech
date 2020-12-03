// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "legacyVtkWriter.h"
#include <util/feElementDefs.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIdList.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnsignedIntArray.h>

rw::writer::LegacyVtkWriter::LegacyVtkWriter(const std::string &filename,
                                             const std::string &compress_type)
    : d_filename(filename), d_compressType(compress_type) {
  std::string f = filename + ".vtu";

  d_myfile.open(f);

  if (!d_myfile.is_open()) {
    std::cerr << "Error: Could not open or generate following file: " << f
              << std::endl;
  }

  d_myfile << "# vtk DataFile Version 3.0" << std::endl;
  d_myfile << "NLMech legacy vtk writer" << std::endl;
  d_myfile << "ASCII" << std::endl;
}

void rw::writer::LegacyVtkWriter::appendNodes(
    const std::vector<util::Point3> *nodes) {
  d_myfile << "DATASET UNSTRUCTURED_GRID" << std::endl;
  d_myfile << "POINTS " << nodes->size() << " double" << std::endl;
  for (auto n : *nodes)
    d_myfile << n[0] << " " << n[1] << " " << n[2] << std::endl;

  d_myfile << "CELLS " << nodes->size() << " " << 2 * nodes->size()
           << std::endl;
  for (size_t i = 0; i < nodes->size(); i++)
    d_myfile << "1"
             << " " << i << std::endl;

  // Write the VTK cell type ( 1 = vtk vertex)
  d_myfile << "CELL_TYPES " << nodes->size() << std::endl;
  for (size_t i = 0; i < nodes->size(); i++) d_myfile << "1" << std::endl;

  d_myfile << "POINT_DATA " << nodes->size() << std::endl;
}

void rw::writer::LegacyVtkWriter::appendNodes(
    const std::vector<util::Point3> *nodes,
    const std::vector<util::Point3> *u) {
  d_myfile << "DATASET UNSTRUCTURED_GRID" << std::endl;
  d_myfile << "POINTS " << nodes->size() << " double" << std::endl;
  for (size_t i = 0; i < nodes->size(); i++) {
    util::Point3 p = (*nodes)[i];
    if (u) p = p + (*u)[i];
    d_myfile << p[0] << " " << p[1] << " " << p[2] << std::endl;
  }

  d_myfile << "CELLS " << nodes->size() << " " << 2 * nodes->size()
           << std::endl;
  for (size_t i = 0; i < nodes->size(); i++)
    d_myfile << "1"
             << " " << i << std::endl;

  // Write the VTK cell type ( 1 = vtk vertex)
  d_myfile << "CELL_TYPES " << nodes->size() << std::endl;
  for (size_t i = 0; i < nodes->size(); i++) d_myfile << "1" << std::endl;

  d_myfile << "POINT_DATA " << nodes->size() << std::endl;
}

void rw::writer::LegacyVtkWriter::appendMesh(
    const std::vector<util::Point3> *nodes, const size_t &element_type,
    const std::vector<size_t> *en_con, const std::vector<util::Point3> *u) {
  std::cerr
      << "Warning: Field data is not implementend in the legacy vtk writer"
      << std::endl;
}

void rw::writer::LegacyVtkWriter::appendPointData(
    const std::string &name, const std::vector<uint8_t> *data) {
  d_myfile << "SCALARS " << name << " int "
           << "1" << std::endl;
  d_myfile << "LOOKUP_TABLE default" << std::endl;
  for (size_t i = 0; i < data->size(); i++) d_myfile << (*data)[i] << std::endl;
}

void rw::writer::LegacyVtkWriter::appendPointData(
    const std::string &name, const std::vector<size_t> *data) {
  d_myfile << "SCALARS " << name << " int "
           << "1" << std::endl;
  d_myfile << "LOOKUP_TABLE default" << std::endl;
  for (size_t i = 0; i < data->size(); i++) d_myfile << (*data)[i] << std::endl;
}

void rw::writer::LegacyVtkWriter::appendPointData(
    const std::string &name, const std::vector<int> *data) {
  d_myfile << "SCALARS " << name << " int "
           << "1" << std::endl;
  d_myfile << "LOOKUP_TABLE default" << std::endl;
  for (size_t i = 0; i < data->size(); i++) d_myfile << (*data)[i] << std::endl;
}

void rw::writer::LegacyVtkWriter::appendPointData(
    const std::string &name, const std::vector<float> *data) {
  d_myfile << "SCALARS " << name << " float "
           << "1" << std::endl;
  d_myfile << "LOOKUP_TABLE default" << std::endl;
  for (size_t i = 0; i < data->size(); i++) d_myfile << (*data)[i] << std::endl;
}

void rw::writer::LegacyVtkWriter::appendPointData(
    const std::string &name, const std::vector<double> *data) {
  d_myfile << "SCALARS " << name << " double "
           << "1" << std::endl;
  d_myfile << "LOOKUP_TABLE default" << std::endl;
  for (size_t i = 0; i < data->size(); i++) d_myfile << (*data)[i] << std::endl;
}

void rw::writer::LegacyVtkWriter::appendPointData(
    const std::string &name, const std::vector<util::Point3> *data) {
  d_myfile << "SCALARS " << name << " double "
           << "3" << std::endl;
  d_myfile << "LOOKUP_TABLE default" << std::endl;
  for (size_t i = 0; i < data->size(); i++)
    d_myfile << (*data)[i][0] << " " << (*data)[i][1] << " " << (*data)[i][2]
             << std::endl;
}

void rw::writer::LegacyVtkWriter::appendPointData(
    const std::string &name, const std::vector<util::SymMatrix3> *data) {
  d_myfile << "TENSORS " << name << " double " << std::endl;

  for (size_t i = 0; i < data->size(); i++) {
    d_myfile << (*data)[i](0, 0) << " " << (*data)[i](0, 1) << " "
             << (*data)[i](0, 2) << std::endl;
    d_myfile << (*data)[i](1, 0) << " " << (*data)[i](1, 1) << " "
             << (*data)[i](1, 2) << std::endl;
    d_myfile << (*data)[i](2, 0) << " " << (*data)[i](2, 1) << " "
             << (*data)[i](2, 2) << std::endl;
  }
}

void rw::writer::LegacyVtkWriter::appendPointData(
    const std::string &name,
    const std::vector<blaze::StaticMatrix<double, 3, 3> > *data) {
  d_myfile << "TENSORS " << name << " double " << std::endl;

  for (size_t i = 0; i < data->size(); i++) {
    d_myfile << (*data)[i](0, 0) << " " << (*data)[i](0, 1) << " "
             << (*data)[i](0, 2) << std::endl;
    d_myfile << (*data)[i](1, 0) << " " << (*data)[i](1, 1) << " "
             << (*data)[i](1, 2) << std::endl;
    d_myfile << (*data)[i](2, 0) << " " << (*data)[i](2, 1) << " "
             << (*data)[i](2, 2) << std::endl;
  }
}

void rw::writer::LegacyVtkWriter::appendCellData(
    const std::string &name, const std::vector<float> *data) {
  std::cerr
      << "Warning: Field data is not implementend in the legacy vtk writer"
      << std::endl;
}

void rw::writer::LegacyVtkWriter::appendCellData(
    const std::string &name, const std::vector<util::SymMatrix3> *data) {
  std::cerr
      << "Warning: Field data is not implementend in the legacy vtk writer"
      << std::endl;
}

void rw::writer::LegacyVtkWriter::addTimeStep(const double &timestep) {}

void rw::writer::LegacyVtkWriter::close() { d_myfile.close(); }

void rw::writer::LegacyVtkWriter::appendFieldData(const std::string &name,
                                                  const double &data) {
  std::cerr
      << "Warning: Field data is not implementend in the legacy vtk writer"
      << std::endl;
}

void rw::writer::LegacyVtkWriter::appendFieldData(const std::string &name,
                                                  const float &data) {
  std::cerr
      << "Warning: Field data is not implementend in the legacy vtk writer"
      << std::endl;
}