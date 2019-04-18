// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "vtkWriter.h"
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIdList.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnsignedIntArray.h>

rw::writer::VtkWriter::VtkWriter(const std::string &filename) {

  std::string f = filename + ".vtu";

  d_writer_p = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  d_writer_p->SetFileName(const_cast<char *>(f.c_str()));
}

void rw::writer::VtkWriter::appendNodes(
    const std::vector<util::Point3> *nodes) {

  auto points = vtkSmartPointer<vtkPoints>::New();

  for (auto p : *nodes)
    points->InsertNextPoint(p.d_x, p.d_y, p.d_z);

  d_grid_p = vtkSmartPointer<vtkUnstructuredGrid>::New();
  d_grid_p->SetPoints(points);
}

void rw::writer::VtkWriter::appendNodes(const std::vector<util::Point3> *nodes,
                                        const std::vector<util::Point3> *u) {

  auto points = vtkSmartPointer<vtkPoints>::New();

  for (size_t i = 0; i < nodes->size(); i++) {

    util::Point3 p = (*nodes)[i] + (*u)[i];
    points->InsertNextPoint(p.d_x, p.d_y, p.d_z);
  }

  d_grid_p = vtkSmartPointer<vtkUnstructuredGrid>::New();
  d_grid_p->SetPoints(points);
}

void rw::writer::VtkWriter::appendMesh(
    const std::vector<util::Point3> *nodes, const size_t &element_type,
    const std::vector<std::vector<size_t>> *en_con,
    const std::vector<util::Point3> *u) {

  // we write following things to the file
  //
  // Node data
  // 1. Coordinates of nodes (current)
  //
  // Element data
  // 1. element node connectivity
  // 2. element type (either triangle or square)

  // add current position of nodes
  this->appendNodes(nodes, u);

  // get the total number of elements
  size_t num_elems = en_con->size();

  //
  // process elements data
  //
  // element node connectivity
  auto cells = vtkSmartPointer<vtkCellArray>::New();

  size_t max_pts = (*en_con)[0].size();
  cells->Allocate(max_pts, num_elems);

  // element type
  int cell_types[num_elems];

  size_t elem_counter = 0;
  for (size_t i = 0; i < num_elems; i++) {

    // if element number is -1 skip the element
    std::vector<size_t> e_nodes = (*en_con)[i];

    size_t n = e_nodes.size();
    vtkIdType ids[n];
    for (size_t k = 0; k < n; k++)
      ids[k] = e_nodes[k];

    cells->InsertNextCell(n, ids);

    cell_types[elem_counter] =
        element_type; // see util::element struct definition

    elem_counter++;
  }

  // element node connectivity
  d_grid_p->SetCells(cell_types, cells);
}

void rw::writer::VtkWriter::appendPointData(const std::string &name,
                                            const std::vector<uint8_t> *data) {

  auto array = vtkSmartPointer<vtkDoubleArray>::New();
  array->SetNumberOfComponents(1);
  array->SetName(name.c_str());

  double value[1];
  for (unsigned char i : *data) {
    value[0] = i;
    array->InsertNextTuple(value);
  }

  d_grid_p->GetPointData()->AddArray(array);
}

void rw::writer::VtkWriter::appendPointData(const std::string &name,
                                            const std::vector<size_t> *data) {

  auto array = vtkSmartPointer<vtkDoubleArray>::New();
  array->SetNumberOfComponents(1);
  array->SetName(name.c_str());

  double value[1];
  for (unsigned long i : *data) {
    value[0] = i;
    array->InsertNextTuple(value);
  }

  d_grid_p->GetPointData()->AddArray(array);
}

void rw::writer::VtkWriter::appendPointData(const std::string &name,
                                            const std::vector<int> *data) {

  auto array = vtkSmartPointer<vtkDoubleArray>::New();
  array->SetNumberOfComponents(1);
  array->SetName(name.c_str());

  double value[1];
  for (int i : *data) {
    value[0] = i;
    array->InsertNextTuple(value);
  }

  d_grid_p->GetPointData()->AddArray(array);
}

void rw::writer::VtkWriter::appendPointData(const std::string &name,
                                            const std::vector<float> *data) {

  auto array = vtkSmartPointer<vtkDoubleArray>::New();
  array->SetNumberOfComponents(1);
  array->SetName(name.c_str());

  double value[1];
  for (float i : *data) {
    value[0] = i;
    array->InsertNextTuple(value);
  }

  d_grid_p->GetPointData()->AddArray(array);
}

void rw::writer::VtkWriter::appendPointData(const std::string &name,
                                            const std::vector<double> *data) {

  auto array = vtkSmartPointer<vtkDoubleArray>::New();
  array->SetNumberOfComponents(1);
  array->SetName(name.c_str());

  double value[1];
  for (double i : *data) {
    value[0] = i;
    array->InsertNextTuple(value);
  }

  d_grid_p->GetPointData()->AddArray(array);
}

void rw::writer::VtkWriter::appendPointData(
    const std::string &name, const std::vector<util::Point3> *data) {

  auto array = vtkSmartPointer<vtkDoubleArray>::New();
  array->SetNumberOfComponents(3);
  array->SetName(name.c_str());

  array->SetComponentName(0, "x");
  array->SetComponentName(1, "y");
  array->SetComponentName(2, "z");

  double value[3];
  for (const auto &i : *data) {
    value[0] = i.d_x;
    value[1] = i.d_y;
    value[2] = i.d_z;
    array->InsertNextTuple(value);
  }

  d_grid_p->GetPointData()->AddArray(array);
}

void rw::writer::VtkWriter::addTimeStep(const double &timestep) {

  auto t = vtkDoubleArray::New();
  t->SetName("TIME");
  t->SetNumberOfTuples(1);
  t->SetTuple1(0, timestep);
  d_grid_p->GetFieldData()->AddArray(t);
}

void rw::writer::VtkWriter::close() {
  d_writer_p->SetInputData(d_grid_p);
  d_writer_p->SetDataModeToAppended();
  d_writer_p->EncodeAppendedDataOn();
  d_writer_p->SetCompressor(0);
  d_writer_p->Write();
}

void rw::writer::VtkWriter::appendFieldData(const std::string &name,
                                            const double &data) {

  auto t = vtkDoubleArray::New();
  t->SetName(name.c_str());
  t->SetNumberOfTuples(1);
  t->SetTuple1(0, data);
  d_grid_p->GetFieldData()->AddArray(t);
}

void rw::writer::VtkWriter::appendFieldData(const std::string &name,
                                            const float &data) {

  auto t = vtkDoubleArray::New();
  t->SetName(name.c_str());
  t->SetNumberOfTuples(1);
  t->SetTuple1(0, data);
  d_grid_p->GetFieldData()->AddArray(t);
}