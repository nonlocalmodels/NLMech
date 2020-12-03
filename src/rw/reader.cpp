////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "reader.h"
#include "csv.h"
#include "mshReader.h"
#include "vtkReader.h"

void rw::reader::readCsvFile(const std::string &filename, size_t dim,
                             std::vector<util::Point3> *nodes,
                             std::vector<double> *volumes) {
  nodes->clear();
  volumes->clear();
  if (dim == 1) {
    io::CSVReader<3> in(filename);
    in.read_header(io::ignore_extra_column, "id", "x", "volume");

    double x;
    double volume;
    int id;
    while (in.read_row(id, x, volume)) {
      volumes->emplace_back(volume);
      nodes->emplace_back(util::Point3(x, 0., 0.));
    }
  }

  if (dim == 2) {
    io::CSVReader<4> in(filename);
    in.read_header(io::ignore_extra_column, "id", "x", "y", "volume");

    double x, y, volume;
    int id;
    while (in.read_row(id, x, y, volume)) {
      volumes->emplace_back(volume);
      nodes->emplace_back(util::Point3(x, y, 0.));
    }
  }

  if (dim == 3) {
    io::CSVReader<5> in(filename);
    in.read_header(io::ignore_extra_column, "id", "x", "y", "z", "volume");

    double x, y, z, volume;
    int id;
    while (in.read_row(id, x, y, z, volume)) {
      volumes->emplace_back(volume);
      nodes->emplace_back(util::Point3(x, y, z));
    }
  }
}

void rw::reader::readVtuFile(const std::string &filename, size_t dim,
                             std::vector<util::Point3> *nodes,
                             size_t &element_type, size_t &num_elem,
                             std::vector<size_t> *enc,
                             std::vector<std::vector<size_t>> *nec,
                             std::vector<double> *volumes, bool is_fd) {
  // call vtk reader
  auto rdr = rw::reader::VtkReader(filename);
  rdr.readMesh(dim, nodes, element_type, num_elem, enc, nec, volumes, is_fd);
  rdr.close();
}

void rw::reader::readVtuFileNodes(const std::string &filename, size_t dim,
                                  std::vector<util::Point3> *nodes,
                                  bool ref_config) {
  // call vtk reader
  auto rdr = rw::reader::VtkReader(filename);

  // below will read the current position of nodes
  rdr.readNodes(nodes);

  // need to subtract the displacement to get reference configuration of nodes
  if (ref_config) {
    std::vector<util::Point3> u;
    if (rdr.readPointData("Displacement", &u)) {
      std::cerr << "Error: Did not find displacement in the vtu file."
                << std::endl;
      exit(1);
    }

    if (u.size() != nodes->size()) {
      std::cerr << "Error: Displacement data and node data size do not match."
                << std::endl;
      exit(1);
    }

    for (size_t i = 0; i < u.size(); i++) (*nodes)[i] -= u[i];
  }

  rdr.close();
}

bool rw::reader::vtuHasPointData(const std::string &filename,
                                 const std::string &tag) {
  // call vtk reader
  auto rdr = rw::reader::VtkReader(filename);
  auto has_data = rdr.vtuHasPointData(tag);
  rdr.close();
  return has_data;
}

bool rw::reader::vtuHasCellData(const std::string &filename,
                                const std::string &tag) {
  // call vtk reader
  auto rdr = rw::reader::VtkReader(filename);
  auto has_data = rdr.vtuHasCellData(tag);
  rdr.close();
  return has_data;
}

std::vector<std::string> rw::reader::readVtuFilePointTags(
    const std::string &filename) {
  // call vtk reader
  auto rdr = rw::reader::VtkReader(filename);
  auto data = rdr.readVtuFilePointTags();
  rdr.close();
  return data;
}

std::vector<std::string> rw::reader::readVtuFileCellTags(
    const std::string &filename) {
  // call vtk reader
  auto rdr = rw::reader::VtkReader(filename);
  auto data = rdr.readVtuFileCellTags();
  rdr.close();
  return data;
}

void rw::reader::readVtuFileCells(const std::string &filename, size_t dim,
                                  size_t &element_type, size_t &num_elem,
                                  std::vector<size_t> *enc,
                                  std::vector<std::vector<size_t>> *nec) {
  // call vtk reader
  auto rdr = rw::reader::VtkReader(filename);

  // below will read the current position of nodes
  rdr.readCells(dim, element_type, num_elem, enc, nec);

  rdr.close();
}

void rw::reader::readVtuFileRestart(const std::string &filename,
                                    std::vector<util::Point3> *u,
                                    std::vector<util::Point3> *v,
                                    const std::vector<util::Point3> *X) {
  // call vtk reader
  auto rdr = rw::reader::VtkReader(filename);
  // if displacement is not in input file, use reference coordinate to get
  // displacement
  if (!rdr.readPointData("Displacement", u)) {
    std::vector<util::Point3> y;
    rdr.readNodes(&y);
    if (y.size() != X->size()) {
      std::cerr << "Error: Number of nodes in input file = " << filename
                << " and number nodes in data X are not same.\n";
      exit(1);
    }

    u->resize(y.size());
    for (size_t i = 0; i < y.size(); i++) (*u)[i] = y[i] - (*X)[i];
  }

  // get velocity
  rdr.readPointData("Velocity", v);
  rdr.close();
}

bool rw::reader::readVtuFilePointData(const std::string &filename,
                                      const std::string &tag,
                                      std::vector<uint8_t> *data) {
  // call vtk reader
  auto rdr = rw::reader::VtkReader(filename);
  // read data
  auto st = rdr.readPointData(tag, data);
  rdr.close();
  return st;
}

bool rw::reader::readVtuFilePointData(const std::string &filename,
                                      const std::string &tag,
                                      std::vector<size_t> *data) {
  // call vtk reader
  auto rdr = rw::reader::VtkReader(filename);
  // read data
  auto st = rdr.readPointData(tag, data);
  rdr.close();
  return st;
}

bool rw::reader::readVtuFilePointData(const std::string &filename,
                                      const std::string &tag,
                                      std::vector<int> *data) {
  // call vtk reader
  auto rdr = rw::reader::VtkReader(filename);
  // read data
  auto st = rdr.readPointData(tag, data);
  rdr.close();
  return st;
}

bool rw::reader::readVtuFilePointData(const std::string &filename,
                                      const std::string &tag,
                                      std::vector<float> *data) {
  // call vtk reader
  auto rdr = rw::reader::VtkReader(filename);
  // read data
  auto st = rdr.readPointData(tag, data);
  rdr.close();
  return st;
}

bool rw::reader::readVtuFilePointData(const std::string &filename,
                                      const std::string &tag,
                                      std::vector<double> *data) {
  // call vtk reader
  auto rdr = rw::reader::VtkReader(filename);
  // read data
  auto st = rdr.readPointData(tag, data);
  rdr.close();
  return st;
}

bool rw::reader::readVtuFilePointData(const std::string &filename,
                                      const std::string &tag,
                                      std::vector<util::Point3> *data) {
  // call vtk reader
  auto rdr = rw::reader::VtkReader(filename);
  // read data
  auto st = rdr.readPointData(tag, data);
  rdr.close();
  return st;
}

bool rw::reader::readVtuFilePointData(const std::string &filename,
                                      const std::string &tag,
                                      std::vector<util::SymMatrix3> *data) {
  // call vtk reader
  auto rdr = rw::reader::VtkReader(filename);
  // read data
  auto st = rdr.readPointData(tag, data);
  rdr.close();
  return st;
}

bool rw::reader::readVtuFilePointData(const std::string &filename,
                                      const std::string &tag,
                                      std::vector<util::Matrix33> *data) {
  // call vtk reader
  auto rdr = rw::reader::VtkReader(filename);
  // read data
  auto st = rdr.readPointData(tag, data);
  rdr.close();
  return st;
}

bool rw::reader::readVtuFileCellData(const std::string &filename,
                                     const std::string &tag,
                                     std::vector<float> *data) {
  // call vtk reader
  auto rdr = rw::reader::VtkReader(filename);
  // read data
  auto st = rdr.readCellData(tag, data);
  rdr.close();
  return st;
}

bool rw::reader::readVtuFileCellData(const std::string &filename,
                                     const std::string &tag,
                                     std::vector<double> *data) {
  // call vtk reader
  auto rdr = rw::reader::VtkReader(filename);
  // read data
  auto st = rdr.readCellData(tag, data);
  rdr.close();
  return st;
}

bool rw::reader::readVtuFileCellData(const std::string &filename,
                                     const std::string &tag,
                                     std::vector<util::Point3> *data) {
  // call vtk reader
  auto rdr = rw::reader::VtkReader(filename);
  // read data
  auto st = rdr.readCellData(tag, data);
  rdr.close();
  return st;
}

bool rw::reader::readVtuFileCellData(const std::string &filename,
                                     const std::string &tag,
                                     std::vector<util::SymMatrix3> *data) {
  // call vtk reader
  auto rdr = rw::reader::VtkReader(filename);
  // read data
  auto st = rdr.readCellData(tag, data);
  rdr.close();
  return st;
}

bool rw::reader::readVtuFileCellData(const std::string &filename,
                                     const std::string &tag,
                                     std::vector<util::Matrix33> *data) {
  // call vtk reader
  auto rdr = rw::reader::VtkReader(filename);
  // read data
  auto st = rdr.readCellData(tag, data);
  rdr.close();
  return st;
}

void rw::reader::readMshFile(const std::string &filename, size_t dim,
                             std::vector<util::Point3> *nodes,
                             size_t &element_type, size_t &num_elem,
                             std::vector<size_t> *enc,
                             std::vector<std::vector<size_t>> *nec,
                             std::vector<double> *volumes, bool is_fd) {
  // call vtk reader
  auto rdr = rw::reader::MshReader(filename);
  rdr.readMesh(dim, nodes, element_type, num_elem, enc, nec, volumes, is_fd);
  rdr.close();
}

void rw::reader::readMshFileRestart(const std::string &filename,
                                    std::vector<util::Point3> *u,
                                    std::vector<util::Point3> *v,
                                    const std::vector<util::Point3> *X) {
  // call vtk reader
  auto rdr = rw::reader::MshReader(filename);
  // if displacement is not in input file, use reference coordinate to get
  // displacement
  if (!rdr.readPointData("Displacement", u)) {
    std::vector<util::Point3> y;
    rdr.readNodes(&y);
    if (y.size() != X->size()) {
      std::cerr << "Error: Number of nodes in input file = " << filename
                << " and number nodes in data X are not same.\n";
      exit(1);
    }

    u->resize(y.size());
    for (size_t i = 0; i < y.size(); i++) (*u)[i] = y[i] - (*X)[i];
  }

  // get velocity
  rdr.readPointData("Velocity", v);
  rdr.close();
}

bool rw::reader::readMshFilePointData(const std::string &filename,
                                      const std::string &tag,
                                      std::vector<double> *data) {
  // call vtk reader
  auto rdr = rw::reader::VtkReader(filename);
  // get velocity
  auto st = rdr.readPointData(tag, data);
  rdr.close();
  return st;
}
