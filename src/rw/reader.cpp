// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "reader.h"
#include "vtkReader.h"
#include "mshReader.h"
#include "../external/csv.h"

void rw::reader::readCsvFile(const std::string& filename, size_t dim,
    std::vector<util::Point3>* nodes,
                 std::vector<double>* volumes) {

  nodes->clear();
  volumes->clear();
  if (dim == 1) {

    io::CSVReader<3> in(filename);
    in.read_header(io::ignore_extra_column, "id", "x", "volume");

    double x; double volume;
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
                 std::vector<util::Point3> *nodes, size_t &element_type,
                 std::vector<size_t> *enc, std::vector<std::vector<size_t>> *nec,
                 std::vector<double> *volumes, bool is_fd) {

  // call vtk reader
  rw::reader::VtkReader rdr = rw::reader::VtkReader(filename);
  rdr.readMesh(dim, nodes, element_type, enc, nec, volumes, is_fd);
  rdr.close();
}

void rw::reader::readMshFile(const std::string &filename, size_t dim,
                             std::vector<util::Point3> *nodes, size_t &element_type,
                             std::vector<size_t> *enc, std::vector<std::vector<size_t>> *nec,
                             std::vector<double> *volumes, bool is_fd) {

  // call vtk reader
  rw::reader::MshReader rdr = rw::reader::MshReader(filename);
  rdr.readMesh(dim, nodes, element_type, enc, nec, volumes, is_fd);
}