#include "QuasiStaticModel.h"

#include <hpx/lcos/when_all.hpp>
#include <vector>

#include "BlazeIterative.hpp"
#include "data/DataManager.h"
#include "fe/mesh.h"
#include "geometry/neighbor.h"
#include "geometry/volumeCorrection.h"
#include "inp/decks/loadingDeck.h"
#include "inp/decks/materialDeck.h"
#include "inp/decks/modelDeck.h"
#include "inp/decks/outputDeck.h"
#include "inp/decks/solverDeck.h"
#include "inp/input.h"
#include "loading/fLoading.h"
#include "loading/initialCondition.h"
#include "loading/uLoading.h"
#include "material/materials.h"
#include "material/pdMaterial.h"
#include "util/parallel.h"
#include "util/stateBasedHelperFunctions.h"
#include "model/util.h"

template <class T>
model::QuasiStaticModel<T>::QuasiStaticModel(inp::Input *deck)
    : d_modelDeck_p(nullptr), d_outputDeck_p(nullptr) {
  d_osThreads = hpx::get_os_thread_count();

  // Generate as many data manager as os threads are avaibale
  for (size_t i = 0; i < d_osThreads; i++)
    d_dataManagers.push_back(new data::DataManager());

  d_dataManager_p = new data::DataManager();

  d_dataManager_p->setModelDeckP(deck->getModelDeck());
  d_dataManager_p->setOutputDeckP(deck->getOutputDeck());

  // Generate as many data manager as os threads are avaibale
  for (size_t i = 0; i < d_osThreads; i++) {
    d_dataManagers[i]->setModelDeckP(deck->getModelDeck());
    d_dataManagers[i]->setOutputDeckP(deck->getOutputDeck());
  }

  d_input_p = deck;

  // d_policy_p = inp::Policy::getInstance(d_input_p->getPolicyDeck());

  // if (d_modelDeck_p->d_isRestartActive)
  // restart(deck);
  // else
  // run(deck);

  initHObjects();

  solver();
}

template <class T>
model::QuasiStaticModel<T>::~QuasiStaticModel() {
  delete d_dataManager_p->getMeshP();
  delete d_dataManager_p->getDisplacementLoadingP();
  delete d_dataManager_p->getForceLoadingP();
  delete d_dataManager_p->getNeighborP();
  delete d_dataManager_p->getVolumeCorrectionP();
  delete d_dataManager_p->getBodyForceP();
  delete d_dataManager_p->getDisplacementP();
  delete d_dataManager_p->getVelocityP();

  for (size_t i = 0; i < d_osThreads; i++) delete d_dataManagers[i];

  delete d_material_p;
  delete d_dataManager_p;
}

template <class T>
void model::QuasiStaticModel<T>::initHObjects() {
  d_name = "QSModel";

  std::cout << d_name << ": Initializing high level objects." << std::endl;
  // read mesh data
  std::cout << d_name << ": Creating mesh." << std::endl;
  d_dataManager_p->setMeshP(new fe::Mesh(d_input_p->getMeshDeck()));
  d_dataManager_p->getMeshP()->clearElementData();

  d_nnodes = d_dataManager_p->getMeshP()->getNumNodes();

  std::cout << "Number of nodes = " << d_nnodes << std::endl;

  // initialize major simulation data
  d_dataManager_p->setBodyForceP(
      new std::vector<util::Point3>(d_nnodes, util::Point3()));

  d_dataManager_p->setForceP(
      new std::vector<util::Point3>(d_nnodes, util::Point3()));

  d_dataManager_p->setDisplacementP(
      new std::vector<util::Point3>(d_nnodes, util::Point3()));

  d_dataManager_p->setVelocityP(
      new std::vector<util::Point3>(d_nnodes, util::Point3()));

  // initialize loading class
  std::cout << d_name << ": Initializing displacement loading object."
            << std::endl;
  d_dataManager_p->setDisplacementLoadingP(new loading::ULoading(
      d_input_p->getLoadingDeck(), d_dataManager_p->getMeshP()));
  std::cout << d_name << ": Initializing force loading object." << std::endl;
  d_dataManager_p->setForceLoadingP(new loading::FLoading(
      d_input_p->getLoadingDeck(), d_dataManager_p->getMeshP()));
  // initialize neighbor class class
  std::cout << d_name << ": Creating neighbor list." << std::endl;

  d_dataManager_p->setNeighborP(new geometry::Neighbor(
      d_dataManager_p->getModelDeckP()->d_horizon, d_input_p->getNeighborDeck(),
      d_dataManager_p->getMeshP()->getNodesP()));

  // initialize the volume correction and weighted volumes
  d_dataManager_p->setVolumeCorrectionP(
      new geometry::VolumeCorrection(d_dataManager_p));

  // initialize material class
  std::cout << d_name << ": Initializing material object." << std::endl;

  d_material_p = new T(d_input_p->getMaterialDeck(), d_dataManager_p);

  // initialize jacobian matrix
  size_t dim = d_dataManager_p->getModelDeckP()->d_dim;
  size_t matrixSize = d_nnodes * dim;
  std::cout << d_name << ": Initializing Jacobian matrix (" << matrixSize << "x"
            << matrixSize << ")." << std::endl;
  jacobian = util::Matrixij(matrixSize, matrixSize, 0.);

  for (size_t i = 0; i < d_osThreads; i++) {
    d_dataManagers[i]->setMeshP(d_dataManager_p->getMeshP());
    d_dataManagers[i]->setBodyForceP(
        new std::vector<util::Point3>(d_nnodes, util::Point3()));
    d_dataManagers[i]->setForceP(
        new std::vector<util::Point3>(d_nnodes, util::Point3()));
    d_dataManagers[i]->setDisplacementP(
        new std::vector<util::Point3>(d_nnodes, util::Point3()));
    d_dataManagers[i]->setVelocityP(
        new std::vector<util::Point3>(d_nnodes, util::Point3()));

    d_dataManagers[i]->setForceLoadingP(d_dataManager_p->getForceLoadingP());
    d_dataManagers[i]->setDisplacementLoadingP(
        d_dataManager_p->getDisplacementLoadingP());
    d_dataManagers[i]->setVolumeCorrectionP(
        d_dataManager_p->getVolumeCorrectionP());
    d_dataManagers[i]->setNeighborP(d_dataManager_p->getNeighborP());
  }

  if (d_dataManager_p->getOutputDeckP()->isTagInOutput("Strain_Energy")) {
    d_dataManager_p->setStrainEnergyP(new std::vector<float>(d_nnodes, 0.));
  }

  if (d_dataManager_p->getOutputDeckP()->isTagInOutput("Stress_Tensor")) {
    d_dataManager_p->setStressTensorP(
        new std::vector<util::Matrix33>(d_nnodes));
  }

  if (d_dataManager_p->getOutputDeckP()->isTagInOutput("Strain_Tensor")) {
    d_dataManager_p->setStrainTensorP(
        new std::vector<util::Matrix33>(d_nnodes));
  }
}

template <class T>
void model::QuasiStaticModel<T>::computeForces(bool full) {
  d_material_p->update();

  // Clear the vector

  hpx::parallel::for_loop(hpx::parallel::execution::par, 0, d_nnodes,
                          [&](boost::uint64_t i) {
                            (*d_dataManager_p->getForceP())[i].d_x = 0.;
                            (*d_dataManager_p->getForceP())[i].d_y = 0.;
                            (*d_dataManager_p->getForceP())[i].d_z = 0.;
                          });

  hpx::lcos::local::mutex m;

  hpx::parallel::for_loop(
      hpx::parallel::execution::par, 0, d_nnodes,
      [&](boost::uint64_t i) {
        util::Point3 force_i = util::Point3();

        // inner loop over neighbors
        std::vector<std::size_t> i_neighs =
            d_dataManager_p->getNeighborP()->getNeighbors(i);

        for (size_t j = 0; j < i_neighs.size(); j++) {
          size_t j_id = i_neighs[j];

          auto res = d_material_p->getBondEF(size_t(i), size_t(j_id));

          force_i += res.first *
                     (*d_dataManager_p->getMeshP()->getNodalVolumesP())[j_id];

          m.lock();
          (*d_dataManager_p->getForceP())[j_id] -=
              res.first * (*d_dataManager_p->getMeshP()->getNodalVolumesP())[i];
          m.unlock();

          if (full) {
            if (d_dataManager_p->getOutputDeckP()->isTagInOutput(
                    "Strain_Energy"))

            {
              (*d_dataManager_p->getStrainEnergyP())[i] += (float)res.second;
            }
          }
        }

        if (full) {
          if (d_dataManager_p->getOutputDeckP()->isTagInOutput("Strain_Tensor"))

          {
            (*d_dataManager_p->getStrainTensorP())[i] =
                d_material_p->getStrain(size_t(i));
          }

          if (d_dataManager_p->getOutputDeckP()->isTagInOutput("Stress_Tensor"))

          {
            (*d_dataManager_p->getStressTensorP())[i] =
                d_material_p->getStress(size_t(i));
          }
        }

        // update force and energy
        m.lock();

        (*d_dataManager_p->getForceP())[i] += force_i;
        m.unlock();
      }  // loop over nodes

  );  // end of parallel for loop
}

template <class T>
inline void model::QuasiStaticModel<T>::computePertubatedForces(size_t thread) {
  material::pd::BaseMaterial *material =
      new T(d_input_p->getMaterialDeck(), d_dataManagers[thread]);

  material->update();

  // Clear the vector

  hpx::parallel::for_loop(hpx::parallel::execution::par, 0, d_nnodes,
                          [&](boost::uint64_t i) {
                            (*d_dataManagers[thread]->getForceP())[i].d_x = 0.;
                            (*d_dataManagers[thread]->getForceP())[i].d_y = 0.;
                            (*d_dataManagers[thread]->getForceP())[i].d_z = 0.;
                          });

  hpx::lcos::local::mutex m;

  hpx::parallel::for_loop(
      hpx::parallel::execution::par, 0, d_nnodes,
      [&](boost::uint64_t i) {
        util::Point3 force_i = util::Point3();

        // inner loop over neighbors
        std::vector<std::size_t> i_neighs =
            d_dataManager_p->getNeighborP()->getNeighbors(i);

        for (size_t j = 0; j < i_neighs.size(); j++) {
          size_t j_id = i_neighs[j];

          auto res = material->getBondEF(size_t(i), size_t(j_id));

          force_i += res.first *
                     (*d_dataManager_p->getMeshP()->getNodalVolumesP())[j_id];

          m.lock();
          (*d_dataManagers[thread]->getForceP())[j_id] -=
              res.first * (*d_dataManager_p->getMeshP()->getNodalVolumesP())[i];
          m.unlock();
        }
        // update force and energy
        m.lock();

        (*d_dataManagers[thread]->getForceP())[i] += force_i;
        m.unlock();
      }  // loop over nodes

  );  // end of parallel for loop

  delete material;
}

template <class T>
void model::QuasiStaticModel<T>::removeRow(util::Matrixij &matrix,
                                           size_t rowToRemove) {
  for (size_t i = 0; i < matrix.rows(); i++)

    for (size_t j = rowToRemove; j < matrix.columns() - 1; j++)

      matrix(j, i) = matrix(j + 1, i);

  matrix.resize(matrix.rows() - 1, matrix.columns(), false);
  matrix.shrinkToFit();
}

template <class T>
void model::QuasiStaticModel<T>::removeCol(util::Matrixij &matrix,
                                           size_t colToRemove) {
  for (size_t i = 0; i < matrix.rows(); i++)

    for (size_t j = colToRemove; j < matrix.columns() - 1; j++)

      matrix(i, j) = matrix(i, j + 1);

  matrix.resize(matrix.rows(), matrix.columns() - 1, true);
  matrix.shrinkToFit();
}

template <class T>
void model::QuasiStaticModel<T>::removeRow(util::VectorXi &vector,
                                           size_t rowToRemove) {
  for (size_t i = rowToRemove; i < vector.size() - 1; i++)

    vector[i] = vector[i + 1];

  vector.resize(vector.size() - 1, false);
  vector.shrinkToFit();
}

template <class T>
void model::QuasiStaticModel<T>::assembly_jacobian_matrix() {
  size_t dim = d_dataManager_p->getModelDeckP()->d_dim;
  size_t matrixSize = d_nnodes * dim;

  jacobian.resize(matrixSize, matrixSize);
  reset(jacobian);

  size_t slice = int(d_nnodes / d_osThreads);

  std::vector<hpx::future<void>> futures;

  for (size_t thread = 0; thread < d_osThreads; thread++) {
    size_t start = thread * slice;
    size_t end = 0;

    if (thread < d_osThreads - 1)
      end = (thread + 1) * slice;
    else
      end = d_nnodes;

    futures.push_back(hpx::async([this, start, end, thread]() {
      this->assembly_jacobian_matrix_part(start, end, thread);
    }));
  }

  for (size_t i = 0; i < futures.size(); i++) futures[i].get();

  // hpx::when_all(futures);
}

template <class T>
inline void model::QuasiStaticModel<T>::assembly_jacobian_matrix_part(
    size_t begin, size_t end, size_t thread) {
  double eps = d_input_p->getSolverDeck()->d_perturbation *
               d_dataManager_p->getMeshP()->getMeshSize();

  size_t dim = d_dataManager_p->getModelDeckP()->d_dim;

  std::vector<util::Point3> backup(d_nnodes, util::Point3());

  util::parallel::copy<std::vector<util::Point3>>(
      *d_dataManager_p->getDisplacementP(), backup);

  std::vector<size_t> removeId;

  // get instance of ULoading class (d_uLoading_p)
  auto bcD = d_dataManager_p->getDisplacementLoadingP()->d_bcData;
  auto bcN = d_dataManager_p->getDisplacementLoadingP()->d_bcNodes;

  // get fixed nodes
  size_t k = 0;
  for (auto bc : bcD) {
    auto direction = bc.d_direction;

    for (auto d : direction) {
      for (auto n : bcN[k]) {
        removeId.push_back(n);
      }
    }
    k++;
  }

  for (size_t i = begin; i < end; i++) {
    if (!(std::find(removeId.begin(), removeId.end(), i) != removeId.end())) {
      std::vector<size_t> *traversal_list = new std::vector<size_t>;

      traversal_list->push_back(i);

      std::vector<size_t> i_neighs =
          d_dataManager_p->getNeighborP()->getNeighbors(i);
      for (auto j : i_neighs) traversal_list->push_back(j);

      for (auto j : *traversal_list) {
        for (size_t r = 0; r < dim; r++) {
          std::vector<util::Point3> eps_vector =
              std::vector<util::Point3>(d_nnodes, util::Point3());

          switch (r) {
            case 0:
              eps_vector[j].d_x = eps;
              break;
            case 1:
              eps_vector[j].d_y = eps;
              break;
            case 2:
              eps_vector[j].d_z = eps;
              break;
          }

          std::vector<util::Point3> *tmp =
              new std::vector<util::Point3>(d_nnodes, util::Point3());
          util::parallel::copy(backup, *tmp);
          util::parallel::addInplace(*tmp, eps_vector);

          d_dataManagers[thread]->setDisplacementP(tmp);

          computePertubatedForces(thread);

          util::Point3 force_p = (*d_dataManagers[thread]->getForceP())[i];

          util::parallel::copy(backup, *tmp);
          util::parallel::subInplace(*tmp, eps_vector);

          d_dataManagers[thread]->setDisplacementP(tmp);

          computePertubatedForces(thread);

          util::Point3 force_m = (*d_dataManagers[thread]->getForceP())[i];

          util::Point3 f_diff = force_p - force_m;

          delete tmp;

          for (size_t s = 0; s < dim; s++) {
            if (r == s)
              jacobian(i * dim + r, j * dim + s) = f_diff[r] / (2. * eps);
          }
        }
      }

      delete traversal_list;
    }
  }
}

template <class T>
util::VectorXi model::QuasiStaticModel<T>::newton_step(util::VectorXi &res) {
  this->assembly_jacobian_matrix();

  std::vector<size_t> removeId;

  // get instance of ULoading class (d_uLoading_p)
  auto bcD = d_dataManager_p->getDisplacementLoadingP()->d_bcData;
  auto bcN = d_dataManager_p->getDisplacementLoadingP()->d_bcNodes;

  // get dimension
  auto dim = d_dataManager_p->getModelDeckP()->d_dim;

  size_t k = 0;
  for (auto bc : bcD) {
    auto direction = bc.d_direction;

    for (auto d : direction) {
      for (auto n : bcN[k]) {
        removeId.push_back(n * dim + d - 1);
      }
    }
    k++;
  }

  std::sort(removeId.begin(), removeId.end());
  std::reverse(removeId.begin(), removeId.end());

  for (auto id : removeId) {
    this->removeRow(jacobian, id);
    this->removeCol(jacobian, id);
    this->removeRow(res, id);
  }

  util::VectorXi x = util::VectorXi(res.size(), 0.);

  blaze::iterative::ConjugateGradientTag tag;

  if (d_input_p->getSolverDeck()->d_solverType == "BiCGSTAB")
    blaze::iterative::BiCGSTABTag tag;

  res *= -1.;

  x = blaze::iterative::solve(jacobian, res, tag);

  std::vector<bool> mask = std::vector<bool>(
      d_nnodes * d_dataManager_p->getModelDeckP()->d_dim, true);

  for (size_t i = 0; i < mask.size(); i++) {
    if (std::find(removeId.begin(), removeId.end(), i) != removeId.end())
      mask[i] = false;
  }

  util::VectorXi new_disp =
      util::VectorXi(d_nnodes * d_dataManager_p->getModelDeckP()->d_dim, 0.);

  size_t i = 0;
  size_t j = 0;

  for (auto m : mask) {
    if (m == true) {
      new_disp[i] = x[j];
      j += 1;
    }
    i += 1;
  }

  return new_disp;
}

template <class T>
void model::QuasiStaticModel<T>::solver() {
  size_t dim = d_dataManager_p->getModelDeckP()->d_dim;
  d_n = 1;
  d_time = 0;

  double delta_t = d_dataManager_p->getModelDeckP()->d_dt;

  // Write the initial data
  model::Output(d_input_p, d_dataManager_p, d_n - 1, d_time);

  for (; d_n < d_input_p->getModelDeck()->d_Nt + 1; d_n++) {
    d_time = d_n * delta_t;

    double residual = std::numeric_limits<double>::max();

    // Apply the force loading
    d_dataManager_p->getForceLoadingP()->apply(
        d_time, d_dataManager_p->getBodyForceP(), d_dataManager_p->getMeshP());

    // Apply the displacement loading
    d_dataManager_p->getDisplacementLoadingP()->apply(
        d_time, d_dataManager_p->getDisplacementP(),
        d_dataManager_p->getVelocityP(), d_dataManager_p->getMeshP());

    auto res = this->computeResidual();

    size_t iteration = 0;

    residual = util::l2Norm(res);

    std::cout << "It: " << iteration << " Res: " << residual << std::endl;

    while (residual >= d_input_p->getSolverDeck()->d_tol and
           iteration < d_input_p->getSolverDeck()->d_maxIters) {
      auto new_disp = this->newton_step(res);

      hpx::parallel::for_loop(
          hpx::parallel::execution::par, 0, d_nnodes, [&](boost::uint64_t i) {
            size_t id = i * dim;

            (*d_dataManager_p->getDisplacementP())[i].d_x =
                (*d_dataManager_p->getDisplacementP())[i].d_x + new_disp[id];
            if (dim >= 2)
              (*d_dataManager_p->getDisplacementP())[i].d_y =
                  (*d_dataManager_p->getDisplacementP())[i].d_y +
                  new_disp[id + 1];
            if (dim == 3)
              (*d_dataManager_p->getDisplacementP())[i].d_z =
                  (*d_dataManager_p->getDisplacementP())[i].d_z +
                  new_disp[id + 2];
          });

      this->computeForces();

      res = this->computeResidual();

      residual = util::l2Norm(res);
      iteration++;

      std::cout << "It: " << iteration << " Res: " << residual << std::endl;
    }

    this->computeForces(true);

    // Do the output after one successful iteration
    model::Output(d_input_p, d_dataManager_p, d_n, d_time);
  }
}

template <class T>
util::VectorXi model::QuasiStaticModel<T>::computeResidual() {
  size_t dim = d_dataManager_p->getModelDeckP()->d_dim;
  util::VectorXi res = util::VectorXi(d_nnodes * dim, 0.);

  auto bcN = d_dataManager_p->getDisplacementLoadingP()->d_bcNodes;
  auto bcD = d_dataManager_p->getDisplacementLoadingP()->d_bcData;

  for (size_t i = 0; i < d_nnodes; i++) {
    std::vector<size_t> dimensions(dim);
    std::iota(std::begin(dimensions), std::end(dimensions), 0);

    size_t k = 0;
    for (auto bc : bcD) {
      auto direction = bc.d_direction;

      for (auto d : direction) {
        for (auto n : bcN[k])

          if (std::find(bcN[k].begin(), bcN[k].end(), i) != bcN[k].end()) {
            // dimensions.erase(dimensions.begin()+d-1);
            dimensions[d - 1] = std::numeric_limits<std::size_t>::max();
          }
      }
      k++;
    }

    size_t id = i * dim;

    for (auto e : dimensions) {
      switch (e) {
        case 0:
          res[id] = (*d_dataManager_p->getForceP())[i][0] +
                    (*d_dataManager_p->getBodyForceP())[i][0]; 
          break;
        case 1:
          res[id + 1] =
              (*d_dataManager_p->getForceP())[i][1] +
              (*d_dataManager_p->getBodyForceP())[i][1];
          break;
        case 2:
          res[id + 2] =
              (*d_dataManager_p->getForceP())[i][2] +
              (*d_dataManager_p->getBodyForceP())[i][2];
          break;
        default:
          break;
      }
    }
  }

  return res;
}

template class model::QuasiStaticModel<material::pd::ElasticState>;
