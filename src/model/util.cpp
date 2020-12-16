#include "util.h"
#include "util/fastMethods.h"

model::Output::Output(inp::Input *d_input_p, data::DataManager *d_dataManager_p,
                      size_t d_n, double d_time) {

  if  (d_input_p->getModelDeck()->d_isRestartActive)
    d_n +=  d_input_p->getRestartDeck()->d_step;

  std::cout << "Output: time step = " << d_n << "\n";

  // write out % completion of simulation at 10% interval
  {
    float p = float(d_n) * 100. / d_input_p->getModelDeck()->d_Nt;
    int m = std::max(1, int(d_input_p->getModelDeck()->d_Nt / 10));
    if (d_n % m == 0 && int(p) > 0)
      std::cout << "Message: Simulation " << int(p) << "% complete.\n";
  }

  // filename
  // use smaller dt_out as the tag for files
  size_t dt_out = d_input_p->getOutputDeck()->d_dtOutCriteria;
  std::string filename = d_input_p->getOutputDeck()->d_path + "output_" +
                         std::to_string(d_n / dt_out);

  // open
  auto writer =
      rw::writer::Writer(filename, d_input_p->getOutputDeck()->d_outFormat,
                         d_input_p->getOutputDeck()->d_compressType);

  // write mesh
  if (d_dataManager_p->getMeshP()->getNumElements() != 0 &&
      d_input_p->getOutputDeck()->d_performFEOut)
    writer.appendMesh(d_dataManager_p->getMeshP()->getNodesP(),
                      d_dataManager_p->getMeshP()->getElementType(),
                      d_dataManager_p->getMeshP()->getElementConnectivitiesP(),
                      d_dataManager_p->getDisplacementP());
  else
    writer.appendNodes(d_dataManager_p->getMeshP()->getNodesP(),
                       d_dataManager_p->getDisplacementP());

  //
  // major simulation data
  //
  std::string tag = "Displacement";
  if (d_input_p->getOutputDeck()->isTagInOutput(tag))
    writer.appendPointData(tag, d_dataManager_p->getDisplacementP());

  tag = "Velocity";
  if (d_input_p->getOutputDeck()->isTagInOutput(tag))
    writer.appendPointData(tag, d_dataManager_p->getVelocityP());

  tag = "Force";
  if (d_input_p->getOutputDeck()->isTagInOutput(tag)) {
    std::vector<util::Point3> force(d_dataManager_p->getMeshP()->getNumNodes(),
                                    util::Point3());

    for (size_t i = 0; i < d_dataManager_p->getForceP()->size(); i++)
      force[i] = (*d_dataManager_p->getForceP())[i] *
                 d_dataManager_p->getMeshP()->getNodalVolume(i);

    writer.appendPointData(tag, &force);
  }

  writer.addTimeStep(d_time);

  //
  // minor simulation data
  //
  if (!d_input_p->getPolicyDeck()->d_enablePostProcessing) {
    writer.close();
    return;
  }

  tag = "Force_Density";
  if (d_input_p->getOutputDeck()->isTagInOutput(tag))
    writer.appendPointData(tag, d_dataManager_p->getForceP());

  tag = "Reaction_Force";
  if (d_input_p->getOutputDeck()->isTagInOutput(tag)) {
    writer.appendPointData(tag, d_dataManager_p->getReactionForceP());
  }

  tag = "Total_Reaction_Force";
  if (d_input_p->getOutputDeck()->isTagInOutput(tag)) {
    double sum =
        std::accumulate((*d_dataManager_p->getTotalReactionForceP()).begin(),
                        (*d_dataManager_p->getTotalReactionForceP()).end(), 0);

    // Computation of the area
    auto delta = d_input_p->getModelDeck()->d_horizon;
    auto min_x = d_dataManager_p->getMeshP()->getBoundingBox().first[0];
    auto max_x = d_dataManager_p->getMeshP()->getBoundingBox().second[0];
    auto max_y = d_dataManager_p->getMeshP()->getBoundingBox().second[1];
    auto min_y = d_dataManager_p->getMeshP()->getBoundingBox().first[1];

    double area =
        (std::abs(max_x - min_x) - delta) * (std::abs(max_y - min_y) - delta);

    writer.appendFieldData("Total_Reaction_Force", sum * area);
  }

  auto d_policy_p = inp::Policy::getInstance(d_input_p->getPolicyDeck());

  tag = "Strain_Energy";
  if (d_input_p->getOutputDeck()->isTagInOutput(tag) &&
      d_policy_p->populateData("Model_d_e"))
    writer.appendPointData(tag, d_dataManager_p->getStrainEnergyP());

  tag = "Work_Done";
  if (d_input_p->getOutputDeck()->isTagInOutput(tag) &&
      d_policy_p->populateData("Model_d_w"))
    writer.appendPointData(tag, d_dataManager_p->getWorkDoneP());

  tag = "Fixity";
  if (d_input_p->getOutputDeck()->isTagInOutput(tag))
    writer.appendPointData(tag, d_dataManager_p->getMeshP()->getFixityP());

  tag = "Node_Volume";
  if (d_input_p->getOutputDeck()->isTagInOutput(tag))
    writer.appendPointData(tag,
                           d_dataManager_p->getMeshP()->getNodalVolumesP());

  tag = "Damage_Phi";
  if (d_input_p->getOutputDeck()->isTagInOutput(tag) &&
      d_policy_p->populateData("Model_d_phi"))
    writer.appendPointData(tag, d_dataManager_p->getPhiP());

  tag = "Damage_Z";
  if (d_input_p->getOutputDeck()->isTagInOutput(tag) &&
      d_policy_p->populateData("Model_d_Z"))
    writer.appendPointData(tag, d_dataManager_p->getDamageFunctionP());

  tag = "Fracture_Perienergy_Bond";
  if (d_input_p->getOutputDeck()->isTagInOutput(tag) &&
      d_policy_p->populateData("Model_d_eFB"))
    writer.appendPointData(tag, d_dataManager_p->getBBFractureEnergyP());

  tag = "Fracture_Perienergy_Total";
  if (d_input_p->getOutputDeck()->isTagInOutput(tag) &&
      d_policy_p->populateData("Model_d_eF"))
    writer.appendPointData(tag, d_dataManager_p->getFractureEnergyP());

  tag = "Total_Energy";

  if (d_input_p->getOutputDeck()->isTagInOutput(tag) &&
      d_policy_p->populateData("Model_d_e")) {
    double d_te = util::methods::add((*d_dataManager_p->getStrainEnergyP()));
    double d_tw = util::methods::add((*d_dataManager_p->getWorkDoneP()));
    double d_tk = util::methods::add((*d_dataManager_p->getKineticEnergyP()));
    double te = d_te - d_tw + d_tk;
    writer.appendFieldData(tag, te);
  }

  tag = "Total_Fracture_Perienergy_Bond";
  if (d_input_p->getOutputDeck()->isTagInOutput(tag) &&
      d_policy_p->populateData("Model_d_eFB")) {
    double d_teFB =
        util::methods::add((*d_dataManager_p->getBBFractureEnergyP()));
    writer.appendFieldData(tag, d_teFB);
  }

  tag = "Total_Fracture_Perienergy_Total";
  if (d_input_p->getOutputDeck()->isTagInOutput(tag) &&
      d_policy_p->populateData("Model_d_eF")) {
    double d_teF = util::methods::add((*d_dataManager_p->getFractureEnergyP()));
    writer.appendFieldData(tag, d_teF);
  }

  tag = "Neighbors";
  if (d_input_p->getOutputDeck()->isTagInOutput(tag)) {
    std::vector<size_t> amountNeighbors;
    size_t nodes = d_dataManager_p->getMeshP()->getNumNodes();
    for (size_t i = 0; i < nodes; i++)
      amountNeighbors.push_back(
          d_dataManager_p->getNeighborP()->getNeighbors(i).size());
    writer.appendPointData(tag, &amountNeighbors);
  }

  tag = "Strain_Energy";
  if (d_input_p->getOutputDeck()->isTagInOutput(tag))
    writer.appendPointData(tag, d_dataManager_p->getStrainEnergyP());

  tag = "Strain_Tensor";
  if (d_input_p->getOutputDeck()->isTagInOutput(tag))
    writer.appendPointData("Strain_Tensor",
                           d_dataManager_p->getStrainTensorP());

  tag = "Stress_Tensor";
  if (d_input_p->getOutputDeck()->isTagInOutput(tag))
    writer.appendPointData("Stress_Tensor",
                           d_dataManager_p->getStressTensorP());

  writer.close();
}
