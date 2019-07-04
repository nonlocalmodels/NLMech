// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "policy.h"
#include "decks/policyDeck.h"
#include <iostream>

static void addTag(std::vector<std::string> list, const std::string &tag) {

  bool found = false;
  for (const auto &s : list)
    if (s == tag)
      found = true;

  if (!found)
    list.emplace_back(tag);
}

inp::Policy *inp::Policy::d_instance_p = nullptr;

inp::Policy *inp::Policy::getInstance(inp::PolicyDeck *deck) {
  if (d_instance_p == nullptr)
    d_instance_p = new inp::Policy(deck);

  return d_instance_p;
}

void inp::Policy::destroyInstance() { delete d_instance_p; }

inp::Policy::~Policy() = default;

inp::Policy::Policy(inp::PolicyDeck *deck): d_enablePostProcessing(true),
d_memControlFlag(0), d_modelTag("Model"), d_maxLevel(0) {
  if (deck != nullptr) {
    d_enablePostProcessing = deck->d_enablePostProcessing;
    d_memControlFlag = deck->d_memControlFlag;
  }
  init();
}

void inp::Policy::init() {

  if (d_modelTag == "Model") {
    d_maxLevel = 3;
    if (d_lTags.empty()) {
      d_lTags.resize(d_maxLevel+1);

      // level 0 tags: none
      d_lTags[0] = std::vector<std::string>(0);

      // level 1 tags
      d_lTags[1].emplace_back("Model_d_e");
      d_lTags[1].emplace_back("Model_d_w");
      d_lTags[1].emplace_back("Model_d_phi");
      d_lTags[1].emplace_back("Model_d_strain");
      d_lTags[1].emplace_back("Model_d_stress");

      // level 2 tags
      d_lTags[2] = d_lTags[1];
      d_lTags[2].emplace_back("Model_d_eF");
      d_lTags[2].emplace_back("Model_d_eFB");

      // level 3 tags
      d_lTags[3] = d_lTags[2];
      d_lTags[3].emplace_back("Model_d_Z");
    }
  }
}

void inp::Policy::addToTags(const size_t &level, const std::string &tag) {

  if (level >= 0 && level <= d_maxLevel)
    for (size_t i=level; i<=d_maxLevel; i++)
      addTag(d_lTags[i], tag);
}

void inp::Policy::removeTag(const std::string &tag) {

  // find tag in d_memControlFlag and remove it from tag list
  std::vector<std::string> temp;
  for (const auto &s : d_lTags[d_memControlFlag]) {
    if (s == tag)
      continue;
    temp.emplace_back(s);
  }
  d_lTags[d_memControlFlag] = temp;
}

bool inp::Policy::populateData(const std::string &tag) {

  // depending on level find if variable_name is in any of the
  // tag list
  for (const auto &s : d_lTags[d_memControlFlag])
    if (s == tag)
      return false;

  return true;
}

int inp::Policy::getMemoryControlFlag() {
  return d_memControlFlag;
}

bool inp::Policy::enablePostProcessing() {
  return d_enablePostProcessing;
}
