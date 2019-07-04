////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

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

inp::Policy *inp::Policy::getInstance() {
  if (d_instance_p == nullptr)
    d_instance_p = new inp::Policy();

  return d_instance_p;
}

inp::Policy *inp::Policy::getInstance(inp::PolicyDeck *deck) {
  if (d_instance_p == nullptr)
    d_instance_p = new inp::Policy(deck);

  return d_instance_p;
}

void inp::Policy::destroyInstance() { delete d_instance_p; }

inp::Policy::Policy() {
  d_policyDeck_p = new inp::PolicyDeck();
  d_policyDeck_p->d_memControlFlag = 0; // default
  init();
}

inp::Policy::~Policy() = default;

inp::Policy::Policy(inp::PolicyDeck *deck) {
  d_policyDeck_p = deck;
  init();
}

void inp::Policy::init() {

  if (d_l1Tags.empty() == 0) {

    // level 0 tags: none
    d_l0Tags = std::vector<std::string>(0);

    // level 1 tags
    d_l1Tags.emplace_back("Model_d_e");
    d_l1Tags.emplace_back("Model_d_w");
    d_l1Tags.emplace_back("Model_d_phi");
    d_l1Tags.emplace_back("Model_d_strain");
    d_l1Tags.emplace_back("Model_d_stress");

    // level 2 tags
    d_l2Tags.emplace_back("Model_d_Z");
    d_l2Tags.emplace_back("Model_d_eF");
    d_l2Tags.emplace_back("Model_d_eFB");
  }
}

void inp::Policy::addToTags(const size_t &level, const std::string &tag) {

  if (level == 0)
    addTag(d_l0Tags, tag);
  else if (level == 1)
    addTag(d_l1Tags, tag);
  else if (level == 2)
    addTag(d_l2Tags, tag);
}

void inp::Policy::removeTag(const std::string &tag) {

  // depending on level find if variable_name is in any of the
  // tag list
  switch (d_policyDeck_p->d_memControlFlag) {

  case 0: {
    std::vector<std::string> temp;
    for (const auto &s : d_l0Tags) {
      if (s == tag)
        continue;
      temp.emplace_back(s);
    }
    d_l0Tags = temp;
  } break;

  case 1: {

    std::vector<std::string> temp;
    for (const auto &s : d_l0Tags) {
      if (s == tag)
        continue;
      temp.emplace_back(s);
    }
    d_l0Tags = temp;

    temp.clear();
    for (const auto &s : d_l1Tags) {
      if (s == tag)
        continue;
    temp.emplace_back(s);
    }
    d_l1Tags = temp;
  } break;

  case 2: {
    std::vector<std::string> temp;
    for (const auto &s : d_l0Tags) {
      if (s == tag)
        continue;
      temp.emplace_back(s);
    }
    d_l0Tags = temp;

    temp.clear();
    for (const auto &s : d_l1Tags) {
      if (s == tag)
        continue;
      temp.emplace_back(s);
    }
    d_l1Tags = temp;

    temp.clear();
    for (const auto &s : d_l2Tags) {
      if (s == tag)
        continue;
      temp.emplace_back(s);
    }
    d_l2Tags = temp;
  } break;

  default: {
    std::cerr << "Error: Check memory consumption flag.\n";
    exit(1);
  } break;
  }
}

bool inp::Policy::populateData(const std::string &tag) {

  // depending on level find if variable_name is in any of the
  // tag list
  switch (d_policyDeck_p->d_memControlFlag) {

  case 0: {
    for (const auto &s : d_l0Tags)
      if (s == tag)
        return false;

    return true;
  } break;

  case 1: {

    for (const auto &s : d_l0Tags)
      if (s == tag)
        return false;

    for (const auto &s : d_l1Tags)
      if (s == tag)
        return false;

    return true;
  } break;

  case 2: {
    for (const auto &s : d_l0Tags)
      if (s == tag)
        return false;

    for (const auto &s : d_l1Tags)
      if (s == tag)
        return false;

    for (const auto &s : d_l2Tags)
      if (s == tag)
        return false;

    return true;
  } break;

  default: {
    std::cerr << "Error: Check memory consumption flag.\n";
    exit(1);
  } break;
  }
}

int inp::Policy::getMemoryControlFlag() {
  return d_policyDeck_p->d_memControlFlag;
}

bool inp::Policy::enablePostProcessing() {
  return d_policyDeck_p->d_enablePostProcessing;
}
