//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LevelSetContactor.h"

InputParameters
LevelSetContactor::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addClassDescription("Base class for a rigid contactor described implicitly by a "
                             "signed-distance (level-set) function.");
  return params;
}

LevelSetContactor::LevelSetContactor(const InputParameters & parameters)
  : GeneralUserObject(parameters)
{
}
