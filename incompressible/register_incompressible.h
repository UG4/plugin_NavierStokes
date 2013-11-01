/*
 * register_incompressible.h
 *
 *  Created on: 31.10.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__REGISTER_INCOMPRESSIBLE__
#define __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__REGISTER_INCOMPRESSIBLE__

#include "registry/registry.h"
#include <string>

namespace ug{

void Init___IncompressibleNavierStokes(ug::bridge::Registry* reg, std::string grp);

}// namespace ug

#endif /* __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__REGISTER_INCOMPRESSIBLE__ */
