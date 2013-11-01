/*
 * register_navier_stokes.h
 *
 *  Created on: 01.11.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UG__PLUGINS__NAVIER_STOKES__REGISTER_NAVIER_STOKES__
#define __H__UG__PLUGINS__NAVIER_STOKES__REGISTER_NAVIER_STOKES__

#include "registry/registry.h"
#include <string>

namespace ug{

void Init___NavierStokes(ug::bridge::Registry* reg, std::string grp);

}// namespace ug

#endif /* __H__UG__PLUGINS__NAVIER_STOKES__REGISTER_NAVIER_STOKES__ */
