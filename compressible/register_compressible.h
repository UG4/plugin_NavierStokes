/*
 * register_compressible.h
 *
 *  Created on: 31.10.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UG__PLUGINS__NAVIER_STOKES__COMPRESSIBLE__REGISTER_COMPRESSIBLE__
#define __H__UG__PLUGINS__NAVIER_STOKES__COMPRESSIBLE__REGISTER_COMPRESSIBLE__

#include "registry/registry.h"
#include <string>

namespace ug{

void Init___CompressibleNavierStokes(ug::bridge::Registry* reg, std::string grp);

}// namespace ug

#endif /* __H__UG__PLUGINS__NAVIER_STOKES__COMPRESSIBLE__REGISTER_COMPRESSIBLE__ */

