/*
 * register_fv1.h
 *
 *  Created on: 31.10.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UG__PLUGINS__NAVIER_STOKES__COMPRESSIBLE__FV1__REGISTER_FV1__
#define __H__UG__PLUGINS__NAVIER_STOKES__COMPRESSIBLE__FV1__REGISTER_FV1__

#include "registry/registry.h"
#include <string>

namespace ug{

void Init___CompressibleNavierStokes___FV1(ug::bridge::Registry* reg, std::string grp);

}// namespace ug

#endif /* __H__UG__PLUGINS__NAVIER_STOKES__COMPRESSIBLE__FV1__REGISTER_FV1__ */


