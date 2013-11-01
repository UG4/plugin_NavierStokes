/*
 * register_fe.h
 *
 *  Created on: 07.03.2013
 *      Author: andreasvogel
 */

#ifndef __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FE__REGISTER_FE__
#define __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FE__REGISTER_FE__

#include "registry/registry.h"
#include <string>

namespace ug{

void Init___NavierStokes___FE(ug::bridge::Registry* reg, std::string grp);

}// namespace ug

#endif /* __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FE__REGISTER_FE__ */
