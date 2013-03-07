/*
 * register_fvcr.h
 *
 *  Created on: 07.03.2013
 *      Author: andreasvogel
 */

#ifndef __H__UG__PLUGINS__NAVIER_STOKES__REGISTER_FVCR__
#define __H__UG__PLUGINS__NAVIER_STOKES__REGISTER_FVCR__

#include "registry/registry.h"
#include <string>

namespace ug{

void Init___NavierStokes___FVCR(ug::bridge::Registry* reg, std::string grp);

}// namespace ug

#endif /* __H__UG__PLUGINS__NAVIER_STOKES__REGISTER_FVCR__ */
