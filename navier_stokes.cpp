/*
 * navier_stokes.cpp
 *
 *  Created on: 04.07.2012
 *      Author: Christian Wehner
 */

#include "navier_stokes_common.h"
#include "navier_stokes_fv1.h"
#include "navier_stokes_fvho.h"
#include "navier_stokes_cr.h"

namespace ug{
namespace NavierStokes{

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class NavierStokes<Domain1d>;
#endif
#ifdef UG_DIM_2
template class NavierStokes<Domain2d>;
#endif
#ifdef UG_DIM_3
template class NavierStokes<Domain3d>;
#endif

} // namespace NavierStokes
} // end namespace ug

