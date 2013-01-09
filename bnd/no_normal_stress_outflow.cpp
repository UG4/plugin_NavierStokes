/*
 * no_normal_stress_outflow.cpp
 *
 *  Created on: 27.03.2012
 *  D. Logashenko, A. Vogel
 */

#include "no_normal_stress_outflow.h"
#include "no_normal_stress_outflow_common.h"
#include "no_normal_stress_outflow_fv1.h"

namespace ug{
namespace NavierStokes{

template class NavierStokesNoNormalStressOutflow<Domain1d>;
template class NavierStokesNoNormalStressOutflow<Domain2d>;
template class NavierStokesNoNormalStressOutflow<Domain3d>;

} // namespace NavierStokes
} // namespace ug
