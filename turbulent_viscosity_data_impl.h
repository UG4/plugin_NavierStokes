/*
 * turbulent_viscosity_data_impl.h
 *
 *  Created on: 01.11.2012
 *      Author: Christian Wehner
 */

#ifndef __H__UG__NAVIER_STOKES_TURBULENT_VISCOSITY_DATA_IMPL_H_
#define __H__UG__NAVIER_STOKES_TURBULENT_VISCOSITY_DATA_IMPL_H_

#include "common/common.h"

#include "turbulent_viscosity_data.h"
#include "lib_disc/common/subset_group.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/quadrature/quadrature.h"
#include "lib_disc/local_finite_element/local_shape_function_set.h"
#include "lib_disc/spatial_disc/user_data/user_data.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"


namespace ug{
namespace NavierStokes{

template<typename TGridFunction>
bool CRSmagorinskyTurbViscData<TGridFunction>::update(const TGridFunction& u){
	return true;
}

template<typename TGridFunction>
bool CRDynamicTurbViscData<TGridFunction>::update(const TGridFunction& u){
	UG_LOG("Dynamic model not implemented yet.\n");
	return true;
}

} // namespace NavierStokes
} // end namespace ug

#endif /* __H__UG__NAVIER_STOKES_TURBULENT_VISCOSITY_DATA_IMPL_H_ */