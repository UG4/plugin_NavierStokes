/*
 * inflow_base_impl.h
 *
 *  Created on: 12.04.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__NAVIER_STOKES__BND__INFLOW_BASE_IMPL__
#define __H__UG__NAVIER_STOKES__BND__INFLOW_BASE_IMPL__

#include "inflow_base.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{
namespace NavierStokes{

template <typename TDomain, typename TAlgebra>
void NavierStokesInflowBase<TDomain,TAlgebra>::
add(const std::vector<number>& vVel, const char* subsetsBND)
{
	if((int)vVel.size() != dim)
		UG_THROW("NavierStokesInflow: Setting velocity vector of dimension "<<
		         vVel.size()<<" to a Discretization for world dim " << dim);

	add(SmartPtr<ConstUserVector<dim> >(new ConstUserVector<dim>(vVel)), subsetsBND);
}


#ifdef UG_FOR_LUA
template <typename TDomain, typename TAlgebra>
void NavierStokesInflowBase<TDomain,TAlgebra>::
add(const char* fctName, const char* subsetsBND)
{
	add(LuaUserDataFactory<MathVector<dim>, dim>::create(fctName), subsetsBND);
}
#endif

} // namespace NavierStokes
} // end namespace ug

#endif /* __H__UG__NAVIER_STOKES__BND__INFLOW_BASE_IMPL__ */
