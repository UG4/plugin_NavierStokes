/*
 * wall_impl.h
 *
 *  Created on: 12.04.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__NAVIER_STOKES__BND__WALL_IMPL__
#define __H__UG__NAVIER_STOKES__BND__WALL_IMPL__

#include "wall.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{
namespace NavierStokes{

template <typename TDomain, typename TAlgebra>
NavierStokesWall<TDomain,TAlgebra>::
NavierStokesWall(SmartPtr< IncompressibleNavierStokesBase<TDomain> > spMaster)
	: m_spDirichletConstraint(new DirichletBoundary<TDomain,TAlgebra>)
{
	m_vFctName = spMaster->symb_fcts();

	if(m_vFctName.size() != TDomain::dim + 1)
		UG_THROW("NavierStokesWall::set_functions: This Boundary "
				"Condition works on exactly dim+1 (velocity+pressure) "
				"components, but "<<m_vFctName.size()<<"components given.");
}

template <typename TDomain, typename TAlgebra>
void NavierStokesWall<TDomain,TAlgebra>::
add(const char* subsetsBND)
{
	if(m_vFctName.empty())
		UG_THROW("NavierStokesWall::add: Symbolic names for"
				" velocity and pressure not set. Please set them first.");

	for(int i = 0; i < TDomain::dim; ++i)
	{
		m_spDirichletConstraint->add(0.0, m_vFctName[i].c_str(), subsetsBND);
	}
}

} // namespace NavierStokes
} // end namespace ug

#endif /* __H__UG__NAVIER_STOKES__BND__WALL_IMPL__ */
