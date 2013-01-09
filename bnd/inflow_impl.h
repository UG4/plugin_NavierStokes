/*
 * navier_stokes_bnd_impl.h
 *
 *  Created on: 12.04.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__NAVIER_STOKES_BND_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__NAVIER_STOKES_BND_IMPL__

#include "inflow.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{
namespace NavierStokes{

////////////////////////////////////////////////////////////////////////////////
//	NavierStokesInflow
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
NavierStokesInflow<TDomain,TAlgebra>::
NavierStokesInflow(SmartPtr< NavierStokes<TDomain> > spMaster)
	: m_spNeumannDisc(new NeumannBoundary<TDomain>(spMaster->symb_subsets())),
	  m_spDirichletConstraint(new DirichletBoundary<TDomain,TAlgebra>),
	  m_spMaster(spMaster)
{
	m_vFctName = spMaster->symb_fcts();

	if(m_vFctName.size() != TDomain::dim + 1)
		UG_THROW("NavierStokesInflow::set_functions: This Boundary "
				"Condition works on exactly dim+1 (velocity+pressure) "
				"components, but "<<m_vFctName.size()<<"components given.");
}

template <typename TDomain, typename TAlgebra>
void NavierStokesInflow<TDomain,TAlgebra>::
add(SmartPtr<UserData<MathVector<dim>, dim> > user, const char* subsetsBND)
{
	if(m_vFctName.empty())
		UG_THROW("NavierStokesInflow::add: Symbolic names for"
				" velocity and pressure not set. Please set them first.");

	if(m_spMaster->disc_scheme() != "staggered")
		m_spNeumannDisc->add(user, m_vFctName[dim].c_str(), subsetsBND);

	std::string velNames;
	for(int i=0;i<TDomain::dim; ++i)
	{
		if(i>0) velNames.append(",");
		velNames.append(m_vFctName[i]);
	}
	m_spDirichletConstraint->add(user, velNames.c_str(), subsetsBND);
}

template <typename TDomain, typename TAlgebra>
void NavierStokesInflow<TDomain,TAlgebra>::
add(const std::vector<number>& vVel, const char* subsetsBND)
{
	if((int)vVel.size() != dim)
		UG_THROW("NavierStokesInflow: Setting velocity vector of dimension "<<
		         vVel.size()<<" to a Discretization for world dim " << dim);

	add(SmartPtr<ConstUserVector<dim> >(new ConstUserVector<dim>(vVel)), subsetsBND);
}


#ifdef UG_FOR_LUA
template <typename TDomain, typename TAlgebra>
void NavierStokesInflow<TDomain,TAlgebra>::
add(const char* fctName, const char* subsetsBND)
{
	add(LuaUserDataFactory<MathVector<dim>, dim>::create(fctName), subsetsBND);
}
#endif


////////////////////////////////////////////////////////////////////////////////
//	NavierStokesWall
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
NavierStokesWall<TDomain,TAlgebra>::
NavierStokesWall(SmartPtr< NavierStokes<TDomain> > spMaster)
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

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__NAVIER_STOKES_BND_IMPL__ */
