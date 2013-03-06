/*
 * inflow_impl.h
 *
 *  Created on: 12.04.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__NAVIER_STOKES__BND__INFLOW_IMPL__
#define __H__UG__NAVIER_STOKES__BND__INFLOW_IMPL__

#include "inflow.h"
#include "lib_disc/spatial_disc/elem_disc/neumann_boundary/fv1/neumann_boundary_fv1.h"
#include "lib_disc/spatial_disc/elem_disc/neumann_boundary/fv/neumann_boundary_fv.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{
namespace NavierStokes{

template <typename TDomain, typename TAlgebra>
NavierStokesInflow<TDomain,TAlgebra>::
NavierStokesInflow(SmartPtr< NavierStokes<TDomain> > spMaster)
	: m_spNeumannDisc(NULL),
	  m_spDirichletConstraint(new DirichletBoundary<TDomain,TAlgebra>),
	  m_spMaster(spMaster)
{
	m_vFctName = m_spMaster->symb_fcts();

	if(m_vFctName.size() != TDomain::dim + 1)
		UG_THROW("NavierStokesInflow::set_functions: This Boundary "
				"Condition works on exactly dim+1 (velocity+pressure) "
				"components, but "<<m_vFctName.size()<<"components given.");

	if(m_spMaster->disc_scheme() == "fv1"){
		m_spNeumannDisc = SmartPtr<NeumannBoundaryFV1<TDomain> >(new NeumannBoundaryFV1<TDomain>(m_vFctName[dim].c_str()));
	}
	else if(m_spMaster->disc_scheme() == "fv"){
		m_spNeumannDisc = SmartPtr<NeumannBoundaryFV<TDomain> >(new NeumannBoundaryFV<TDomain>(m_vFctName[dim].c_str()));
	}
	else if(m_spMaster->disc_scheme() != "staggered"){
		// nothing
	}
	else{
		UG_THROW("NavierStokesInflow: Master Disc scheme not recognized.");
	}
}

template <typename TDomain, typename TAlgebra>
void NavierStokesInflow<TDomain,TAlgebra>::
add(SmartPtr<UserData<MathVector<dim>, dim> > user, const char* subsetsBND)
{
	if(m_vFctName.empty())
		UG_THROW("NavierStokesInflow::add: Symbolic names for"
				" velocity and pressure not set. Please set them first.");

	std::string innerSubsets;
	for(size_t s = 0; s < m_spMaster->symb_subsets().size(); ++s){
		if(s > 0) innerSubsets.append(",");
		innerSubsets.append(m_spMaster->symb_subsets()[s]);
	}


	if(m_spMaster->disc_scheme() == "fv1"){
		m_spNeumannDisc->add(user, subsetsBND, innerSubsets.c_str());
	}
	else if(m_spMaster->disc_scheme() == "fv"){
		m_spNeumannDisc->add(user, subsetsBND, innerSubsets.c_str());
	}
	else if(m_spMaster->disc_scheme() != "staggered"){
		// nothing
	}
	else{
		UG_THROW("NavierStokesInflow: Master Disc scheme not recognized.");
	}

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

} // namespace NavierStokes
} // end namespace ug

#endif /* __H__UG__NAVIER_STOKES__BND__INFLOW_IMPL__ */
