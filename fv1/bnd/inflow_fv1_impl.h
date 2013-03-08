/*
 * inflow_fv1_impl.h
 *
 *  Created on: 12.04.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__NAVIER_STOKES__BND__INFLOW_FV1_IMPL__
#define __H__UG__NAVIER_STOKES__BND__INFLOW_FV1_IMPL__

#include "inflow_fv1.h"
#include "lib_disc/spatial_disc/elem_disc/neumann_boundary/fv1/neumann_boundary_fv1.h"

namespace ug{
namespace NavierStokes{

template <typename TDomain, typename TAlgebra>
NavierStokesInflowFV1<TDomain,TAlgebra>::
NavierStokesInflowFV1(SmartPtr< NavierStokesFV1<TDomain> > spMaster)
	: m_spNeumannDisc(new NeumannBoundaryFV1<TDomain>(spMaster->symb_fcts()[dim].c_str())),
	  m_spDirichletConstraint(new DirichletBoundary<TDomain,TAlgebra>),
	  m_spMaster(spMaster)
{
	const std::vector<std::string>& vFctName = m_spMaster->symb_fcts();

	if(vFctName.size() != TDomain::dim + 1)
		UG_THROW("NavierStokesInflow::set_functions: This Boundary "
				"Condition works on exactly dim+1 (velocity+pressure) "
				"components, but "<<vFctName.size()<<"components given.");
}

template <typename TDomain, typename TAlgebra>
void NavierStokesInflowFV1<TDomain,TAlgebra>::
add(SmartPtr<UserData<MathVector<dim>, dim> > user, const char* subsetsBND)
{
	const std::vector<std::string>& vFctName = m_spMaster->symb_fcts();
	if(vFctName.empty())
		UG_THROW("NavierStokesInflow::add: Symbolic names for"
				" velocity and pressure not set. Please set them first.");

	std::string innerSubsets;
	for(size_t s = 0; s < m_spMaster->symb_subsets().size(); ++s){
		if(s > 0) innerSubsets.append(",");
		innerSubsets.append(m_spMaster->symb_subsets()[s]);
	}


	m_spNeumannDisc->add(user, subsetsBND, innerSubsets.c_str());

	std::string velNames;
	for(int i=0;i<TDomain::dim; ++i)
	{
		if(i>0) velNames.append(",");
		velNames.append(vFctName[i]);
	}
	m_spDirichletConstraint->add(user, velNames.c_str(), subsetsBND);
}

} // namespace NavierStokes
} // end namespace ug

#endif /* __H__UG__NAVIER_STOKES__BND__INFLOW_FV1_IMPL__ */