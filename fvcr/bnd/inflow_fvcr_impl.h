/*
 * inflow_fvcr_impl.h
 *
 *  Created on: 12.04.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__NAVIER_STOKES__BND__INFLOW_FVCR_IMPL__
#define __H__UG__NAVIER_STOKES__BND__INFLOW_FVCR_IMPL__

#include "inflow_fvcr.h"

namespace ug{
namespace NavierStokes{

template <typename TDomain, typename TAlgebra>
NavierStokesInflowFVCR<TDomain,TAlgebra>::
NavierStokesInflowFVCR(SmartPtr< NavierStokesFVCR<TDomain> > spMaster)
	: m_spDirichletConstraint(new DirichletBoundary<TDomain,TAlgebra>),
	  m_spMaster(spMaster)
{
	const std::vector<std::string>& vFctName = m_spMaster->symb_fcts();

	if(vFctName.size() != TDomain::dim + 1)
		UG_THROW("NavierStokesInflow::set_functions: This Boundary "
				"Condition works on exactly dim+1 (velocity+pressure) "
				"components, but "<<vFctName.size()<<"components given.");
}

template <typename TDomain, typename TAlgebra>
void NavierStokesInflowFVCR<TDomain,TAlgebra>::
add(SmartPtr<CplUserData<MathVector<dim>, dim> > user, const char* subsetsBND)
{
	const std::vector<std::string>& vFctName = m_spMaster->symb_fcts();
	if(vFctName.empty())
		UG_THROW("NavierStokesInflow::add: Symbolic names for"
				" velocity and pressure not set. Please set them first.");

	std::string velNames;
	for(int i=0;i<TDomain::dim; ++i)
	{
		if(i>0) velNames.append(",");
		velNames.append(vFctName[i]);
	}
	m_spDirichletConstraint->add(user.template cast_dynamic<UserData<MathVector<dim>, dim> >(), velNames.c_str(), subsetsBND);
}

} // namespace NavierStokes
} // end namespace ug

#endif /* __H__UG__NAVIER_STOKES__BND__INFLOW_FVCR_IMPL__ */
