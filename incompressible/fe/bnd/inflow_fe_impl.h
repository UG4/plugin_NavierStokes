/*
 * Copyright (c) 2013:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG__NAVIER_STOKES__BND__INFLOW_FE_IMPL__
#define __H__UG__NAVIER_STOKES__BND__INFLOW_FE_IMPL__

#include "inflow_fe.h"
#include "lib_disc/spatial_disc/elem_disc/neumann_boundary/fe/neumann_boundary_fe.h"

namespace ug{
namespace NavierStokes{

template <typename TDomain, typename TAlgebra>
NavierStokesInflowFE<TDomain,TAlgebra>::
NavierStokesInflowFE(SmartPtr< NavierStokesFE<TDomain> > spMaster)
	: m_spNeumannDisc(new NeumannBoundaryFE<TDomain>(spMaster->symb_fcts()[dim].c_str())),
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
void NavierStokesInflowFE<TDomain,TAlgebra>::
add(SmartPtr<CplUserData<MathVector<dim>, dim> > user, const char* subsetsBND)
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
	m_spDirichletConstraint->add(user.template cast_dynamic<UserData<MathVector<dim>, dim> >(), velNames.c_str(), subsetsBND);
}

} // namespace NavierStokes
} // end namespace ug

#endif /* __H__UG__NAVIER_STOKES__BND__INFLOW_FE_IMPL__ */
