/*
 * no_normal_stress_outflow_common.h
 *
 *  Created on: 27.03.2012
 *  D. Logashenko, A. Vogel
 */

#include "no_normal_stress_outflow.h"

#include "common/util/provider.h"
#include "lib_disc/spatial_disc/disc_util/finite_volume_geometry.h"

namespace ug{
namespace NavierStokes{

template<typename TDomain>
bool NavierStokesNoNormalStressOutflow<TDomain>::
request_finite_element_id(const std::vector<LFEID>& vLfeID)
{
//	check number
	if(vLfeID.size() != dim+1) return false;

//	check trial spaces
	if(m_spMaster->disc_scheme() == "stab"){
		for(size_t i = 0; i < vLfeID.size(); ++i)
			if(vLfeID[i] != LFEID(LFEID::LAGRANGE, 1)) return false;
	}
	else if(m_spMaster->disc_scheme() == "staggered"){
		for(int i = 0; i < dim; ++i)
			if(vLfeID[i] != LFEID(LFEID::CROUZEIX_RAVIART, 1)) return false;

		if(vLfeID[dim] != LFEID(LFEID::PIECEWISE_CONSTANT, 0)) return false;
	}
	return true;
}

template<typename TDomain>
bool NavierStokesNoNormalStressOutflow<TDomain>::
request_non_regular_grid(bool bNonRegular)
{
//	switch, which assemble functions to use.
	if(bNonRegular)
	{
		UG_LOG("ERROR in 'NavierStokes::request_non_regular_grid':"
				" Non-regular grid not implemented.\n");
		return false;
	}

//	this disc supports regular grids
	return true;
}

/**
 * converts the subset names where the BC is imposed to the corresponding subset
 * indices (i.e. m_vScheduledBndSubSets -> m_vBndSubSetIndex):
 */
template<typename TDomain>
void NavierStokesNoNormalStressOutflow<TDomain>::extract_scheduled_data()
{
//	clear all extracted data
	m_vBndSubSetIndex.clear();

//	loop all scheduled subsets
	for(size_t i = 0; i < m_vScheduledBndSubSets.size(); ++i)
	{
	//	create Subset Group
		SubsetGroup subsetGroup;

	//	convert strings
		try{
			subsetGroup = this->approx_space()->subset_grp_by_name(m_vScheduledBndSubSets[i].c_str());
		}UG_CATCH_THROW("'NavierStokesNoNormalStressOutflow:extract_scheduled_data':"
						" Subsets '" <<m_vScheduledBndSubSets[i].c_str() <<"' not"
						" all contained in ApproximationSpace.");
	
	//	get subsethandler
		const ISubsetHandler& rSH = *this->function_pattern().subset_handler();

	// 	loop subsets
		for(size_t si = 0; si < subsetGroup.size(); ++si)
		{
		//	get subset index
			const int subsetIndex = subsetGroup[si];
		
		//	check that subsetIndex is valid
			if(subsetIndex < 0 || subsetIndex >= rSH.num_subsets())
			{
				UG_LOG("ERROR in 'NavierStokesNoNormalStressOutflow:extract_scheduled_data':"
						" Invalid subset Index " << subsetIndex <<
						". (Valid is 0, .. , " << rSH.num_subsets() <<").\n");
				return;
			}
		
		// save the index
			m_vBndSubSetIndex.push_back(subsetIndex);
		}
	}
}

/**
 * The add method for the boundary subsets:
 */
template<typename TDomain>
void NavierStokesNoNormalStressOutflow<TDomain>::add
(
	const char* subsets // string with the ','-separated names of the subsets
)
{
	m_vScheduledBndSubSets.push_back(subsets);

	if(this->fct_pattern_set()) extract_scheduled_data();
}

////////////////////////////////////////////////////////////////////////////////
//	Constructor
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
NavierStokesNoNormalStressOutflow<TDomain>::
NavierStokesNoNormalStressOutflow(SmartPtr< NavierStokes<TDomain> > spMaster)
: IDomainElemDisc<TDomain>(spMaster->symb_fcts(), spMaster->symb_subsets()), m_spMaster (spMaster)
{
//	check number of functions
	if(this->num_fct() != dim+1)
		UG_THROW("Wrong number of functions: The ElemDisc 'NavierStokes'"
					   " needs exactly "<<dim+1<<" symbolic function.");
	
//	yet no boundary subsets
	m_vBndSubSetIndex.clear ();

//	register imports
	register_import(m_imKinViscosity);
	register_import(m_imDensity);

//	initialize the imports from the master discretization
	m_imKinViscosity.set_data(spMaster->get_kinematic_viscosity_data ());
	m_imDensity.set_data(spMaster->get_density ());

//	register assemble functions
	if(spMaster->disc_scheme() == "stab") this->register_all_fv1_funcs(false);
	else if(spMaster->disc_scheme() == "staggered") this->register_all_cr_funcs(false);
	else UG_THROW("NavierStokesNoNormalStressOutflow: Disc scheme not supported.");
}

} // namespace NavierStokes
} // namespace ug
