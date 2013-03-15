/*
 * no_normal_stress_outflow_base.cpp
 *
 *  Created on: 27.03.2012
 *  D. Logashenko, A. Vogel
 */

#include "no_normal_stress_outflow_base.h"

#include "common/util/provider.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"

namespace ug{
namespace NavierStokes{

/**
 * converts the subset names where the BC is imposed to the corresponding subset
 * indices (i.e. m_vScheduledBndSubSets -> m_vBndSubSetIndex):
 */
template<typename TDomain>
void NavierStokesNoNormalStressOutflowBase<TDomain>::extract_scheduled_data()
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
void NavierStokesNoNormalStressOutflowBase<TDomain>::add
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
NavierStokesNoNormalStressOutflowBase<TDomain>::
NavierStokesNoNormalStressOutflowBase(SmartPtr< NavierStokesBase<TDomain> > spMaster)
	: IDomainElemDisc<TDomain>(spMaster->symb_fcts(), spMaster->symb_subsets()),
	  m_spMaster (spMaster)
{
//	check number of functions
	if(this->num_fct() != dim+1)
		UG_THROW("Wrong number of functions: The ElemDisc 'NavierStokes'"
					   " needs exactly "<<dim+1<<" symbolic function.");
	
//	yet no boundary subsets
	m_vBndSubSetIndex.clear ();
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class NavierStokesNoNormalStressOutflowBase<Domain1d>;
#endif
#ifdef UG_DIM_2
template class NavierStokesNoNormalStressOutflowBase<Domain2d>;
#endif
#ifdef UG_DIM_3
template class NavierStokesNoNormalStressOutflowBase<Domain3d>;
#endif

} // namespace NavierStokes
} // namespace ug
