/*
 * no_normal_stress_outflow.cpp
 *
 *  Created on: 27.03.2011
 */

#include "no_normal_stress_outflow.h"

#include "common/util/provider.h"
#include "lib_disc/spatial_disc/disc_util/finite_volume_geometry.h"

namespace ug{


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Provide a generic implementation for all elements
//  (since this discretization can be implemented in a generic way)
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/**
 * converts the subset names where the BC is imposed to the corresponding subset
 * indices (i.e. m_vScheduledBndSubSets -> m_vBndSubSetIndex):
 */
template<typename TDomain>
void
FVNavierStokesNoNormalStressOutflow<TDomain>::
extract_scheduled_data()
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
		}UG_CATCH_THROW("'FVNavierStokesNoNormalStressOutflow:extract_scheduled_data':"
						" Subsets '" <<m_vScheduledBndSubSets[i].c_str() <<"' not"
						" all contained in ApproximationSpace.");
	
	//	get subsethandler
		const ISubsetHandler& rSH = *this->function_pattern().subset_handler();

	// 	loop subsets
		for(size_t si = 0; si < subsetGroup.num_subsets(); ++si)
		{
		//	get subset index
			const int subsetIndex = subsetGroup[si];
		
		//	check that subsetIndex is valid
			if(subsetIndex < 0 || subsetIndex >= rSH.num_subsets())
			{
				UG_LOG("ERROR in 'FVNavierStokesNoNormalStressOutflow:extract_scheduled_data':"
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
void
FVNavierStokesNoNormalStressOutflow<TDomain>::
add
(
	const char* subsets // string with the ','-separated names of the subsets
)
{
	m_vScheduledBndSubSets.push_back(subsets);

	if(this->fct_pattern_set()) extract_scheduled_data();
}

/**
 * Prepares the element loop for a given element type: computes the FV-geo, ...
 */
template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void FVNavierStokesNoNormalStressOutflow<TDomain>::
prepare_element_loop()
{
//	register subsetIndex at Geometry
	static TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();

// 	Only first order implementation
	if(!(TFVGeom<TElem, dim>::order == 1))
		UG_THROW_FATAL("Only first order implementation, but other Finite Volume"
						" Geometry set.");

//	check, that kinematic Viscosity has been set
	if(!m_imKinViscosity.data_given())
		UG_THROW_FATAL("FVNavierStokesNoNormalStressOutflow::prepare_element_loop:"
						" Kinematic Viscosity has not been set, but is required.\n");

//	set local positions for imports
	typedef typename reference_element_traits<TElem>::reference_element_type
																ref_elem_type;
	static const int refDim = ref_elem_type::dim;

	if(!TFVGeom<TElem, dim>::usesHangingNodes)
	{
		TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();
		m_imKinViscosity.template set_local_ips<refDim>(geo.scvf_local_ips(),
		                                                geo.num_scvf_ips());
	}

//	request the subset indices as boundary subset. This will force the
//	creation of boundary subsets when calling geo.update
	typename std::vector<int>::const_iterator subsetIter;
	for(subsetIter = m_vBndSubSetIndex.begin();
			subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
		geo.add_boundary_subset(*subsetIter);
}

template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void FVNavierStokesNoNormalStressOutflow<TDomain>::
finish_element_loop()
{
	static TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();

//	remove the bnd subsets
	typename std::vector<int>::const_iterator subsetIter;
	for(subsetIter = m_vBndSubSetIndex.begin();
			subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
		geo.remove_boundary_subset(*subsetIter);
}


template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void FVNavierStokesNoNormalStressOutflow<TDomain>::
prepare_element(TElem* elem, const LocalVector& u)
{
//	get corners
	m_vCornerCoords = this->template element_corners<TElem>(elem);

// 	Update Geometry for this element
	TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();
	if(!geo.update(elem, &m_vCornerCoords[0], &(this->subset_handler())))
		UG_THROW_FATAL("FVNavierStokesNoNormalStressOutflow::prepare_element:"
						" Cannot update Finite Volume Geometry.\n");

//	set local positions for imports
	if(TFVGeom<TElem, dim>::usesHangingNodes)
	{
	//	set local positions for imports
		typedef typename reference_element_traits<TElem>::reference_element_type
																	ref_elem_type;
		static const int refDim = ref_elem_type::dim;

	//	request ip series
		m_imKinViscosity.template set_local_ips<refDim>(geo.scvf_local_ips(),
		                                                geo.num_scvf_ips());
	}

//	set global positions for imports
	m_imKinViscosity.set_global_ips(geo.scvf_global_ips(), geo.num_scvf_ips());
}

template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void FVNavierStokesNoNormalStressOutflow<TDomain>::
ass_JA_elem(LocalMatrix& J, const LocalVector& u)
{
// 	Only first order implementation
	UG_ASSERT((TFVGeom<TElem, dim>::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	static const TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();

// 	loop registered boundary segments
	typename std::vector<int>::const_iterator subsetIter;
	for(subsetIter = m_vBndSubSetIndex.begin(); subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
	{
	//	get subset index corresponding to boundary
		const int bndSubset = *subsetIter;
		
	// 	loop Boundary Faces
		for(size_t i = 0; i < geo.num_bf(bndSubset); ++i)
		{
		// get current BF
			const typename TFVGeom<TElem, dim>::BF& bf = geo.bf(bndSubset, i);

		// 	loop shape functions
			for(size_t sh = 0; sh < bf.num_sh(); ++sh)
			{
			//	1. Compute the total flux
				MathMatrix<dim,dim> diffFlux, tang_diffFlux;
		
			//	Add \nabla u
				MatSet (diffFlux, 0);
				MatDiagSet (diffFlux, VecDot (bf.global_grad (sh), bf.normal ()));
		
			//	Add (\nabla u)^T
				if(!m_spMaster->get_laplace())
					for (size_t d1 = 0; d1 < (size_t)dim; ++d1)
						for(size_t d2 = 0; d2 < (size_t)dim; ++d2)
							diffFlux(d1,d2) += bf.global_grad (sh) [d1] * bf.normal () [d2];
			
			//	2. Subtract the normal part:
				tang_diffFlux = diffFlux;
				for (size_t d1 = 0; d1 < (size_t)dim; ++d1)
					for(size_t d2 = 0; d2 < (size_t)dim; ++d2)
						for (size_t k = 0; k < (size_t)dim; ++k)
							tang_diffFlux(d1,d2) -= bf.normal () [d1]
								* diffFlux(d2,k) * bf.normal () [k];
		
			//	3. Scale by viscosity
				tang_diffFlux *= - m_imKinViscosity[i];
		
			//	4. Add flux to local defect
				for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
					for(size_t d2 = 0; d2 < (size_t)dim; ++d2)
						J(d1, bf.node_id(), d2, sh) += tang_diffFlux (d1, d2);
			}
		}
	}
}

template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void FVNavierStokesNoNormalStressOutflow<TDomain>::
ass_dA_elem(LocalVector& d, const LocalVector& u)
{
// 	Only first order implemented
	UG_ASSERT((TFVGeom<TElem, dim>::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	static const TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();

// 	loop registered boundary segments
	typename std::vector<int>::const_iterator subsetIter;
	for(subsetIter = m_vBndSubSetIndex.begin(); subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
	{
	//	get subset index corresponding to boundary
		const int bndSubset = *subsetIter;

	// 	loop Boundary Faces
		for(size_t i = 0; i < geo.num_bf(bndSubset); ++i)
		{
		// get current BF
			const typename TFVGeom<TElem, dim>::BF& bf = geo.bf(bndSubset, i);

		// 	1. Interpolate Functional Matrix of velocity at ip
			MathMatrix<dim, dim> gradVel;
			for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
				for(size_t d2 = 0; d2 < (size_t)dim; ++d2)
				{
				//	sum up contributions of each shape
					gradVel(d1, d2) = 0.0;
					for(size_t sh = 0; sh < bf.num_sh(); ++sh)
						gradVel(d1, d2) += bf.global_grad(sh)[d2] * u(d1, sh);
				}
	
		//	2. Compute the total flux
			MathVector<dim> diffFlux;
	
		//	Add (\nabla u) \cdot \vec{n}
			MatVecMult(diffFlux, gradVel, bf.normal());
	
		//	Add (\nabla u)^T \cdot \vec{n}
			if(!m_spMaster->get_laplace())
				TransposedMatVecMultAdd(diffFlux, gradVel, bf.normal());
		
		//	3. Subtract the normal part:
			VecScaleAppend (diffFlux, - VecDot (diffFlux, bf.normal ()), bf.normal ());
	
		//	4. Scale by viscosity
			diffFlux *= - m_imKinViscosity[i];
	
		//	5. Add flux to local defect
			for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
				d(d1, bf.node_id()) += diffFlux[d1];
		}
	}
}


template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void FVNavierStokesNoNormalStressOutflow<TDomain>::
ass_JM_elem(LocalMatrix& J, const LocalVector& u)
{}


template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void FVNavierStokesNoNormalStressOutflow<TDomain>::
ass_dM_elem(LocalVector& d, const LocalVector& u)
{}


template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void FVNavierStokesNoNormalStressOutflow<TDomain>::
ass_rhs_elem(LocalVector& d)
{}


////////////////////////////////////////////////////////////////////////////////
//	Constructor
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
FVNavierStokesNoNormalStressOutflow<TDomain>::
FVNavierStokesNoNormalStressOutflow(SmartPtr< NavierStokes<TDomain> > spMaster)
: IDomainElemDisc<TDomain>(spMaster->symb_fcts(), spMaster->symb_subsets()), m_spMaster (spMaster)
{
//	check number of functions
	if(this->num_fct() != dim+1)
		UG_THROW_FATAL("Wrong number of functions: The ElemDisc 'NavierStokes'"
					   " needs exactly "<<dim+1<<" symbolic function.");
	
//	yet no boundary subsets
	m_vBndSubSetIndex.clear ();

//	register imports
	register_import(m_imKinViscosity);

//	initialize the imports from the master discretization
	m_imKinViscosity.set_data(spMaster->get_kinematic_viscosity_data ());

//	register assemble functions
	register_all_fv1_funcs(false);
}


////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

// register for 1D
template<typename TDomain>
void
FVNavierStokesNoNormalStressOutflow<TDomain>::
register_all_fv1_funcs(bool bHang)
{
//	get all grid element types in this dimension and below
	typedef typename domain_traits<dim>::AllElemList ElemList;

//	switch assemble functions
	if(!bHang) boost::mpl::for_each<ElemList>( RegisterFV1<FV1Geometry>(this) );
	else throw(UGFatalError("Not implemented."));
}

template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void
FVNavierStokesNoNormalStressOutflow<TDomain>::
register_fv1_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->enable_fast_ass_elem(true);
	set_prep_elem_loop_fct(id, &T::template prepare_element_loop<TElem, TFVGeom>);
	set_prep_elem_fct(	 id, &T::template prepare_element<TElem, TFVGeom>);
	set_fsh_elem_loop_fct( id, &T::template finish_element_loop<TElem, TFVGeom>);
	set_ass_JA_elem_fct(		 id, &T::template ass_JA_elem<TElem, TFVGeom>);
	set_ass_JM_elem_fct(		 id, &T::template ass_JM_elem<TElem, TFVGeom>);
	set_ass_dA_elem_fct(		 id, &T::template ass_dA_elem<TElem, TFVGeom>);
	set_ass_dM_elem_fct(		 id, &T::template ass_dM_elem<TElem, TFVGeom>);
	set_ass_rhs_elem_fct(	 id, &T::template ass_rhs_elem<TElem, TFVGeom>);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

template class FVNavierStokesNoNormalStressOutflow<Domain1d>;
template class FVNavierStokesNoNormalStressOutflow<Domain2d>;
template class FVNavierStokesNoNormalStressOutflow<Domain3d>;


} // namespace ug
