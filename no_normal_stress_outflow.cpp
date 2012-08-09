/*
 * no_normal_stress_outflow.cpp
 *
 *  Created on: 27.03.2012
 *  D. Logashenko, A. Vogel
 */

#include "no_normal_stress_outflow.h"

#include "common/util/provider.h"
#include "lib_disc/spatial_disc/disc_util/finite_volume_geometry.h"

namespace ug{
namespace NavierStokes{


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
 * Note that there are separate loops for every type of the grid elements.
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
		UG_THROW("Only first order implementation, but other Finite Volume"
						" Geometry set.");

//	check, that kinematic Viscosity has been set
	if(!m_imKinViscosity.data_given())
		UG_THROW("FVNavierStokesNoNormalStressOutflow::prepare_element_loop:"
						" Kinematic Viscosity has not been set, but is required.\n");

//	request the subset indices as boundary subset. This will force the
//	creation of boundary subsets when calling geo.update
	typename std::vector<int>::const_iterator subsetIter;
	for(subsetIter = m_vBndSubSetIndex.begin();
			subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
		geo.add_boundary_subset(*subsetIter);
}

/**
 * Finalizes the element loop for a given element type.
 */
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


/**
 * General initializations of a given grid element for the assembling.
 */
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
		UG_THROW("FVNavierStokesNoNormalStressOutflow::prepare_element:"
						" Cannot update Finite Volume Geometry.\n");

//	find and set the local and the global positions of the IPs for imports
	typedef typename TFVGeom<TElem, dim>::BF BF;
	typename std::vector<int>::const_iterator subsetIter;
	
	m_vLocIP.clear(); m_vGloIP.clear();
	for(subsetIter = m_vBndSubSetIndex.begin();
			subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
	{
		const int bndSubset = *subsetIter;
		if(geo.num_bf(bndSubset) == 0) continue;
		const std::vector<BF>& vBF = geo.bf(bndSubset);
		for(size_t i = 0; i < vBF.size(); ++i)
		{
			m_vLocIP.push_back(vBF[i].local_ip());
			m_vGloIP.push_back(vBF[i].global_ip());
		}
	}
	// REMARK: The loop above determines the ordering of the integration points:
	// The "outer ordering" corresponds to the ordering of the subsets in
	// m_vBndSubSetIndex, and "inside" of this ordering, the ip's are ordered
	// according to the order of the boundary faces in the FV geometry structure.

	m_imKinViscosity.set_local_ips(&m_vLocIP[0], m_vLocIP.size());
	m_imKinViscosity.set_global_ips(&m_vGloIP[0], m_vGloIP.size());
}

/// Assembling of the diffusive flux (due to the viscosity) in the Jacobian of the momentum eq.
template<typename TDomain>
template<typename BF>
void FVNavierStokesNoNormalStressOutflow<TDomain>::
ass_diffusive_flux_Jac
(
	const size_t ip, // index of the integration point (for the viscosity)
	const BF& bf, // boundary face to assemble
	LocalMatrix& J, // local Jacobian to update
	const LocalVector& u // local solution
)
{
	MathMatrix<dim,dim> diffFlux, tang_diffFlux;
	MathVector<dim> normalStress;
	
	for(size_t sh = 0; sh < bf.num_sh(); ++sh) // loop shape functions
	{
	//	1. Compute the total flux
	//	- add \nabla u
		MatSet (diffFlux, 0);
		MatDiagSet (diffFlux, VecDot (bf.global_grad(sh), bf.normal()));
	
	//	- add (\nabla u)^T
		if(!m_spMaster->get_laplace())
			for (size_t d1 = 0; d1 < (size_t)dim; ++d1)
				for (size_t d2 = 0; d2 < (size_t)dim; ++d2)
					diffFlux(d1,d2) += bf.global_grad(sh)[d1] * bf.normal()[d2];
	
	//	2. Subtract the normal part:
		tang_diffFlux = diffFlux;
		TransposedMatVecMult(normalStress, diffFlux, bf.normal ());
		for (size_t d2 = 0; d2 < (size_t)dim; ++d2)
			for (size_t d1 = 0; d1 < (size_t)dim; ++d1)
				tang_diffFlux(d1,d2) -= bf.normal()[d1] * normalStress[d2];
	
	//	3. Scale by viscosity
		tang_diffFlux *= - m_imKinViscosity[ip];
	
	//	4. Add flux to local Jacobian
		for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
			for(size_t d2 = 0; d2 < (size_t)dim; ++d2)
				J(d1, bf.node_id(), d2, sh) += tang_diffFlux (d1, d2);
	}
}

/// Assembling of the diffusive flux (due to the viscosity) in the defect of the momentum eq.
template<typename TDomain>
template<typename BF>
void FVNavierStokesNoNormalStressOutflow<TDomain>::
ass_diffusive_flux_defect
(
	const size_t ip, // index of the integration point (for the viscosity)
	const BF& bf, // boundary face to assemble
	LocalVector& d, // local defect to update
	const LocalVector& u // local solution
)
{
	MathMatrix<dim, dim> gradVel;
	MathVector<dim> diffFlux;
	
// 	1. Get the gradient of the velocity at ip
	for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
		for(size_t d2 = 0; d2 < (size_t)dim; ++d2)
		{
		//	sum up contributions of each shape
			gradVel(d1, d2) = 0.0;
			for(size_t sh = 0; sh < bf.num_sh(); ++sh)
				gradVel(d1, d2) += bf.global_grad(sh)[d2] * u(d1, sh);
		}

//	2. Compute the total flux

//	- add (\nabla u) \cdot \vec{n}
	MatVecMult(diffFlux, gradVel, bf.normal());

//	- add (\nabla u)^T \cdot \vec{n}
	if(!m_spMaster->get_laplace())
		TransposedMatVecMultAdd(diffFlux, gradVel, bf.normal());

//	3. Subtract the normal part:
	VecScaleAppend (diffFlux, - VecDot (diffFlux, bf.normal()), bf.normal());

//	A4. Scale by viscosity
	diffFlux *= - m_imKinViscosity[ip];

//	5. Add flux to local defect
	for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
		d(d1, bf.node_id()) += diffFlux[d1];
}

/// Assembling of the convective flux (due to the quadratic inertial term) in the Jacobian of the momentum eq.
template<typename TDomain>
template<typename BF>
void FVNavierStokesNoNormalStressOutflow<TDomain>::
ass_convective_flux_Jac
(
	const size_t ip, // index of the integration point (for the density)
	const BF& bf, // boundary face to assemble
	LocalMatrix& J, // local Jacobian to update
	const LocalVector& u // local solution
)
{
	MathVector<dim> StdVel;
	number old_momentum_flux, t;
	
// The convection velocity according to the current approximation:
	for(size_t sh = 0; sh < bf.num_sh(); ++sh)
		for(size_t d1 = 0; d1 < (size_t) dim; ++d1)
			StdVel[d1] += u(d1, sh) * bf.shape(sh);
	old_momentum_flux = VecDot (StdVel, bf.normal ());
	// Multiply old_momentum_flux by the density here!
	
// We assume that there should be no inflow through the outflow boundary:
	if (old_momentum_flux < 0)
		old_momentum_flux = 0;
	
//	Add flux to local Jacobian
	for(size_t sh = 0; sh < bf.num_sh(); ++sh)
	{
		t = old_momentum_flux * bf.shape(sh);
		for(size_t d1 = 0; d1 < (size_t) dim; ++d1)
			J(d1, bf.node_id(), d1, sh) += t;
	}
}

/// Assembling of the convective flux (due to the quadratic inertial term) in the defect of the momentum eq.
template<typename TDomain>
template<typename BF>
void FVNavierStokesNoNormalStressOutflow<TDomain>::
ass_convective_flux_defect
(
	const size_t ip, // index of the integration point (for the density)
	const BF& bf, // boundary face to assemble
	LocalVector& d, // local defect to update
	const LocalVector& u // local solution
)
{
	MathVector<dim> StdVel;
	number old_momentum_flux;
	
// The convection velocity according to the current approximation:
	for(size_t sh = 0; sh < bf.num_sh(); ++sh)
		for(size_t d1 = 0; d1 < (size_t) dim; ++d1)
			StdVel[d1] += u(d1, sh) * bf.shape(sh);
	old_momentum_flux = VecDot (StdVel, bf.normal ());
	// Multiply old_momentum_flux by the density here!
	
// We assume that there should be no inflow through the outflow boundary:
	if (old_momentum_flux < 0)
		old_momentum_flux = 0;
	
// Add the flux to the defect:
	for(size_t d1 = 0; d1 < (size_t) dim; ++d1)
		d(d1, bf.node_id()) += old_momentum_flux * StdVel[d1];
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
	typedef typename TFVGeom<TElem, dim>::BF BF;

// 	loop registered boundary segments
	typename std::vector<int>::const_iterator subsetIter;
	size_t ip = 0;
	for(subsetIter = m_vBndSubSetIndex.begin();
		subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
	{
	//	get subset index corresponding to boundary
		const int bndSubset = *subsetIter;
		
	//	get the list of the ip's:
		if(geo.num_bf(bndSubset) == 0) continue;
		const std::vector<BF>& vBF = geo.bf(bndSubset);

	// 	loop the boundary faces
		typename std::vector<BF>::const_iterator bf;
		for(bf = vBF.begin(); bf != vBF.end(); ++bf)
		{
		//	A. The momentum equation:
			ass_diffusive_flux_Jac<BF> (ip, *bf, J, u);
			if (!m_spMaster->get_stokes ())
				ass_convective_flux_Jac<BF> (ip, *bf, J, u);
			
		//	B. The continuity equation
			for(size_t sh = 0; sh < bf->num_sh(); ++sh) // loop shape functions
				for (size_t d2 = 0; d2 < (size_t)dim; ++d2)
					J(_P_, bf->node_id (), d2, sh) += bf->shape(sh) * bf->normal()[d2];
		
		// Next IP:
			ip++;
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
	typedef typename TFVGeom<TElem, dim>::BF BF;

// 	loop registered boundary segments
	typename std::vector<int>::const_iterator subsetIter;
	size_t ip = 0;
	for(subsetIter = m_vBndSubSetIndex.begin();
		subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
	{
	//	get subset index corresponding to boundary
		const int bndSubset = *subsetIter;
		
	//	get the list of the ip's:
		if(geo.num_bf(bndSubset) == 0) continue;
		const std::vector<BF>& vBF = geo.bf(bndSubset);

	// 	loop the boundary faces
		typename std::vector<BF>::const_iterator bf;
		for(bf = vBF.begin(); bf != vBF.end(); ++bf)
		{
		// A. Momentum equation:
			ass_diffusive_flux_defect<BF> (ip, *bf, d, u);
			if (!m_spMaster->get_stokes ())
				ass_convective_flux_defect<BF> (ip, *bf, d, u);
		
		// B. Continuity equation:
			{
				MathVector<dim> stdVel;
				VecSet (stdVel, 0);
				for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
					for(size_t sh = 0; sh < bf->num_sh(); ++sh)
						stdVel[d1] += u(d1, sh) * bf->shape(sh);
				d(_P_, bf->node_id()) += VecDot (stdVel, bf->normal());
			}
		
		// Next IP:
			ip++;
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
		UG_THROW("Wrong number of functions: The ElemDisc 'NavierStokes'"
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
	typedef typename domain_traits<dim>::DimElemList ElemList;
	// REMARK: Note that we register this boundary condition only
	// for the full-dimensional elements (DimElemList instead of AllElemList).

//	switch assemble functions
	if(!bHang) boost::mpl::for_each<ElemList>( RegisterFV1<FV1Geometry>(this) );
	else throw(UGError("Not implemented."));
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


} // namespace NavierStokes
} // namespace ug
