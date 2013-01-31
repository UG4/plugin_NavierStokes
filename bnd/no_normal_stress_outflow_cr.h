/*
 * no_normal_stress_outflow_cr.h
 *
 *  Created on: 23.07.2012
 *      Author: Christian Wehner
 */

#ifndef __H__UG__NAVIER_STOKES__BND__NO_NORMAL_STRESS_OUTFLOW_CR__
#define __H__UG__NAVIER_STOKES__BND__NO_NORMAL_STRESS_OUTFLOW_CR__

#include "no_normal_stress_outflow.h"

namespace ug{
namespace NavierStokes{

/**
 * Prepares the element loop for a given element type: computes the FV-geo, ...
 * Note that there are separate loops for every type of the grid elements.
 */
template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void NavierStokesNoNormalStressOutflow<TDomain>::
prep_elem_loop_cr()
{
//	register subsetIndex at Geometry
	static TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();

// 	Only first order implementation
	if(!(TFVGeom<TElem, dim>::order == 1))
		UG_THROW("Only first order implementation, but other Finite Volume"
						" Geometry set.");

//	check if kinematic Viscosity has been set
	if(!m_imKinViscosity.data_given())
		UG_THROW("NavierStokesNoNormalStressOutflow::prep_elem_loop_cr:"
						" Kinematic Viscosity has not been set, but is required.\n");

//	check if Density has been set
	if(!m_imDensity.data_given())
		UG_THROW("NavierStokesNoNormalStressOutflow::prep_elem_loop_cr:"
						" Density has not been set, but is required.\n");

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
void NavierStokesNoNormalStressOutflow<TDomain>::
fsh_elem_loop_cr()
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
void NavierStokesNoNormalStressOutflow<TDomain>::
prep_elem_cr(TElem* elem, const LocalVector& u)
{
// 	Update Geometry for this element
	TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();
	try{
		geo.update(elem,
	               this->template element_corners<TElem>(elem),
	               &(this->subset_handler()));
	}
	UG_CATCH_THROW("NavierStokesNoNormalStressOutflow::prep_elem_cr:"
						" Cannot update Finite Volume Geometry.");

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

	m_imDensity.set_local_ips(&m_vLocIP[0], m_vLocIP.size());
	m_imDensity.set_global_ips(&m_vGloIP[0], m_vGloIP.size());
}

/// Assembling of the diffusive flux (due to the viscosity) in the Jacobian of the momentum eq.
template<typename TDomain>
template<typename BF>
void NavierStokesNoNormalStressOutflow<TDomain>::
diffusive_flux_Jac_cr
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
void NavierStokesNoNormalStressOutflow<TDomain>::
diffusive_flux_defect_cr
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
void NavierStokesNoNormalStressOutflow<TDomain>::
convective_flux_Jac_cr
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
	old_momentum_flux = VecDot (StdVel, bf.normal ()) * m_imDensity [ip];

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
void NavierStokesNoNormalStressOutflow<TDomain>::
convective_flux_defect_cr
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
	old_momentum_flux = VecDot (StdVel, bf.normal ()) * m_imDensity [ip];

// We assume that there should be no inflow through the outflow boundary:
	if (old_momentum_flux < 0)
		old_momentum_flux = 0;

// Add the flux to the defect:
	for(size_t d1 = 0; d1 < (size_t) dim; ++d1)
		d(d1, bf.node_id()) += old_momentum_flux * StdVel[d1];
}

template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void NavierStokesNoNormalStressOutflow<TDomain>::
add_JA_elem_cr(LocalMatrix& J, const LocalVector& u)
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
			diffusive_flux_Jac_cr<BF> (ip, *bf, J, u);
			if (!m_spMaster->get_stokes ())
				convective_flux_Jac_cr<BF> (ip, *bf, J, u);

		//	B. The continuity equation
		//	for(size_t sh = 0; sh < bf->num_sh(); ++sh) // loop shape functions
		//		for (size_t d2 = 0; d2 < (size_t)dim; ++d2)
		//			J(_P_, bf->node_id (), d2, sh) += bf->shape(sh) * bf->normal()[d2]
		//				* m_imDensity [ip];

		// Next IP:
			ip++;
		}
	}
}

template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void NavierStokesNoNormalStressOutflow<TDomain>::
add_dA_elem_cr(LocalVector& d, const LocalVector& u)
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
			diffusive_flux_defect_cr<BF> (ip, *bf, d, u);
			if (!m_spMaster->get_stokes ())
				convective_flux_defect_cr<BF> (ip, *bf, d, u);

		// Next IP:
			ip++;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

// register for 1D
template<typename TDomain>
void
NavierStokesNoNormalStressOutflow<TDomain>::
register_all_cr_funcs(bool bHang)
{
//	get all grid element types in this dimension and below
	typedef typename domain_traits<dim>::DimElemList ElemList;
	// REMARK: Note that we register this boundary condition only
	// for the full-dimensional elements (DimElemList instead of AllElemList).

//	switch assemble functions
	if(!bHang) boost::mpl::for_each<ElemList>( RegisterCR<CRFVGeometry>(this) );
	else UG_THROW("Not implemented.");
}

template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void
NavierStokesNoNormalStressOutflow<TDomain>::
register_cr_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->enable_fast_add_elem(true);
	this->set_prep_elem_loop_fct(	id, &T::template prep_elem_loop_cr<TElem, TFVGeom>);
	this->set_prep_elem_fct(	 	id, &T::template prep_elem_cr<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct( 	id, &T::template fsh_elem_loop_cr<TElem, TFVGeom>);
	this->set_add_jac_A_elem_fct(		id, &T::template add_JA_elem_cr<TElem, TFVGeom>);
	this->set_add_jac_M_elem_fct(		id, &T::template add_JM_elem<TElem, TFVGeom>);
	this->set_add_def_A_elem_fct(		id, &T::template add_dA_elem_cr<TElem, TFVGeom>);
	this->set_add_def_M_elem_fct(		id, &T::template add_dM_elem<TElem, TFVGeom>);
	this->set_add_rhs_elem_fct(		id, &T::template add_rhs_elem<TElem, TFVGeom>);
}

} // namespace NavierStokes
} // end namespace ug

#endif /* __H__UG__NAVIER_STOKES__BND__NO_NORMAL_STRESS_OUTFLOW_CR__ */
