/*
 * no_normal_stress_outflow_fv.cpp
 *
 *  Created on: 27.03.2012
 *  D. Logashenko, A. Vogel
 */

#include "no_normal_stress_outflow_fv.h"

#include "lib_disc/spatial_disc/disc_util/fvho_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"

namespace ug{
namespace NavierStokes{

////////////////////////////////////////////////////////////////////////////////
//	Constructor - set default values
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
NavierStokesNoNormalStressOutflowFV<TDomain>::
NavierStokesNoNormalStressOutflowFV(SmartPtr< IncompressibleNavierStokesBase<TDomain> > spMaster)
: NavierStokesNoNormalStressOutflowBase<TDomain>(spMaster)
{
//	register imports
	this->register_import(m_imKinViscosity);
	this->register_import(m_imDensity);
	this->register_import(m_imDensityP);

//	initialize the imports from the master discretization
	set_kinematic_viscosity(spMaster->kinematic_viscosity ());
	set_density(spMaster->density ());

	//	update assemble functions
	this->enable_fast_add_elem(true);
};


template<typename TDomain>
void NavierStokesNoNormalStressOutflowFV<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	if(bNonRegularGrid)
		UG_THROW("NavierStokesNoNormalStressOutflow: only regular grid implemented.");

	//	check number
	if(vLfeID.size() != dim+1)
		UG_THROW("NavierStokesNoNormalStressOutflow: Needs exactly "<<dim+1<<" functions");

	for(int d = 1; d < dim; ++d)
		if(vLfeID[0] != vLfeID[d])
			UG_THROW("NavierStokesNoNormalStressOutflow: trial spaces for velocity"
					" expected to be identical for all velocity components.");

//	remember lfeID;
	m_vLFEID = vLfeID[0];
	m_pLFEID = vLfeID[dim];
	m_quadOrder = std::max(m_vLFEID.order(), m_pLFEID.order()) + 1;

	//	update assemble functions
	register_all_funcs(m_vLFEID, m_pLFEID);
}

////////////////////////////////////////////////////////////////////////////////
//	assembling functions
////////////////////////////////////////////////////////////////////////////////

/**
 * Prepares the element loop for a given element type: computes the FV-geo, ...
 * Note that there are separate loops for every type of the grid elements.
 */
template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesNoNormalStressOutflowFV<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
//	register subsetIndex at Geometry
	VGeom& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);
	PGeom& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_quadOrder);

	try{
		vgeo.update_local(roid, m_vLFEID);
		pgeo.update_local(roid, m_pLFEID);
	}
	UG_CATCH_THROW("NavierStokesNoNormalStressOutflowFV: "
					"Cannot update Finite Volume Geometry.");

//	check if kinematic Viscosity has been set
	if(!m_imKinViscosity.data_given())
		UG_THROW("NavierStokesNoNormalStressOutflowFV::prep_elem_loop:"
						" Kinematic Viscosity has not been set, but is required.\n");

//	check if Density has been set
	if(!m_imDensity.data_given())
		UG_THROW("NavierStokesNoNormalStressOutflowFV::prep_elem_loop:"
						" Density has not been set, but is required.\n");

//	extract indices of boundary
	this->extract_scheduled_data();

//	request the subset indices as boundary subset. This will force the
//	creation of boundary subsets when calling geo.update
	typename std::vector<int>::const_iterator subsetIter;
	for(subsetIter = m_vBndSubSetIndex.begin();
			subsetIter != m_vBndSubSetIndex.end(); ++subsetIter){
		vgeo.add_boundary_subset(*subsetIter);
		pgeo.add_boundary_subset(*subsetIter);
	}
}

/**
 * Finalizes the element loop for a given element type.
 */
template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesNoNormalStressOutflowFV<TDomain>::
fsh_elem_loop()
{
	VGeom& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);
	PGeom& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_quadOrder);

//	remove the bnd subsets
	typename std::vector<int>::const_iterator subsetIter;
	for(subsetIter = m_vBndSubSetIndex.begin();
			subsetIter != m_vBndSubSetIndex.end(); ++subsetIter){
		vgeo.remove_boundary_subset(*subsetIter);
		pgeo.remove_boundary_subset(*subsetIter);
	}
}


/**
 * General initializations of a given grid element for the assembling.
 */
template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesNoNormalStressOutflowFV<TDomain>::
prep_elem(const LocalVector& u, GeometricObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	Update Geometry for this element
	VGeom& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);
	PGeom& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_quadOrder);
	try{
		vgeo.update(elem, vCornerCoords, &(this->subset_handler()));
	    pgeo.update(elem, vCornerCoords, &(this->subset_handler()));
	}
	UG_CATCH_THROW("NavierStokesNoNormalStressOutflowFV::prep_elem:"
						" Cannot update Finite Volume Geometry.");

//	find and set the local and the global positions of the IPs for imports
	typename std::vector<int>::const_iterator subsetIter;
	
	m_vLocIPv.clear(); m_vGloIPv.clear();
	for(subsetIter = m_vBndSubSetIndex.begin();
			subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
	{
		const int bndSubset = *subsetIter;
		if(vgeo.num_bf(bndSubset) == 0) continue;
		typedef typename VGeom::BF BF;
		const std::vector<BF>& vBF = vgeo.bf(bndSubset);
		for(size_t i = 0; i < vBF.size(); ++i)
			for(size_t ip = 0; ip < vBF[i].num_ip(); ++ip){
				m_vLocIPv.push_back(vBF[i].local_ip(ip));
				m_vGloIPv.push_back(vBF[i].global_ip(ip));
			}
	}

	m_vLocIPp.clear(); m_vGloIPp.clear();
	for(subsetIter = m_vBndSubSetIndex.begin();
			subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
	{
		const int bndSubset = *subsetIter;
		if(vgeo.num_bf(bndSubset) == 0) continue;
		typedef typename PGeom::BF BF;
		const std::vector<BF>& vBF = pgeo.bf(bndSubset);
		for(size_t i = 0; i < vBF.size(); ++i)
			for(size_t ip = 0; ip < vBF[i].num_ip(); ++ip){
				m_vLocIPp.push_back(vBF[i].local_ip(ip));
				m_vGloIPp.push_back(vBF[i].global_ip(ip));
			}
	}

	const LocalShapeFunctionSet<dim>& rVTrialSpace =
		LocalFiniteElementProvider::get<dim>(elem->reference_object_id(),m_vLFEID);
	m_vvVShape.resize(m_vLocIPp.size());
	for(size_t ip = 0; ip < m_vvVShape.size(); ++ip){
		m_vvVShape[ip].resize(rVTrialSpace.num_sh());
		for(size_t sh = 0; sh < rVTrialSpace.num_sh(); ++sh){
			m_vvVShape[ip][sh] = rVTrialSpace.shape(sh, m_vLocIPp[ip]);
		}
	}

	// REMARK: The loop above determines the ordering of the integration points:
	// The "outer ordering" corresponds to the ordering of the subsets in
	// m_vBndSubSetIndex, and "inside" of this ordering, the ip's are ordered
	// according to the order of the boundary faces in the FV geometry structure.

	m_imKinViscosity.set_local_ips(&m_vLocIPv[0], m_vLocIPv.size());
	m_imKinViscosity.set_global_ips(&m_vGloIPv[0], m_vGloIPv.size());
	
	m_imDensity.set_local_ips(&m_vLocIPv[0], m_vLocIPv.size());
	m_imDensity.set_global_ips(&m_vGloIPv[0], m_vGloIPv.size());

	m_imDensityP.set_local_ips(&m_vLocIPp[0], m_vLocIPp.size());
	m_imDensityP.set_global_ips(&m_vGloIPp[0], m_vGloIPp.size());
}

/// Assembling of the diffusive flux (due to the viscosity) in the Jacobian of the momentum eq.
template<typename TDomain>
template<typename BF>
void NavierStokesNoNormalStressOutflowFV<TDomain>::
diffusive_flux_Jac
(
	const size_t i, // index of the integration point with in boundary face
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
		MatDiagSet (diffFlux, VecDot (bf.global_grad(i, sh), bf.normal()));
	
	//	- add (\nabla u)^T
		if(!m_spMaster->laplace())
			for (size_t d1 = 0; d1 < (size_t)dim; ++d1)
				for (size_t d2 = 0; d2 < (size_t)dim; ++d2)
					diffFlux(d1,d2) += bf.global_grad(i, sh)[d1] * bf.normal()[d2];
	
	//	2. Subtract the normal part:
		tang_diffFlux = diffFlux;
		TransposedMatVecMult(normalStress, diffFlux, bf.normal ());
		for (size_t d2 = 0; d2 < (size_t)dim; ++d2)
			for (size_t d1 = 0; d1 < (size_t)dim; ++d1)
				tang_diffFlux(d1,d2) -= bf.normal()[d1] * normalStress[d2];
	
	//	3. Scale by viscosity
		tang_diffFlux *= - m_imKinViscosity[ip] * bf.weight(i);
	
	//	4. Add flux to local Jacobian
		for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
			for(size_t d2 = 0; d2 < (size_t)dim; ++d2)
				J(d1, bf.node_id(), d2, sh) += tang_diffFlux (d1, d2);
	}
}

/// Assembling of the diffusive flux (due to the viscosity) in the defect of the momentum eq.
template<typename TDomain>
template<typename BF>
void NavierStokesNoNormalStressOutflowFV<TDomain>::
diffusive_flux_defect
(
	const size_t i, // index of the integration point with in boundary face
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
				gradVel(d1, d2) += bf.global_grad(i, sh)[d2] * u(d1, sh);
		}

//	2. Compute the total flux

//	- add (\nabla u) \cdot \vec{n}
	MatVecMult(diffFlux, gradVel, bf.normal());

//	- add (\nabla u)^T \cdot \vec{n}
	if(!m_spMaster->laplace())
		TransposedMatVecMultAdd(diffFlux, gradVel, bf.normal());

//	3. Subtract the normal part:
	VecScaleAppend (diffFlux, - VecDot (diffFlux, bf.normal()), bf.normal());

//	A4. Scale by viscosity
	diffFlux *= - m_imKinViscosity[ip] * bf.weight(i);

//	5. Add flux to local defect
	for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
		d(d1, bf.node_id()) += diffFlux[d1];
}

/// Assembling of the convective flux (due to the quadratic inertial term) in the Jacobian of the momentum eq.
template<typename TDomain>
template<typename BF>
void NavierStokesNoNormalStressOutflowFV<TDomain>::
convective_flux_Jac
(
	const size_t i, // index of the integration point with in boundary face
	const size_t ip, // index of the integration point (for the density)
	const BF& bf, // boundary face to assemble
	LocalMatrix& J, // local Jacobian to update
	const LocalVector& u // local solution
)
{
// The convection velocity according to the current approximation:
	MathVector<dim> StdVel(0.0);
	for(size_t sh = 0; sh < bf.num_sh(); ++sh)
		for(size_t d1 = 0; d1 < (size_t) dim; ++d1)
			StdVel[d1] += u(d1, sh) * bf.shape(i, sh);
	
	number old_momentum_flux = VecDot (StdVel, bf.normal ())
								* m_imDensity [ip] * bf.weight(i);

// We assume that there should be no inflow through the outflow boundary:
	if (old_momentum_flux < 0)
		old_momentum_flux = 0;
	
//	Add flux to local Jacobian
	for(size_t sh = 0; sh < bf.num_sh(); ++sh)
	{
		const number t = old_momentum_flux * bf.shape(i, sh);
		for(size_t d1 = 0; d1 < (size_t) dim; ++d1)
			J(d1, bf.node_id(), d1, sh) += t;
	}
}

/// Assembling of the convective flux (due to the quadratic inertial term) in the defect of the momentum eq.
template<typename TDomain>
template<typename BF>
void NavierStokesNoNormalStressOutflowFV<TDomain>::
convective_flux_defect
(
	const size_t i, // index of the integration point with in boundary face
	const size_t ip, // index of the integration point (for the density)
	const BF& bf, // boundary face to assemble
	LocalVector& d, // local defect to update
	const LocalVector& u // local solution
)
{
// Compute Velocity at ip
	MathVector<dim> StdVel(0.0);
	for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
		for(size_t sh = 0; sh < bf.num_sh(); ++sh)
			StdVel[d1] += u(d1, sh) * bf.shape(i, sh);

	number old_momentum_flux = VecDot (StdVel, bf.normal ())
								* m_imDensity [ip] * bf.weight(i);
	
// We assume that there should be no inflow through the outflow boundary:
	if (old_momentum_flux < 0)
		old_momentum_flux = 0;
	
// Add the flux to the defect:
	for(size_t d1 = 0; d1 < (size_t) dim; ++d1)
		d(d1, bf.node_id()) += old_momentum_flux * StdVel[d1];
}

template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesNoNormalStressOutflowFV<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GeometricObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	get finite volume geometry
	VGeom& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);
	PGeom& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_quadOrder);

// 	loop registered boundary segments
	typename std::vector<int>::const_iterator subsetIter;
	size_t ipV = 0, ipP = 0;
	for(subsetIter = m_vBndSubSetIndex.begin();
		subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
	{
	//	get subset index corresponding to boundary
		const int bndSubset = *subsetIter;
		
	//	A. The momentum equation:
		{
			typedef typename VGeom::BF BF;
			const std::vector<BF>& vBF = vgeo.bf(bndSubset);
			typename std::vector<BF>::const_iterator bf;
			for(bf = vBF.begin(); bf != vBF.end(); ++bf){
				for(size_t i = 0; i < bf->num_ip(); ++i, ++ipV){
					diffusive_flux_Jac<BF> (i, ipV, *bf, J, u);
					if (!m_spMaster->stokes ())
						convective_flux_Jac<BF> (i, ipV, *bf, J, u);
				}
			}
		}

		//	B. The continuity equation
		{
			typedef typename PGeom::BF BF;
			const std::vector<BF>& vBF = pgeo.bf(bndSubset);
			typename std::vector<BF>::const_iterator bf;
			for(bf = vBF.begin(); bf != vBF.end(); ++bf){
				for(size_t i = 0; i < bf->num_ip(); ++i, ++ipP){
					for(size_t sh = 0; sh < bf->num_sh(); ++sh) // loop shape functions
						for (size_t d2 = 0; d2 < (size_t)dim; ++d2)
							J(_P_, bf->node_id (), d2, sh) += m_vvVShape[ipP][sh] * bf->normal()[d2]
															   * m_imDensityP [ipP] * bf->weight(i);
				}
			}
		}
	}
}

template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesNoNormalStressOutflowFV<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u, GeometricObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	get finite volume geometry
	VGeom& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);
	PGeom& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_quadOrder);

// 	loop registered boundary segments
	typename std::vector<int>::const_iterator subsetIter;
	size_t ipV = 0, ipP = 0;
	for(subsetIter = m_vBndSubSetIndex.begin();
		subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
	{
	//	get subset index corresponding to boundary
		const int bndSubset = *subsetIter;
		
	//	A. The momentum equation:
		{
			typedef typename VGeom::BF BF;
			const std::vector<BF>& vBF = vgeo.bf(bndSubset);
			typename std::vector<BF>::const_iterator bf;
			for(bf = vBF.begin(); bf != vBF.end(); ++bf){
				for(size_t i = 0; i < bf->num_ip(); ++i, ++ipV){
					diffusive_flux_defect<BF> (i, ipV, *bf, d, u);
					if (!m_spMaster->stokes ())
						convective_flux_defect<BF> (i, ipV, *bf, d, u);
				}
			}
		}

	//	B. The continuity equation
		{
			typedef typename PGeom::BF BF;
			const std::vector<BF>& vBF = pgeo.bf(bndSubset);
			typename std::vector<BF>::const_iterator bf;
			for(bf = vBF.begin(); bf != vBF.end(); ++bf){
				for(size_t i = 0; i < bf->num_ip(); ++i, ++ipP){

				// Compute Velocity at ip
					MathVector<dim> stdVel(0.0);
					for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
						for(size_t sh = 0; sh < bf->num_sh(); ++sh)
							stdVel[d1] += u(d1, sh) * m_vvVShape[ipP][sh];

					d(_P_, bf->node_id()) += VecDot (stdVel, bf->normal())
											 * m_imDensityP[ipP] * bf->weight(i);
				}
			}
		}
	}
}


////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template<>
void NavierStokesNoNormalStressOutflowFV<Domain1d>::
register_all_funcs(const LFEID& vLfeID, const LFEID& pLfeID)
{
	UG_THROW("Not implemented.");
}
#endif

#ifdef UG_DIM_2
template<>
void NavierStokesNoNormalStressOutflowFV<Domain2d>::
register_all_funcs(const LFEID& vLfeID, const LFEID& pLfeID)
{
	typedef DimFVGeometry<dim> FVGeom;
	register_func<Triangle, FVGeom, FVGeom >();
	register_func<Quadrilateral, FVGeom, FVGeom >();
}
#endif

#ifdef UG_DIM_3
template<>
void NavierStokesNoNormalStressOutflowFV<Domain3d>::
register_all_funcs(const LFEID& vLfeID, const LFEID& pLfeID)
{
	typedef DimFVGeometry<dim> FVGeom;
	register_func<Tetrahedron, FVGeom, FVGeom >();
	register_func<Prism, FVGeom, FVGeom >();
	register_func<Hexahedron, FVGeom, FVGeom >();
}
#endif

template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void
NavierStokesNoNormalStressOutflowFV<TDomain>::
register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->enable_fast_add_elem(true);
	this->set_prep_elem_loop_fct(	id, &T::template prep_elem_loop<TElem, VGeom, PGeom>);
	this->set_prep_elem_fct(	 	id, &T::template prep_elem<TElem, VGeom, PGeom>);
	this->set_fsh_elem_loop_fct( 	id, &T::template fsh_elem_loop<TElem, VGeom, PGeom>);
	this->set_add_jac_A_elem_fct(	id, &T::template add_jac_A_elem<TElem, VGeom, PGeom>);
	this->set_add_jac_M_elem_fct(	id, &T::template add_jac_M_elem<TElem, VGeom, PGeom>);
	this->set_add_def_A_elem_fct(	id, &T::template add_def_A_elem<TElem, VGeom, PGeom>);
	this->set_add_def_M_elem_fct(	id, &T::template add_def_M_elem<TElem, VGeom, PGeom>);
	this->set_add_rhs_elem_fct(	id, &T::template add_rhs_elem<TElem, VGeom, PGeom>);
}


////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class NavierStokesNoNormalStressOutflowFV<Domain1d>;
#endif
#ifdef UG_DIM_2
template class NavierStokesNoNormalStressOutflowFV<Domain2d>;
#endif
#ifdef UG_DIM_3
template class NavierStokesNoNormalStressOutflowFV<Domain3d>;
#endif

} // namespace NavierStokes
} // namespace ug
