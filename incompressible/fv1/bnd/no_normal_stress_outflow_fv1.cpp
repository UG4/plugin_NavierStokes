/*
 * Copyright (c) 2012-2014:  G-CSC, Goethe University Frankfurt
 * Authors: Dmitry Logashenko, Andreas Vogel
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


#include "no_normal_stress_outflow_fv1.h"

#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"

namespace ug{
namespace NavierStokes{

////////////////////////////////////////////////////////////////////////////////
//	Constructor - set default values
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
NavierStokesNoNormalStressOutflowFV1<TDomain>::
NavierStokesNoNormalStressOutflowFV1(SmartPtr< IncompressibleNavierStokesBase<TDomain> > spMaster)
: NavierStokesNoNormalStressOutflowBase<TDomain>(spMaster)
{
//	register imports
	this->register_import(m_imKinViscosity);
	this->register_import(m_imDensity);
	this->register_import(m_imBinghamViscosity);
	this->register_import(m_imYieldStress);

//	initialize the imports from the master discretization
	set_kinematic_viscosity(spMaster->kinematic_viscosity ());
	set_density(spMaster->density ());
	set_bingham_viscosity(spMaster->bingham_viscosity ());
	set_yield_stress(spMaster->yield_stress ());

	//	update assemble functions
	register_all_funcs(false);
};


template<typename TDomain>
void NavierStokesNoNormalStressOutflowFV1<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	if(bNonRegularGrid)
		UG_THROW("NavierStokes: only regular grid implemented.");

//	check number
	if(vLfeID.size() != dim+1)
		UG_THROW("NavierStokes: Need exactly "<<dim+1<<" functions");

	for(int d = 0; d <= dim; ++d)
		if(vLfeID[d].type() != LFEID::LAGRANGE || vLfeID[d].order() != 1)
			UG_THROW("NavierStokes: 'fv1' expects Lagrange P1 trial space "
					"for velocity and pressure.");

	//	update assemble functions
	register_all_funcs(false);
}

////////////////////////////////////////////////////////////////////////////////
//	assembling functions
////////////////////////////////////////////////////////////////////////////////

/**
 * Prepares the element loop for a given element type: computes the FV-geo, ...
 * Note that there are separate loops for every type of the grid elements.
 */
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesNoNormalStressOutflowFV1<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
//	register subsetIndex at Geometry
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();

// 	Only first order implementation
	if(!(TFVGeom::order == 1))
		UG_THROW("Only first order implementation, but other Finite Volume"
						" Geometry set.");

//	check if kinematic Viscosity has been set
	if(!m_imKinViscosity.data_given())
		UG_THROW("NavierStokesNoNormalStressOutflowFV1::prep_elem_loop:"
						" Kinematic Viscosity has not been set, but is required.\n");

//	check if Density has been set
	if(!m_imDensity.data_given())
		UG_THROW("NavierStokesNoNormalStressOutflowFV1::prep_elem_loop:"
						" Density has not been set, but is required.\n");

//	extract indices of boundary
	this->extract_scheduled_data();

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
template<typename TElem, typename TFVGeom>
void NavierStokesNoNormalStressOutflowFV1<TDomain>::
fsh_elem_loop()
{
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();

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
template<typename TElem, typename TFVGeom>
void NavierStokesNoNormalStressOutflowFV1<TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
{
// 	Update Geometry for this element
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();
	try{
		geo.update(elem, vCornerCoords, &(this->subset_handler()));
	}
	UG_CATCH_THROW("NavierStokesNoNormalStressOutflowFV1::prep_elem:"
						" Cannot update Finite Volume Geometry.");

//	find and set the local and the global positions of the IPs for imports
	typedef typename TFVGeom::BF BF;
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

	m_imBinghamViscosity.set_local_ips(&m_vLocIP[0], m_vLocIP.size());
	m_imBinghamViscosity.set_global_ips(&m_vGloIP[0], m_vGloIP.size());
	
	m_imYieldStress.set_local_ips(&m_vLocIP[0], m_vLocIP.size());
	m_imYieldStress.set_global_ips(&m_vGloIP[0], m_vGloIP.size());

	if(m_spMaster->bingham())
	{
		// get old solution
		const LocalVector *pOldSol = NULL;
		const LocalVectorTimeSeries* vLocSol = this->local_time_solutions();
		pOldSol = &vLocSol->solution(1);

		// get const point of viscosity at integration points
		const number* pVisco = m_imKinViscosity.values();

		// cast constness away
		number* vVisco = const_cast<number*>(pVisco);

	// 	get finite volume geometry
		static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
		typedef typename TFVGeom::BF BF;

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
			for(bf = vBF.begin(); bf != vBF.end(); ++bf, ++ip)
			{
				number innerSum = 0.0;
				for(int d1 = 0; d1 < dim; ++d1)
				{
					for(int d2 = 0; d2 < dim; ++d2)
					{
						for(size_t sh = 0; sh < bf->num_sh(); ++sh)
						{
							if(m_spMaster->laplace())
								innerSum += (bf->global_grad(sh)[d1] * (*pOldSol)(d2, sh) * bf->global_grad(sh)[d2] * (*pOldSol)(d1, sh));
							else
							{
								innerSum += ((bf->global_grad(sh)[d1] * (*pOldSol)(d2, sh) + bf->global_grad(sh)[d2] * (*pOldSol)(d1, sh))
									*(bf->global_grad(sh)[d1] * (*pOldSol)(d2, sh) + bf->global_grad(sh)[d2] * (*pOldSol)(d1, sh)));
							}
						}
					}
				}
				number secondInvariant = 1.0/(pow(2, dim))*innerSum;
		
				vVisco[ip] = (m_imBinghamViscosity[ip] + m_imYieldStress[ip]/sqrt(0.1+secondInvariant))/m_imDensity[ip];
			}
		}
	}
}

/// Assembling of the diffusive flux (due to the viscosity) in the Jacobian of the momentum eq.
template<typename TDomain>
template<typename BF>
void NavierStokesNoNormalStressOutflowFV1<TDomain>::
diffusive_flux_Jac
(
	const size_t ip, // index of the integration point (for the viscosity)
	const BF& bf, // boundary face to assemble
	LocalMatrix& J, // local Jacobian to update
	const LocalVector& u // local solution
)
{
	MathMatrix<dim,dim> diffFlux, tang_diffFlux;
	MathVector<dim> normalStress;

	//UG_LOG("ip: " << ip << "\n");

	for(size_t sh = 0; sh < bf.num_sh(); ++sh) // loop shape functions
	{
		//UG_LOG("sh: " << sh << "\n");
	//	1. Compute the total flux
	//	- add \nabla u
		MatSet (diffFlux, 0);
		MatDiagSet (diffFlux, VecDot (bf.global_grad(sh), bf.normal()));
	
	//	- add (\nabla u)^T
		if(!m_spMaster->laplace())
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
		tang_diffFlux *= - m_imKinViscosity[ip] * m_imDensity [ip];
	
	//	4. Add flux to local Jacobian
		for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
			for(size_t d2 = 0; d2 < (size_t)dim; ++d2)
				J(d1, bf.node_id(), d2, sh) += tang_diffFlux (d1, d2);
	}
}

/// Assembling of the diffusive flux (due to the viscosity) in the defect of the momentum eq.
template<typename TDomain>
template<typename BF>
void NavierStokesNoNormalStressOutflowFV1<TDomain>::
diffusive_flux_defect
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
	if(!m_spMaster->laplace())
		TransposedMatVecMultAdd(diffFlux, gradVel, bf.normal());

//	3. Subtract the normal part:
	VecScaleAppend (diffFlux, - VecDot (diffFlux, bf.normal()), bf.normal());

//	A4. Scale by viscosity
	diffFlux *= - m_imKinViscosity[ip] * m_imDensity [ip];

//	5. Add flux to local defect
	for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
		d(d1, bf.node_id()) += diffFlux[d1];
}

/// Assembling of the convective flux (due to the quadratic inertial term) in the Jacobian of the momentum eq.
template<typename TDomain>
template<typename BF>
void NavierStokesNoNormalStressOutflowFV1<TDomain>::
convective_flux_Jac
(
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
			StdVel[d1] += u(d1, sh) * bf.shape(sh);
	number old_momentum_flux = VecDot (StdVel, bf.normal ()) * m_imDensity [ip];
	
// We assume that there should be no inflow through the outflow boundary:
	if (old_momentum_flux < 0)
		old_momentum_flux = 0;
	
//	Add flux to local Jacobian
	for(size_t sh = 0; sh < bf.num_sh(); ++sh)
	{
		number t = old_momentum_flux * bf.shape(sh);
		for(size_t d1 = 0; d1 < (size_t) dim; ++d1)
			J(d1, bf.node_id(), d1, sh) += t;
	}
}

/// Assembling of the convective flux (due to the quadratic inertial term) in the defect of the momentum eq.
template<typename TDomain>
template<typename BF>
void NavierStokesNoNormalStressOutflowFV1<TDomain>::
convective_flux_defect
(
	const size_t ip, // index of the integration point (for the density)
	const BF& bf, // boundary face to assemble
	LocalVector& d, // local defect to update
	const LocalVector& u, // local solution
	const MathVector<dim>& StdVel // velocity at ip
)
{
// The convection velocity according to the current approximation:
	number old_momentum_flux = VecDot (StdVel, bf.normal ()) * m_imDensity [ip];
	
// We assume that there should be no inflow through the outflow boundary:
	if (old_momentum_flux < 0)
		old_momentum_flux = 0;
	
// Add the flux to the defect:
	for(size_t d1 = 0; d1 < (size_t) dim; ++d1)
		d(d1, bf.node_id()) += old_momentum_flux * StdVel[d1];
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesNoNormalStressOutflowFV1<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	Only first order implementation
	UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
	typedef typename TFVGeom::BF BF;

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
		for(bf = vBF.begin(); bf != vBF.end(); ++bf, ++ip)
		{
		//	A. The momentum equation:
			diffusive_flux_Jac<BF> (ip, *bf, J, u);
			if (!m_spMaster->stokes ())
				convective_flux_Jac<BF> (ip, *bf, J, u);
			
		//	B. The continuity equation
			for(size_t sh = 0; sh < bf->num_sh(); ++sh) // loop shape functions
				for (size_t d2 = 0; d2 < (size_t)dim; ++d2)
					J(_P_, bf->node_id (), d2, sh) += bf->shape(sh) * bf->normal()[d2]
						* m_imDensity [ip];
		}
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesNoNormalStressOutflowFV1<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	Only first order implemented
	UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
	typedef typename TFVGeom::BF BF;

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
		for(bf = vBF.begin(); bf != vBF.end(); ++bf, ++ip)
		{
		// A. Compute Velocity at ip
			MathVector<dim> stdVel(0.0);
			for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
				for(size_t sh = 0; sh < bf->num_sh(); ++sh)
					stdVel[d1] += u(d1, sh) * bf->shape(sh);

		// B. Momentum equation:
			diffusive_flux_defect<BF> (ip, *bf, d, u);
			if (!m_spMaster->stokes ())
				convective_flux_defect<BF> (ip, *bf, d, u, stdVel);
		
		// c. Continuity equation:
			d(_P_, bf->node_id()) += VecDot (stdVel, bf->normal()) * m_imDensity[ip];
		}
	}
}


////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template<>
void NavierStokesNoNormalStressOutflowFV1<Domain1d>::
register_all_funcs(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
	}
	else
	{
		UG_THROW("NavierStokesNoNormalStressOutflowFV1: Hanging Nodes not implemented.")
	}
}
#endif

#ifdef UG_DIM_2
template<>
void NavierStokesNoNormalStressOutflowFV1<Domain2d>::
register_all_funcs(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_func<Triangle, FV1Geometry<Triangle, dim> >();
		register_func<Quadrilateral, FV1Geometry<Quadrilateral, dim> >();
	}
	else
	{
		UG_THROW("NavierStokesNoNormalStressOutflowFV1: Hanging Nodes not implemented.")
	}
}
#endif

#ifdef UG_DIM_3
template<>
void NavierStokesNoNormalStressOutflowFV1<Domain3d>::
register_all_funcs(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_func<Tetrahedron, FV1Geometry<Tetrahedron, dim> >();
		register_func<Prism, FV1Geometry<Prism, dim> >();
		register_func<Pyramid, FV1Geometry<Pyramid, dim> >();
		register_func<Hexahedron, FV1Geometry<Hexahedron, dim> >();
	}
	else
	{
		UG_THROW("NavierStokesNoNormalStressOutflowFV1: Hanging Nodes not implemented.")
	}
}
#endif

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void
NavierStokesNoNormalStressOutflowFV1<TDomain>::
register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->clear_add_fct(id);
	this->set_prep_elem_loop_fct(	id, &T::template prep_elem_loop<TElem, TFVGeom>);
	this->set_prep_elem_fct(	 	id, &T::template prep_elem<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct( 	id, &T::template fsh_elem_loop<TElem, TFVGeom>);
	this->set_add_jac_A_elem_fct(	id, &T::template add_jac_A_elem<TElem, TFVGeom>);
	this->set_add_jac_M_elem_fct(	id, &T::template add_jac_M_elem<TElem, TFVGeom>);
	this->set_add_def_A_elem_fct(	id, &T::template add_def_A_elem<TElem, TFVGeom>);
	this->set_add_def_M_elem_fct(	id, &T::template add_def_M_elem<TElem, TFVGeom>);
	this->set_add_rhs_elem_fct(	id, &T::template add_rhs_elem<TElem, TFVGeom>);
}


////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class NavierStokesNoNormalStressOutflowFV1<Domain1d>;
#endif
#ifdef UG_DIM_2
template class NavierStokesNoNormalStressOutflowFV1<Domain2d>;
#endif
#ifdef UG_DIM_3
template class NavierStokesNoNormalStressOutflowFV1<Domain3d>;
#endif

} // namespace NavierStokes
} // namespace ug
