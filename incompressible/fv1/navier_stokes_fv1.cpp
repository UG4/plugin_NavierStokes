/*
 * Copyright (c) 2010-2014:  G-CSC, Goethe University Frankfurt
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

#include "navier_stokes_fv1.h"

#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"

namespace ug{
namespace NavierStokes{

////////////////////////////////////////////////////////////////////////////////
//	Constructor - set default values
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
NavierStokesFV1<TDomain>::NavierStokesFV1(const char* functions,
                                          const char* subsets)
: IncompressibleNavierStokesBase<TDomain>(functions, subsets)
{
	init();
};

template<typename TDomain>
NavierStokesFV1<TDomain>::NavierStokesFV1(const std::vector<std::string>& vFct,
                                          const std::vector<std::string>& vSubset)
: IncompressibleNavierStokesBase<TDomain>(vFct, vSubset)
{
	init();
};


template<typename TDomain>
void NavierStokesFV1<TDomain>::init()
{
//	check number of functions
	if(this->num_fct() != dim+1)
		UG_THROW("Wrong number of functions: The ElemDisc 'NavierStokes'"
					   " needs exactly "<<dim+1<<" symbolic function.");

//	register imports
	this->register_import(m_imSourceSCV);
	this->register_import(m_imSourceSCVF);
	this->register_import(m_imKinViscosity);
	this->register_import(m_imDensitySCVF);
	this->register_import(m_imDensitySCV);
	this->register_import(m_imBinghamViscosity);
	this->register_import(m_imYieldStress);

	m_imSourceSCV.set_rhs_part();
	m_imSourceSCVF.set_rhs_part();
	m_imDensitySCV.set_mass_part();

	//	default value for density
	base_type::set_density(1.0);

	//	update assemble functions
	register_all_funcs(false);
}

template<typename TDomain>
void NavierStokesFV1<TDomain>::
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

template<typename TDomain>
void NavierStokesFV1<TDomain>::
set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > data)
{
	m_imKinViscosity.set_data(data);
}

template<typename TDomain>
void NavierStokesFV1<TDomain>::
set_density(SmartPtr<CplUserData<number, dim> > data)
{
	m_imDensitySCVF.set_data(data);
	m_imDensitySCV.set_data(data);
}

template<typename TDomain>
void NavierStokesFV1<TDomain>::
set_source(SmartPtr<CplUserData<MathVector<dim>, dim> > data)
{
	m_imSourceSCV.set_data(data);
	m_imSourceSCVF.set_data(data);
}


template<typename TDomain>
void NavierStokesFV1<TDomain>::
set_bingham_viscosity(SmartPtr<CplUserData<number, dim> > data)
{
	m_imBinghamViscosity.set_data(data);
}

template<typename TDomain>
void NavierStokesFV1<TDomain>::
set_yield_stress(SmartPtr<CplUserData<number, dim> > data)
{
	m_imYieldStress.set_data(data);
}

//template<typename TDomain>
//void NavierStokesFV1<TDomain>::
//set_regularize_delta(number data)
//{
//	m_imRegularizeDelta.set_data(make_sp(new ConstUserNumber<dim>(data)));
//}


////////////////////////////////////////////////////////////////////////////////
//	assembling functions
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesFV1<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
// 	Only first order implementation
	if(!(TFVGeom::order == 1))
		UG_THROW("Only first order implementation, but other Finite Volume"
						" Geometry set.");

//	check, that stabilization has been set
	if(m_spStab.invalid())
		UG_THROW("Stabilization has not been set.");

//	init stabilization for element type
	m_spStab->template set_geometry_type<TFVGeom >();

	if (! m_bStokes) // no convective terms in the Stokes eq. => no upwinding
	{
	//	check, that convective upwinding has been set
		if(m_spConvStab.invalid()  && m_spConvUpwind.invalid())
			UG_THROW("Upwinding for convective Term in Momentum eq. not set.");
	
	//	init convection stabilization for element type
		if(m_spConvStab.valid())
			m_spConvStab->template set_geometry_type<TFVGeom >();
	
	//	init convection stabilization for element type
		if(m_spConvUpwind.valid())
			m_spConvUpwind->template set_geometry_type<TFVGeom >();
	}

//	check, that kinematic Viscosity has been set
	if(!m_bBingham)
		if(!m_imKinViscosity.data_given())
			UG_THROW("NavierStokes::prep_elem_loop:"
							" Kinematic Viscosity has not been set, but is required.");

//	check, that Density has been set
	if(!m_imDensitySCVF.data_given())
		UG_THROW("NavierStokes::prep_elem_loop:"
						" Density has not been set, but is required.");

//	check, that Density has been set
	if(!m_imDensitySCV.data_given())
		UG_THROW("NavierStokes::prep_elem_loop:"
						" Density has not been set, but is required.");

	if(m_bBingham){
		if(!m_imBinghamViscosity.data_given())
			UG_THROW("NavierStokes::prep_elem_loop:"
							" Bingham Viscosity has not been set, but is required.");

		if(!m_imYieldStress.data_given())
			UG_THROW("NavierStokes::prep_elem_loop:"
							" Yield Stress has not been set, but is required.");

		//if(!m_imRegularizeDelta.data_given())
		//	UG_THROW("NavierStokes::prep_elem_loop:"
		//					" regularizing delta has not been set, but is required.");

		if(!m_bStokes)
			UG_THROW("NavierStokes::prep_elem_loop:"
							" Bingham behaviour is only available for Stokes.");

		if(!this->is_time_dependent()){
			UG_THROW("NavierStokes::prep_elem_loop:"
							" Bingham behaviour is only available for time-dependent problems.");
		}
	}

//	set local positions for imports

	if(!TFVGeom::usesHangingNodes)
	{
		static const int refDim = TElem::dim;
		TFVGeom& geo = GeomProvider<TFVGeom>::get();
		const MathVector<refDim>* vSCVFip = geo.scvf_local_ips();
		const size_t numSCVFip = geo.num_scvf_ips();
		const MathVector<refDim>* vSCVip = geo.scv_local_ips();
		const size_t numSCVip = geo.num_scv_ips();

		m_imKinViscosity.template set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imDensitySCVF.template set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imDensitySCV.template set_local_ips<refDim>(vSCVip,numSCVip);
		m_imSourceSCV.template set_local_ips<refDim>(vSCVip,numSCVip);
		m_imSourceSCVF.template set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imBinghamViscosity.template set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imYieldStress.template set_local_ips<refDim>(vSCVFip,numSCVFip);
		//m_imRegularizeDelta.template set_local_ips<refDim>(vSCVFip,numSCVFip);
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesFV1<TDomain>::
fsh_elem_loop()
{}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesFV1<TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
{
// 	Update Geometry for this element
	TFVGeom& geo = GeomProvider<TFVGeom>::get();
	try{
		geo.update(elem, vCornerCoords, &(this->subset_handler()));
	}
	UG_CATCH_THROW("NavierStokes::prep_elem:"
						" Cannot update Finite Volume Geometry.");

//	set local positions for imports
	if(TFVGeom::usesHangingNodes)
	{
	//	request ip series
		static const int refDim = TElem::dim;
		const MathVector<refDim>* vSCVFip = geo.scvf_local_ips();
		const size_t numSCVFip = geo.num_scvf_ips();
		const MathVector<refDim>* vSCVip = geo.scv_local_ips();
		const size_t numSCVip = geo.num_scv_ips();
		m_imKinViscosity.template set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imDensitySCVF.template set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imDensitySCV.template set_local_ips<refDim>(vSCVip,numSCVip);
		m_imSourceSCV.template set_local_ips<refDim>(vSCVip,numSCVip);
		m_imSourceSCVF.template set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imBinghamViscosity.template set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imYieldStress.template set_local_ips<refDim>(vSCVFip,numSCVFip);
		//m_imRegularizeDelta.template set_local_ips<refDim>(vSCVFip,numSCVFip);

	}

//	set global positions for imports
	const MathVector<dim>* vSCVFip = geo.scvf_global_ips();
	const size_t numSCVFip = geo.num_scvf_ips();
	const MathVector<dim>* vSCVip = geo.scv_global_ips();
	const size_t numSCVip = geo.num_scv_ips();
	m_imKinViscosity.set_global_ips(vSCVFip, numSCVFip);
	m_imDensitySCVF.set_global_ips(vSCVFip, numSCVFip);
	m_imDensitySCV.set_global_ips(vSCVip, numSCVip);
	m_imSourceSCV.set_global_ips(vSCVip, numSCVip);
	m_imSourceSCVF.set_global_ips(vSCVFip, numSCVFip);
	m_imBinghamViscosity.set_global_ips(vSCVFip, numSCVFip);
	m_imYieldStress.set_global_ips(vSCVFip, numSCVFip);
	//m_imRegularizeDelta.set_global_ips(vSCVFip, numSCVFip);

//	bingham behaviour
	if(m_bBingham)
	{
		//get old solution
		const LocalVector *pOldSol = NULL;
		const LocalVectorTimeSeries* vLocSol = this->local_time_solutions();
		pOldSol = &vLocSol->solution(1);
		
		// get const point of viscosity at integration points
		const number* pVisco = m_imKinViscosity.values();

		// cast constness away
		number* vVisco = const_cast<number*>(pVisco);

		number innerSum = 0.0;
		number secondInvariant = 0.0;

		for (size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
			// 	get current SCVF
			const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
			for(int d1 = 0; d1 < dim; ++d1)
			{
				for(int d2 = 0; d2 < dim; ++d2)
				{
					for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					{
						if(m_bLaplace)
						{
							innerSum += (scvf.global_grad(sh)[d1] * (*pOldSol)(d2, sh)* scvf.global_grad(sh)[d2] * (*pOldSol)(d1, sh));
						}
						else
						{
							innerSum += ((scvf.global_grad(sh)[d1] * (*pOldSol)(d2, sh) + scvf.global_grad(sh)[d2] * (*pOldSol)(d1, sh))
								*(scvf.global_grad(sh)[d1] * (*pOldSol)(d2, sh) + scvf.global_grad(sh)[d2] * (*pOldSol)(d1, sh)));
						}
					}
				}
			}
			// overwrite viscosity
			secondInvariant = 1.0/(pow(2, dim))*innerSum;
			vVisco[ip] = (m_imBinghamViscosity[ip] + m_imYieldStress[ip]/sqrt(0.1+secondInvariant))/m_imDensitySCVF[ip];
		}
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesFV1<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	Only first order implementation
	UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	check for source term to pass to the stabilization
	const DataImport<MathVector<dim>, dim>* pSource = NULL;
	if(m_imSourceSCVF.data_given())	pSource = &m_imSourceSCVF;

//	check for solutions to pass to stabilization in time-dependent case
	const LocalVector *pSol = &u, *pOldSol = NULL;
	number dt = 0.0;
	if(this->is_time_dependent())
	{
	//	get and check current and old solution
		const LocalVectorTimeSeries* vLocSol = this->local_time_solutions();
		if(vLocSol->size() != 2)
			UG_THROW("NavierStokes::add_jac_A_elem: "
							" Stabilization needs exactly two time points.");

	//	remember local solutions
		pSol = &vLocSol->solution(0);
		pOldSol = &vLocSol->solution(1);
		dt = vLocSol->time(0) - vLocSol->time(1);
	}

//	interpolate velocity at ip with standard lagrange interpolation
	static const size_t numSCVF = TFVGeom::numSCVF;
	MathVector<dim> StdVel[numSCVF];
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
		VecSet(StdVel[ip], 0.0);

		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			for(int d1 = 0; d1 < dim; ++d1)
				StdVel[ip][d1] += u(d1, sh) * scvf.shape(sh);
	}

//	compute stabilized velocities and shapes for continuity equation
	m_spStab->update(&geo, *pSol, StdVel, m_bStokes, m_imKinViscosity, m_imDensitySCVF, pSource, pOldSol, dt);

	if (! m_bStokes) // no convective terms in the Stokes eq. => no upwinding
	{
	//	compute stabilized velocities and shapes for convection upwind
		if(m_spConvStab.valid())
			if(m_spConvStab != m_spStab)
				m_spConvStab->update(&geo, *pSol, StdVel, false, m_imKinViscosity, m_imDensitySCVF, pSource, pOldSol, dt);
	
	//	compute upwind shapes
		if(m_spConvUpwind.valid())
			if(m_spStab->upwind() != m_spConvUpwind)
				m_spConvUpwind->update(&geo, StdVel);
	}

//	get a const (!!) reference to the stabilization
	const INavierStokesFV1Stabilization<dim>& stab = *m_spStab;
	const INavierStokesFV1Stabilization<dim>& convStab = *m_spConvStab;
	const INavierStokesUpwind<dim>& upwind = *m_spConvUpwind;

// 	loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// 	get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

	// 	loop shape functions
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
			////////////////////////////////////////////////////
			////////////////////////////////////////////////////
			// Momentum Equation (conservation of momentum)
			////////////////////////////////////////////////////
			////////////////////////////////////////////////////

			////////////////////////////////////////////////////
			// Diffusive Term (Momentum Equation)
			////////////////////////////////////////////////////

		// 	Compute flux derivative at IP
			const number flux_sh =  -1.0 * m_imKinViscosity[ip] * m_imDensitySCVF[ip]
									* VecDot(scvf.global_grad(sh), scvf.normal());

		// 	Add flux derivative  to local matrix
			for(int d1 = 0; d1 < dim; ++d1)
			{
				J(d1, scvf.from(), d1, sh) += flux_sh;
				J(d1, scvf.to()  , d1, sh) -= flux_sh;
			}
			
			if(!m_bLaplace)
			{
				for(int d1 = 0; d1 < dim; ++d1)
					for(int d2 = 0; d2 < dim; ++d2)
					{
						const number flux2_sh = -1.0 * m_imKinViscosity[ip] * m_imDensitySCVF[ip]
												* scvf.global_grad(sh)[d1] * scvf.normal()[d2];
						J(d1, scvf.from(), d2, sh) += flux2_sh;
						J(d1, scvf.to()  , d2, sh) -= flux2_sh;
					}
			}

			////////////////////////////////////////////////////
			// Pressure Term (Momentum Equation)
			////////////////////////////////////////////////////

		//	Add flux derivative for local matrix
			for(int d1 = 0; d1 < dim; ++d1)
			{
				const number flux_sh = scvf.shape(sh) * scvf.normal()[d1];
				J(d1, scvf.from(), _P_, sh) += flux_sh;
				J(d1, scvf.to()  , _P_, sh) -= flux_sh;
			}

			////////////////////////////////////////////////////
			// Convective Term (Momentum Equation)
			////////////////////////////////////////////////////

			// note: StdVel will be the advecting velocity (not upwinded)
			//       UpwindVel/StabVel will be the transported velocity

			if (! m_bStokes) // no convective terms in the Stokes equation
			{
			//	compute upwind velocity
				MathVector<dim> UpwindVel;
	
			//	switch PAC
				if(m_spConvUpwind.valid())  UpwindVel = upwind.upwind_vel(ip, u, StdVel);
				else if (m_spConvStab.valid()) UpwindVel = convStab.stab_vel(ip);
				else UG_THROW("Cannot find upwind for convective term.");
	
			//	peclet blend
				number w = 1.0;
				if(m_bPecletBlend)
					w = peclet_blend(UpwindVel, geo, ip, StdVel[ip], m_imKinViscosity[ip]);

			//	compute product of stabilized vel and normal
				const number prod = VecProd(StdVel[ip], scvf.normal()) * m_imDensitySCVF[ip];
	
			///////////////////////////////////
			//	Add fixpoint linearization
			///////////////////////////////////
	
			//	Stabilization used as upwind
				if(m_spConvStab.valid())
				{
				//	velocity derivatives
					if(stab.vel_comp_connected())
						for(int d1 = 0; d1 < dim; ++d1)
							for(int d2 = 0; d2 < dim; ++d2)
							{
								const number convFlux_vel = prod * w * convStab.stab_shape_vel(ip, d1, d2, sh);
								J(d1, scvf.from(), d2, sh) += convFlux_vel;
								J(d1, scvf.to()  , d2, sh) -= convFlux_vel;
							}
					else
						for(int d1 = 0; d1 < dim; ++d1)
						{
							const number convFlux_vel = prod * w * convStab.stab_shape_vel(ip, d1, d1, sh);
							J(d1, scvf.from(), d1, sh) += convFlux_vel;
							J(d1, scvf.to()  , d1, sh) -= convFlux_vel;
						}
	
				//	pressure derivative
					for(int d1 = 0; d1 < dim; ++d1)
					{
						const number convFlux_p = prod * w * convStab.stab_shape_p(ip, d1, sh);
	
						J(d1, scvf.from(), _P_, sh) += convFlux_p;
						J(d1, scvf.to()  , _P_, sh) -= convFlux_p;
					}
				}
	
			//	Upwind used as upwind
				if(m_spConvUpwind.valid())
				{
					number convFlux_vel = upwind.upwind_shape_sh(ip, sh);

				//	in some cases (e.g. PositiveUpwind, RegularUpwind) the upwind
				//	velocity in an ip depends also on the upwind velocity in
				//	other ips. This is reflected by the fact, that the ip
				//	shapes are non-zero. In that case, we can interpolate an
				//	approximate upwind only from the corner velocities by using
				//	u_up = \sum shape_co U_co + \sum shape_ip \tilde{u}_ip
				//	     = \sum shape_co U_co + \sum \sum shape_ip norm_shape_co|_ip * U_co
					if(m_spConvUpwind->non_zero_shape_ip())
					{
						for(size_t ip2 = 0; ip2 < geo.num_scvf(); ++ip2)
						{
							const typename TFVGeom::SCVF& scvf2 = geo.scvf(ip2);
							convFlux_vel += scvf2.shape(sh) * upwind.upwind_shape_ip(ip, ip2);
						}
					}

					convFlux_vel *= prod * w;

					for(int d1 = 0; d1 < dim; ++d1)
					{
						J(d1, scvf.from(), d1, sh) += convFlux_vel;
						J(d1, scvf.to()  , d1, sh) -= convFlux_vel;
					}
				}
	
			//	derivative due to peclet blending
				if(m_bPecletBlend)
				{
					const number convFluxPe = prod * (1.0-w) * scvf.shape(sh);
					for(int d1 = 0; d1 < dim; ++d1)
					{
						J(d1, scvf.from(), d1, sh) += convFluxPe;
						J(d1, scvf.to()  , d1, sh) -= convFluxPe;
					}
				}
	
			/////////////////////////////////////////
			//	Add full jacobian (remaining part)
			/////////////////////////////////////////
	
			//	Add remaining term for exact jacobian
				if(m_bFullNewtonFactor)
				{
				//	Stabilization used as upwind
					if(m_spConvStab.valid())
					{
					//	loop defect components
						for(int d1 = 0; d1 < dim; ++d1)
						{
							for(int d2 = 0; d2 < dim; ++d2)
							{
						//	derivatives w.r.t. velocity
						//	Compute n * derivs
							number prod_vel = 0.0;
	
						//	Compute sum_j n_j * \partial_{u_i^sh} u_j
							if(stab.vel_comp_connected())
								for(size_t k = 0; k < (size_t)dim; ++k)
									prod_vel += w * convStab.stab_shape_vel(ip, k, d2, sh)
													* scvf.normal()[k];
							else
								prod_vel = convStab.stab_shape_vel(ip, d1, d1, sh)
													* scvf.normal()[d1];
	
							prod_vel *= m_bFullNewtonFactor * m_imDensitySCVF[ip];

							J(d1, scvf.from(), d2, sh) += prod_vel * UpwindVel[d1];
							J(d1, scvf.to()  , d2, sh) -= prod_vel * UpwindVel[d1];
							}
	
						//	derivative w.r.t pressure
						//	Compute n * derivs
							number prod_p = 0.0;
	
						//	Compute sum_j n_j * \parial_{u_i^sh} u_j
							for(int k = 0; k < dim; ++k)
								prod_p += convStab.stab_shape_p(ip, k, sh)
													* scvf.normal()[k];
	
							prod_p *= m_bFullNewtonFactor * m_imDensitySCVF[ip];

							J(d1, scvf.from(), _P_, sh) += prod_p * UpwindVel[d1];
							J(d1, scvf.to()  , _P_, sh) -= prod_p * UpwindVel[d1];
						}
					}
	
				//	Upwind used as upwind
					if(m_spConvUpwind.valid())
					{
					//	loop defect components
						for(int d1 = 0; d1 < dim; ++d1)
							for(int d2 = 0; d2 < dim; ++d2)
							{
							//	derivatives w.r.t. velocity
								number prod_vel = w * upwind.upwind_shape_sh(ip,sh)
													* scvf.normal()[d2] * m_imDensitySCVF[ip];

								J(d1, scvf.from(), d2, sh) += prod_vel * UpwindVel[d1];
								J(d1, scvf.to()  , d2, sh) -= prod_vel * UpwindVel[d1];
							}
					}
	
				//	derivative due to peclet blending
					if(m_bPecletBlend)
					{
						for(int d1 = 0; d1 < dim; ++d1)
							for(int d2 = 0; d2 < dim; ++d2)
							{
								const number convFluxPe = UpwindVel[d1] * (1.0-w)
														  * scvf.shape(sh)
														  * scvf.normal()[d2]
														  * m_imDensitySCVF[ip];
								J(d1, scvf.from(), d2, sh) += convFluxPe;
								J(d1, scvf.to()  , d2, sh) -= convFluxPe;
							}
					}
				} // end exact jacobian part
			
			} // end of if (! m_bStokes) for the convective terms

			////////////////////////////////////////////////////
			////////////////////////////////////////////////////
			// Continuity Equation (conservation of mass)
			////////////////////////////////////////////////////
			////////////////////////////////////////////////////

		//	Add derivative of stabilized flux w.r.t velocity comp to local matrix
			if(stab.vel_comp_connected())
			{
				for(int d1 = 0; d1 < dim; ++d1)
				{
					number contFlux_vel = 0.0;
					for(int d2 = 0; d2 < dim; ++d2)
						contFlux_vel += stab.stab_shape_vel(ip, d2, d1, sh)
										* scvf.normal()[d2] * m_imDensitySCVF[ip];

					J(_P_, scvf.from(), d1, sh) += contFlux_vel;
					J(_P_, scvf.to()  , d1, sh) -= contFlux_vel;
				}
			}
			else
			{
				for(int d1 = 0; d1 < dim; ++d1)
				{
					const number contFlux_vel = stab.stab_shape_vel(ip, d1, d1, sh)
													* scvf.normal()[d1] * m_imDensitySCVF[ip];

					J(_P_, scvf.from(), d1, sh) += contFlux_vel;
					J(_P_, scvf.to()  , d1, sh) -= contFlux_vel;
				}
			}

		//	Add derivative of stabilized flux w.r.t pressure to local matrix
			number contFlux_p = 0.0;
			for(int d1 = 0; d1 < dim; ++d1)
				contFlux_p += stab.stab_shape_p(ip, d1, sh) * scvf.normal()[d1] * m_imDensitySCVF[ip];

			J(_P_, scvf.from(), _P_, sh) += contFlux_p;
			J(_P_, scvf.to()  , _P_, sh) -= contFlux_p;
		}
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesFV1<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	Only first order implemented
	UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	check for source term to pass to the stabilization
	const DataImport<MathVector<dim>, dim>* pSource = NULL;
	if(m_imSourceSCVF.data_given())	pSource = &m_imSourceSCVF;

//	check for solutions to pass to stabilization in time-dependent case
	const LocalVector *pSol = &u, *pOldSol = NULL;
	number dt = 0.0;
	if(this->is_time_dependent())
	{
	//	get and check current and old solution 
		const LocalVectorTimeSeries* vLocSol = this->local_time_solutions();
		if(vLocSol->size() != 2)
			UG_THROW("NavierStokes::add_def_A_elem: "
							" Stabilization needs exactly two time points.");

	//	remember local solutions
		pSol = &vLocSol->solution(0);
		pOldSol = &vLocSol->solution(1);
		dt = vLocSol->time(0) - vLocSol->time(1);
	}

//	interpolate velocity at ip with standard lagrange interpolation
	static const size_t numSCVF = TFVGeom::numSCVF;
	MathVector<dim> StdVel[numSCVF];
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
		VecSet(StdVel[ip], 0.0);

		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			for(int d1 = 0; d1 < dim; ++d1)
				StdVel[ip][d1] += u(d1, sh) * scvf.shape(sh);
	}

//	compute stabilized velocities and shapes for continuity equation
	// \todo: (optional) Here we can skip the computation of shapes, implement?
	m_spStab->update(&geo, *pSol, StdVel, m_bStokes, m_imKinViscosity, m_imDensitySCVF, pSource, pOldSol, dt);

	if (! m_bStokes) // no convective terms in the Stokes eq. => no upwinding
	{
	//	compute stabilized velocities and shapes for convection upwind
		if(m_spConvStab.valid())
			if(m_spConvStab != m_spStab)
				m_spConvStab->update(&geo, *pSol, StdVel, false, m_imKinViscosity, m_imDensitySCVF, pSource, pOldSol, dt);
	
	//	compute upwind shapes
		if(m_spConvUpwind.valid())
			if(m_spStab->upwind() != m_spConvUpwind)
				m_spConvUpwind->update(&geo, StdVel);
	}

//	get a const (!!) reference to the stabilization
	const INavierStokesFV1Stabilization<dim>& stab = *m_spStab;
	const INavierStokesFV1Stabilization<dim>& convStab = *m_spConvStab;
	const INavierStokesUpwind<dim>& upwind = *m_spConvUpwind;

// 	loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// 	get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

		////////////////////////////////////////////////////
		////////////////////////////////////////////////////
		// Momentum Equation (conservation of momentum)
		////////////////////////////////////////////////////
		////////////////////////////////////////////////////

	//	\todo: we could also add all fluxes at once in order to save several
	//			accesses to the local defect, implement and loose clarity?

		////////////////////////////////////////////////////
		// Diffusive Term (Momentum Equation)
		////////////////////////////////////////////////////

	// 	1. Interpolate Functional Matrix of velocity at ip
		MathMatrix<dim, dim> gradVel;
		for(int d1 = 0; d1 < dim; ++d1)
			for(int d2 = 0; d2 <dim; ++d2)
			{
			//	sum up contributions of each shape
				gradVel(d1, d2) = 0.0;
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					gradVel(d1, d2) += scvf.global_grad(sh)[d2]
					                    * u(d1, sh);
			}

	//	2. Compute flux
		MathVector<dim> diffFlux;

	//	Add (\nabla u) \cdot \vec{n}
		MatVecMult(diffFlux, gradVel, scvf.normal());

	//	Add (\nabla u)^T \cdot \vec{n}
		if(!m_bLaplace)
			TransposedMatVecMultAdd(diffFlux, gradVel, scvf.normal());

	//	scale by viscosity
		VecScale(diffFlux, diffFlux, (-1.0) * m_imKinViscosity[ip] * m_imDensitySCVF[ip]);

	//	3. Add flux to local defect
		for(int d1 = 0; d1 < dim; ++d1)
		{
			d(d1, scvf.from()) += diffFlux[d1];
			d(d1, scvf.to()  ) -= diffFlux[d1];
		}

		////////////////////////////////////////////////////
		// Convective Term (Momentum Equation)
		////////////////////////////////////////////////////

		// note: StdVel will be the advecting velocity (not upwinded)
		//       UpwindVel/StabVel will be the transported velocity

		if (! m_bStokes) // no convective terms in the Stokes equation
		{
		//	find the upwind velocity at ip
			MathVector<dim> UpwindVel;
	
		//	switch PAC
			if(m_spConvUpwind.valid())  UpwindVel = upwind.upwind_vel(ip, u, StdVel);
			else if (m_spConvStab.valid()) UpwindVel = convStab.stab_vel(ip);
			else UG_THROW("Cannot find upwind for convective term.");
	
		//	Peclet Blend
			if(m_bPecletBlend)
				peclet_blend(UpwindVel, geo, ip, StdVel[ip], m_imKinViscosity[ip]);
	
		//	compute product of standard velocity and normal
			const number prod = VecProd(StdVel[ip], scvf.normal()) * m_imDensitySCVF[ip];
	
		//	Add contributions to local velocity components
			for(int d1 = 0; d1 < dim; ++d1)
			{
				d(d1, scvf.from()) += UpwindVel[d1] * prod;
				d(d1, scvf.to()  ) -= UpwindVel[d1] * prod;
			}
		}

		////////////////////////////////////////////////////
		// Pressure Term (Momentum Equation)
		////////////////////////////////////////////////////

	//	1. Interpolate pressure at ip
		number pressure = 0.0;
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			pressure += scvf.shape(sh) * u(_P_, sh);

	//	2. Add contributions to local defect
		for(int d1 = 0; d1 < dim; ++d1)
		{
			d(d1, scvf.from()) += pressure * scvf.normal()[d1];
			d(d1, scvf.to()  ) -= pressure * scvf.normal()[d1];
		}

		////////////////////////////////////////////////////
		////////////////////////////////////////////////////
		// Continuity Equation (conservation of mass)
		////////////////////////////////////////////////////
		////////////////////////////////////////////////////

	//	compute flux at ip
		const number contFlux = VecProd(stab.stab_vel(ip), scvf.normal()) * m_imDensitySCVF[ip];

	//	Add contributions to local defect
		d(_P_, scvf.from()) += contFlux;
		d(_P_, scvf.to()  ) -= contFlux;
	}
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesFV1<TDomain>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	Only first order implementation
	UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int sh = scv.node_id();

	// 	loop velocity components
		for(int d1 = 0; d1 < dim; ++d1)
		{
		// 	Add to local matrix
			J(d1, sh, d1, sh) += scv.volume() * m_imDensitySCV[ip];
		}
	}
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesFV1<TDomain>::
add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	Only first order implementation
	UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int sh = scv.node_id();

	// 	loop velocity components
		for(int d1 = 0; d1 < dim; ++d1)
		{
		// 	Add to local matrix
			d(d1, sh) += u(d1, sh) * scv.volume() * m_imDensitySCV[ip];
		}
	}
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesFV1<TDomain>::
add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	Only first order implementation
	UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

//	if zero data given, return
	if(!m_imSourceSCV.data_given()) return;

// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int sh = scv.node_id();

	// 	Add to local rhs
		for(int d1 = 0; d1 < dim; ++d1){
			d(d1, sh) += m_imSourceSCV[ip][d1] * scv.volume() * m_imDensitySCV[ip];
		}
	}
}

template<typename TDomain>
template<typename TFVGeom>
inline
number
NavierStokesFV1<TDomain>::
peclet_blend(MathVector<dim>& UpwindVel, const TFVGeom& geo, size_t ip,
             const MathVector<dim>& StdVel, number kinVisco)
{
	const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
//	compute peclet number
	number Pe = VecProd(StdVel, scvf.normal())/VecTwoNormSq(scvf.normal())
	 * VecDistance(geo.corners() [scvf.to()], geo.corners() [scvf.from()]) / kinVisco;

//	compute weight
	const number Pe2 = Pe * Pe;
	const number w = Pe2 / (5.0 + Pe2);

//	compute upwind vel
	VecScaleAdd(UpwindVel, w, UpwindVel, (1.0-w), StdVel);

	return w;
}


////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template<>
void NavierStokesFV1<Domain1d>::
register_all_funcs(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
	}
	else
	{
		UG_THROW("NavierStokesFV1: Hanging Nodes not implemented.")
	}
}
#endif

#ifdef UG_DIM_2
template<>
void NavierStokesFV1<Domain2d>::
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
		UG_THROW("NavierStokesFV1: Hanging Nodes not implemented.")
	}
}
#endif

#ifdef UG_DIM_3
template<>
void NavierStokesFV1<Domain3d>::
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
		UG_THROW("NavierStokesFV1: Hanging Nodes not implemented.")
	}
}
#endif

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void
NavierStokesFV1<TDomain>::
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

#ifdef UG_DIM_2
template class NavierStokesFV1<Domain2d>;
#endif
#ifdef UG_DIM_3
template class NavierStokesFV1<Domain3d>;
#endif

} // namespace NavierStokes
} // namespace ug
