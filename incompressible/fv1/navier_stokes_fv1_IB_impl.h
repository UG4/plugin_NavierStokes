/*
 * navier_stokes_fv1_IB_impl.h
 *
 *  Created on: 06.02.2014
 *      Author: suze
 */

#ifndef NAVIER_STOKES_FV1_IB_IMPL_H_
#define NAVIER_STOKES_FV1_IB_IMPL_H_

#include "navier_stokes_fv1_IB.h"

#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/fv1ib_geom.h" //-> in register_substitutions: FVGeometryIB :-)
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"


namespace ug{
namespace NavierStokes{

////////////////////////////////////////////////////////////////////////////////
//	Constructor - set default values
////////////////////////////////////////////////////////////////////////////////
/*
template<typename TDomain>
void register_substitution_funcs(bool bHang)
{
#ifdef UG_DIM_1
//	switch assemble functions
	if(!bHang)
	{
		//register_substitutes<Edge, FV1Geometry<Edge, dim> >();
		register_substitutes<Edge, FV1Geometry<Edge, 1> >();
	}
	else
	{
		UG_THROW("NavierStokesFV1: Hanging Nodes not implemented.")
	}
#endif
#ifdef UG_DIM_2
	if(!bHang)
	{
		register_substitutes<Triangle, FV1Geometry<Triangle, 2> >();
		//register_substitutes<Triangle, FV1IBGeometry<Triangle, dim> >();
		register_substitutes<Quadrilateral, FV1Geometry<Quadrilateral, 2> >();
		//register_substitutes<Quadrilateral, FV1IBGeometry<Quadrilateral, dim> >();
	}
	else
	{
		UG_THROW("NavierStokesFV1IB: Hanging Nodes not implemented.")
	}
#endif
#ifdef UG_DIM_3
	if(!bHang)
	{
		register_substitutes<Tetrahedron, FV1Geometry<Tetrahedron, 3> >();
		//register_substitutes<Tetrahedron, FV1IBGeometry<Tetrahedron, dim> >();

		register_substitutes<Prism, FV1Geometry<Prism, 3> >();
		//register_substitutes<Prism, FV1IBGeometry<Prism, dim> >();

		register_substitutes<Pyramid, FV1Geometry<Pyramid, 3> >();
		//register_substitutes<Pyramid, FV1IBGeometry<Pyramid, dim> >();

		register_substitutes<Hexahedron, FV1Geometry<Hexahedron, 3> >();
		//register_substitutes<Hexahedron, FV1IBGeometry<Hexahedron, dim> >();

 	}
	else
	{
		UG_THROW("NavierStokesFV1IB: Hanging Nodes not implemented.")
	}
#endif
}
*/
template<typename TDomain>
NavierStokesFV1IB<TDomain>::NavierStokesFV1IB(const char* functions,
                                          const char* subsets)
: NavierStokesFV1<TDomain>(functions, subsets)
{
	UG_LOG("call init()...\n");

	init();


};


template<typename TDomain>
void NavierStokesFV1IB<TDomain>::init()
{
	UG_LOG("IN init() 1...\n");
	//this->NavierStokesFV1<TDomain>::init(); // passiert schon beim constructor von NavierStokesFV1

	register_substitution_funcs(false);

}




template<typename TDomain>
void NavierStokesFV1IB<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{

	this->NavierStokesFV1<TDomain>::template prepare_setting(vLfeID, bNonRegularGrid);

	//	update assemble functions
	register_substitution_funcs(false);
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesFV1IB<TDomain>::
prep_elem_IB(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{


// 	Update Geometry for this element
  	FV1IBGeometry<TElem, dim>& geo = GeomProvider<FV1IBGeometry<TElem, dim> >::get();
  	//DimFV1IBGeometry<dim>& geo = GeomProvider<DimFV1IBGeometry<dim> >::get(LFEID(LFEID::LAGRANGE, dim, 1), 1);
	try{
		geo.update(elem, vCornerCoords, &(this->subset_handler()));
	}
	UG_CATCH_THROW("NavierStokes_IB::prep_elem:"
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
		this->NavierStokesFV1<TDomain>::m_imKinViscosity.template set_local_ips<refDim>(vSCVFip,numSCVFip);
		this->NavierStokesFV1<TDomain>::m_imDensitySCVF.template set_local_ips<refDim>(vSCVFip,numSCVFip);
		this->NavierStokesFV1<TDomain>::m_imDensitySCV.template set_local_ips<refDim>(vSCVip,numSCVip);
		this->NavierStokesFV1<TDomain>::m_imSourceSCV.template set_local_ips<refDim>(vSCVip,numSCVip);
		this->NavierStokesFV1<TDomain>::m_imSourceSCVF.template set_local_ips<refDim>(vSCVFip,numSCVFip);
	}

//	set global positions for imports
	const MathVector<dim>* vSCVFip = geo.scvf_global_ips();
	const size_t numSCVFip = geo.num_scvf_ips();
	const MathVector<dim>* vSCVip = geo.scv_global_ips();
	const size_t numSCVip = geo.num_scv_ips();
	this->NavierStokesFV1<TDomain>::m_imKinViscosity.set_global_ips(vSCVFip, numSCVFip);
	this->NavierStokesFV1<TDomain>::m_imDensitySCVF.set_global_ips(vSCVFip, numSCVFip);
	this->NavierStokesFV1<TDomain>::m_imDensitySCV.set_global_ips(vSCVip, numSCVip);
	this->NavierStokesFV1<TDomain>::m_imSourceSCV.set_global_ips(vSCVip, numSCVip);
	this->NavierStokesFV1<TDomain>::m_imSourceSCVF.set_global_ips(vSCVFip, numSCVFip);

// first get elem_data for adaptions from ParticleProvider:
	//MathVector<dim>* IBNormals = m_myParticle->get_normals(elem);
	//MathVector<dim>* IBIntegrationPoints = m_myParticle->get_ips(elem);

	UG_LOG("in prep_elem...\n");
	m_myParticle->compute_normals(geo, elem);

 //  finally call the methods of the 'FV1IBGeometry'-class to adapt the element geometry data
//  for elements which are cut by the inner boundary:
	geo.adapt(elem, vCornerCoords);
	//geo.adapt_normals(elem, vCornerCoords);
	//geo.adapt_integration_points(elem, vCornerCoords);


}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesFV1IB<TDomain>::
add_jac_A_elem_IB(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{

	size_t _P_ = dim;
	// 	Only first order implementation
		UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

	// 	get finite volume geometry
		static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	//	check for source term to pass to the stabilization
		const DataImport<MathVector<dim>, dim>* pSource = NULL;
		if(this->m_imSourceSCVF.data_given())	pSource = &this->m_imSourceSCVF;

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
		this->m_spStab->update(&geo, *pSol, StdVel, this->m_bStokes, this->m_imKinViscosity, pSource, pOldSol, dt);

		if (! this->m_bStokes) // no convective terms in the Stokes eq. => no upwinding
		{
		//	compute stabilized velocities and shapes for convection upwind
			if(this->m_spConvStab.valid())
				if(this->m_spConvStab != this->m_spStab)
					this->m_spConvStab->update(&geo, *pSol, StdVel, false, this->m_imKinViscosity, pSource, pOldSol, dt);

		//	compute upwind shapes
			if(this->m_spConvUpwind.valid())
				if(this->m_spStab->upwind() != this->m_spConvUpwind)
					this->m_spConvUpwind->update(&geo, StdVel);
		}

	//	get a const (!!) reference to the stabilization
		const INavierStokesFV1Stabilization<dim>& stab = *this->m_spStab;
		const INavierStokesFV1Stabilization<dim>& convStab = *this->m_spConvStab;
		const INavierStokesUpwind<dim>& upwind = *this->m_spConvUpwind;

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
				const number flux_sh =  -1.0 * this->m_imKinViscosity[ip] * this->m_imDensitySCVF[ip]
										* VecDot(scvf.global_grad(sh), scvf.normal());

			// 	Add flux derivative  to local matrix
				for(int d1 = 0; d1 < dim; ++d1)
				{
					J(d1, scvf.from(), d1, sh) += flux_sh;
					J(d1, scvf.to()  , d1, sh) -= flux_sh;
				}

				if(!this->m_bLaplace)
				{
					for(int d1 = 0; d1 < dim; ++d1)
						for(int d2 = 0; d2 < dim; ++d2)
						{
							const number flux2_sh = -1.0 * this->m_imKinViscosity[ip] * this->m_imDensitySCVF[ip]
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

				if (! this->m_bStokes) // no convective terms in the Stokes equation
				{
				//	compute upwind velocity
					MathVector<dim> UpwindVel;

				//	switch PAC
					if(this->m_spConvUpwind.valid())  UpwindVel = upwind.upwind_vel(ip, u, StdVel);
					else if (this->m_spConvStab.valid()) UpwindVel = convStab.stab_vel(ip);
					else UG_THROW("Cannot find upwind for convective term.");

				//	peclet blend
					number w = 1.0;

				//	compute product of stabilized vel and normal
					const number prod = VecProd(StdVel[ip], scvf.normal()) * this->m_imDensitySCVF[ip];

				///////////////////////////////////
				//	Add fixpoint linearization
				///////////////////////////////////

				//	Stabilization used as upwind
					if(this->m_spConvStab.valid())
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
					if(this->m_spConvUpwind.valid())
					{
						number convFlux_vel = upwind.upwind_shape_sh(ip, sh);

					//	in some cases (e.g. PositiveUpwind, RegularUpwind) the upwind
					//	velocity in an ip depends also on the upwind velocity in
					//	other ips. This is reflected by the fact, that the ip
					//	shapes are non-zero. In that case, we can interpolate an
					//	approximate upwind only from the corner velocities by using
					//	u_up = \sum shape_co U_co + \sum shape_ip \tilde{u}_ip
					//	     = \sum shape_co U_co + \sum \sum shape_ip norm_shape_co|_ip * U_co
						if(this->m_spConvUpwind->non_zero_shape_ip())
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
					if(this->m_bPecletBlend)
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
					if(this->m_bFullNewtonFactor)
					{
					//	Stabilization used as upwind
						if(this->m_spConvStab.valid())
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

								prod_vel *= this->m_bFullNewtonFactor * this->m_imDensitySCVF[ip];

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

								prod_p *= this->m_bFullNewtonFactor * this->m_imDensitySCVF[ip];

								J(d1, scvf.from(), _P_, sh) += prod_p * UpwindVel[d1];
								J(d1, scvf.to()  , _P_, sh) -= prod_p * UpwindVel[d1];
							}
						}

					//	Upwind used as upwind
						if(this->m_spConvUpwind.valid())
						{
						//	loop defect components
							for(int d1 = 0; d1 < dim; ++d1)
								for(int d2 = 0; d2 < dim; ++d2)
								{
								//	derivatives w.r.t. velocity
									number prod_vel = w * upwind.upwind_shape_sh(ip,sh)
														* scvf.normal()[d2] * this->m_imDensitySCVF[ip];

									J(d1, scvf.from(), d2, sh) += prod_vel * UpwindVel[d1];
									J(d1, scvf.to()  , d2, sh) -= prod_vel * UpwindVel[d1];
								}
						}

					//	derivative due to peclet blending
						if(this->m_bPecletBlend)
						{
							for(int d1 = 0; d1 < dim; ++d1)
								for(int d2 = 0; d2 < dim; ++d2)
								{
									const number convFluxPe = UpwindVel[d1] * (1.0-w)
															  * scvf.shape(sh)
															  * scvf.normal()[d2]
															  * this->m_imDensitySCVF[ip];
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
											* scvf.normal()[d2] * this->m_imDensitySCVF[ip];

						J(_P_, scvf.from(), d1, sh) += contFlux_vel;
						J(_P_, scvf.to()  , d1, sh) -= contFlux_vel;
					}
				}
				else
				{
					for(int d1 = 0; d1 < dim; ++d1)
					{
						const number contFlux_vel = stab.stab_shape_vel(ip, d1, d1, sh)
														* scvf.normal()[d1] * this->m_imDensitySCVF[ip];

						J(_P_, scvf.from(), d1, sh) += contFlux_vel;
						J(_P_, scvf.to()  , d1, sh) -= contFlux_vel;
					}
				}

			//	Add derivative of stabilized flux w.r.t pressure to local matrix
				number contFlux_p = 0.0;
				for(int d1 = 0; d1 < dim; ++d1)
					contFlux_p += stab.stab_shape_p(ip, d1, sh) * scvf.normal()[d1] * this->m_imDensitySCVF[ip];

				J(_P_, scvf.from(), _P_, sh) += contFlux_p;
				J(_P_, scvf.to()  , _P_, sh) -= contFlux_p;
			}
			UG_LOG("scvf.normal() = " << scvf.normal() << "\n");

		}

	UG_LOG("in 'add_jac_A_elem_IB': jeppa!...\n");

}



template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesFV1IB<TDomain>::
add_def_A_elem_IB(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	size_t _P_ = dim;

	// 	Only first order implemented
		UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

	// 	get finite volume geometry
		static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	//	check for source term to pass to the stabilization
		const DataImport<MathVector<dim>, dim>* pSource = NULL;
		if(this->m_imSourceSCVF.data_given())	pSource = &this->m_imSourceSCVF;

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
		this->m_spStab->update(&geo, *pSol, StdVel, this->m_bStokes, this->m_imKinViscosity, pSource, pOldSol, dt);

		if (! this->m_bStokes) // no convective terms in the Stokes eq. => no upwinding
		{
		//	compute stabilized velocities and shapes for convection upwind
			if(this->m_spConvStab.valid())
				if(this->m_spConvStab != this->m_spStab)
					this->m_spConvStab->update(&geo, *pSol, StdVel, false, this->m_imKinViscosity, pSource, pOldSol, dt);

		//	compute upwind shapes
			if(this->m_spConvUpwind.valid())
				if(this->m_spStab->upwind() != this->m_spConvUpwind)
					this->m_spConvUpwind->update(&geo, StdVel);
		}

	//	get a const (!!) reference to the stabilization
		const INavierStokesFV1Stabilization<dim>& stab = *this->m_spStab;
		const INavierStokesFV1Stabilization<dim>& convStab = *this->m_spConvStab;
		const INavierStokesUpwind<dim>& upwind = *this->m_spConvUpwind;

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
			if(!this->m_bLaplace)
				TransposedMatVecMultAdd(diffFlux, gradVel, scvf.normal());

		//	scale by viscosity
			VecScale(diffFlux, diffFlux, (-1.0) * this->m_imKinViscosity[ip] * this->m_imDensitySCVF[ip]);

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

			if (! this->m_bStokes) // no convective terms in the Stokes equation
			{
			//	find the upwind velocity at ip
				MathVector<dim> UpwindVel;

			//	switch PAC
				if(this->m_spConvUpwind.valid())  UpwindVel = upwind.upwind_vel(ip, u, StdVel);
				else if (this->m_spConvStab.valid()) UpwindVel = convStab.stab_vel(ip);
				else UG_THROW("Cannot find upwind for convective term.");

			//	compute product of standard velocity and normal
				const number prod = VecProd(StdVel[ip], scvf.normal()) * this->m_imDensitySCVF[ip];

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
			const number contFlux = VecProd(stab.stab_vel(ip), scvf.normal()) * this->m_imDensitySCVF[ip];

		//	Add contributions to local defect
			d(_P_, scvf.from()) += contFlux;
			d(_P_, scvf.to()  ) -= contFlux;
		}

}

////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template<>
void NavierStokesFV1IB<Domain1d>::
register_substitution_funcs(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
<<<<<<< .mine
		register_substitutes<Edge, FV1Geometry<Edge, dim> >();
		//register_substitutes<Edge, FV1IBGeometry<Edge, dim> >();
=======
		//register_substitutes<RegularEdge, FV1Geometry<RegularEdge, dim> >();
		register_substitutes<RegularEdge, FV1IBGeometry<RegularEdge, dim> >();
>>>>>>> .r13732
	}
	else
	{
		UG_THROW("NavierStokesFV1: Hanging Nodes not implemented.")
	}
}
#endif

#ifdef UG_DIM_2
template<>
void NavierStokesFV1IB<Domain2d>::
register_substitution_funcs(bool bHang)
{
	UG_LOG("in register_substitution_funcs() IB...\n");
//	switch assemble functions
	if(!bHang)
	{
		register_substitutes<Triangle, FV1Geometry<Triangle, dim> >();
		//register_substitutes<Triangle, FV1IBGeometry<Triangle, dim> >();
		register_substitutes<Quadrilateral, FV1Geometry<Quadrilateral, dim> >();
		//register_substitutes<Quadrilateral, FV1IBGeometry<Quadrilateral, dim> >();
	}
	else
	{
		UG_THROW("NavierStokesFV1IB: Hanging Nodes not implemented.")
	}
}
#endif

#ifdef UG_DIM_3
template<>
void NavierStokesFV1IB<Domain3d>::
register_substitution_funcs(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_substitutes<Tetrahedron, FV1Geometry<Tetrahedron, dim> >();
		//register_substitutes<Tetrahedron, FV1IBGeometry<Tetrahedron, dim> >();

		register_substitutes<Prism, FV1Geometry<Prism, dim> >();
		//register_substitutes<Prism, FV1IBGeometry<Prism, dim> >();

		register_substitutes<Pyramid, FV1Geometry<Pyramid, dim> >();
		//register_substitutes<Pyramid, FV1IBGeometry<Pyramid, dim> >();

		register_substitutes<Hexahedron, FV1Geometry<Hexahedron, dim> >();
		//register_substitutes<Hexahedron, FV1IBGeometry<Hexahedron, dim> >();

 	}
	else
	{
		UG_THROW("NavierStokesFV1IB: Hanging Nodes not implemented.")
	}
}
#endif


/*
template<typename TAlgebra>
class NavierStokesFV1IB<Domain1d,TAlgebra>
{
	void register_substitution_funcs(bool bHang)
	{
		if(!bHang)
		{
			//register_substitutes<Edge, FV1Geometry<Edge, dim> >();
			register_substitutes<Edge, FV1Geometry<Edge, dim> >();
		}
		else
		{
			UG_THROW("NavierStokesFV1: Hanging Nodes not implemented.")
		}
	}
};
template<typename TAlgebra>
class NavierStokesFV1IB<Domain2d,TAlgebra>
{
	void register_substitution_funcs(bool bHang)
	{
		//	switch assemble functions
		if(!bHang)
		{
			register_substitutes<Triangle, FV1Geometry<Triangle, dim> >();
			//register_substitutes<Triangle, FV1IBGeometry<Triangle, dim> >();
			register_substitutes<Quadrilateral, FV1Geometry<Quadrilateral, dim> >();
			//register_substitutes<Quadrilateral, FV1IBGeometry<Quadrilateral, dim> >();
		}
		else
		{
			UG_THROW("NavierStokesFV1IB: Hanging Nodes not implemented.")
		}
	}
};
template<typename TAlgebra>
class NavierStokesFV1IB<Domain3d,TAlgebra>
{
	void register_substitution_funcs(bool bHang)
	{
		//	switch assemble functions
		if(!bHang)
		{
			register_substitutes<Tetrahedron, FV1Geometry<Tetrahedron, dim> >();
			//register_substitutes<Tetrahedron, FV1IBGeometry<Tetrahedron, dim> >();

			register_substitutes<Prism, FV1Geometry<Prism, dim> >();
			//register_substitutes<Prism, FV1IBGeometry<Prism, dim> >();

			register_substitutes<Pyramid, FV1Geometry<Pyramid, dim> >();
			//register_substitutes<Pyramid, FV1IBGeometry<Pyramid, dim> >();

			register_substitutes<Hexahedron, FV1Geometry<Hexahedron, dim> >();
			//register_substitutes<Hexahedron, FV1IBGeometry<Hexahedron, dim> >();

		}
		else
		{
			UG_THROW("NavierStokesFV1IB: Hanging Nodes not implemented.")
		}
	}
};
*/
/*


template<typename TDomain>
void NavierStokesFV1IB<TDomain>::
register_substitution_funcs(bool bHang)
{
#ifdef UG_DIM_1
//	switch assemble functions
	if(!bHang)
	{
		//register_substitutes<Edge, FV1Geometry<Edge, dim> >();
		register_substitutes<Edge, FV1Geometry<Edge, 1> >();
	}
	else
	{
		UG_THROW("NavierStokesFV1: Hanging Nodes not implemented.")
	}
#endif
#ifdef UG_DIM_2
	if(!bHang)
	{
		register_substitutes<Triangle, FV1Geometry<Triangle, 2> >();
		//register_substitutes<Triangle, FV1IBGeometry<Triangle, dim> >();
		register_substitutes<Quadrilateral, FV1Geometry<Quadrilateral, 2> >();
		//register_substitutes<Quadrilateral, FV1IBGeometry<Quadrilateral, dim> >();
	}
	else
	{
		UG_THROW("NavierStokesFV1IB: Hanging Nodes not implemented.")
	}
#endif
#ifdef UG_DIM_3
	if(!bHang)
	{
		register_substitutes<Tetrahedron, FV1Geometry<Tetrahedron, 3> >();
		//register_substitutes<Tetrahedron, FV1IBGeometry<Tetrahedron, dim> >();

		register_substitutes<Prism, FV1Geometry<Prism, 3> >();
		//register_substitutes<Prism, FV1IBGeometry<Prism, dim> >();

		register_substitutes<Pyramid, FV1Geometry<Pyramid, 3> >();
		//register_substitutes<Pyramid, FV1IBGeometry<Pyramid, dim> >();

		register_substitutes<Hexahedron, FV1Geometry<Hexahedron, 3> >();
		//register_substitutes<Hexahedron, FV1IBGeometry<Hexahedron, dim> >();

 	}
	else
	{
		UG_THROW("NavierStokesFV1IB: Hanging Nodes not implemented.")
	}
#endif
}

*/

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesFV1IB<TDomain>::
register_substitutes()
{

	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->enable_fast_add_elem(true);

	this->set_prep_elem_fct(	 	id, &T::template prep_elem<TElem, TFVGeom>);
   	this->set_add_jac_A_elem_fct(	id, &T::template add_jac_A_elem<TElem, TFVGeom>);
 	this->set_add_def_A_elem_fct(	id, &T::template add_def_A_elem<TElem, TFVGeom>);


}

} // namespace NavierStokes
} // namespace ug

#endif /* NAVIER_STOKES_FV1_IB_IMPL_H_ */
