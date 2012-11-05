/*
 * navier_stokes_fv1.h
 *
 *  Created on: 20.09.2010
 *      Author: andreasvogel
 */

#include "navier_stokes.h"

#include "common/util/provider.h"
#include "lib_disc/spatial_disc/disc_util/finite_volume_geometry.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{
namespace NavierStokes{

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Provide a generic implementation for all elements
//  (since this discretization can be implemented in a generic way)
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void NavierStokes<TDomain>::
prepare_element_loop_fv1()
{
// 	Only first order implementation
	if(!(TFVGeom<TElem, dim>::order == 1))
		UG_THROW("Only first order implementation, but other Finite Volume"
						" Geometry set.");

//	check, that stabilization has been set
	if(m_spStab.invalid())
		UG_THROW("Stabilization has not been set.");

//	init stabilization for element type
	m_spStab->template set_geometry_type<TFVGeom<TElem, dim> >();

	if (! m_bStokes) // no convective terms in the Stokes eq. => no upwinding
	{
	//	check, that convective upwinding has been set
		if(m_spConvStab.invalid()  && m_spConvUpwind.invalid())
			UG_THROW("Upwinding for convective Term in Momentum eq. not set.");
	
	//	init convection stabilization for element type
		if(m_spConvStab.valid())
			m_spConvStab->template set_geometry_type<TFVGeom<TElem, dim> >();
	
	//	init convection stabilization for element type
		if(m_spConvUpwind.valid())
			m_spConvUpwind->template set_geometry_type<TFVGeom<TElem, dim> >();
	}

//	check, that kinematic Viscosity has been set
	if(!m_imKinViscosity.data_given())
		UG_THROW("NavierStokes::prepare_element_loop:"
						" Kinematic Viscosity has not been set, but is required.");

//	check, that Density has been set
	if(!m_imDensitySCVF.data_given())
		UG_THROW("NavierStokes::prepare_element_loop:"
						" Density has not been set, but is required.");

//	check, that Density has been set
	if(!m_imDensitySCV.data_given())
		UG_THROW("NavierStokes::prepare_element_loop:"
						" Density has not been set, but is required.");

//	set local positions for imports
	typedef typename reference_element_traits<TElem>::reference_element_type
																ref_elem_type;
	static const int refDim = ref_elem_type::dim;

	if(!TFVGeom<TElem, dim>::usesHangingNodes)
	{
		TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();
		m_imKinViscosity.template set_local_ips<refDim>(geo.scvf_local_ips(),
		                                                geo.num_scvf_ips());
		m_imDensitySCVF.template set_local_ips<refDim>(geo.scvf_local_ips(),
		                                                geo.num_scvf_ips());
		m_imDensitySCV.template set_local_ips<refDim>(geo.scv_local_ips(),
		                                          geo.num_scv_ips());
		m_imSource.template set_local_ips<refDim>(geo.scv_local_ips(),
		                                          geo.num_scv_ips());
	}
}

template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void NavierStokes<TDomain>::
finish_element_loop_fv1()
{}


template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void NavierStokes<TDomain>::
prepare_element_fv1(TElem* elem, const LocalVector& u)
{
//	get corners
	m_vCornerCoords = this->template element_corners<TElem>(elem);

// 	Update Geometry for this element
	TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();
	if(!geo.update(elem, &m_vCornerCoords[0], &(this->subset_handler())))
		UG_THROW("NavierStokes::prepare_element:"
						" Cannot update Finite Volume Geometry.");

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
		m_imDensitySCVF.template set_local_ips<refDim>(geo.scvf_local_ips(),
		                                                geo.num_scvf_ips());
		m_imDensitySCV.template set_local_ips<refDim>(geo.scv_local_ips(),
		                                          geo.num_scv_ips());
		m_imSource.template set_local_ips<refDim>(geo.scv_local_ips(),
		                                          geo.num_scv_ips());
	}

//	set global positions for imports
	m_imKinViscosity.set_global_ips(geo.scvf_global_ips(), geo.num_scvf_ips());
	m_imDensitySCVF.set_global_ips(geo.scvf_global_ips(), geo.num_scvf_ips());
	m_imDensitySCV.set_global_ips(geo.scv_global_ips(), geo.num_scv_ips());
	m_imSource.set_global_ips(geo.scv_global_ips(), geo.num_scv_ips());
}

template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void NavierStokes<TDomain>::
ass_JA_elem_fv1(LocalMatrix& J, const LocalVector& u)
{
// 	Only first order implementation
	UG_ASSERT((TFVGeom<TElem, dim>::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	static const TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();

//	check for source term to pass to the stabilization
	const DataImport<MathVector<dim>, dim>* pSource = NULL;
	if(m_imSource.data_given())	pSource = &m_imSource;

//	check for solutions to pass to stabilization in time-dependent case
	const LocalVector *pSol = &u, *pOldSol = NULL;
	number dt = 0.0;
	if(this->is_time_dependent())
	{
	//	get and check current and old solution
		const LocalVectorTimeSeries* vLocSol = this->local_time_solutions();
		if(vLocSol->size() != 2)
			UG_THROW("NavierStokes::ass_JA_elem: "
							" Stabilization needs exactly two time points.");

	//	remember local solutions
		pSol = &vLocSol->solution(0);
		pOldSol = &vLocSol->solution(1);
		dt = vLocSol->time(0) - vLocSol->time(1);
	}

//	interpolate velocity at ip with standard lagrange interpolation
	static const size_t numSCVF = TFVGeom<TElem, dim>::numSCVF;
	MathVector<dim> StdVel[numSCVF];
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
		const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(ip);
		VecSet(StdVel[ip], 0.0);

		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			for(int d1 = 0; d1 < dim; ++d1)
				StdVel[ip][d1] += u(d1, sh) * scvf.shape(sh);
	}

//	compute stabilized velocities and shapes for continuity equation
	m_spStab->update(&geo, *pSol, StdVel, m_bStokes, m_imKinViscosity, pSource, pOldSol, dt);

	if (! m_bStokes) // no convective terms in the Stokes eq. => no upwinding
	{
	//	compute stabilized velocities and shapes for convection upwind
		if(m_spConvStab.valid())
			if(m_spConvStab != m_spStab)
				m_spConvStab->update(&geo, *pSol, StdVel, false, m_imKinViscosity, pSource, pOldSol, dt);
	
	//	compute upwind shapes
		if(m_spConvUpwind.valid())
			if(m_spStab->upwind() != m_spConvUpwind)
				m_spConvUpwind->update(&geo, StdVel);
	}

//	get a const (!!) reference to the stabilization
	const INavierStokesFV1Stabilization<dim>& stab = *m_spStab;
	const INavierStokesFV1Stabilization<dim>& convStab = *m_spConvStab;
	const INavierStokesFV1Upwind<dim>& upwind = *m_spConvUpwind;

// 	loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// 	get current SCVF
		const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(ip);

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
					w = peclet_blend_fv1(UpwindVel, geo, ip, StdVel[ip], m_imKinViscosity[ip]);

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
							const typename TFVGeom<TElem, dim>::SCVF& scvf2 = geo.scvf(ip2);
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
				if(m_bExactJacobian)
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
	
							prod_vel *= m_imDensitySCVF[ip];

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
	
							prod_p *= m_imDensitySCVF[ip];

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
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void NavierStokes<TDomain>::
ass_dA_elem_fv1(LocalVector& d, const LocalVector& u)
{
// 	Only first order implemented
	UG_ASSERT((TFVGeom<TElem, dim>::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	static const TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();

//	check for source term to pass to the stabilization
	const DataImport<MathVector<dim>, dim>* pSource = NULL;
	if(m_imSource.data_given())	pSource = &m_imSource;

//	check for solutions to pass to stabilization in time-dependent case
	const LocalVector *pSol = &u, *pOldSol = NULL;
	number dt = 0.0;
	if(this->is_time_dependent())
	{
	//	get and check current and old solution 
		const LocalVectorTimeSeries* vLocSol = this->local_time_solutions();
		if(vLocSol->size() != 2)
			UG_THROW("NavierStokes::ass_dA_elem: "
							" Stabilization needs exactly two time points.");

	//	remember local solutions
		pSol = &vLocSol->solution(0);
		pOldSol = &vLocSol->solution(1);
		dt = vLocSol->time(0) - vLocSol->time(1);
	}

//	interpolate velocity at ip with standard lagrange interpolation
	static const size_t numSCVF = TFVGeom<TElem, dim>::numSCVF;
	MathVector<dim> StdVel[numSCVF];
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
		const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(ip);
		VecSet(StdVel[ip], 0.0);

		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			for(int d1 = 0; d1 < dim; ++d1)
				StdVel[ip][d1] += u(d1, sh) * scvf.shape(sh);
	}

//	compute stabilized velocities and shapes for continuity equation
	// \todo: (optional) Here we can skip the computation of shapes, implement?
	m_spStab->update(&geo, *pSol, StdVel, m_bStokes, m_imKinViscosity, pSource, pOldSol, dt);

	if (! m_bStokes) // no convective terms in the Stokes eq. => no upwinding
	{
	//	compute stabilized velocities and shapes for convection upwind
		if(m_spConvStab.valid())
			if(m_spConvStab != m_spStab)
				m_spConvStab->update(&geo, *pSol, StdVel, false, m_imKinViscosity, pSource, pOldSol, dt);
	
	//	compute upwind shapes
		if(m_spConvUpwind.valid())
			if(m_spStab->upwind() != m_spConvUpwind)
				m_spConvUpwind->update(&geo, StdVel);
	}

//	get a const (!!) reference to the stabilization
	const INavierStokesFV1Stabilization<dim>& stab = *m_spStab;
	const INavierStokesFV1Stabilization<dim>& convStab = *m_spConvStab;
	const INavierStokesFV1Upwind<dim>& upwind = *m_spConvUpwind;

// 	loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// 	get current SCVF
		const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(ip);

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
				peclet_blend_fv1(UpwindVel, geo, ip, StdVel[ip], m_imKinViscosity[ip]);
	
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
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void NavierStokes<TDomain>::
ass_JM_elem_fv1(LocalMatrix& J, const LocalVector& u)
{
// 	Only first order implementation
	UG_ASSERT((TFVGeom<TElem, dim>::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	const static TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(ip);

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
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void NavierStokes<TDomain>::
ass_dM_elem_fv1(LocalVector& d, const LocalVector& u)
{
// 	Only first order implementation
	UG_ASSERT((TFVGeom<TElem, dim>::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	const static TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(ip);

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
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void NavierStokes<TDomain>::
ass_rhs_elem_fv1(LocalVector& d)
{
// 	Only first order implementation
	UG_ASSERT((TFVGeom<TElem, dim>::order == 1), "Only first order implemented.");

//	if zero data given, return
	if(!m_imSource.data_given()) return;

// 	get finite volume geometry
	const static TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int sh = scv.node_id();

	// 	Add to local rhs
		for(int d1 = 0; d1 < dim; ++d1)
			d(d1, sh) += m_imSource[ip][d1] * scv.volume();
	}
}

template<typename TDomain>
template<typename TFVGeom>
inline
number
NavierStokes<TDomain>::
peclet_blend_fv1(MathVector<dim>& UpwindVel, const TFVGeom& geo, size_t ip,
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

// register for 1D
template<typename TDomain>
void
NavierStokes<TDomain>::
register_all_fv1_funcs(bool bHang)
{
//	get all grid element types in this dimension and below
	typedef typename domain_traits<dim>::AllElemList ElemList;

//	switch assemble functions
	if(!bHang) boost::mpl::for_each<ElemList>( RegisterFV1<FV1Geometry>(this) );
	else throw(UGError("Not implemented."));
}

template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void
NavierStokes<TDomain>::
register_fv1_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->enable_fast_ass_elem(true);
	set_prep_elem_loop_fct(	id, &T::template prepare_element_loop_fv1<TElem, TFVGeom>);
	set_prep_elem_fct(	 	id, &T::template prepare_element_fv1<TElem, TFVGeom>);
	set_fsh_elem_loop_fct( 	id, &T::template finish_element_loop_fv1<TElem, TFVGeom>);
	set_ass_JA_elem_fct(	id, &T::template ass_JA_elem_fv1<TElem, TFVGeom>);
	set_ass_JM_elem_fct(	id, &T::template ass_JM_elem_fv1<TElem, TFVGeom>);
	set_ass_dA_elem_fct(	id, &T::template ass_dA_elem_fv1<TElem, TFVGeom>);
	set_ass_dM_elem_fct(	id, &T::template ass_dM_elem_fv1<TElem, TFVGeom>);
	set_ass_rhs_elem_fct(	id, &T::template ass_rhs_elem_fv1<TElem, TFVGeom>);
}

} // namespace NavierStokes
} // namespace ug