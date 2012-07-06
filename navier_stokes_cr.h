/*
 * navier_stokes_cr.cpp
 *
 *  Created on: 04.07.2012
 *      Author: Christian Wehner
 */

#include "navier_stokes.h"

#include "common/util/provider.h"
#include "lib_disc/spatial_disc/disc_util/cr_finite_volume_geometry.h"
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
prepare_element_loop_cr()
{
// 	Only first order implementation
	if(!(TFVGeom<TElem, dim>::order == 1))
		UG_THROW("Only first order implementation, but other Finite Volume"
						" Geometry set.");

	if (! m_bStokes) // no convective terms in the Stokes eq. => no upwinding
	{
	//	check, that convective upwinding has been set
		if(m_spConvCRUpwind.invalid()){
			UG_THROW("Upwinding for convective Term in Momentum eq. not set.");
		}else
			m_spConvCRUpwind->template set_geometry_type<TFVGeom<TElem, dim> >();
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
finish_element_loop_cr()
{}


template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void NavierStokes<TDomain>::
prepare_element_cr(TElem* elem, const LocalVector& u)
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
ass_JA_elem_cr(LocalMatrix& J, const LocalVector& u)
{
	// 	Only first order implementation
		UG_ASSERT((TFVGeom<TElem, dim>::order == 1), "Only first order implemented.");

	// 	get finite volume geometry
		static const TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();

	/*	check for source term to pass to the stabilization
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
		}*/

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
	//	m_spStab->update(&geo, *pSol, StdVel, m_bStokes, m_imKinViscosity, pSource, pOldSol, dt);
		const INavierStokesCRUpwind<dim>& upwind = *m_spConvCRUpwind;

		if (! m_bStokes) // no convective terms in the Stokes eq. => no upwinding
		{
		//	compute upwind shapes
			m_spConvCRUpwind->update(&geo, StdVel);
		}

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
				// Convective Term (Momentum Equation)
				////////////////////////////////////////////////////

				if (! m_bStokes) // no convective terms in the Stokes equation
				{
				//	compute upwind velocity
					MathVector<dim> UpwindVel;

				//	switch PAC
					UpwindVel = upwind.upwind_vel(ip, u, StdVel);

				//	peclet blend
					number w = 1.0;
				//	if(m_bPecletBlend)
				//		w = peclet_blend(UpwindVel, scvf, StdVel[ip], m_imKinViscosity[ip]);

				//	compute product of stabilized vel and normal todo which is better upwindVel or StdVel?
					const number prod = VecProd(StdVel[ip], scvf.normal()) * m_imDensitySCVF[ip];

				///////////////////////////////////
				//	Add fixpoint linearization
				///////////////////////////////////

					number convFlux_vel = upwind.upwind_shape_sh(ip, sh);

				//	in some cases (e.g. PositiveUpwind, RegularUpwind) the upwind
				//	velocity in an ip depends also on the upwind velocity in
				//	other ips. This is reflected by the fact, that the ip
				//	shapes are non-zero. In that case, we can interpolate an
				//	approximate upwind only from the corner velocities by using
				//	u_up = \sum shape_co U_co + \sum shape_ip \tilde{u}_ip
				//	     = \sum shape_co U_co + \sum \sum shape_ip norm_shape_co|_ip * U_co
					if(upwind.non_zero_shape_ip())
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

				//	derivative due to peclet blending
				//	if(m_bPecletBlend)
				//	{
				//		const number convFluxPe = prod * (1.0-w) * scvf.shape(sh);
				//		for(int d1 = 0; d1 < dim; ++d1)
				//		{
				//			J(d1, scvf.from(), d1, sh) += convFluxPe;
				//			J(d1, scvf.to()  , d1, sh) -= convFluxPe;
				//		}
				//	}

				/////////////////////////////////////////
				//	Add full jacobian (remaining part)
				/////////////////////////////////////////

				//	Add remaining term for exact jacobian
					if(m_bExactJacobian)
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
			}// end of loop shape functions

			////////////////////////////////////////////////////
			// Pressure Term (Momentum Equation)
			////////////////////////////////////////////////////

			//	Add flux derivative for local matrix
			for(int d1 = 0; d1 < dim; ++d1)
			{
				// pressure is constant over element
				J(d1, scvf.from(), _P_, 0) += scvf.normal()[d1];
				J(d1, scvf.to()  , _P_, 0) -= scvf.normal()[d1];
			}
		}// end of loop ips

		////////////////////////////////////////////////////
		////////////////////////////////////////////////////
		// Continuity Equation (conservation of mass)
		////////////////////////////////////////////////////
		////////////////////////////////////////////////////
		for(size_t sh = 0; sh < geo.num_scv(); ++sh){
			const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(sh);
			for(int d1 = 0; d1 < dim; ++d1)
			{
				J(_P_, 0 , d1, sh) += scv.normal()[d1];
			}
		}
};

template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void NavierStokes<TDomain>::
ass_dA_elem_cr(LocalVector& d, const LocalVector& u)
{
	// 	Only first order implemented
		UG_ASSERT((TFVGeom<TElem, dim>::order == 1), "Only first order implemented.");

	// 	get finite volume geometry
		static const TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();

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
		//	m_spStab->update(&geo, *pSol, StdVel, m_bStokes, m_imKinViscosity, pSource, pOldSol, dt);
		const INavierStokesCRUpwind<dim>& upwind = *m_spConvCRUpwind;

		if (! m_bStokes) // no convective terms in the Stokes eq. => no upwinding
		{
			//	compute upwind shapes
			m_spConvCRUpwind->update(&geo, StdVel);
		}

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
				UpwindVel = upwind.upwind_vel(ip, u, StdVel);

			//	Peclet Blend
			//	if(m_bPecletBlend)
			//		peclet_blend(UpwindVel, scvf, StdVel[ip], m_imKinViscosity[ip]);

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
			number pressure = u(_P_, 0);
			for(int d1 = 0; d1 < dim; ++d1)
			{
				d(d1, scvf.from()) += pressure * scvf.normal()[d1];
				d(d1, scvf.to()  ) -= pressure * scvf.normal()[d1];
			}
		}
		////////////////////////////////////////////////////
		////////////////////////////////////////////////////
		// Continuity Equation (conservation of mass)
		////////////////////////////////////////////////////
		////////////////////////////////////////////////////
		for(size_t sh = 0; sh < geo.num_scv(); ++sh){
			const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(sh);
			for(int d1 = 0; d1 < dim; ++d1)
			{
				d(_P_, 0 ) += scv.normal()[d1] * u(d1,sh);
			}
		}
}


template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void NavierStokes<TDomain>::
ass_JM_elem_cr(LocalMatrix& J, const LocalVector& u)
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
ass_dM_elem_cr(LocalVector& d, const LocalVector& u)
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
ass_rhs_elem_cr(LocalVector& d)
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


////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

// register for 1D
template<typename TDomain>
void
NavierStokes<TDomain>::
register_all_cr_funcs(bool bHang)
{
//	get all grid element types in this dimension and below
	typedef typename domain_traits<dim>::AllElemList ElemList;

//	switch assemble functions
	if(!bHang) boost::mpl::for_each<ElemList>( RegisterCR<CRFVGeometry>(this) );
	else throw(UGError("Not implemented."));
}

template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void
NavierStokes<TDomain>::
register_cr_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->enable_fast_ass_elem(true);
	set_prep_elem_loop_fct(	id, &T::template prepare_element_loop_cr<TElem, TFVGeom>);
	set_prep_elem_fct(	 	id, &T::template prepare_element_cr<TElem, TFVGeom>);
	set_fsh_elem_loop_fct( 	id, &T::template finish_element_loop_cr<TElem, TFVGeom>);
	set_ass_JA_elem_fct(	id, &T::template ass_JA_elem_cr<TElem, TFVGeom>);
	set_ass_JM_elem_fct(	id, &T::template ass_JM_elem_cr<TElem, TFVGeom>);
	set_ass_dA_elem_fct(	id, &T::template ass_dA_elem_cr<TElem, TFVGeom>);
	set_ass_dM_elem_fct(	id, &T::template ass_dM_elem_cr<TElem, TFVGeom>);
	set_ass_rhs_elem_fct(	id, &T::template ass_rhs_elem_cr<TElem, TFVGeom>);
}

} // namespace NavierStokes
} // namespace ug
