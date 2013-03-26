/*
 * navier_stokes_fvcr.cpp
 *
 *  Created on: 04.07.2012
 *      Author: Christian Wehner
 */

#include "navier_stokes_fvcr.h"

#include "lib_disc/spatial_disc/disc_util/fvcr_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"

namespace ug{
namespace NavierStokes{

////////////////////////////////////////////////////////////////////////////////
//	Constructor - set default values
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
NavierStokesFVCR<TDomain>::NavierStokesFVCR(const char* functions,
                                          const char* subsets)
: NavierStokesBase<TDomain>(functions, subsets)
{
	init();
};

template<typename TDomain>
NavierStokesFVCR<TDomain>::NavierStokesFVCR(const std::vector<std::string>& vFct,
                                          const std::vector<std::string>& vSubset)
: NavierStokesBase<TDomain>(vFct, vSubset)
{
	init();
};


template<typename TDomain>
void NavierStokesFVCR<TDomain>::init()
{
//	check number of functions
	if(this->num_fct() != dim+1)
		UG_THROW("Wrong number of functions: The ElemDisc 'NavierStokes'"
					   " needs exactly "<<dim+1<<" symbolic function.");
//	register imports
	this->register_import(m_imSource);
	this->register_import(m_imKinViscosity);
	this->register_import(m_imDensitySCVF);
	this->register_import(m_imDensitySCV);

	m_imSource.set_rhs_part();
	m_imDensitySCV.set_mass_part();

	//	default value for density
	base_type::set_density(1.0);

	m_bDefectUpwind = true;

	//	update assemble functions
	register_all_funcs(false);
}

template<typename TDomain>
bool NavierStokesFVCR<TDomain>::request_non_regular_grid(bool bNonRegular)
{
//	switch, which assemble functions to use.
	if(bNonRegular)
	{
		UG_LOG("ERROR in 'NavierStokes::request_non_regular_grid':"
				" Non-regular grid not implemented.\n");
		return false;
	}

//	this disc supports regular grids
	return true;
}

template<typename TDomain>
bool NavierStokesFVCR<TDomain>::
request_finite_element_id(const std::vector<LFEID>& vLfeID)
{
//	check number
	if(vLfeID.size() != dim+1)
	{
		UG_LOG("NavierStokes:"
				" Wrong number of functions given. Need exactly "<<dim+1<<"\n");
		return false;
	}

	for(int d = 0; d < dim; ++d)
		if(vLfeID[d].type() != LFEID::CROUZEIX_RAVIART)
		{
			UG_LOG("NavierStokes: 'fvcr' expects Crouzeix-Raviart trial"
					" space for velocity.\n");
			return false;
		}
	if(vLfeID[dim].type() != LFEID::PIECEWISE_CONSTANT)
	{
		UG_LOG("NavierStokes: 'fvcr' expects piecewise constant trial"
				" space for pressure.\n");
		return false;
	}

	//	update assemble functions
	register_all_funcs(false);

	//	is supported
	return true;
}

template<typename TDomain>
void NavierStokesFVCR<TDomain>::
set_kinematic_viscosity(SmartPtr<UserData<number, dim> > data)
{
	m_imKinViscosity.set_data(data);
}

template<typename TDomain>
void NavierStokesFVCR<TDomain>::
set_density(SmartPtr<UserData<number, dim> > data)
{
	m_imDensitySCVF.set_data(data);
	m_imDensitySCV.set_data(data);
}

template<typename TDomain>
void NavierStokesFVCR<TDomain>::
set_source(SmartPtr<UserData<MathVector<dim>, dim> > data)
{
	m_imSource.set_data(data);
}

////////////////////////////////////////////////////////////////////////////////
//	assembling functions
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesFVCR<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
// 	Only first order implementation
	if(!(TFVGeom::order == 1))
		UG_THROW("Only first order implementation, but other Finite Volume"
						" Geometry set.");

	if (! m_bStokes) // no convective terms in the Stokes eq. => no upwinding
	{
	//	check, that convective upwinding has been set
		if(m_spConvUpwind.invalid()){
			UG_THROW("Upwinding for convective Term in Momentum eq. not set.");
		}else
			m_spConvUpwind->template set_geometry_type<TFVGeom >();
	}

//	check, that kinematic Viscosity has been set
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

//	set local positions for imports
	typedef typename reference_element_traits<TElem>::reference_element_type
																ref_elem_type;
	static const int refDim = ref_elem_type::dim;

	if(!TFVGeom::usesHangingNodes)
	{
		static TFVGeom& geo = GeomProvider<TFVGeom>::get();
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
template<typename TElem, typename TFVGeom>
void NavierStokesFVCR<TDomain>::
fsh_elem_loop()
{}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesFVCR<TDomain>::
prep_elem(TElem* elem, LocalVector& u)
{
//	get corners
	m_vCornerCoords = this->template element_corners<TElem>(elem);

// 	Update Geometry for this element
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();
	try{
		geo.update(elem, &m_vCornerCoords[0], &(this->subset_handler()));
	}
	UG_CATCH_THROW("NavierStokes::prep_elem:"
					" Cannot update Finite Volume Geometry.");

//	set local positions for imports
	if(TFVGeom::usesHangingNodes)
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
template<typename TFVGeom>
inline
number
NavierStokesFVCR<TDomain>::
peclet_blend(MathVector<dim>& UpwindVel, const TFVGeom& geo, size_t ip,
             const MathVector<dim>& StdVel, number kinVisco)
{
	const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
//	compute peclet number
	number Pe = VecProd(StdVel, scvf.normal())/VecTwoNormSq(scvf.normal())
	 * VecDistance( geo.scv(scvf.to()).global_ip(), geo.scv(scvf.from()).global_ip() ) / kinVisco;

//	compute weight
	const number Pe2 = Pe * Pe;
	number w = Pe2 / (5.0 + Pe2);

//	compute upwind vel
	VecScaleAdd(UpwindVel, w, UpwindVel, (1.0-w), StdVel);

	return w;
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesFVCR<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u)
{
	// 	Only first order implementation
		UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

	// 	get finite volume geometry
		static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

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

		const INavierStokesUpwind<dim>& upwind = *m_spConvUpwind;

		if (! m_bStokes) // no convective terms in the Stokes eq. => no upwinding
		{
		//	compute upwind shapes
			m_spConvUpwind->update(&geo, StdVel);
		}

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
//				UG_LOG("---------------------------------------------------\n");
				if(!m_bLaplace)
				{
					for(int d1 = 0; d1 < dim; ++d1)
						for(int d2 = 0; d2 < dim; ++d2)
						{
							const number flux2_sh = -1.0 * m_imKinViscosity[ip] * m_imDensitySCVF[ip]
													* scvf.global_grad(sh)[d1] * scvf.normal()[d2];
//	UG_LOG("d1 = " << d1  << " scvf.global_grad(sh)[d1] = " << scvf.global_grad(sh)[d1] << " d2 = " << d2 << " scvf.normal()[d2] = " << scvf.normal()[d2] << "\n");
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
					if(m_bPecletBlend)
						w = peclet_blend(UpwindVel, geo, ip, StdVel[ip], m_imKinViscosity[ip]);

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
			const typename TFVGeom::SCV& scv = geo.scv(sh);
			for(int d1 = 0; d1 < dim; ++d1)
			{
				J(_P_, 0 , d1, sh) += scv.normal()[d1];
			}
		}
};

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesFVCR<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u)
{
	// 	Only first order implemented
		UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

	// 	get finite volume geometry
		static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

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

		const INavierStokesUpwind<dim>& upwind = *m_spConvUpwind;

		if ((! m_bStokes) && (m_bDefectUpwind == true))
		{
			//	compute upwind shapes
			m_spConvUpwind->update(&geo, StdVel);
		}

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

			if (! m_bStokes) // no convective terms in the Stokes equation
			{
				if (m_bDefectUpwind == true){
				//	find the upwind velocity at ip
					MathVector<dim> UpwindVel;

				//	switch PAC
					UpwindVel = upwind.upwind_vel(ip, u, StdVel);

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
				} else {
					for(int d1 = 0; d1 < dim; ++d1)
					{
						//	compute product of standard velocity and normal
						const number prod = VecProd(StdVel[ip], scvf.normal()) * m_imDensitySCVF[ip];
						d(d1, scvf.from()) += StdVel[ip][d1] * prod;
						d(d1, scvf.to()  ) -= StdVel[ip][d1] * prod;
					}
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
			const typename TFVGeom::SCV& scv = geo.scv(sh);
			for(int d1 = 0; d1 < dim; ++d1)
			{
				d(_P_, 0 ) += scv.normal()[d1] * u(d1,sh);
			}
		}
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesFVCR<TDomain>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u)
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
void NavierStokesFVCR<TDomain>::
add_def_M_elem(LocalVector& d, const LocalVector& u)
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
void NavierStokesFVCR<TDomain>::
add_rhs_elem(LocalVector& d)
{
// 	Only first order implementation
	UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

//	if zero data given, return
	if(!m_imSource.data_given()) return;

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
		for(int d1 = 0; d1 < dim; ++d1)
			d(d1, sh) += m_imSource[ip][d1] * scv.volume();
	}
}

////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

template<>
void NavierStokesFVCR<Domain1d>::
register_all_funcs(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_func<Edge, CRFVGeometry<Edge, dim> >();
	}
	else
	{
		UG_THROW("ConvectionDiffusion"
						" Hanging nodes not supported for CRFV discretization.");
	}
}

template<>
void NavierStokesFVCR<Domain2d>::
register_all_funcs(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_func<Triangle, CRFVGeometry<Triangle, dim> >();
		register_func<Quadrilateral, CRFVGeometry<Quadrilateral, dim> >();
	}
	else
	{
		UG_THROW("ConvectionDiffusion"
						" Hanging nodes not supported for CRFV discretization.");
	}
}

template<>
void NavierStokesFVCR<Domain3d>::
register_all_funcs(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_func<Tetrahedron, CRFVGeometry<Tetrahedron, dim> >();
		register_func<Prism, CRFVGeometry<Prism, dim> >();
		register_func<Pyramid, CRFVGeometry<Pyramid, dim> >();
		register_func<Hexahedron, CRFVGeometry<Hexahedron, dim> >();
	}
	else
	{
		UG_THROW("ConvectionDiffusion"
						" Hanging nodes not supported for CRFV discretization.");
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void
NavierStokesFVCR<TDomain>::
register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->enable_fast_add_elem(true);
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
template class NavierStokesFVCR<Domain2d>;
#endif
#ifdef UG_DIM_3
template class NavierStokesFVCR<Domain3d>;
#endif

} // namespace NavierStokes
} // namespace ug
