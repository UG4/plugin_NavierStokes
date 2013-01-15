/*
 * navier_stokes_fvho.h
 *
 *  Created on: 20.09.2010
 *      Author: andreasvogel
 */

#include "navier_stokes.h"

#include "common/util/provider.h"
#include "lib_disc/spatial_disc/disc_util/finite_volume_geometry.h"

namespace ug{
namespace NavierStokes{

template<typename TDomain>
template<typename TElem, typename VGeomProvider, typename PGeomProvider>
void NavierStokes<TDomain>::
prep_elem_loop_fvho()
{
//	check, that kinematic Viscosity has been set
	if(!m_imKinViscosity.data_given())
		UG_THROW("NavierStokes: Kinematic Viscosity has not been set, but is required.");

//	check, that Density has been set
	if(!m_imDensitySCVF.data_given())
		UG_THROW("NavierStokes: Density has not been set, but is required.");

//	check, that Density has been set
	if(!m_imDensitySCV.data_given())
		UG_THROW("NavierStokes: Density has not been set, but is required.");

//	request geometry
	typedef typename reference_element_traits<TElem>::reference_element_type reference_element_type;
	static const ReferenceObjectID roid = reference_element_type::REFERENCE_OBJECT_ID;

	static typename VGeomProvider::Type& vgeo = VGeomProvider::get();
	static typename PGeomProvider::Type& pgeo = PGeomProvider::get();
	try{
		vgeo.update_local(roid, m_order);
		pgeo.update_local(roid, m_order-1);
	}
	UG_CATCH_THROW("NavierStokes: Cannot update Finite Volume Geometry.");

//	set local positions for imports
	typedef typename reference_element_traits<TElem>::reference_element_type
																ref_elem_type;
	static const int refDim = ref_elem_type::dim;

	if(!VGeomProvider::Type::usesHangingNodes)
	{
		m_imKinViscosity.template set_local_ips<refDim>(vgeo.scvf_local_ips(),
		                                                vgeo.num_scvf_ips());
		m_imDensitySCVF.template set_local_ips<refDim>(vgeo.scvf_local_ips(),
		                                               vgeo.num_scvf_ips());
		m_imDensitySCV.template set_local_ips<refDim>(vgeo.scv_local_ips(),
		                                              vgeo.num_scv_ips());
		m_imSource.template set_local_ips<refDim>(vgeo.scv_local_ips(),
		                                          vgeo.num_scv_ips());
	}

	if(!PGeomProvider::Type::usesHangingNodes)
	{
		m_imDensitySCVFp.template set_local_ips<refDim>(pgeo.scvf_local_ips(),
		                                               pgeo.num_scvf_ips());

		const LocalShapeFunctionSet<dim>& rVTrialSpace =
			LocalShapeFunctionSetProvider::get<dim>(roid, LFEID(LFEID::LAGRANGE, m_order));
		const MathVector<dim>* PLocIP = pgeo.scv_local_ips();

		m_vvVShape.resize(pgeo.num_scvf_ips());
		for(size_t ip = 0; ip < m_vvVShape.size(); ++ip){
			m_vvVShape[ip].resize(rVTrialSpace.num_sh());
			for(size_t sh = 0; sh < rVTrialSpace.num_sh(); ++sh){
				m_vvVShape[ip][sh] = rVTrialSpace.shape(sh, PLocIP[ip]);
			}
		}

		const LocalShapeFunctionSet<dim>& rPTrialSpace =
			LocalShapeFunctionSetProvider::get<dim>(roid, LFEID(LFEID::LAGRANGE, m_order-1));
		const MathVector<dim>* VLocIP = vgeo.scv_local_ips();

		m_vvPShape.resize(vgeo.num_scvf_ips());
		for(size_t ip = 0; ip < m_vvPShape.size(); ++ip){
			m_vvPShape[ip].resize(rPTrialSpace.num_sh());
			for(size_t sh = 0; sh < rPTrialSpace.num_sh(); ++sh){
				m_vvPShape[ip][sh] = rPTrialSpace.shape(sh, VLocIP[ip]);
			}
		}
	}
}

template<typename TDomain>
template<typename TElem, typename VGeomProvider, typename PGeomProvider>
void NavierStokes<TDomain>::
fsh_elem_loop_fvho()
{}


template<typename TDomain>
template<typename TElem, typename VGeomProvider, typename PGeomProvider>
void NavierStokes<TDomain>::
prep_elem_fvho(TElem* elem, const LocalVector& u)
{
// 	Update Geometry for this element
	static typename VGeomProvider::Type& vgeo = VGeomProvider::get();
	static typename PGeomProvider::Type& pgeo = PGeomProvider::get();
	try{
		vgeo.update(elem, this->template element_corners<TElem>(elem),
	               &(this->subset_handler()));
	    pgeo.update(elem, this->template element_corners<TElem>(elem),
	               &(this->subset_handler()));
	}
	UG_CATCH_THROW("NavierStokes: Cannot update Finite Volume Geometry.");

//	set local positions for imports
	if(VGeomProvider::Type::usesHangingNodes)
	{
	//	set local positions for imports
		typedef typename reference_element_traits<TElem>::reference_element_type
																	ref_elem_type;
		static const int refDim = ref_elem_type::dim;

	//	request ip series
		m_imKinViscosity.template set_local_ips<refDim>(vgeo.scvf_local_ips(),
		                                                vgeo.num_scvf_ips());
		m_imDensitySCVF.template set_local_ips<refDim>(vgeo.scvf_local_ips(),
		                                                vgeo.num_scvf_ips());
		m_imDensitySCV.template set_local_ips<refDim>(vgeo.scv_local_ips(),
		                                          vgeo.num_scv_ips());
		m_imSource.template set_local_ips<refDim>(vgeo.scv_local_ips(),
		                                          vgeo.num_scv_ips());
	}

//	set local positions for imports
	if(PGeomProvider::Type::usesHangingNodes)
	{
	//	set local positions for imports
		typedef typename reference_element_traits<TElem>::reference_element_type
																	ref_elem_type;
		static const int refDim = ref_elem_type::dim;

	//	request ip series
		m_imDensitySCVFp.template set_local_ips<refDim>(pgeo.scvf_local_ips(),
														pgeo.num_scvf_ips());
	}

//	set global positions for imports
	m_imKinViscosity.set_global_ips(vgeo.scvf_global_ips(), vgeo.num_scvf_ips());
	m_imDensitySCVF.set_global_ips(vgeo.scvf_global_ips(), vgeo.num_scvf_ips());
	m_imDensitySCV.set_global_ips(vgeo.scv_global_ips(), vgeo.num_scv_ips());
	m_imSource.set_global_ips(vgeo.scv_global_ips(), vgeo.num_scv_ips());

	m_imDensitySCVFp.set_global_ips(pgeo.scvf_global_ips(), pgeo.num_scvf_ips());
}

template<typename TDomain>
template<typename TElem, typename VGeomProvider, typename PGeomProvider>
void NavierStokes<TDomain>::
add_jac_A_elem_fvho(LocalMatrix& J, const LocalVector& u)
{
//	request geometry
	static const typename VGeomProvider::Type& vgeo = VGeomProvider::get();
	static const typename PGeomProvider::Type& pgeo = PGeomProvider::get();

// 	loop Sub Control Volume Faces (SCVF)
	for(size_t i = 0, ipCnt = 0; i < vgeo.num_scvf(); ++i)
	{
	// 	get current SCVF
		const typename VGeomProvider::Type::SCVF& scvf = vgeo.scvf(i);

	//	loop integration points
		for(size_t ip = 0; ip < scvf.num_ip(); ++ip)
		{
			////////////////////////////////////////////////////
			////////////////////////////////////////////////////
			// Momentum Equation (conservation of momentum)
			////////////////////////////////////////////////////
			////////////////////////////////////////////////////

			////////////////////////////////////////////////////
			// Diffusive Term (Momentum Equation)
			////////////////////////////////////////////////////
		// 	loop shape functions
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			{

			// 	Compute flux derivative at IP
				const number flux_sh =  -1.0 * m_imKinViscosity[ipCnt]
				                        * m_imDensitySCVF[ipCnt]
										* VecDot(scvf.global_grad(ip, sh), scvf.normal())
										* scvf.weight(ip);

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
							const number flux2_sh = -1.0 * m_imKinViscosity[ipCnt]
							                        * m_imDensitySCVF[ipCnt]
													* scvf.global_grad(ip, sh)[d1]
													* scvf.normal()[d2]
									                * scvf.weight(ip);
							J(d1, scvf.from(), d2, sh) += flux2_sh;
							J(d1, scvf.to()  , d2, sh) -= flux2_sh;
						}
				}
			}

			////////////////////////////////////////////////////
			// Pressure Term (Momentum Equation)
			////////////////////////////////////////////////////

		//	Add flux derivative for local matrix
			for(size_t sh = 0; sh < m_vvPShape[ipCnt].size(); ++sh)
				for(int d1 = 0; d1 < dim; ++d1)
				{
					const number flux_sh = m_vvPShape[ipCnt][sh]
					                       * scvf.normal()[d1]
	                                       * scvf.weight(ip);
					J(d1, scvf.from(), _P_, sh) += flux_sh;
					J(d1, scvf.to()  , _P_, sh) -= flux_sh;
				}

			if(!m_bStokes)
				UG_THROW("Only Stokes implemented");

			ipCnt++;
		} // end ip
	} // end scvf


	////////////////////////////////////////////////////
	////////////////////////////////////////////////////
	// Continuity Equation (conservation of mass)
	////////////////////////////////////////////////////
	////////////////////////////////////////////////////

// 	loop Sub Control Volume Faces (SCVF)
	for(size_t i = 0, ipCnt = 0; i < pgeo.num_scvf(); ++i)
	{
// 	get current SCVF
	const typename PGeomProvider::Type::SCVF& scvf = pgeo.scvf(i);

//	loop integration points
	for(size_t ip = 0; ip < scvf.num_ip(); ++ip)
	{
		for(size_t sh = 0; sh < m_vvVShape[ipCnt].size(); ++sh)
		{
		//	compute flux at ip
			const number contFlux = m_vvVShape[ipCnt][sh]
			                        * m_imDensitySCVFp[ipCnt]
			                        * scvf.weight(ip);

		//	Add contributions to local defect
			for(int d1 = 0; d1 < dim; ++d1){
				J(_P_, scvf.from(), d1, sh) += contFlux * scvf.normal()[d1];
				J(_P_, scvf.to()  , d1, sh) -= contFlux * scvf.normal()[d1];
			}
		}

		ipCnt++;
	}
	}
}

template<typename TDomain>
template<typename TElem, typename VGeomProvider, typename PGeomProvider>
void NavierStokes<TDomain>::
add_def_A_elem_fvho(LocalVector& d, const LocalVector& u)
{
//	request geometry
	static const typename VGeomProvider::Type& vgeo = VGeomProvider::get();
	static const typename PGeomProvider::Type& pgeo = PGeomProvider::get();

// 	loop Sub Control Volume Faces (SCVF)
	for(size_t i = 0, ipCnt = 0; i < vgeo.num_scvf(); ++i)
	{
// 	get current SCVF
	const typename VGeomProvider::Type::SCVF& scvf = vgeo.scvf(i);

//	loop integration points
	for(size_t ip = 0; ip < scvf.num_ip(); ++ip)
	{

		////////////////////////////////////////////////////
		////////////////////////////////////////////////////
		// Momentum Equation (conservation of momentum)
		////////////////////////////////////////////////////
		////////////////////////////////////////////////////

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
					gradVel(d1, d2) += scvf.global_grad(ip, sh)[d2]
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
		VecScale(diffFlux, diffFlux, (-1.0) * m_imKinViscosity[ipCnt] * m_imDensitySCVF[ipCnt]);

	//	3. Add flux to local defect
		for(int d1 = 0; d1 < dim; ++d1)
		{
			d(d1, scvf.from()) += diffFlux[d1] * scvf.weight(ip);
			d(d1, scvf.to()  ) -= diffFlux[d1] * scvf.weight(ip);
		}

		////////////////////////////////////////////////////
		// Convective Term (Momentum Equation)
		////////////////////////////////////////////////////

		if (! m_bStokes) // no convective terms in the Stokes equation
			UG_THROW("Only Stokes implemented");

		////////////////////////////////////////////////////
		// Pressure Term (Momentum Equation)
		////////////////////////////////////////////////////

	//	1. Interpolate pressure at ip
		number pressure = 0.0;
		for(size_t sh = 0; sh < m_vvPShape[ipCnt].size(); ++sh)
			pressure += m_vvPShape[ipCnt][sh] * u(_P_, sh);

	//	2. Add contributions to local defect
		for(int d1 = 0; d1 < dim; ++d1)
		{
			d(d1, scvf.from()) += pressure * scvf.normal()[d1] * scvf.weight(ip);
			d(d1, scvf.to()  ) -= pressure * scvf.normal()[d1] * scvf.weight(ip);
		}
		ipCnt++;
	}
	}

//	interpolate velocity at ip with standard lagrange interpolation
	std::vector<MathVector<dim> > PStdVel(m_vvVShape.size());
	for(size_t ip = 0; ip < m_vvVShape.size(); ++ip){
		VecSet(PStdVel[ip], 0.0);
		for(size_t sh = 0; sh < m_vvVShape[ip].size(); ++sh)
			for(int d1 = 0; d1 < dim; ++d1)
				PStdVel[ip][d1] += u(d1, sh) * m_vvVShape[ip][sh];
	}

// 	loop Sub Control Volume Faces (SCVF)
	for(size_t i = 0, ipCnt = 0; i < pgeo.num_scvf(); ++i)
	{
// 	get current SCVF
	const typename PGeomProvider::Type::SCVF& scvf = pgeo.scvf(i);

//	loop integration points
	for(size_t ip = 0; ip < scvf.num_ip(); ++ip)
	{
		////////////////////////////////////////////////////
		////////////////////////////////////////////////////
		// Continuity Equation (conservation of mass)
		////////////////////////////////////////////////////
		////////////////////////////////////////////////////

	//	compute flux at ip
		const number contFlux = VecProd(PStdVel[ipCnt], scvf.normal()) * m_imDensitySCVFp[ipCnt];

	//	Add contributions to local defect
		d(_P_, scvf.from()) += contFlux * scvf.weight(ip);
		d(_P_, scvf.to()  ) -= contFlux * scvf.weight(ip);

		ipCnt++;
	}
	}
}


template<typename TDomain>
template<typename TElem, typename VGeomProvider, typename PGeomProvider>
void NavierStokes<TDomain>::
add_jac_M_elem_fvho(LocalMatrix& J, const LocalVector& u)
{
//	request geometry
	static const typename VGeomProvider::Type& vgeo = VGeomProvider::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0, ipOffset = 0; ip < vgeo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename VGeomProvider::Type::SCV& scv = vgeo.scv(ip);

	//	loop shapes
		for(size_t sh = 0; sh < scv.num_sh(); ++sh)
		{
		//	reset integral
			number integral = 0;

		//	loop integration points
			for(size_t ip = 0; ip < scv.num_ip(); ++ip)
			{
				integral += scv.shape(ip, sh) * scv.weight(ip) * m_imDensitySCV[ipOffset+ip];
			}

		// 	loop velocity components
			for(int d1 = 0; d1 < dim; ++d1)
			{
			// 	Add to local matrix
				J(d1, sh, d1, sh) += integral * scv.volume();
			}
		}

	//	increase ip offset
		ipOffset += scv.num_ip();
	}
}


template<typename TDomain>
template<typename TElem, typename VGeomProvider, typename PGeomProvider>
void NavierStokes<TDomain>::
add_def_M_elem_fvho(LocalVector& d, const LocalVector& u)
{
//	request geometry
	static const typename VGeomProvider::Type& vgeo = VGeomProvider::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t i = 0, ipCnt = 0; i < vgeo.num_scv(); ++i)
	{
	// 	get current SCV
		const typename VGeomProvider::Type::SCV& scv = vgeo.scv(i);

	// 	get associated node
		const int co = scv.node_id();

	//	loop integration points
		for(size_t ip = 0; ip < scv.num_ip(); ++ip)
		{
			MathVector<dim> Vel; VecSet(Vel, 0.0);
			for(size_t sh = 0; sh < scv.num_sh(); ++sh)
				for(int d1 = 0; d1 < dim; ++d1)
					Vel[d1] += u(d1, sh) * scv.shape(ip, sh);

		// 	Add to local defect
			for(int d1 = 0; d1 < dim; ++d1)
				d(d1, co) +=  Vel[d1] * m_imDensitySCV[ipCnt]
                              * scv.weight(ip) * scv.volume();

		//	next ip
			ipCnt++;
		}
	}
}


template<typename TDomain>
template<typename TElem, typename VGeomProvider, typename PGeomProvider>
void NavierStokes<TDomain>::
add_rhs_elem_fvho(LocalVector& d)
{
//	if zero data given, return
	if(!m_imSource.data_given()) return;

	UG_THROW("Not implemented.")
}

////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

// register for all dim
template<>
void NavierStokes<Domain1d>::
register_all_fvho_funcs(int order)
{
	UG_THROW("Not implemented.");
}

// register for all dim
template<>
void NavierStokes<Domain2d>::
register_all_fvho_funcs(int order)
{
	typedef DimFVGeometry<2, dim> FVGeom;
	register_fvho_func<Triangle, VFlexGeomProvider<FVGeom>, PFlexGeomProvider<FVGeom> >();
	register_fvho_func<Quadrilateral, VFlexGeomProvider<FVGeom>, PFlexGeomProvider<FVGeom> >();
}

// register for all dim
template<>
void NavierStokes<Domain3d>::
register_all_fvho_funcs(int order)
{
	typedef DimFVGeometry<3, dim> FVGeom;
	register_fvho_func<Tetrahedron, VFlexGeomProvider<FVGeom>, PFlexGeomProvider<FVGeom> >();
	register_fvho_func<Prism, VFlexGeomProvider<FVGeom>, PFlexGeomProvider<FVGeom> >();
	register_fvho_func<Hexahedron, VFlexGeomProvider<FVGeom>, PFlexGeomProvider<FVGeom> >();
}

template<typename TDomain>
template<typename TElem, typename VGeomProvider, typename PGeomProvider>
void NavierStokes<TDomain>::register_fvho_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->enable_fast_add_elem(true);
	set_prep_elem_loop_fct(	id, &T::template prep_elem_loop_fvho<TElem, VGeomProvider, PGeomProvider>);
	set_prep_elem_fct(	 	id, &T::template prep_elem_fvho<TElem, VGeomProvider, PGeomProvider>);
	set_fsh_elem_loop_fct( 	id, &T::template fsh_elem_loop_fvho<TElem, VGeomProvider, PGeomProvider>);
	set_add_jac_A_elem_fct(	id, &T::template add_jac_A_elem_fvho<TElem, VGeomProvider, PGeomProvider>);
	set_add_jac_M_elem_fct(	id, &T::template add_jac_M_elem_fvho<TElem, VGeomProvider, PGeomProvider>);
	set_add_def_A_elem_fct(	id, &T::template add_def_A_elem_fvho<TElem, VGeomProvider, PGeomProvider>);
	set_add_def_M_elem_fct(	id, &T::template add_def_M_elem_fvho<TElem, VGeomProvider, PGeomProvider>);
	set_add_rhs_elem_fct(	id, &T::template add_rhs_elem_fvho<TElem, VGeomProvider, PGeomProvider>);
}

} // namespace NavierStokes
} // namespace ug
