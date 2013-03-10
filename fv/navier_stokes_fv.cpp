/*
 * navier_stokes_fv.cpp
 *
 *  Created on: 20.09.2010
 *      Author: andreasvogel
 */

#include "navier_stokes_fv.h"

#include "common/util/provider.h"
#include "lib_disc/spatial_disc/disc_util/fvho_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"

namespace ug{
namespace NavierStokes{

////////////////////////////////////////////////////////////////////////////////
//	Constructor - set default values
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
NavierStokesFV<TDomain>::NavierStokesFV(const char* functions,
                                          const char* subsets)
: NavierStokesBase<TDomain>(functions, subsets)
{
	init();
};

template<typename TDomain>
NavierStokesFV<TDomain>::NavierStokesFV(const std::vector<std::string>& vFct,
                                          const std::vector<std::string>& vSubset)
: NavierStokesBase<TDomain>(vFct, vSubset)
{
	init();
};


template<typename TDomain>
void NavierStokesFV<TDomain>::init()
{
//	check number of functions
	if(this->num_fct() != dim+1)
		UG_THROW("Wrong number of functions: The ElemDisc 'NavierStokes'"
					   " needs exactly "<<dim+1<<" symbolic function.");
//	register imports
	this->register_import(m_imSource);
	this->register_import(m_imKinViscosity);
	this->register_import(m_imDensitySCVF);
	this->register_import(m_imDensitySCVFp);
	this->register_import(m_imDensitySCV);

	m_imSource.set_rhs_part();
	m_imDensitySCV.set_mass_part();

	//	default value for density
	base_type::set_density(1.0);

	//	update assemble functions
	this->enable_fast_add_elem(true);
}

template<typename TDomain>
bool NavierStokesFV<TDomain>::request_non_regular_grid(bool bNonRegular)
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
bool NavierStokesFV<TDomain>::
request_finite_element_id(const std::vector<LFEID>& vLfeID)
{
//	check number
	if(vLfeID.size() != dim+1)
	{
		UG_LOG("NavierStokes:"
				" Wrong number of functions given. Need exactly "<<dim+1<<"\n");
		return false;
	}

	for(int d = 1; d < dim; ++d)
		if(vLfeID[0] != vLfeID[d])
		{
			UG_LOG("NavierStokes: trial spaces for velocity expected to be"
					" identical for all velocity components.\n");
			return false;
		}

	for(int d = 0; d <= dim; ++d)
		if(vLfeID[d].type() != LFEID::LAGRANGE)
		{
			UG_LOG("NavierStokes: 'fv' expects Lagrange trial space "
					"for velocity and pressure.\n");
			return false;
		}
	for(int d = 0; d < dim; ++d)
		if(vLfeID[d].order() != vLfeID[dim].order() + 1)
		{
			UG_LOG("NavierStokes: 'fv' expects Lagrange trial space "
					"P_k for velocity and P_{k-1} pressure.\n");
			return false;
		}

//	remember lfeID;
	m_vLFEID = vLfeID[0];
	m_pLFEID = vLfeID[dim];

	//	update assemble functions
	register_all_funcs(m_vLFEID, m_pLFEID);

	//	is supported
	return true;
}

template<typename TDomain>
void NavierStokesFV<TDomain>::
set_kinematic_viscosity(SmartPtr<UserData<number, dim> > data)
{
	m_imKinViscosity.set_data(data);
}

template<typename TDomain>
void NavierStokesFV<TDomain>::
set_density(SmartPtr<UserData<number, dim> > data)
{
	m_imDensitySCVF.set_data(data);
	m_imDensitySCVFp.set_data(data);
	m_imDensitySCV.set_data(data);
}

template<typename TDomain>
void NavierStokesFV<TDomain>::
set_source(SmartPtr<UserData<MathVector<dim>, dim> > data)
{
	m_imSource.set_data(data);
}

////////////////////////////////////////////////////////////////////////////////
//	assembling functions
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFV<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
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

	VGeom& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_vLFEID.order()+1);
	PGeom& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_pLFEID.order()+1);
	try{
		vgeo.update_local(roid, m_vLFEID.order());
		pgeo.update_local(roid, m_pLFEID.order());
	}
	UG_CATCH_THROW("NavierStokes: Cannot update Finite Volume Geometry.");

	static const int refDim = TElem::dim;
	if(!VGeom::usesHangingNodes)
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

	if(!PGeom::usesHangingNodes)
	{
		m_imDensitySCVFp.template set_local_ips<refDim>(pgeo.scvf_local_ips(),
		                                               pgeo.num_scvf_ips());

		const LocalShapeFunctionSet<dim>& rVTrialSpace =
			LocalShapeFunctionSetProvider::get<dim>(roid, LFEID(LFEID::LAGRANGE, m_vLFEID.order()));
		const MathVector<dim>* PLocIP = pgeo.scvf_local_ips();

		m_vvVShape.resize(pgeo.num_scvf_ips());
		for(size_t ip = 0; ip < m_vvVShape.size(); ++ip){
			m_vvVShape[ip].resize(rVTrialSpace.num_sh());
			for(size_t sh = 0; sh < rVTrialSpace.num_sh(); ++sh){
				m_vvVShape[ip][sh] = rVTrialSpace.shape(sh, PLocIP[ip]);
			}
		}

		const LocalShapeFunctionSet<dim>& rPTrialSpace =
			LocalShapeFunctionSetProvider::get<dim>(roid, LFEID(LFEID::LAGRANGE, m_pLFEID.order()));
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
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFV<TDomain>::
fsh_elem_loop()
{}


template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFV<TDomain>::
prep_elem(TElem* elem, const LocalVector& u)
{
// 	Update Geometry for this element
	VGeom& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_vLFEID.order()+1);
	PGeom& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_pLFEID.order()+1);
	try{
		vgeo.update(elem, this->template element_corners<TElem>(elem),
	               &(this->subset_handler()));
	    pgeo.update(elem, this->template element_corners<TElem>(elem),
	               &(this->subset_handler()));
	}
	UG_CATCH_THROW("NavierStokes: Cannot update Finite Volume Geometry.");

//	set local positions for imports
	if(VGeom::usesHangingNodes)
	{
		static const int refDim = TElem::dim;

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
	if(PGeom::usesHangingNodes)
	{
		static const int refDim = TElem::dim;

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
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFV<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u)
{
//	request geometry
	const VGeom& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_vLFEID.order()+1);
	const PGeom& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_pLFEID.order()+1);

// 	loop Sub Control Volume Faces (SCVF)
	for(size_t i = 0, ipCnt = 0; i < vgeo.num_scvf(); ++i)
	{
	// 	get current SCVF
		const typename VGeom::SCVF& scvf = vgeo.scvf(i);

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
	const typename PGeom::SCVF& scvf = pgeo.scvf(i);

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
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFV<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u)
{
//	request geometry
	const VGeom& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_vLFEID.order()+1);
	const PGeom& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_pLFEID.order()+1);

// 	loop Sub Control Volume Faces (SCVF)
	for(size_t i = 0, ipCnt = 0; i < vgeo.num_scvf(); ++i)
	{
// 	get current SCVF
	const typename VGeom::SCVF& scvf = vgeo.scvf(i);

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
		const typename PGeom::SCVF& scvf = pgeo.scvf(i);

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
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFV<TDomain>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u)
{
//	request geometry
	const VGeom& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_vLFEID.order()+1);

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0, ipOffset = 0; ip < vgeo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename VGeom::SCV& scv = vgeo.scv(ip);

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
				J(d1, sh, d1, sh) += integral;
			}
		}

	//	increase ip offset
		ipOffset += scv.num_ip();
	}
}


template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFV<TDomain>::
add_def_M_elem(LocalVector& d, const LocalVector& u)
{
//	request geometry
	const VGeom& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_vLFEID.order()+1);

// 	loop Sub Control Volumes (SCV)
	for(size_t i = 0, ipCnt = 0; i < vgeo.num_scv(); ++i)
	{
	// 	get current SCV
		const typename VGeom::SCV& scv = vgeo.scv(i);

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
                              * scv.weight(ip);

		//	next ip
			ipCnt++;
		}
	}
}


template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFV<TDomain>::
add_rhs_elem(LocalVector& d)
{
//	if zero data given, return
	if(!m_imSource.data_given()) return;

//	UG_THROW("Not implemented.")
}

////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

// register for all dim
template<>
void NavierStokesFV<Domain1d>::
register_all_funcs(const LFEID& vLfeID, const LFEID& pLfeID)
{
	UG_THROW("Not implemented.");
}

// register for all dim
template<>
void NavierStokesFV<Domain2d>::
register_all_funcs(const LFEID& vLfeID, const LFEID& pLfeID)
{
	typedef DimFVGeometry<dim, 2> FVGeom;
	register_func<Triangle, FVGeom, FVGeom >();
	register_func<Quadrilateral, FVGeom, FVGeom >();
}

// register for all dim
template<>
void NavierStokesFV<Domain3d>::
register_all_funcs(const LFEID& vLfeID, const LFEID& pLfeID)
{
	typedef DimFVGeometry<dim, 3> FVGeom;
	register_func<Tetrahedron, FVGeom, FVGeom >();
	register_func<Prism, FVGeom, FVGeom >();
	register_func<Hexahedron, FVGeom, FVGeom >();
}

template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFV<TDomain>::register_func()
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
template class NavierStokesFV<Domain1d>;
#endif
#ifdef UG_DIM_2
template class NavierStokesFV<Domain2d>;
#endif
#ifdef UG_DIM_3
template class NavierStokesFV<Domain3d>;
#endif

} // namespace NavierStokes
} // namespace ug
