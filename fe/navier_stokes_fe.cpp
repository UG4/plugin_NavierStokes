/*
 * navier_stokes_fe.cpp
 *
 *  Created on: 20.09.2010
 *      Author: andreasvogel
 */

#include "navier_stokes_fe.h"

#include "common/util/provider.h"
#include "lib_disc/spatial_disc/disc_util/fe_geom.h"

namespace ug{
namespace NavierStokes{

////////////////////////////////////////////////////////////////////////////////
//	Constructor - set default values
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
NavierStokesFE<TDomain>::NavierStokesFE(const char* functions,
                                          const char* subsets)
: NavierStokesBase<TDomain>(functions, subsets)
{
	init();
};

template<typename TDomain>
NavierStokesFE<TDomain>::NavierStokesFE(const std::vector<std::string>& vFct,
                                          const std::vector<std::string>& vSubset)
: NavierStokesBase<TDomain>(vFct, vSubset)
{
	init();
};


template<typename TDomain>
void NavierStokesFE<TDomain>::init()
{
//	check number of functions
	if(this->num_fct() != dim+1)
		UG_THROW("Wrong number of functions: The ElemDisc 'NavierStokes'"
					   " needs exactly "<<dim+1<<" symbolic function.");

//	register imports
	this->register_import(m_imSource);
	this->register_import(m_imKinViscosity);
	this->register_import(m_imDensity);

	m_imSource.set_rhs_part();
	m_imDensity.set_mass_part();

	//	default value for density
	base_type::set_density(1.0);

	//	update assemble functions
	register_all_funcs(1, 1);
}

template<typename TDomain>
bool NavierStokesFE<TDomain>::request_non_regular_grid(bool bNonRegular)
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
bool NavierStokesFE<TDomain>::
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

//	remember lfeID;
	m_vLFEID = vLfeID[0];
	m_pLFEID = vLfeID[dim];

//	set order
	m_vorder = vLfeID[0].order();
	m_porder = vLfeID[dim].order();

	//	update assemble functions
	register_all_funcs(m_vorder, m_porder);

	//	is supported
	return true;
}

template<typename TDomain>
void NavierStokesFE<TDomain>::
set_kinematic_viscosity(SmartPtr<UserData<number, dim> > data)
{
	m_imKinViscosity.set_data(data);
}

template<typename TDomain>
void NavierStokesFE<TDomain>::
set_density(SmartPtr<UserData<number, dim> > data)
{
	m_imDensity.set_data(data);
}

template<typename TDomain>
void NavierStokesFE<TDomain>::
set_source(SmartPtr<UserData<MathVector<dim>, dim> > data)
{
	m_imSource.set_data(data);
}

////////////////////////////////////////////////////////////////////////////////
//	assembling functions
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFE<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
//	check, that kinematic Viscosity has been set
	if(!m_imKinViscosity.data_given())
		UG_THROW("NavierStokes: Kinematic Viscosity has not been set, but is required.");

//	check, that Density has been set
	if(!m_imDensity.data_given())
		UG_THROW("NavierStokes: Density has not been set, but is required.");

	static DimFEGeometry<dim, dim>& vgeo = Provider<VGeom>::get(m_vorder);
	static DimFEGeometry<dim, dim>& pgeo = Provider<PGeom>::get(m_porder);
	try{
		vgeo.update_local(roid, m_vLFEID, 2*m_vorder+1);
		pgeo.update_local(roid, m_pLFEID, 2*m_porder+1);
	}
	UG_CATCH_THROW("NavierStokes: Cannot update Finite Element Geometry.");

//	set local positions for imports
	m_imKinViscosity.set_local_ips(vgeo.local_ips(), vgeo.num_ip());
	m_imDensity. set_local_ips(vgeo.local_ips(), vgeo.num_ip());
	m_imSource. set_local_ips(vgeo.local_ips(), vgeo.num_ip());
}

template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFE<TDomain>::
fsh_elem_loop()
{}


template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFE<TDomain>::
prep_elem(TElem* elem, const LocalVector& u)
{
// 	Update Geometry for this element
	static DimFEGeometry<dim, dim>& vgeo = Provider<VGeom>::get(m_vorder);
	static DimFEGeometry<dim, dim>& pgeo = Provider<PGeom>::get(m_porder);
	try{
		vgeo.update(elem, this->template element_corners<TElem>(elem), m_vLFEID, 2*m_vorder+1);
	    pgeo.update(elem, this->template element_corners<TElem>(elem), m_pLFEID, 2*m_porder+1);
	}
	UG_CATCH_THROW("NavierStokes: Cannot update Finite Element Geometry.");

//	set global positions for imports
	m_imKinViscosity.set_global_ips(vgeo.global_ips(), vgeo.num_ip());
	m_imDensity.set_global_ips(vgeo.global_ips(), vgeo.num_ip());
	m_imSource.set_global_ips(vgeo.global_ips(), vgeo.num_ip());
}

template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFE<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u)
{
//	request geometry
	static const DimFEGeometry<dim, dim>& vgeo = Provider<VGeom>::get(m_vorder);
	static const DimFEGeometry<dim, dim>& pgeo = Provider<PGeom>::get(m_porder);

	for (size_t ip = 0; ip < vgeo.num_ip(); ++ip){

		for (int j = 0; j < dim; ++j){
			for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){

				for (size_t ush = 0; ush < vgeo.num_sh(); ++ush){
					for (int i = 0; i < dim; ++i) {

						J(j, vsh, j, ush) +=  m_imKinViscosity[ip]
						               * vgeo.global_grad(ip, ush)[i]
						               * vgeo.global_grad(ip, vsh)[i]
						               * vgeo.weight(ip);

						if(!m_bLaplace) UG_THROW("Not implemented.");
					}
				}
			}
		}

		for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){
			for (int vdim = 0; vdim < dim; ++vdim){
				for (size_t psh = 0; psh < pgeo.num_sh(); ++psh)
					J(vdim, vsh, _P_, psh) -= pgeo.shape(ip, psh)
								   * vgeo.global_grad(ip, vsh)[vdim]
								   * vgeo.weight(ip);
			}
		}

		for (size_t psh = 0; psh < pgeo.num_sh(); ++psh){
			for (size_t ush = 0; ush < vgeo.num_sh(); ++ush){
				for (int udim = 0; udim < dim; ++udim) {
						J(_P_, psh, udim, ush) +=   vgeo.global_grad(ip, ush)[udim]
						               * pgeo.shape(ip, psh)
						               * vgeo.weight(ip);
				}
			}
		}
	}

}

template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFE<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u)
{
//	request geometry
	static const DimFEGeometry<dim, dim>& vgeo = Provider<VGeom>::get(m_vorder);
	static const DimFEGeometry<dim, dim>& pgeo = Provider<PGeom>::get(m_porder);

	for (size_t ip = 0; ip < vgeo.num_ip(); ++ip){

		for (int j = 0; j < dim; ++j){
			for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){

				for (size_t ush = 0; ush < vgeo.num_sh(); ++ush){
					for (int i = 0; i < dim; ++i) {

						d(j, vsh) +=  m_imKinViscosity[ip]
						               * u(j, ush)
						               * vgeo.global_grad(ip, ush)[i]
						               * vgeo.global_grad(ip, vsh)[i]
						               * vgeo.weight(ip);

						if(!m_bLaplace) UG_THROW("Not implemented.");
					}
				}
			}
		}

		number pressure = 0.0;
		for (size_t psh = 0; psh < pgeo.num_sh(); ++psh)
			pressure += u(_P_, psh) * pgeo.shape(ip, psh);

		for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){
			for (int vdim = 0; vdim < dim; ++vdim){
				d(vdim, vsh) -= pressure
							   * vgeo.global_grad(ip, vsh)[vdim]
							   * vgeo.weight(ip);
			}
		}

		number divu = 0.0;
		for (size_t ush = 0; ush < vgeo.num_sh(); ++ush){
			for (int udim = 0; udim < dim; ++udim) {
				divu += u(udim, ush) * vgeo.global_grad(ip, ush)[udim];
			}
		}

		for (size_t psh = 0; psh < pgeo.num_sh(); ++psh){
						d(_P_, psh) += divu
						               * pgeo.shape(ip, psh)
						               * vgeo.weight(ip);
		}
	}

}


template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFE<TDomain>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u)
{
//	request geometry
//	static const DimFEGeometry<dim, dim>& vgeo = Provider<VGeom>::get(m_vorder);

	UG_THROW("Not implemented.");
}


template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFE<TDomain>::
add_def_M_elem(LocalVector& d, const LocalVector& u)
{
//	request geometry
//	static const DimFEGeometry<dim, dim>& vgeo = Provider<VGeom>::get(m_vorder);

	UG_THROW("Not implemented.");
}


template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFE<TDomain>::
add_rhs_elem(LocalVector& d)
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
void NavierStokesFE<Domain1d>::
register_all_funcs(int vorder, int porder)
{
	UG_THROW("Not implemented.");
}

// register for all dim
template<>
void NavierStokesFE<Domain2d>::
register_all_funcs(int vorder, int porder)
{
	typedef DimFEGeometry<dim, dim> FVGeom;
	register_func<Triangle, FVGeom, FVGeom >();
	register_func<Quadrilateral, FVGeom, FVGeom >();
}

// register for all dim
template<>
void NavierStokesFE<Domain3d>::
register_all_funcs(int vorder, int porder)
{
	typedef DimFEGeometry<dim, dim> FVGeom;
	register_func<Tetrahedron, FVGeom, FVGeom >();
	register_func<Prism, FVGeom, FVGeom >();
	register_func<Hexahedron, FVGeom, FVGeom >();
}

template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFE<TDomain>::register_func()
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
template class NavierStokesFE<Domain1d>;
#endif
#ifdef UG_DIM_2
template class NavierStokesFE<Domain2d>;
#endif
#ifdef UG_DIM_3
template class NavierStokesFE<Domain3d>;
#endif

} // namespace NavierStokes
} // namespace ug
