/*
 * navier_stokes_fe.h
 *
 *  Created on: 20.09.2010
 *      Author: andreasvogel
 */

#include "navier_stokes.h"

#include "common/util/provider.h"
#include "lib_disc/spatial_disc/disc_util/finite_element_geometry.h"

namespace ug{
namespace NavierStokes{

template<typename TDomain>
template<typename TElem, typename VGeomProvider, typename PGeomProvider>
void NavierStokes<TDomain>::
prep_elem_loop_fe()
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

	static DimFEGeometry<dim, dim>& vgeo = VGeomProvider::get();
	static DimFEGeometry<dim, dim>& pgeo = PGeomProvider::get();
	try{
		vgeo.update_local(roid, m_vLFEID, m_quadOrder);
		pgeo.update_local(roid, m_pLFEID, m_quadOrder);
	}
	UG_CATCH_THROW("NavierStokes: Cannot update Finite Element Geometry.");

//	set local positions for imports
	typedef typename reference_element_traits<TElem>::reference_element_type
																ref_elem_type;
	m_imKinViscosity.set_local_ips(vgeo.local_ips(), vgeo.num_ip());
	m_imDensitySCVF. set_local_ips(vgeo.local_ips(), vgeo.num_ip());
	m_imDensitySCV.  set_local_ips(vgeo.local_ips(), vgeo.num_ip());
	m_imSource. set_local_ips(vgeo.local_ips(), vgeo.num_ip());
}

template<typename TDomain>
template<typename TElem, typename VGeomProvider, typename PGeomProvider>
void NavierStokes<TDomain>::
fsh_elem_loop_fe()
{}


template<typename TDomain>
template<typename TElem, typename VGeomProvider, typename PGeomProvider>
void NavierStokes<TDomain>::
prep_elem_fe(TElem* elem, const LocalVector& u)
{
// 	Update Geometry for this element
	static DimFEGeometry<dim, dim>& vgeo = VGeomProvider::get();
	static DimFEGeometry<dim, dim>& pgeo = PGeomProvider::get();
	try{
		vgeo.update(elem, this->template element_corners<TElem>(elem), m_vLFEID, m_quadOrder);
	    pgeo.update(elem, this->template element_corners<TElem>(elem), m_pLFEID, m_quadOrder);
	}
	UG_CATCH_THROW("NavierStokes: Cannot update Finite Element Geometry.");

//	set global positions for imports
	m_imKinViscosity.set_global_ips(vgeo.global_ips(), vgeo.num_ip());
	m_imDensitySCVF.set_global_ips(vgeo.global_ips(), vgeo.num_ip());
	m_imDensitySCV.set_global_ips(vgeo.global_ips(), vgeo.num_ip());
	m_imSource.set_global_ips(vgeo.global_ips(), vgeo.num_ip());
}

template<typename TDomain>
template<typename TElem, typename VGeomProvider, typename PGeomProvider>
void NavierStokes<TDomain>::
add_jac_A_elem_fe(LocalMatrix& J, const LocalVector& u)
{
//	request geometry
	static const DimFEGeometry<dim, dim>& vgeo = VGeomProvider::get();
	static const DimFEGeometry<dim, dim>& pgeo = PGeomProvider::get();

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
template<typename TElem, typename VGeomProvider, typename PGeomProvider>
void NavierStokes<TDomain>::
add_def_A_elem_fe(LocalVector& d, const LocalVector& u)
{
//	request geometry
	static const DimFEGeometry<dim, dim>& vgeo = VGeomProvider::get();
	static const DimFEGeometry<dim, dim>& pgeo = PGeomProvider::get();

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
template<typename TElem, typename VGeomProvider, typename PGeomProvider>
void NavierStokes<TDomain>::
add_jac_M_elem_fe(LocalMatrix& J, const LocalVector& u)
{
//	request geometry
//	static const DimFEGeometry<dim, dim>& vgeo = VGeomProvider::get();

	UG_THROW("Not implemented.");
}


template<typename TDomain>
template<typename TElem, typename VGeomProvider, typename PGeomProvider>
void NavierStokes<TDomain>::
add_def_M_elem_fe(LocalVector& d, const LocalVector& u)
{
//	request geometry
//	static const DimFEGeometry<dim, dim>& vgeo = VGeomProvider::get();

	UG_THROW("Not implemented.");
}


template<typename TDomain>
template<typename TElem, typename VGeomProvider, typename PGeomProvider>
void NavierStokes<TDomain>::
add_rhs_elem_fe(LocalVector& d)
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
register_all_fe_funcs(int order)
{
	UG_THROW("Not implemented.");
}

// register for all dim
template<>
void NavierStokes<Domain2d>::
register_all_fe_funcs(int order)
{
	typedef DimFEGeometry<dim, dim> FVGeom;
	register_fe_func<Triangle, VFlexGeomProvider<FVGeom>, PFlexGeomProvider<FVGeom> >();
	register_fe_func<Quadrilateral, VFlexGeomProvider<FVGeom>, PFlexGeomProvider<FVGeom> >();
}

// register for all dim
template<>
void NavierStokes<Domain3d>::
register_all_fe_funcs(int order)
{
	typedef DimFEGeometry<dim, dim> FVGeom;
	register_fe_func<Tetrahedron, VFlexGeomProvider<FVGeom>, PFlexGeomProvider<FVGeom> >();
	register_fe_func<Prism, VFlexGeomProvider<FVGeom>, PFlexGeomProvider<FVGeom> >();
	register_fe_func<Hexahedron, VFlexGeomProvider<FVGeom>, PFlexGeomProvider<FVGeom> >();
}

template<typename TDomain>
template<typename TElem, typename VGeomProvider, typename PGeomProvider>
void NavierStokes<TDomain>::register_fe_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->enable_fast_add_elem(true);
	set_prep_elem_loop_fct(	id, &T::template prep_elem_loop_fe<TElem, VGeomProvider, PGeomProvider>);
	set_prep_elem_fct(	 	id, &T::template prep_elem_fe<TElem, VGeomProvider, PGeomProvider>);
	set_fsh_elem_loop_fct( 	id, &T::template fsh_elem_loop_fe<TElem, VGeomProvider, PGeomProvider>);
	set_add_jac_A_elem_fct(	id, &T::template add_jac_A_elem_fe<TElem, VGeomProvider, PGeomProvider>);
	set_add_jac_M_elem_fct(	id, &T::template add_jac_M_elem_fe<TElem, VGeomProvider, PGeomProvider>);
	set_add_def_A_elem_fct(	id, &T::template add_def_A_elem_fe<TElem, VGeomProvider, PGeomProvider>);
	set_add_def_M_elem_fct(	id, &T::template add_def_M_elem_fe<TElem, VGeomProvider, PGeomProvider>);
	set_add_rhs_elem_fct(	id, &T::template add_rhs_elem_fe<TElem, VGeomProvider, PGeomProvider>);
}

} // namespace NavierStokes
} // namespace ug
