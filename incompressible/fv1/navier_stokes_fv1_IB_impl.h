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
prep_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	UG_LOG("hier Aufruf:  'prep_elem()...\n");
	UG_LOG("juhuuuu!...prep_elem()\n");


// 	Update Geometry for this element
	TFVGeom& geo = GeomProvider<TFVGeom>::get();
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


//  finally call the methods of the 'FV1IBGeometry'-class to adapt the element geometry data
//  for elements which are cut by the inner boundary:
	FV1IBGeometry<TElem, dim>& geo_adapt = GeomProvider<FV1IBGeometry<TElem, dim> >::get();
/*	geo_adapt = static_cast<FV1IBGeometry<TElem, dim>& >(geo);

	geo_adapt.adapt_normals(elem, vCornerCoords);
	geo_adapt.adapt_integration_points(elem, vCornerCoords);

	//static_cast<ug::FV1IBGeometry<TElem, dim> >(geo).adapt_integration_points(elem, vCornerCoords);
*/
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesFV1IB<TDomain>::
add_jac_A_elem_IB(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{

}



template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesFV1IB<TDomain>::
add_def_A_elem_IB(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{

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
		//register_substitutes<Edge, FV1Geometry<Edge, dim> >();
		register_substitutes<Edge, FV1IBGeometry<Edge, dim> >();
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


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NavierStokesFV1IB<TDomain>::
register_substitutes()
{
	UG_LOG("in register_substitutes() IB...\n");

	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->enable_fast_add_elem(true);

	this->set_prep_elem_fct(	 	id, &T::template prep_elem<TElem, TFVGeom>);
   	this->set_add_jac_A_elem_fct(	id, &T::template add_jac_A_elem<TElem, TFVGeom>);
 	this->set_add_def_A_elem_fct(	id, &T::template add_def_A_elem<TElem, TFVGeom>);

	UG_LOG("in register_substitutes() IB...jeppaaa......\n");

}

} // namespace NavierStokes
} // namespace ug

#endif /* NAVIER_STOKES_FV1_IB_IMPL_H_ */
