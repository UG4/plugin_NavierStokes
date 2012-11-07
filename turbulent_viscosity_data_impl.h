/*
 * turbulent_viscosity_data_impl.h
 *
 *  Created on: 01.11.2012
 *      Author: Christian Wehner
 */

#ifndef __H__UG__NAVIER_STOKES_TURBULENT_VISCOSITY_DATA_IMPL_H_
#define __H__UG__NAVIER_STOKES_TURBULENT_VISCOSITY_DATA_IMPL_H_

#include "common/common.h"

#include "turbulent_viscosity_data.h"
#include "lib_disc/common/subset_group.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/quadrature/quadrature.h"
#include "lib_disc/local_finite_element/local_shape_function_set.h"
#include "lib_disc/spatial_disc/user_data/user_data.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"


namespace ug{
namespace NavierStokes{

template<typename TGridFunction>
void CRSmagorinskyTurbViscData<TGridFunction>::calculate_deformation_vol(const TGridFunction& u){
	//	get domain of grid function
	domain_type& domain = *u.domain().get();

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();

	DimCRFVGeometry<dim> geo;

	//	create Multiindex
	std::vector<MultiIndex<2> > multInd;

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	// assemble deformation tensor fluxes
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		//	get iterators
		ElemIterator iter = u.template begin<elem_type>(si);
		ElemIterator iterEnd = u.template end<elem_type>(si);

		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			elem_type* elem = *iter;
			//	get vertices and extract corner coordinates
			const size_t numVertices = elem->num_vertices();
			for(size_t i = 0; i < numVertices; ++i){
				vVrt[i] = elem->vertex(i);
				coCoord[i] = posAcc[vVrt[i]];
			};

			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

			static const size_t MaxNumSidesOfElem = 10;

			typedef MathVector<dim> MVD;
			std::vector<MVD> uValue(MaxNumSidesOfElem);
			std::vector<MVD> ipVelocity(MaxNumSidesOfElem);

			typename grid_type::template traits<side_type>::secure_container sides;

			UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Only elements of type elem_type are currently supported");

			m_grid->associated_elements_sorted(sides, static_cast<elem_type*>(elem) );

			size_t nofsides = geo.num_scv();

			size_t nip = geo.num_scvf();

			for (size_t s=0;s < nofsides;s++)
			{
				const typename DimCRFVGeometry<dim>::SCV& scv = geo.scv(s);
				for (size_t d=0;d<dim;d++){
				//	get indices of function fct on vertex
					u.inner_multi_indices(sides[s], d, multInd);

					//	read value of index from vector
					uValue[s][d]=DoFRef(u,multInd[0]);
				}
				m_acVolume += scv.volume();
			}

			for (size_t ip=0;ip<nip;ip++){
				// 	get current SCVF
				ipVelocity[ip] = 0;
				const typename DimCRFVGeometry<dim>::SCVF& scvf = geo.scvf(ip);
				for (size_t s=0;s < nofsides;s++){
					for (size_t d=0;d<dim;d++){
					    ipVelocity[ip][d] += scvf.shape(s)*uValue[s][d];
					};
				};
				dimMat ipDefTensorFlux;
				ipDefTensorFlux = 0;
				for (size_t d=0;d<dim;d++){
					for (size_t j=0;j<d;j++){
						ipDefTensorFlux[d][j]= 0.5 * (ipVelocity[d] * scvf.normal()[j] + ipVelocity[j] * scvf.normal()[d]);
					}
				}
				m_acDeformation[sides[scvf.from()]]+=ipDefTensorFlux;
				m_acDeformation[sides[scvf.to()]]-=ipDefTensorFlux;
			}
		}
		// average and compute turbulent viscosity , loop over sides
		iter = u.template begin<side_type>(si);
		iterEnd = u.template end<side_type>(si);
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			number tensorNorm=0;
			elem_type* elem = *iter;
			dimMat defTensor = m_acDeformation[elem];
			number delta = m_acVolume[elem];
			// complete deformation tensor computation by averaging
			defTensor/=delta;
			// compute norm of tensor
			for (size_t d1=0;d1<dim;d1++)
				for (size_t d2=0;d2<dim;d2++){
					tensorNorm += 2 * defTensor[d1][d2] * defTensor[d1][d2];
				}
			tensorNorm = sqrt(tensorNorm);
			// for possible other choices of delta see Fršhlich p 160
			delta = pow(delta,(number)1.0/(number)dim);
			// calculate viscosity constant
			m_acTurbulentViscosity[elem] = m_c * delta*delta * tensorNorm;
		}
		// transfer to lower levels

	}
}

template<typename TGridFunction>
bool CRSmagorinskyTurbViscData<TGridFunction>::update(const TGridFunction& u){
	return true;
}

template<typename TGridFunction>
bool CRDynamicTurbViscData<TGridFunction>::update(const TGridFunction& u){
	UG_LOG("Dynamic model not implemented yet.\n");
	return true;
}

} // namespace NavierStokes
} // end namespace ug

#endif /* __H__UG__NAVIER_STOKES_TURBULENT_VISCOSITY_DATA_IMPL_H_ */
