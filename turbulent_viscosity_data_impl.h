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
#include "lib_grid/algorithms/attachment_util.h"

namespace ug{
namespace NavierStokes{

template<typename TGridFunction>
void CRSmagorinskyTurbViscData<TGridFunction>::update(){
	//	get domain of grid function
	domain_type& domain = *m_u->domain().get();

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();

	DimCRFVGeometry<dim> geo;

	// initialize attachment values with 0
	SetAttachmentValues(m_acDeformation , m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);
	SetAttachmentValues(m_acTurbulentViscosity, m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);
	SetAttachmentValues(m_acVolume,m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);

	//	create Multiindex
	std::vector<MultiIndex<2> > multInd;

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	// assemble deformation tensor fluxes
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		//	get iterators
		ElemIterator iter = m_u->template begin<elem_type>(si);
		ElemIterator iterEnd = m_u->template end<elem_type>(si);

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
				// UG_LOG("co_coord(" << i<< "+1,:)=" << coCoord[i] << "\n");
			};

			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

			static const size_t MaxNumSidesOfElem = 10;

			typedef MathVector<dim> MVD;
			std::vector<MVD> uValue(MaxNumSidesOfElem);
			MVD ipVelocity;

			typename grid_type::template traits<side_type>::secure_container sides;

			UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Only elements of type elem_type are currently supported");

			domain.grid()->associated_elements_sorted(sides, static_cast<elem_type*>(elem) );

			size_t nofsides = geo.num_scv();

			size_t nip = geo.num_scvf();

			for (size_t s=0;s < nofsides;s++)
			{
				const typename DimCRFVGeometry<dim>::SCV& scv = geo.scv(s);
				for (int d=0;d<dim;d++){
				//	get indices of function fct on vertex
					m_u->multi_indices(sides[s], d, multInd);

					//	read value of index from vector
					uValue[s][d]=DoFRef(*m_u,multInd[0]);
				}
				// UG_LOG("scv.volume(" << s << ")=" << scv.volume() << "\n");
				m_acVolume[sides[s]] += scv.volume();
			}

			for (size_t ip=0;ip<nip;ip++){
				// 	get current SCVF
				ipVelocity = 0;
				const typename DimCRFVGeometry<dim>::SCVF& scvf = geo.scvf(ip);
				for (size_t s=0;s < nofsides;s++){
					for (int d=0;d<dim;d++){
					    ipVelocity[d] += scvf.shape(s)*uValue[s][d];
					};
				};
				dimMat ipDefTensorFlux;
				ipDefTensorFlux = 0;
				for (int d=0;d<dim;d++){
					for (int j=0;j<d;j++){
						ipDefTensorFlux[d][j]= 0.5 * (ipVelocity[d] * scvf.normal()[j] + ipVelocity[j] * scvf.normal()[d]);
						//UG_LOG("ipTensorFlux(" << d << "," << j << ")=" << ipDefTensorFlux[d][j] << "\n");
					}
				}
				//UG_LOG("Tensor(" << scvf.from() << ")=" << m_acDeformation[sides[scvf.from()]] << "\n");
				//UG_LOG("Tensor(" << scvf.to() << ")=" << m_acDeformation[sides[scvf.to()]] << "\n");
				m_acDeformation[sides[scvf.from()]]+=ipDefTensorFlux;
				m_acDeformation[sides[scvf.to()]]-=ipDefTensorFlux;
				//UG_LOG("&side_from = " << sides[scvf.from()] << "\n");
				//UG_LOG("Tensor(" << scvf.from() << ")=" << m_acDeformation[sides[scvf.from()]] << "\n");
				//UG_LOG("&side_to = " << sides[scvf.to()] << "\n");
				//UG_LOG("Tensor(" << scvf.to() << ")=" << m_acDeformation[sides[scvf.to()]] << "\n");
			}
		}
		// average and compute turbulent viscosity , loop over sides
		SideIterator sideIter = m_u->template begin<side_type>(si);
		SideIterator sideIterEnd = m_u->template end<side_type>(si);
		for(  ;sideIter !=sideIterEnd; sideIter++)
		{
			//	get Elem
			number tensorNorm=0;
			side_type* elem = *sideIter;
			dimMat defTensor = m_acDeformation[elem];
			number delta = m_acVolume[elem];
			// UG_LOG("delta=" << delta <<  "\n");
			// UG_LOG("&side = " << elem << "\n");
			// complete deformation tensor computation by averaging
			defTensor/=delta;
			// compute norm of tensor
			for (int d1=0;d1<dim;d1++)
				for (int d2=0;d2<dim;d2++){
					// UG_LOG("tensor(" << d1 << "," << d2 << ")=" << defTensor[d1][d2] << "\n");
					tensorNorm += 2 * defTensor[d1][d2] * defTensor[d1][d2];
				}
			tensorNorm = sqrt(tensorNorm);
			// for possible other choices of delta see Fršhlich p 160
			delta = pow(delta,(number)1.0/(number)dim);
			// calculate viscosity constant
			// UG_LOG("c=" << m_c << " delta=" << delta << " tensor=" << tensorNorm << "\n");
			m_acTurbulentViscosity[elem] = m_c * delta*delta * tensorNorm;
			// UG_LOG(m_acTurbulentViscosity[elem] << "\n");
		}
		// transfer to lower levels, averaging over child edges (2d) / child faces (3d)
		for (size_t lev=m_spApproxSpace->num_levels()-2;(int)lev>=0;lev--){
			// UG_LOG("level=" << lev << "\n");
			const LevelDoFDistribution& lDD = *m_spApproxSpace->level_dof_distribution(lev);
			const MultiGrid& grid = lDD.multi_grid();
			typedef typename LevelDoFDistribution::template traits<side_type>::const_iterator coarseLevelSideIter;
			coarseLevelSideIter clsIter, clsIterEnd;
			clsIter = lDD.template begin<side_type>(si);
			clsIterEnd = lDD.template end<side_type>(si);
			for (;clsIter != clsIterEnd;clsIter++){
				side_type* elem = *clsIter;
				size_t numChildren = grid.num_children<side_type>(elem);
				number avgValue=0;
				for (size_t i=0;i<numChildren;i++){
					avgValue += m_acTurbulentViscosity[grid.get_child<side_type>(elem, i)];
				}
				avgValue/=numChildren;
				m_acTurbulentViscosity[elem] = avgValue;
				// UG_LOG(" " << m_acTurbulentViscosity[elem] << "\n");
			}
			if (lev==0) break;
		}
	}
}

template<typename TGridFunction>
void CRDynamicTurbViscData<TGridFunction>::update(){
	UG_LOG("Dynamic model not implemented yet.\n");
}

} // namespace NavierStokes
} // end namespace ug

#endif /* __H__UG__NAVIER_STOKES_TURBULENT_VISCOSITY_DATA_IMPL_H_ */
