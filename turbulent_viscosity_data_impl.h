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

template <typename TData, int dim, typename TImpl,typename TGridFunction>
void StdTurbulentViscosityData<TData,dim,TImpl,TGridFunction>::transferToLowerLevels(aSideNumber& aaData,ApproximationSpace<domain_type>& approximationSpace){
	for(int si = 0; si < approximationSpace.domain()->subset_handler()->num_subsets(); ++si){
		// transfer to lower levels, averaging over child edges (2d) / child faces (3d)
		for (size_t lev=approximationSpace.num_levels()-2;(int)lev>=0;lev--){
			const LevelDoFDistribution& lDD = *approximationSpace.level_dof_distribution(lev);
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
					avgValue += aaData[grid.get_child<side_type>(elem, i)];
				}
				avgValue/=numChildren;
				aaData[elem] = avgValue;
			}
			if (lev==0) break;
		}
	}
}

template <typename TData, int dim, typename TImpl,typename TGridFunction>
void StdTurbulentViscosityData<TData,dim,TImpl,TGridFunction>::fillAttachment(aSideDimVector& aaU,SmartPtr<TGridFunction> u){
	//	get domain
	domain_type& domain = *u->domain().get();
	//	create Multiindex
	std::vector<MultiIndex<2> > multInd;
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
		ElemIterator iter = u->template begin<side_type>(si);
		ElemIterator iterEnd = u->template end<side_type>(si);
		for(  ;iter !=iterEnd; ++iter)
		{
			side_type* side = iter;
			for (int d=0;d<dim;d++){
				u->multi_indices(side, d, multInd);
				aaU[side][d]=DoFRef(*u,multInd[0]);
			}
		}
	}
}

// go over all elements, interpolate data to barycenter, average by multiplying with corresponding element volume and deviding by complete adjacent element volume
template <typename TData, int dim, typename TImpl,typename TGridFunction>
template <typename VType>
void StdTurbulentViscosityData<TData,dim,TImpl,TGridFunction>::elementFilter(PeriodicAttachmentAccessor<side_type,Attachment<VType> >& aaUHat,aSideNumber& aaVol,const PeriodicAttachmentAccessor<side_type,Attachment<VType> >& aaU){
	//	get domain of grid function
	domain_type& domain = *m_uInfo->domain().get();
	DimCRFVGeometry<dim> geo;

	// set attachment values to zero
	SetAttachmentValues(aaUHat , m_uInfo->template begin<side_type>(), m_uInfo->template end<side_type>(), 0);
	SetAttachmentValues(aaVol , m_uInfo->template begin<side_type>(), m_uInfo->template end<side_type>(), 0);

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();

	// assemble deformation tensor fluxes
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		if (si>0) continue;
		//	get iterators
		ElemIterator iter = m_uInfo->template begin<elem_type>(si);
		ElemIterator iterEnd = m_uInfo->template end<elem_type>(si);

		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			elem_type* elem = *iter;
			//	get vertices and extract corner coordinates
			const size_t numVertices = elem->num_vertices();
			MathVector<dim> bary,localBary;
			bary = 0;

			for(size_t i = 0; i < numVertices; ++i){
				vVrt[i] = elem->vertex(i);
				coCoord[i] = posAcc[vVrt[i]];
				bary += coCoord[i];
			};
			bary /= numVertices;

			//	reference object id
			ReferenceObjectID roid = elem->reference_object_id();

			//	get trial space
			const LocalShapeFunctionSet<dim>& rTrialSpace =
			LocalShapeFunctionSetProvider::get<dim>(roid, LFEID(LFEID::CROUZEIX_RAVIART, 1));

			//	get Reference Mapping
			DimReferenceMapping<dim, dim>& map = ReferenceMappingProvider::get<dim, dim>(roid, coCoord);

			map.global_to_local(localBary,bary);

			typename grid_type::template traits<side_type>::secure_container sides;

			UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Only elements of type elem_type are currently supported");

			domain.grid()->associated_elements_sorted(sides, static_cast<elem_type*>(elem) );

			//	memory for shapes
			std::vector<number> vShape;

			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

			rTrialSpace.shapes(vShape, localBary);

			size_t nofsides = geo.num_scv();

			VType value;
			value = 0;
			number elementVolume = 0;
			for (size_t s=0;s<nofsides;s++){
				const typename DimCRFVGeometry<dim>::SCV& scv = geo.scv(s);
				VType localValue = aaU[sides[s]];
				localValue *= vShape[s];
				value += localValue;
				elementVolume += scv.volume();
			}
			value *= elementVolume;
			for (size_t s=0;s<nofsides;s++){
				aaUHat[sides[s]] += value;
			}
		}
	}
	PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
	// average
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		SideIterator sideIter = m_uInfo->template begin<side_type>(si);
		SideIterator sideIterEnd = m_uInfo->template end<side_type>(si);
		for(  ;sideIter !=sideIterEnd; sideIter++)
		{
			side_type* elem = *sideIter;
			if (pbm && pbm->is_slave(elem)) continue;
			aaUHat[elem]/=(number)aaVol[elem];
		}
	}
}

// go over all elements, interpolate data to barycenter, average by multiplying with corresponding element volume and deviding by complete adjacent element volume
template <typename TData, int dim, typename TImpl,typename TGridFunction>
void StdTurbulentViscosityData<TData,dim,TImpl,TGridFunction>::elementFilter(aSideDimVector& aaUHat,aSideNumber& aaVol,SmartPtr<TGridFunction> u){
	//	get domain of grid function
	domain_type& domain = *m_uInfo->domain().get();
	DimCRFVGeometry<dim> geo;

	 std::vector<MultiIndex<2> > multInd;

	// set attachment values to zero
	SetAttachmentValues(aaUHat , m_uInfo->template begin<side_type>(), m_uInfo->template end<side_type>(), 0);
	SetAttachmentValues(aaVol , m_uInfo->template begin<side_type>(), m_uInfo->template end<side_type>(), 0);

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();

	// assemble deformation tensor fluxes
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		if (si>0) continue;
		//	get iterators
		ElemIterator iter = m_uInfo->template begin<elem_type>(si);
		ElemIterator iterEnd = m_uInfo->template end<elem_type>(si);

		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			elem_type* elem = *iter;
			//	get vertices and extract corner coordinates
			const size_t numVertices = elem->num_vertices();
			MathVector<dim> bary,localBary;
			bary = 0;

			for(size_t i = 0; i < numVertices; ++i){
				vVrt[i] = elem->vertex(i);
				coCoord[i] = posAcc[vVrt[i]];
				bary += coCoord[i];
			};
			bary /= numVertices;

			//	reference object id
			ReferenceObjectID roid = elem->reference_object_id();

			//	get trial space
			const LocalShapeFunctionSet<dim>& rTrialSpace =
			LocalShapeFunctionSetProvider::get<dim>(roid, LFEID(LFEID::CROUZEIX_RAVIART, 1));

			//	get Reference Mapping
			DimReferenceMapping<dim, dim>& map = ReferenceMappingProvider::get<dim, dim>(roid, coCoord);

			map.global_to_local(localBary,bary);

			typename grid_type::template traits<side_type>::secure_container sides;

			UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Only elements of type elem_type are currently supported");

			domain.grid()->associated_elements_sorted(sides, static_cast<elem_type*>(elem) );

			//	memory for shapes
			std::vector<number> vShape;

			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

			rTrialSpace.shapes(vShape, localBary);

			size_t nofsides = geo.num_scv();

			MathVector<dim> value;
			value = 0;
			number elementVolume = 0;
			for (size_t s=0;s<nofsides;s++){
				const typename DimCRFVGeometry<dim>::SCV& scv = geo.scv(s);
				MathVector<dim> localValue;
				for (int d=0;d<dim;d++){
					 localValue[d]=DoFRef(*u,multInd[0]);
					 u->multi_indices(sides[s], d, multInd);
				}
				localValue *= vShape[s];
				value += localValue;
				elementVolume += scv.volume();
			}
			value *= elementVolume;
			for (size_t s=0;s<nofsides;s++){
				aaUHat[sides[s]] += value;
			}
		}
	}
	PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
	// average
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		SideIterator sideIter = m_uInfo->template begin<side_type>(si);
		SideIterator sideIterEnd = m_uInfo->template end<side_type>(si);
		for(  ;sideIter !=sideIterEnd; sideIter++)
		{
			side_type* elem = *sideIter;
			if (pbm && pbm->is_slave(elem)) continue;
			aaUHat[elem]/=(number)aaVol[elem];
		}
	}
}


template <typename TData, int dim, typename TImpl,typename TGridFunction>
void StdTurbulentViscosityData<TData,dim,TImpl,TGridFunction>::assembleDeformationTensor(aSideTensor& aaDefTensor,aSideNumber& aaVol,SmartPtr<TGridFunction> u){
	//	get domain
	domain_type& domain = *u->domain().get();
	// get periodic boundary manager
	PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();

	//	create Multiindex
	std::vector<MultiIndex<2> > multInd;

	DimCRFVGeometry<dim> geo;

	// set attachment values to zero
	SetAttachmentValues(aaDefTensor , u->template begin<side_type>(), u->template end<side_type>(), 0);
	SetAttachmentValues(aaVol , u->template begin<side_type>(), u->template end<side_type>(), 0);

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();

	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
		if (si>0) continue;
			//	get iterators
			ElemIterator iter = u->template begin<elem_type>(si);
			ElemIterator iterEnd = u->template end<elem_type>(si);

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
						u->multi_indices(sides[s], d, multInd);
						//	read value of index from vector
						uValue[s][d]=DoFRef(*u,multInd[0]);
					}
					aaVol[sides[s]] += scv.volume();
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
						for (int j=0;j<dim;j++){
							ipDefTensorFlux[d][j]= 0.5 * (ipVelocity[d] * scvf.normal()[j] + ipVelocity[j] * scvf.normal()[d]);
						}
					}
					aaDefTensor[sides[scvf.from()]]+=ipDefTensorFlux;
					aaDefTensor[sides[scvf.to()]]-=ipDefTensorFlux;
				}
			}
	}
	// average
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
		SideIterator sideIter = u->template begin<side_type>(si);
		SideIterator sideIterEnd = u->template end<side_type>(si);
		for(  ;sideIter !=sideIterEnd; sideIter++)
		{
			side_type* elem = *sideIter;
			if (pbm && pbm->is_slave(elem)) continue;
			aaDefTensor[elem]/=(number)aaVol[elem];
/*			MathVector<dim> posCo;
						for (int d=0;d<dim;d++)
										posCo[d] = 0.5*posAcc[elem->vertex(0)][d] + 0.5*posAcc[elem->vertex(1)][d];
						UG_LOG(" c=" << posCo << "\n");
						for (int d1=0;d1<dim;d1++)
							for (int d2=0;d2<dim;d2++)
									UG_LOG(" tensor(" << d1 << "," << d2 << ")=" << aaDefTensor[elem][d1][d2] << "\n");*/
		}
	}
}

template <typename TData, int dim, typename TImpl,typename TGridFunction>
void StdTurbulentViscosityData<TData,dim,TImpl,TGridFunction>::normalizeTensor(aSideTensor& aaTensor){
	//	get domain of grid function
	domain_type& domain = *m_uInfo->domain().get();
	// average
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
		if (si>0) continue;
		SideIterator sideIter = m_uInfo->template begin<side_type>(si);
		SideIterator sideIterEnd = m_uInfo->template end<side_type>(si);
		for(  ;sideIter !=sideIterEnd; sideIter++)
		{
			side_type* side = *sideIter;
			aaTensor[side]/=(number)FNorm(aaTensor[side]);
		}
	}
}

template <typename TData, int dim, typename TImpl,typename TGridFunction>
void StdTurbulentViscosityData<TData,dim,TImpl,TGridFunction>::assembleDeformationTensor(aSideTensor& aaDefTensor,aSideNumber& aaVol,aSideDimVector aaU){
	//	get domain of grid function
	domain_type& domain = *m_uInfo->domain().get();

	//	create Multiindex
	std::vector<MultiIndex<2> > multInd;

	DimCRFVGeometry<dim> geo;

	SetAttachmentValues(aaDefTensor , m_uInfo->template begin<side_type>(), m_uInfo->template end<side_type>(), 0);
	SetAttachmentValues(aaVol , m_uInfo->template begin<side_type>(), m_uInfo->template end<side_type>(), 0);

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();

	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
		if (si>0) continue;
			//	get iterators
			ElemIterator iter = m_uInfo->template begin<elem_type>(si);
			ElemIterator iterEnd = m_uInfo->template end<elem_type>(si);

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
					uValue[s]=aaU[sides[s]];
					aaVol[sides[s]] += scv.volume();
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
						}
					}
					aaDefTensor[sides[scvf.from()]]+=ipDefTensorFlux;
					aaDefTensor[sides[scvf.to()]]-=ipDefTensorFlux;
				}
			}
	}
	// average
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
		if (si>0) continue;
		SideIterator sideIter = m_uInfo->template begin<side_type>(si);
		SideIterator sideIterEnd = m_uInfo->template end<side_type>(si);
		for(  ;sideIter !=sideIterEnd; sideIter++)
		{
			side_type* elem = *sideIter;
			MathVector<dim> posCo;
			for (int d=0;d<dim;d++)
							posCo[d] = 0.5*posAcc[elem->vertex(0)][d] + 0.5*posAcc[elem->vertex(1)][d];
			UG_LOG(" c=" << posCo << "\n");
			for (int d1=0;d1<dim;d1++)
				for (int d2=0;d2<dim;d2++)
						UG_LOG(" tensor(" << d1 << "," << d2 << ")=" << aaDefTensor[elem][d1][d2] << "\n");
			aaDefTensor[elem]/=(number)aaVol[elem];
		}
	}
}

// Frobenius norm of dim x dim matrix
template <typename TData, int dim, typename TImpl,typename TGridFunction>
number StdTurbulentViscosityData<TData,dim,TImpl,TGridFunction>::FNorm(MathMatrix<dim,dim> M){
	number norm=0;
	for (int d1=0;d1<dim;d1++)
		for (int d2=0;d2<dim;d2++){
			norm += 2 * M[d1][d2] * M[d1][d2];
		}
	return sqrt(norm);
}

template <typename TData, int dim, typename TImpl,typename TGridFunction>
void StdTurbulentViscosityData<TData,dim,TImpl,TGridFunction>::addUiUjTerm(aSideTensor& aaResult,const number factor,aSideDimVector aaU){
	//	get domain of grid function
	domain_type& domain = *m_uInfo->domain().get();

	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
		if (si>0) continue;
		SideIterator sideIter = m_uInfo->template begin<side_type>(si);
		SideIterator sideIterEnd = m_uInfo->template end<side_type>(si);
		for(  ;sideIter !=sideIterEnd; sideIter++)
		{
			side_type* elem = *sideIter;
			dimMat Tij;
			for (int d1=0;d1 < dim;d1++)
				for (int d2=0;d2 < dim;d2++)
					Tij[d1][d2] = aaU[elem][d1]*aaU[elem][d2];
			Tij*=factor;
			aaResult[elem]+=Tij;
		}
	}
}

template <typename TData, int dim, typename TImpl,typename TGridFunction>
void StdTurbulentViscosityData<TData,dim,TImpl,TGridFunction>::addUiUjTerm(aSideTensor& aaResult,const number factor,SmartPtr<TGridFunction> u){
	//	get domain of grid function
	domain_type& domain = *m_uInfo->domain().get();
	//	create Multiindex
	std::vector<MultiIndex<2> > multInd;

	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
		if (si>0) continue;
		SideIterator sideIter = m_uInfo->template begin<side_type>(si);
		SideIterator sideIterEnd = m_uInfo->template end<side_type>(si);
		for(  ;sideIter !=sideIterEnd; sideIter++)
		{
			side_type* elem = *sideIter;
			dimMat Tij;
			MathVector<dim> uValue;
			for (int d=0;d<dim;d++){
				//	get indices of function fct on vertex
				u->multi_indices(elem, d, multInd);
				//	read value of index from vector
				uValue[d]=DoFRef(*u,multInd[0]);
			}
			for (int d1=0;d1 < dim;d1++)
				for (int d2=0;d2 < dim;d2++)
					Tij[d1][d2] = uValue[d1]*uValue[d2];
			Tij*=factor;
			aaResult[elem]+=Tij;
		}
	}
}

template<typename TGridFunction>
void CRSmagorinskyTurbViscData<TGridFunction>::update(){
	//	get domain of grid function
	domain_type& domain = *m_u->domain().get();

	//	coord and vertex array
//	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
//	VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	// assemble deformation tensor fluxes
	assembleDeformationTensor(m_acDeformation,m_acVolume,m_u);
	// compute turbulent viscosity , loop over sides
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		SideIterator sideIter = m_u->template begin<side_type>(si);
		SideIterator sideIterEnd = m_u->template end<side_type>(si);
		for(  ;sideIter !=sideIterEnd; sideIter++)
		{
			//	get Elem
			side_type* elem = *sideIter;
			if (m_pbm && m_pbm->is_slave(elem)) continue;
			number delta = m_acVolume[elem];
			// for possible other choices of delta see Fršhlich p 160
			delta = pow(delta,(number)1.0/(number)dim);
			number tensorNorm = this->FNorm(m_acDeformation[elem]);
			m_acTurbulentViscosity[elem] = m_c * delta*delta * tensorNorm;
		}
	}
	// transfer attachment data to lower levels
	transferToLowerLevels(m_acTurbulentViscosity,*m_spApproxSpace);
}

template<typename TGridFunction>
void CRDynamicTurbViscData<TGridFunction>::update(){
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
	SetAttachmentValues(m_acTurbulentC,m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);
	SetAttachmentValues(m_acVolumeHat,m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);
	SetAttachmentValues(m_acUHat,m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);
	SetAttachmentValues(m_acDeformationNorm,m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);
	SetAttachmentValues(m_acDeformationHat,m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);
	SetAttachmentValues(m_acLij,m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);
	SetAttachmentValues(m_acMij,m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);

	// compute Lij term \hat{u_i u_j} - \hat{u_i} \hat{u_j}
	// \hat{u}
	elementFilter(m_acUHat,m_acVolumeHat,m_u);
	// use Mij attachment to store first Lij part
	// u_i u_j
	addUiUjTerm(m_acMij,1.0,m_u);
	// \hat{u_i u_j}
	elementFilter(m_acLij,m_acVolumeHat,m_acMij);
	// \hat{u_i u_j} - \hat{u_i} \hat{u_j}
	addUiUjTerm(m_acLij,-1.0,m_acUHat);

	// Mij term
	// first term |\hat{S}| \hat{S}
	// assemble \hat{S} using \hat{u}
	assembleDeformationTensor(m_acDeformationHat,m_acVolume,m_acUHat);
	// normalize \hat{S}
	normalizeTensor(m_acDeformationHat);
	// Mij second term \hat{|S|S}
	// compute S
	assembleDeformationTensor(m_acDeformation,m_acVolumeHat,m_u);
	// normalize S
	normalizeTensor(m_acDeformation);
	// filter |S| S
	elementFilter(m_acMij,m_acVolumeHat,m_acDeformation);

	//	create Multiindex
	std::vector<MultiIndex<2> > multInd;

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	// assemble deformation tensor fluxes
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		if (si>0) continue;
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
			MathVector<dim> bary,localBary;
			bary = 0;

			for(size_t i = 0; i < numVertices; ++i){
				vVrt[i] = elem->vertex(i);
				coCoord[i] = posAcc[vVrt[i]];
				bary += coCoord[i];
			};
			bary /= numVertices;

			//	reference object id
			ReferenceObjectID roid = elem->reference_object_id();

			//	get trial space
			const LocalShapeFunctionSet<dim>& rTrialSpace =
				LocalShapeFunctionSetProvider::get<dim>(roid, LFEID(LFEID::CROUZEIX_RAVIART, 1));

			//	get Reference Mapping
			DimReferenceMapping<dim, dim>& map = ReferenceMappingProvider::get<dim, dim>(roid, coCoord);

			map.global_to_local(localBary,bary);

			//	memory for shapes
			std::vector<number> vShape;

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

			number elementVolume = 0;
			MathVector<dim> baryValue;
			baryValue = 0;

			rTrialSpace.shapes(vShape, localBary);

			for (size_t s=0;s < nofsides;s++)
			{
				const typename DimCRFVGeometry<dim>::SCV& scv = geo.scv(s);
				for (int d=0;d<dim;d++){
					//	get indices of function fct on vertex
					m_u->multi_indices(sides[s], d, multInd);
					//	read value of index from vector
					uValue[s][d]=DoFRef(*m_u,multInd[0]);
					baryValue[d] += vShape[s] * uValue[s][d];
				}
				// UG_LOG("scv.volume(" << s << ")=" << scv.volume() << "\n");
				m_acVolume[sides[s]] += scv.volume();
				elementVolume += scv.volume();
			}
			dimMat Lij;
			// first part of Lij term: \hat{u_i u_j}
			for (int d1=0;d1 < dim;d1++)
				for (int d2=0;d2 < dim;d2++)
					Lij = baryValue[d1] * baryValue[d2];
			baryValue*=elementVolume;
			Lij*=elementVolume;
			// second filter volume, second filter function uHat and first term of Lij tensor
			for (size_t s=0;s < nofsides;s++)
			{
				m_acVolumeHat[sides[s]] += elementVolume;
				m_acUHat[sides[s]] += baryValue;
				m_acLij[sides[s]] += Lij;
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
					}
				}
				m_acDeformation[sides[scvf.from()]]+=ipDefTensorFlux;
				m_acDeformation[sides[scvf.to()]]-=ipDefTensorFlux;
			}
		}
		// average uHat and Lij first term part, multiply deformation tensor with its norm
		SideIterator sideIter = m_u->template begin<side_type>(si);
		SideIterator sideIterEnd = m_u->template end<side_type>(si);
		for(  ;sideIter !=sideIterEnd; sideIter++)
		{
			//	get Elem
			number tensorNorm=0;
			side_type* elem = *sideIter;
			dimMat defTensor = m_acDeformation[elem];
			number delta = m_acVolume[elem];
			// complete deformation tensor computation by averaging
			defTensor/=delta;
			// compute norm of tensor
			for (int d1=0;d1<dim;d1++)
				for (int d2=0;d2<dim;d2++){
					// UG_LOG("tensor(" << d1 << "," << d2 << ")=" << defTensor[d1][d2] << "\n");
					tensorNorm += 2 * defTensor[d1][d2] * defTensor[d1][d2];
				}
			tensorNorm = sqrt(tensorNorm);
			defTensor *= tensorNorm;
			m_acDeformationNorm[elem] = tensorNorm;
			m_acDeformation[elem] = defTensor;
			number deltaHat = m_acVolumeHat[elem];
			m_acUHat[elem]/=deltaHat;
			m_acLij[elem]/=deltaHat;
		}
		// iterate over elements again and complete computation of Lij tensor, assemble hat deformation tensor
		iter = m_u->template begin<elem_type>(si);
		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			elem_type* elem = *iter;
			//	get vertices and extract corner coordinates
			const size_t numVertices = elem->num_vertices();
			MathVector<dim> bary,localBary;
			bary = 0;

			for(size_t i = 0; i < numVertices; ++i){
				vVrt[i] = elem->vertex(i);
				coCoord[i] = posAcc[vVrt[i]];
				bary += coCoord[i];
			};
			bary /= numVertices;

			//	reference object id
			ReferenceObjectID roid = elem->reference_object_id();

			//	get trial space
			const LocalShapeFunctionSet<dim>& rTrialSpace =
				LocalShapeFunctionSetProvider::get<dim>(roid, LFEID(LFEID::CROUZEIX_RAVIART, 1));

			//	get Reference Mapping
			DimReferenceMapping<dim, dim>& map = ReferenceMappingProvider::get<dim, dim>(roid, coCoord);

			map.global_to_local(localBary,bary);

			//	memory for shapes
			std::vector<number> vShape;

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

			number elementVolume = 0;
			dimMat baryTensorValue;
			baryTensorValue = 0;

			rTrialSpace.shapes(vShape, localBary);

			for (size_t s=0;s < nofsides;s++)
			{
				const typename DimCRFVGeometry<dim>::SCV& scv = geo.scv(s);
				for (int d1=0;d1<dim;d1++){
					for (int d2=0;d2<dim;d2++)
						baryTensorValue[d1][d2] += vShape[s] * m_acDeformation[sides[s]][d1][d2];
				}
				// UG_LOG("scv.volume(" << s << ")=" << scv.volume() << "\n");
				elementVolume += scv.volume();
			}

			baryTensorValue*=elementVolume;

			// complete Lij computation
			for (size_t s=0;s < nofsides;s++)
			{
				dimMat Lij;
				for (int d1=0;d1 < dim;d1++)
					for (int d2=0;d2 < dim;d2++)
						Lij[d1][d2] = m_acUHat[sides[s]][d1]*m_acUHat[sides[s]][d2];
				m_acLij[sides[s]] -= Lij;
				m_acMij[sides[s]] += baryTensorValue;
			}
			// compute second filter deformation tensor from \hat{u} values
			for (size_t ip=0;ip<nip;ip++){
				// 	get current SCVF
				ipVelocity = 0;
				const typename DimCRFVGeometry<dim>::SCVF& scvf = geo.scvf(ip);
				for (size_t s=0;s < nofsides;s++){
					for (int d=0;d<dim;d++){
					    ipVelocity[d] += scvf.shape(s)*m_acUHat[sides[s]][d];
					};
				};
				dimMat ipDefTensorFlux;
				ipDefTensorFlux = 0;
				for (int d=0;d<dim;d++){
					for (int j=0;j<d;j++){
						ipDefTensorFlux[d][j]= 0.5 * (ipVelocity[d] * scvf.normal()[j] + ipVelocity[j] * scvf.normal()[d]);
					}
				}
				m_acDeformationHat[sides[scvf.from()]]+=ipDefTensorFlux;
				m_acDeformationHat[sides[scvf.to()]]-=ipDefTensorFlux;
			}
		}
		sideIter = m_u->template begin<side_type>(si);
		for(  ;sideIter !=sideIterEnd; sideIter++)
		{
			//	get Elem
			number tensorNorm=0;
			side_type* elem = *sideIter;
			dimMat defTensor = m_acDeformationHat[elem];
			number deltaHat = m_acVolumeHat[elem];
			// compute second Mij term by averaging
			m_acMij[elem] /= 2*deltaHat;
			// complete deformation tensor computation by averaging
			defTensor/=deltaHat;
			// compute norm of tensor
			for (int d1=0;d1<dim;d1++)
				for (int d2=0;d2<dim;d2++){
					// UG_LOG("tensor(" << d1 << "," << d2 << ")=" << defTensor[d1][d2] << "\n");
					tensorNorm += 2 * defTensor[d1][d2] * defTensor[d1][d2];
				}
			tensorNorm = sqrt(tensorNorm);
			// scale tensor with norm
			defTensor *= 2*tensorNorm;
			m_acDeformationHat[elem] = defTensor;
			// add first Mij term
			m_acMij[elem] -= m_acDeformationHat[elem];
			// compute local c
			number c = 0;
			for (int d1=0;d1<dim;d1++)
				for (int d2=0;d2<dim;d2++)
					c += m_acLij[elem][d1][d2]*m_acMij[elem][d1][d2];
			number denom=0;
			for (int d1=0;d1<dim;d1++)
				for (int d2=0;d2<dim;d2++)
					denom += m_acMij[elem][d1][d2]*m_acMij[elem][d1][d2];
			if (denom>1e-15)
				c/=(number)denom;
			else c=0;
			number delta = m_acVolume[elem];
			// for possible other choices of delta see Fršhlich p 160
			delta = pow(delta,(number)1.0/(number)dim);
			// UG_LOG("c = " << c << " delta = " << delta << " deformNorm = " << m_acDeformationNorm[elem] << " denom = " << denom << "\n");
			m_acTurbulentViscosity[elem] = c * delta*delta * m_acDeformationNorm[elem];
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
		side_type* s;
		m_acTurbulentViscosity[s] = 5;
	}
}

} // namespace NavierStokes
} // end namespace ug

#endif /* __H__UG__NAVIER_STOKES_TURBULENT_VISCOSITY_DATA_IMPL_H_ */
