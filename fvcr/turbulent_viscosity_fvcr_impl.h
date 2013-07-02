/*
 * turbulent_viscosity_data_impl.h
 *
 *  Created on: 01.11.2012
 *      Author: Christian Wehner
 */

#ifndef __H__UG__NAVIER_STOKES_TURBULENT_VISCOSITY_DATA_IMPL_H_
#define __H__UG__NAVIER_STOKES_TURBULENT_VISCOSITY_DATA_IMPL_H_

namespace ug{
namespace NavierStokes{

template <typename TData, int dim, typename TImpl,typename TGridFunction>
void StdTurbulentViscosityData<TData,dim,TImpl,TGridFunction>::transferToLowerLevels(aSideNumber& aaData,ApproximationSpace<domain_type>& approximationSpace){
	for(int si = 0; si < approximationSpace.domain()->subset_handler()->num_subsets(); ++si){
		// transfer to lower levels, averaging over child edges (2d) / child faces (3d)
		for (size_t lev=approximationSpace.num_levels()-2;(int)lev>=0;lev--){
			const DoFDistribution& lDD = *approximationSpace.level_dof_distribution(lev);
			const MultiGrid& grid = lDD.multi_grid();
			typedef typename DoFDistribution::traits<side_type>::const_iterator coarseLevelSideIter;
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
// parameters: filtered values, filter volume (computed in function), original values
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
			LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::CROUZEIX_RAVIART, dim, 1));

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
				//for debug UG_LOG(localValue << "\n");
				localValue *= vShape[s];
				value += localValue;
				elementVolume += scv.volume();
			}
			value *= elementVolume;
			for (size_t s=0;s<nofsides;s++){
				aaUHat[sides[s]] += value;
				aaVol[sides[s]] += elementVolume;
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
			side_type* side = *sideIter;
			if (pbm && pbm->is_slave(side)) continue;
			aaUHat[side]/=(number)aaVol[side];
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
			LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::CROUZEIX_RAVIART, dim, 1));

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
					u->multi_indices(sides[s], d, multInd);
					localValue[d]=DoFRef(*u,multInd[0]);
				}
				//for debug UG_LOG("localValue=" << localValue << "\n");
				//for debug UG_LOG("vShape=" << vShape[s] << "\n");
				localValue *= vShape[s];
				value += localValue;
				elementVolume += scv.volume();
			}
			//for debug UG_LOG("value=" << value << " vol=" << elementVolume << "\n");
			value *= elementVolume;
			for (size_t s=0;s<nofsides;s++){
				aaVol[sides[s]]+=elementVolume;
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
			side_type* side = *sideIter;
			if (pbm && pbm->is_slave(side)) continue;
			//for debug UG_LOG("# " << aaUHat[side] << " " << aaVol[side] << "\n");
			aaUHat[side]/=(number)aaVol[side];
		}
	}
}

// go over all elements, interpolate data to scv barycenter, average by multiplying with corresponding scv volume and deviding by volume of complete control volume
template <typename TData, int dim, typename TImpl,typename TGridFunction>
template <typename VType>
void StdTurbulentViscosityData<TData,dim,TImpl,TGridFunction>::scvFilter(PeriodicAttachmentAccessor<side_type,Attachment<VType> >& aaUHat,aSideNumber& aaVol,const PeriodicAttachmentAccessor<side_type,Attachment<VType> >& aaU){
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
			LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::CROUZEIX_RAVIART, dim, 1));

			//	get Reference Mapping
			DimReferenceMapping<dim, dim>& map = ReferenceMappingProvider::get<dim, dim>(roid, coCoord);

			map.global_to_local(localBary,bary);

			typename grid_type::template traits<side_type>::secure_container sides;

			UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Only elements of type elem_type are currently supported");

			domain.grid()->associated_elements_sorted(sides, static_cast<elem_type*>(elem) );

			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

			size_t nofsides = geo.num_scv();

			MathVector<dim> scvBary,scvLocalBary;
			for (size_t s=0;s<nofsides;s++){
				const typename DimCRFVGeometry<dim>::SCV& scv = geo.scv(s);
				scvBary = 0;
				size_t numScvVertices = sides[s]->num_vertices();
				// compute barycenter of scv (average of side corner nodes + element bary)
				for (size_t i=0;i<numScvVertices;i++){
					scvBary += posAcc[sides[s]->vertex(i)];
				}
				scvBary += bary;
				scvBary/=(numScvVertices+1);
				map.global_to_local(scvLocalBary,scvBary);
				//	memory for shapes
				std::vector<number> vShape;
				rTrialSpace.shapes(vShape, scvLocalBary);
				VType localValue = 0;

				for (size_t j=0;j<nofsides;j++){
					localValue += vShape[j]*aaU[sides[j]];
				}

				localValue *= scv.volume();
				aaVol[sides[s]]  += scv.volume();
				aaUHat[sides[s]] += localValue;
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
			side_type* side = *sideIter;
			if (pbm && pbm->is_slave(side)) continue;
			aaUHat[side]/=(number)aaVol[side];
		}
	}
}

// go over all elements, interpolate data to scv barycenter, average by multiplying with corresponding scv volume and deviding by volume of complete control volume
template <typename TData, int dim, typename TImpl,typename TGridFunction>
void StdTurbulentViscosityData<TData,dim,TImpl,TGridFunction>::scvFilter(aSideDimVector& aaUHat,aSideNumber& aaVol,SmartPtr<TGridFunction> u){
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

	//	create Multiindex
	std::vector<MultiIndex<2> > multInd;

	// assemble deformation tensor fluxes
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
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
			LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::CROUZEIX_RAVIART, dim, 1));

			//	get Reference Mapping
			DimReferenceMapping<dim, dim>& map = ReferenceMappingProvider::get<dim, dim>(roid, coCoord);

			map.global_to_local(localBary,bary);

			typename grid_type::template traits<side_type>::secure_container sides;

			UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Only elements of type elem_type are currently supported");

			domain.grid()->associated_elements_sorted(sides, static_cast<elem_type*>(elem) );

			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

			size_t nofsides = geo.num_scv();

			static const size_t MaxNumSidesOfElem = 10;

			typedef MathVector<dim> MVD;
			std::vector<MVD> uValue(MaxNumSidesOfElem);

			for (size_t s=0;s<nofsides;s++){
				for (int d=0;d<dim;d++){
					u->multi_indices(sides[s], d, multInd);
					uValue[s][d]=DoFRef(*u,multInd[0]);
				}
			};

			MathVector<dim> scvBary,scvLocalBary;
			for (size_t s=0;s<nofsides;s++){
				const typename DimCRFVGeometry<dim>::SCV& scv = geo.scv(s);
				scvBary = 0;
				size_t numScvVertices = sides[s]->num_vertices();
				// compute barycenter of scv (average of side corner nodes + element bary)
				for (size_t i=0;i<numScvVertices;i++){
					scvBary += posAcc[sides[s]->vertex(i)];
				}
				scvBary += bary;
				scvBary/=(numScvVertices+1);
				map.global_to_local(scvLocalBary,scvBary);
				//	memory for shapes
				std::vector<number> vShape;
				rTrialSpace.shapes(vShape, scvLocalBary);
				MathVector<dim> localValue = 0;

				for (size_t j=0;j<nofsides;j++)
					for (int d=0;d<dim;d++)
//						localValue[d] += vShape[j]*uValue[j][d];

				localValue *= scv.volume();
				aaVol[sides[s]]  += scv.volume();
				aaUHat[sides[s]] += localValue;
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
			side_type* side = *sideIter;
			if (pbm && pbm->is_slave(side)) continue;
			aaUHat[side]/=(number)aaVol[side];
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

	// add boundary subsets to enforce boundary subset computations in geo.update()
	for(size_t i = 0; i < this->m_turbZeroSg.size(); ++i){
		geo.add_boundary_subset(this->m_turbZeroSg[i]);
	}

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
						u->multi_indices(sides[s], d, multInd);
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
					for (int i=0;i<dim;i++){
						for (int j=0;j<dim;j++){
							ipDefTensorFlux[i][j]= 0.5 * (ipVelocity[i] * scvf.normal()[j] + ipVelocity[j] * scvf.normal()[i]);
						}
					}
					aaDefTensor[sides[scvf.from()]]+=ipDefTensorFlux;
					aaDefTensor[sides[scvf.to()]]-=ipDefTensorFlux;
				}
				for(size_t sgi = 0; sgi < this->m_turbZeroSg.size(); ++sgi){
					const size_t sgsi=this->m_turbZeroSg[sgi];
					if (geo.num_bf(sgsi) == 0) continue;
					for(size_t bfi = 0; bfi < geo.num_bf(sgsi); ++bfi){
						const typename DimCRFVGeometry<dim>::BF& bf = geo.bf(sgsi, bfi);
						const size_t sideID = bf.node_id();
						for (int i=0;i<dim;i++)
							for (int j=0;j<dim;j++){
								// bf ip and u position are identical for CR-FV-Geometry
								//for debug UG_LOG("[" << i << "," << j << "]" <<  0.5 * (uValue[sideID][i] * bf.normal()[j] + uValue[sideID][j] * bf.normal()[i]) << "\n");
								aaDefTensor[sides[sideID]][i][j] += 0.5 * (uValue[sideID][i] * bf.normal()[j] + uValue[sideID][j] * bf.normal()[i]);
							}
					}
				}
			}
	}
	// average
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
		SideIterator sideIter = u->template begin<side_type>(si);
		SideIterator sideIterEnd = u->template end<side_type>(si);
		//for debug UG_LOG("|||||||||||||||||||||||||||||||||||\n");
		//for debug UG_LOG("si=" << si << "\n");
		for(  ;sideIter !=sideIterEnd; sideIter++)
		{
			side_type* side = *sideIter;
			if (pbm && pbm->is_slave(side)) continue;
			aaDefTensor[side]/=(number)aaVol[side];
			//for debug UG_LOG("$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
			MathVector<dim> posCo;
						for (int d=0;d<dim;d++)
										posCo[d] = 0.5*posAcc[side->vertex(0)][d] + 0.5*posAcc[side->vertex(1)][d];
						//for debug UG_LOG(" c=" << posCo << "\n");
		//for debug				for (int d1=0;d1<dim;d1++)
		//for debug					for (int d2=0;d2<dim;d2++)
									//for debug UG_LOG(" tensor(" << d1 << "," << d2 << ")=" << aaDefTensor[side][d1][d2] << "\n");
			//for debug UG_LOG(" norm=" << FNorm(aaDefTensor[side]) << "\n");
		}
	}
}

template <typename TData, int dim, typename TImpl,typename TGridFunction>
void StdTurbulentViscosityData<TData,dim,TImpl,TGridFunction>::scaleTensorByNorm(aSideTensor& aaTensor){
	//	get domain of grid function
	domain_type& domain = *m_uInfo->domain().get();
	// get periodic boundary manager
	PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
	// average
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
		SideIterator sideIter = m_uInfo->template begin<side_type>(si);
		SideIterator sideIterEnd = m_uInfo->template end<side_type>(si);
		for(  ;sideIter !=sideIterEnd; sideIter++)
		{
			side_type* side = *sideIter;
			if (pbm && pbm->is_slave(side)) continue;
			//for debug UG_LOG("&&&&&&&& norm = " << FNorm(aaTensor[side]) << "\n");
			aaTensor[side]*=(number)FNorm(aaTensor[side]);
		}
	}
}

template <typename TData, int dim, typename TImpl,typename TGridFunction>
void StdTurbulentViscosityData<TData,dim,TImpl,TGridFunction>::assembleDeformationTensor(aSideTensor& aaDefTensor,aSideNumber& aaVol,aSideDimVector aaU){
	//	get domain of grid function
	domain_type& domain = *m_uInfo->domain().get();
	// get periodic boundary manager
	PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
	//	create Multiindex
	std::vector<MultiIndex<2> > multInd;

	DimCRFVGeometry<dim> geo;

	// add boundary subsets to enforce boundary subset computations in geo.update()
	for(size_t i = 0; i < this->m_turbZeroSg.size(); ++i){
		geo.add_boundary_subset(this->m_turbZeroSg[i]);
	}

	SetAttachmentValues(aaDefTensor , m_uInfo->template begin<side_type>(), m_uInfo->template end<side_type>(), 0);
	SetAttachmentValues(aaVol , m_uInfo->template begin<side_type>(), m_uInfo->template end<side_type>(), 0);

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();

	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
			//	get iterators
			ElemIterator iter = m_uInfo->template begin<elem_type>(si);
			ElemIterator iterEnd = m_uInfo->template end<elem_type>(si);

			//for debug UG_LOG("|||||||||||||||||||| si = " << si << "\n");

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
					// //for debug UG_LOG("co_coord(" << i<< "+1,:)=" << coCoord[i] << "\n");
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
					//for debug UG_LOG("uvalue(" << s << ")=" << uValue[s] << "\n");
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
				for(size_t sgi = 0; sgi < this->m_turbZeroSg.size(); ++sgi){
					const size_t sgsi=this->m_turbZeroSg[sgi];
					if (geo.num_bf(sgsi) == 0) continue;
					for(size_t bfi = 0; bfi < geo.num_bf(sgsi); ++bfi){
						const typename DimCRFVGeometry<dim>::BF& bf = geo.bf(sgsi, bfi);
						const size_t sideID = bf.node_id();
						for (int i=0;i<dim;i++)
							for (int j=0;j<dim;j++){
								// bf ip and u position are identical for CR-FV-Geometry
								//for debug UG_LOG("[" << i << "," << j << "]" <<  0.5 * (uValue[sideID][i] * bf.normal()[j] + uValue[sideID][j] * bf.normal()[i]) << "\n");
								aaDefTensor[sides[sideID]][i][j] += 0.5 * (uValue[sideID][i] * bf.normal()[j] + uValue[sideID][j] * bf.normal()[i]);
						}
					}
				}
			}
	}
	// average
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
		SideIterator sideIter = m_uInfo->template begin<side_type>(si);
		SideIterator sideIterEnd = m_uInfo->template end<side_type>(si);
		for(  ;sideIter !=sideIterEnd; sideIter++)
		{
			side_type* side = *sideIter;
			if (pbm && pbm->is_slave(side)) continue;
			aaDefTensor[side]/=(number)aaVol[side];
/*					//for debug UG_LOG("$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
						MathVector<dim> posCo;
									for (int d=0;d<dim;d++)
													posCo[d] = 0.5*posAcc[side->vertex(0)][d] + 0.5*posAcc[side->vertex(1)][d];
									//for debug UG_LOG(" c=" << posCo << "\n");
									for (int d1=0;d1<dim;d1++)
										for (int d2=0;d2<dim;d2++)
												//for debug UG_LOG(" tensor(" << d1 << "," << d2 << ")=" << aaDefTensor[side][d1][d2] << "\n");
						//for debug UG_LOG(" norm=" << FNorm(aaDefTensor[side]) << "\n");*/
		}
	}
}

// Frobenius norm of dim x dim matrix
template <typename TData, int dim, typename TImpl,typename TGridFunction>
number StdTurbulentViscosityData<TData,dim,TImpl,TGridFunction>::FNorm(MathMatrix<dim,dim> M){
	number norm=0;
	for (int d1=0;d1<dim;d1++)
		for (int d2=0;d2<dim;d2++){
			norm +=  M[d1][d2] * M[d1][d2];
		}
	return sqrt(2.0*norm);
}

template <typename TData, int dim, typename TImpl,typename TGridFunction>
void StdTurbulentViscosityData<TData,dim,TImpl,TGridFunction>::addUiUjTerm(aSideTensor& aaResult,const number factor,aSideDimVector aaU){
	//	get domain of grid function
	domain_type& domain = *m_uInfo->domain().get();
	// get periodic boundary manager
	PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
		SideIterator sideIter = m_uInfo->template begin<side_type>(si);
		SideIterator sideIterEnd = m_uInfo->template end<side_type>(si);
		for(  ;sideIter !=sideIterEnd; sideIter++)
		{
			side_type* side = *sideIter;
			if (pbm && pbm->is_slave(side)) continue;
			dimMat Tij;
			for (int d1=0;d1 < dim;d1++)
				for (int d2=0;d2 < dim;d2++)
					Tij[d1][d2] = aaU[side][d1]*aaU[side][d2];
			Tij*=factor;
			aaResult[side]+=Tij;
		}
	}
}

template <typename TData, int dim, typename TImpl,typename TGridFunction>
void StdTurbulentViscosityData<TData,dim,TImpl,TGridFunction>::addUiUjTerm(aSideTensor& aaResult,const number factor,SmartPtr<TGridFunction> u){
	//	get domain of grid function
	domain_type& domain = *m_uInfo->domain().get();
	// get periodic boundary manager
	PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
	//	create Multiindex
	std::vector<MultiIndex<2> > multInd;

	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
		SideIterator sideIter = m_uInfo->template begin<side_type>(si);
		SideIterator sideIterEnd = m_uInfo->template end<side_type>(si);
		for(  ;sideIter !=sideIterEnd; sideIter++)
		{
			side_type* side = *sideIter;
			if (pbm && pbm->is_slave(side)) continue;
			dimMat Tij;
			MathVector<dim> uValue;
			for (int d=0;d<dim;d++){
				u->multi_indices(side, d, multInd);
				uValue[d]=DoFRef(*u,multInd[0]);
			}
			for (int d1=0;d1 < dim;d1++)
				for (int d2=0;d2 < dim;d2++)
					Tij[d1][d2] = uValue[d1]*uValue[d2];
			Tij*=factor;
			aaResult[side]+=Tij;
		}
	}
}

template<typename TGridFunction>
void CRSmagorinskyTurbViscData<TGridFunction>::update(){
	//	get domain of grid function
	domain_type& domain = *m_u->domain().get();
	SetAttachmentValues(m_acTurbulentViscosity, m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);

	//	coord and vertex array
//	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
//	VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	// assemble deformation tensor fluxes
	this->assembleDeformationTensor(m_acDeformation,m_acVolume,m_u);
	// compute turbulent viscosity , loop over sides
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		if ((this->m_turbZeroSg.size()!=0) && (this->m_turbZeroSg.contains(si)==true)) continue;
		SideIterator sideIter = m_u->template begin<side_type>(si);
		SideIterator sideIterEnd = m_u->template end<side_type>(si);
		bool periodic_subset = false;
		for(  ;sideIter !=sideIterEnd; sideIter++)
		{
			//	get Elem
			side_type* side = *sideIter;
			if (m_pbm && m_pbm->is_slave(side)){
				periodic_subset = true;
				continue;
			}
			number delta = m_acVolume[side];
			// for possible other choices of delta see Fr�hlich p 160
			delta = pow(delta,(number)1.0/(number)dim);
			number tensorNorm = this->FNorm(m_acDeformation[side]);
			m_acTurbulentViscosity[side] = m_c * delta*delta * tensorNorm;
		}
		// fill slave elements
		if (periodic_subset==true){
			for(  ;sideIter !=sideIterEnd; sideIter++)
			{
				//	get Elem
				side_type* side = *sideIter;
				if (m_pbm && m_pbm->is_slave(side)){

				}
			}
		}
	}
	// transfer attachment data to lower levels
	this->transferToLowerLevels(m_acTurbulentViscosity,*m_spApproxSpace);
}

template<typename TGridFunction>
void CRDynamicTurbViscData<TGridFunction>::update(){
	//	get domain of grid function
	domain_type& domain = *m_u->domain().get();

	//	get position accessor
	// for debug typedef typename domain_type::position_accessor_type position_accessor_type;
	// for debug const position_accessor_type& posAcc = domain.position_accessor();

	// initialize attachment values with 0
//	SetAttachmentValues(m_acDeformation , m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);
	SetAttachmentValues(m_acTurbulentViscosity, m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);
//	SetAttachmentValues(m_acVolume,m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);
	SetAttachmentValues(m_acTurbulentCNew,m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);
//	SetAttachmentValues(m_acVolumeHat,m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);
//	SetAttachmentValues(m_acUHat,m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);
//	SetAttachmentValues(m_acDeformationHat,m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);
//	SetAttachmentValues(m_acLij,m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);
	SetAttachmentValues(m_acMij,m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);

	// compute Lij term \hat{u_i u_j} - \hat{u_i} \hat{u_j}
	// \hat{u}
	this->elementFilter(m_acUHat,m_acVolumeHat,m_u);
	// use Mij attachment to store first Lij part
	// u_i u_j
	this->addUiUjTerm(m_acMij,1.0,m_u);
	// \hat{u_i u_j}
	this->elementFilter(m_acLij,m_acVolumeHat,m_acMij);
	// \hat{u_i u_j} - \hat{u_i} \hat{u_j}
	this->addUiUjTerm(m_acLij,-1.0,m_acUHat);

	// Mij term
	// first term |\hat{S}| \hat{S}
	// assemble \hat{S} using \hat{u}
	this->assembleDeformationTensor(m_acDeformationHat,m_acVolume,m_acUHat);
	// normalize \hat{S}
	this->scaleTensorByNorm(m_acDeformationHat);
	// Mij second term \hat{|S|S}
	// compute S
	this->assembleDeformationTensor(m_acDeformation,m_acVolume,m_u);
	// compute |S| S
	this->scaleTensorByNorm(m_acDeformation);
	// filter |S| S
	//for debug UG_LOG("------------------------------------------------------\n");
	this->elementFilter(m_acMij,m_acVolumeHat,m_acDeformation);

	bool use_filter = false;

	//	create Multiindex
	std::vector<MultiIndex<2> > multInd;

	// complete Mij term computation by scaling and adding the two terms,
	// solve the local least squares problem and compute local c and local turbulent viscosity
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		if (use_filter==false)
			if ((this->m_turbZeroSg.size()!=0) && (this->m_turbZeroSg.contains(si)==true)) continue;
		SideIterator sideIter = m_u->template begin<side_type>(si);
		SideIterator sideIterEnd = m_u->template end<side_type>(si);
		for(  ;sideIter !=sideIterEnd; sideIter++){
			side_type* side = *sideIter;
			if (m_pbm && m_pbm->is_slave(side)) continue;
			// use c to compute turbulent viscosity
			number delta = m_acVolume[side];
			delta = pow(delta,(number)1.0/(number)dim);
			number deltaHat = m_acVolumeHat[side];
			deltaHat = pow(deltaHat,(number)1.0/(number)dim);
			m_u->multi_indices(side, 0, multInd);
			//for debug UG_LOG("--------------------------------------------------");
			//for debug UG_LOG("co=[" << 0.5*(posAcc[side->vertex(0)][0] + posAcc[side->vertex(1)][0]) << "," << 0.5*(posAcc[side->vertex(0)][1] + posAcc[side->vertex(1)][1]) << "]\n");
			//for debug UG_LOG("u = " << DoFRef(*m_u,multInd[0]) << "\n");
			m_u->multi_indices(side, 1, multInd);
			//for debug UG_LOG("v = " << DoFRef(*m_u,multInd[0]) << "\n");
			//for debug UG_LOG("uHat = " << m_acUHat[side] << "\n");
			//for debug UG_LOG("deformNorm = " << FNorm(m_acDeformation[side]) << "\n");
			//UG_LOG("deformHatNorm = " << FNorm( m_acDeformationHat[side]) << "\n");
			//for debug UG_LOG("hatDeform = " << m_acMij[side] << "\n");
			// complete Mij term computation
			m_acDeformationHat[side] *= -2*deltaHat*deltaHat;
			m_acMij[side] *= 2*delta*delta;
			m_acMij[side] += m_acDeformationHat[side];
			//for debug UG_LOG("Mij " << FNorm(m_acMij[side]) << "\n");
			for (int d1=0;d1<dim;d1++){
				for (int d2=0;d2<dim;d2++){
					//for debug UG_LOG(m_acMij[side][d1][d2] << " ");
				}
				//for debug UG_LOG("\n");
			}
			//for debug UG_LOG("Lij " << FNorm(m_acLij[side]) << "\n");
			for (int d1=0;d1<dim;d1++){
				for (int d2=0;d2<dim;d2++){
					//for debug UG_LOG(m_acLij[side][d1][d2] << " ");
				}
				//for debug UG_LOG("\n");
			}
			// compute local c
			// solve least squares problem
			number c = 0;
			for (int d1=0;d1<dim;d1++)
				for (int d2=0;d2<dim;d2++)
					c += m_acLij[side][d1][d2]*m_acMij[side][d1][d2];
			number denom=0;
			//for debug UG_LOG("c=" << c << "\n");
			for (int d1=0;d1<dim;d1++)
				for (int d2=0;d2<dim;d2++)
					denom += m_acMij[side][d1][d2]*m_acMij[side][d1][d2];
			if (denom>1e-15)
				c/=(number)denom;
			else c=0;

			if (m_spaceFilter==false){
				if (m_timeFilter==false){
					m_acTurbulentViscosity[side] = c * delta*delta * this->FNorm(m_acDeformation[side]);
				} else {
					m_acTurbulentC[side]= (m_timeFilterEps * c + (1-m_timeFilterEps)*m_acTurbulentC[side]);
					m_acTurbulentViscosity[side] = m_acTurbulentC[side] * delta*delta * this->FNorm(m_acDeformation[side]);
				}
				if (m_acTurbulentViscosity[side]+m_viscosityNumber<m_small) m_acTurbulentViscosity[side] = m_viscosityNumber + m_small;			}
			else{
				// store c in viscosity attachment
				m_acTurbulentViscosity[side] = c;
			}
		}
	}
	if (m_spaceFilter==true){
		// filter c
		if (m_timeFilter==false)
			// c has been stored in viscosity attachment
			this->elementFilter(m_acTurbulentC,m_acVolumeHat,m_acTurbulentViscosity);
		else
			this->elementFilter(m_acTurbulentCNew,m_acVolumeHat,m_acTurbulentViscosity);
		// compute turbulent viscosity
		for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
		{
			//for debug UG_LOG("si = " << si << "\n");
			if ((this->m_turbZeroSg.size()!=0) && (this->m_turbZeroSg.contains(si)==true)) continue;
			SideIterator sideIter = m_u->template begin<side_type>(si);
			SideIterator sideIterEnd = m_u->template end<side_type>(si);
			for(  ;sideIter !=sideIterEnd; sideIter++){
				side_type* side = *sideIter;
				if (m_pbm && m_pbm->is_slave(side)) continue;
				number delta = m_acVolume[side];
				delta = pow(delta,(number)1.0/(number)dim);
				if (m_timeFilter==true)
					// time averaging, note that c has been stored in m_acVolumeHat
					m_acTurbulentC[side]= (m_timeFilterEps * m_acTurbulentCNew[side] + (1-m_timeFilterEps)*m_acTurbulentC[side]);
				m_acTurbulentViscosity[side] = m_acTurbulentC[side] * delta*delta * this->FNorm(m_acDeformation[side]);
				if (m_acTurbulentViscosity[side]+m_viscosityNumber<m_small) m_acTurbulentViscosity[side] = m_viscosityNumber+m_small;
				//for debug UG_LOG("nu_t = " << m_acTurbulentViscosity[side]  << " c = " << m_acTurbulentC[side] << " delta = " << delta << " co=[" << 0.5*(posAcc[side->vertex(0)][0] + posAcc[side->vertex(1)][0]) << "," << 0.5*(posAcc[side->vertex(0)][1] + posAcc[side->vertex(1)][1]) << "]\n");
			}
		}
	}
}

} // namespace NavierStokes
} // end namespace ug

#endif /* __H__UG__NAVIER_STOKES_TURBULENT_VISCOSITY_DATA_IMPL_H_ */
