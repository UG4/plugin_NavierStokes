/*
 * filter.h
 *
 *  Created on: 8.04.2013
 *      Author: Christian Wehner
 */

#ifndef __H__UG__NAVIER_STOKES_FILTER__
#define __H__UG__NAVIER_STOKES_FILTER__

#include <vector>
#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_grid/tools/periodic_boundary_manager.h"
#include "lib_grid/lg_base.h"
#include "common/profiler/profiler.h"
#include "lib_disc/spatial_disc/disc_util/fvcr_geom.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_grid/algorithms/attachment_util.h"

namespace ug{

template <int dim,typename elem_type,typename TGridFunction>
void copyAttachmentToGridFunction(SmartPtr<TGridFunction> u,PeriodicAttachmentAccessor<elem_type,Attachment<MathVector<dim> > >& aaU){
	/// element iterator
	typedef typename TGridFunction::template traits<elem_type>::const_iterator ElemIterator;
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;
	//	get domain
	domain_type& domain = *u->domain().get();
	//	create Multiindex
	std::vector<DoFIndex> multInd;
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
		ElemIterator iter = u->template begin<elem_type>(si);
		ElemIterator iterEnd = u->template end<elem_type>(si);
		for(  ;iter !=iterEnd; ++iter)
		{
			elem_type* elem = *iter;
			for (int d=0;d<dim;d++){
				u->dof_indices(elem, d, multInd);
				DoFRef(*u,multInd[0])=aaU[elem][d];
			}
		}
	}
}

template <int dim,typename elem_type,typename TGridFunction>
void copyGridFunctionToAttachment(PeriodicAttachmentAccessor<elem_type,Attachment<MathVector<dim> > >& aaU,SmartPtr<TGridFunction> u){
	/// element iterator
	typedef typename TGridFunction::template traits<elem_type>::const_iterator ElemIterator;
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;
	//	get domain
	domain_type& domain = *u->domain().get();
	//	create Multiindex
	std::vector<DoFIndex> multInd;
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
		ElemIterator iter = u->template begin<elem_type>(si);
		ElemIterator iterEnd = u->template end<elem_type>(si);
		for(  ;iter !=iterEnd; ++iter)
		{
			elem_type* elem = *iter;
			for (int d=0;d<dim;d++){
				u->dof_indices(elem, d, multInd);
				aaU[elem][d]=DoFRef(*u,multInd[0]);
			}
		}
	}
}

template <int dim,typename TGridFunction>
void elementFilterFVCR(SmartPtr<TGridFunction> u){
	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object elem_type;
	/// side type
	typedef typename elem_type::side side_type;
	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;
	/// side iterator
	typedef typename TGridFunction::template traits<side_type>::const_iterator SideIterator;
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;
	///	grid type
	typedef typename domain_type::grid_type grid_type;
	/// position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;

	domain_type& domain = *u->domain().get();
	grid_type& grid = *domain.grid();

	// define attachment stuff
	typedef MathVector<dim> vecDim;
	typedef Attachment<vecDim> AMathVectorDim;
	typedef PeriodicAttachmentAccessor<side_type,ANumber > aSideNumber;
	typedef PeriodicAttachmentAccessor<side_type,AMathVectorDim > aSideDimVector;

	//  filtered u attachment
	aSideDimVector acUHat;
	AMathVectorDim aUHat;
	//  volume attachment
	aSideNumber acVolume;
	ANumber aVolume;

	// attach
	grid.template attach_to<side_type>(aUHat);
	grid.template attach_to<side_type>(aVolume);

	// accessors
	acUHat.access(grid,aUHat);
	acVolume.access(grid,aVolume);

	// initial values for attachments
	SetAttachmentValues(acUHat , u->template begin<side_type>(), u->template end<side_type>(), 0);
	SetAttachmentValues(acVolume , u->template begin<side_type>(), u->template end<side_type>(), 0);

	DimCRFVGeometry<dim> geo;
	std::vector<DoFIndex> multInd;
	const position_accessor_type& posAcc = domain.position_accessor();
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	// assemble fluxes
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
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
					u->dof_indices(sides[s], d, multInd);
					localValue[d]=DoFRef(*u,multInd[0]);
				}
				localValue *= vShape[s];
				value += localValue;
				elementVolume += scv.volume();
			}
			value *= elementVolume;
			for (size_t s=0;s<nofsides;s++){
				acVolume[sides[s]]+=elementVolume;
				acUHat[sides[s]] += value;
			}
		}
	}
	PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
	// average
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		SideIterator sideIter = u->template begin<side_type>(si);
		SideIterator sideIterEnd = u->template end<side_type>(si);
		for(  ;sideIter !=sideIterEnd; sideIter++)
		{
			side_type* side = *sideIter;
			if (pbm && pbm->is_slave(side)) continue;
			acUHat[side]/=(number)acVolume[side];
		}
	}
	// set grid function to filtered values
	copyAttachmentToGridFunction<dim,side_type,TGridFunction>(u,acUHat);
	// PeriodicAttachmentAccessor<side_type,Attachment<MathVector<dim> > >& aaUHat,PeriodicAttachmentAccessor<side_type,ANumber >& aaVol,
	grid.template detach_from<side_type>(aUHat);
	grid.template detach_from<side_type>(aVolume);
}

template <int dim,typename TGridFunction>
void scvFilterFVCR(SmartPtr<TGridFunction> u){
	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object elem_type;
	/// side type
	typedef typename elem_type::side side_type;
	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;
	/// side iterator
	typedef typename TGridFunction::template traits<side_type>::const_iterator SideIterator;
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;
	///	grid type
	typedef typename domain_type::grid_type grid_type;
	/// position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;

	domain_type& domain = *u->domain().get();
	grid_type& grid = *domain.grid();

	// define attachment stuff
	typedef MathVector<dim> vecDim;
	typedef Attachment<vecDim> AMathVectorDim;
	typedef PeriodicAttachmentAccessor<side_type,ANumber > aSideNumber;
	typedef PeriodicAttachmentAccessor<side_type,AMathVectorDim > aSideDimVector;

	//  filtered u attachment
	aSideDimVector acUHat;
	AMathVectorDim aUHat;
	//  volume attachment
	aSideNumber acVolume;
	ANumber aVolume;

	grid.template attach_to<side_type>(aUHat);
	grid.template attach_to<side_type>(aVolume);

	// accessors
	acUHat.access(grid,aUHat);
	acVolume.access(grid,aVolume);

	// initial values for attachments
	SetAttachmentValues(acUHat , u->template begin<side_type>(), u->template end<side_type>(), 0);
	SetAttachmentValues(acVolume , u->template begin<side_type>(), u->template end<side_type>(), 0);

	DimCRFVGeometry<dim> geo;
	std::vector<DoFIndex> multInd;
	const position_accessor_type& posAcc = domain.position_accessor();
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	// assemble fluxes
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
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
					u->dof_indices(sides[s], d, multInd);
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
						localValue[d] += vShape[j]*uValue[j][d];


				localValue *= scv.volume();
				acVolume[sides[s]]  += scv.volume();
				acUHat[sides[s]] += localValue;
			}
		}
	 }
	 PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
	// average
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		SideIterator sideIter = u->template begin<side_type>(si);
		SideIterator sideIterEnd = u->template end<side_type>(si);
		for(  ;sideIter !=sideIterEnd; sideIter++)
		{
			side_type* side = *sideIter;
			if (pbm && pbm->is_slave(side)) continue;
			acUHat[side]/=(number)acVolume[side];
		}
	}
	// set grid function to filtered values
	copyAttachmentToGridFunction<dim,side_type,TGridFunction>(u,acUHat);
	// PeriodicAttachmentAccessor<side_type,Attachment<MathVector<dim> > >& aaUHat,PeriodicAttachmentAccessor<side_type,ANumber >& aaVol,
	grid.template detach_from<side_type>(aUHat);
	grid.template detach_from<side_type>(aVolume);
}

template <int dim,typename TGridFunction>
void elementFilterFV1(SmartPtr<TGridFunction> u){
	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object elem_type;
	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;
	/// side iterator
	typedef typename TGridFunction::template traits<VertexBase>::const_iterator VertexIterator;
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;
	///	grid type
	typedef typename domain_type::grid_type grid_type;
	/// position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;

	domain_type& domain = *u->domain().get();
	grid_type& grid = *domain.grid();

	// define attachment stuff
	typedef MathVector<dim> vecDim;
	typedef Attachment<vecDim> AMathVectorDim;
	typedef PeriodicAttachmentAccessor<VertexBase,ANumber > aVertexNumber;
	typedef PeriodicAttachmentAccessor<VertexBase,AMathVectorDim > aVertexDimVector;

	//  filtered u attachment
	aVertexDimVector acUHat;
	AMathVectorDim aUHat;
	//  volume attachment
	aVertexNumber acVolume;
	ANumber aVolume;

	grid.template attach_to<VertexBase>(aUHat);
	grid.template attach_to<VertexBase>(aVolume);

	// accessors
	acUHat.access(grid,aUHat);
	acVolume.access(grid,aVolume);

	// initial values for attachments
	SetAttachmentValues(acUHat , u->template begin<VertexBase>(), u->template end<VertexBase>(), 0);
	SetAttachmentValues(acVolume , u->template begin<VertexBase>(), u->template end<VertexBase>(), 0);

	DimFV1Geometry<dim> geo;
	std::vector<DoFIndex> multInd;
	const position_accessor_type& posAcc = domain.position_accessor();
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
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
			LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::LAGRANGE, dim, 1));

			//	get Reference Mapping
			DimReferenceMapping<dim, dim>& map = ReferenceMappingProvider::get<dim, dim>(roid, coCoord);

			map.global_to_local(localBary,bary);

			//	memory for shapes
			std::vector<number> vShape;

			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

			rTrialSpace.shapes(vShape, localBary);

			size_t noc = elem->num_vertices();

			MathVector<dim> value;
			value = 0;
			number elementVolume = 0;
			for (size_t co=0;co<noc;co++){
				const typename DimFV1Geometry<dim>::SCV& scv = geo.scv(co);
				MathVector<dim> localValue;
				for (int d=0;d<dim;d++){
					u->dof_indices(elem->vertex(co), d, multInd);
					localValue[d]=DoFRef(*u,multInd[0]);
				}
				//for debug UG_LOG("localValue=" << localValue << "\n");
				//for debug UG_LOG("vShape=" << vShape[s] << "\n");
				localValue *= vShape[co];
				value += localValue;
				elementVolume += scv.volume();
			}
			//for debug UG_LOG("value=" << value << " vol=" << elementVolume << "\n");
			value *= elementVolume;
			for (size_t co=0;co<noc;co++){
				acVolume[elem->vertex(co)]+=elementVolume;
				acUHat[elem->vertex(co)] += value;
			}
		}
	}
	PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
	// average
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		VertexIterator vertexIter = u->template begin<VertexBase>(si);
		VertexIterator vertexIterEnd = u->template end<VertexBase>(si);
		for(  ;vertexIter !=vertexIterEnd; vertexIter++)
		{
			VertexBase* vert = *vertexIter;
			if (pbm && pbm->is_slave(vert)) continue;
			acUHat[vert]/=(number)acVolume[vert];
		}
	}
	// set grid function to filtered values
	copyAttachmentToGridFunction<dim,VertexBase,TGridFunction>(u,acUHat);
	grid.template detach_from<VertexBase>(aUHat);
	grid.template detach_from<VertexBase>(aVolume);
}

template <int dim,typename TGridFunction>
void scvFilterFV1(SmartPtr<TGridFunction> u){
	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object elem_type;
	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;
	/// side iterator
	typedef typename TGridFunction::template traits<VertexBase>::const_iterator VertexIterator;
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;
	///	grid type
	typedef typename domain_type::grid_type grid_type;
	/// position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;

	domain_type& domain = *u->domain().get();
	grid_type& grid = *domain.grid();

	// define attachment stuff
	typedef MathVector<dim> vecDim;
	typedef Attachment<vecDim> AMathVectorDim;
	typedef PeriodicAttachmentAccessor<VertexBase,ANumber > aVertexNumber;
	typedef PeriodicAttachmentAccessor<VertexBase,AMathVectorDim > aVertexDimVector;

	//  filtered u attachment
	aVertexDimVector acUHat;
	AMathVectorDim aUHat;
	//  volume attachment
	aVertexNumber acVolume;
	ANumber aVolume;

	grid.template attach_to<VertexBase>(aUHat);
	grid.template attach_to<VertexBase>(aVolume);

	// accessors
	acUHat.access(grid,aUHat);
	acVolume.access(grid,aVolume);

	// initial values for attachments
	SetAttachmentValues(acUHat , u->template begin<VertexBase>(), u->template end<VertexBase>(), 0);
	SetAttachmentValues(acVolume , u->template begin<VertexBase>(), u->template end<VertexBase>(), 0);

	DimFV1Geometry<dim> geo;
	std::vector<DoFIndex> multInd;
	const position_accessor_type& posAcc = domain.position_accessor();
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
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
			LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::LAGRANGE, dim, 1));

			//	get Reference Mapping
			DimReferenceMapping<dim, dim>& map = ReferenceMappingProvider::get<dim, dim>(roid, coCoord);

			map.global_to_local(localBary,bary);

			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

			size_t noc = elem->num_vertices();

			static const size_t MaxNumVerticesOfElem = 10;

			typedef MathVector<dim> MVD;
			std::vector<MVD> uValue(MaxNumVerticesOfElem);

			for (size_t co=0;co<noc;co++){
				for (int d=0;d<dim;d++){
					u->dof_indices(elem->vertex(co), d, multInd);
					uValue[co][d]=DoFRef(*u,multInd[0]);
				}
			};

			MathVector<dim> scvLocalBary;
			for (size_t co=0;co<noc;co++){
				const typename DimFV1Geometry<dim>::SCV& scv = geo.scv(co);
				scvLocalBary = 0;
				// compute barycenter of scv
				for (size_t i=0;i<scv.num_corners();i++){
					scvLocalBary += scv.local_corner(i);
				}
				scvLocalBary/=(number)(scv.num_corners());
				//	memory for shapes
				std::vector<number> vShape;
				rTrialSpace.shapes(vShape, scvLocalBary);
				MathVector<dim> localValue = 0;

				for (size_t j=0;j<noc;j++)
					for (int d=0;d<dim;d++)
						localValue[d] += vShape[j]*uValue[j][d];

				localValue *= scv.volume();
				acVolume[elem->vertex(co)]  += scv.volume();
				acUHat[elem->vertex(co)] += localValue;
			}
		}
    }
	PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
	// average
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		VertexIterator vertexIter = u->template begin<VertexBase>(si);
		VertexIterator vertexIterEnd = u->template end<VertexBase>(si);
		for(  ;vertexIter !=vertexIterEnd; vertexIter++)
		{
			VertexBase* vert = *vertexIter;
			if (pbm && pbm->is_slave(vert)) continue;
			acUHat[vert]/=(number)acVolume[vert];
		}
	}
	// set grid function to filtered values
	copyAttachmentToGridFunction<dim,VertexBase,TGridFunction>(u,acUHat);
	grid.template detach_from<VertexBase>(aUHat);
	grid.template detach_from<VertexBase>(aVolume);
}

template <typename TGridFunction>
void filter(SmartPtr<TGridFunction> u,const std::string& name){
	// world dimension
	static const int dim = TGridFunction::domain_type::dim;

	std::string n = TrimString(name);
	std::transform(n.begin(), n.end(), n.begin(), ::tolower);

	enum filtertype
	{
		ELEM_FILTER = 0,
		SCV_FILTER = 1
	};

	enum disctype
	{
		FVCR = 0,
		FV1  = 1
	};

	filtertype filter;
	disctype disc;

	if (n=="elem") filter=ELEM_FILTER;
	if (n=="scv") filter=SCV_FILTER;

	if (u->local_finite_element_id(0) == LFEID(LFEID::CROUZEIX_RAVIART, dim, 1))
		disc=FVCR; else
	if  (u->local_finite_element_id(0) == LFEID(LFEID::LAGRANGE, dim, 1))
		disc=FV1; else
	UG_THROW("Unsupported local finite element type in filter command.");

	if ((disc==FVCR)&&(filter==ELEM_FILTER)) elementFilterFVCR<dim,TGridFunction>(u);
	if ((disc==FVCR)&&(filter==SCV_FILTER)) scvFilterFVCR<dim,TGridFunction>(u);
	if ((disc==FV1)&&(filter==ELEM_FILTER)) elementFilterFV1<dim,TGridFunction>(u);
	if ((disc==FV1)&&(filter==SCV_FILTER)) scvFilterFV1<dim,TGridFunction>(u);
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__NAVIER_STOKES_FILTER__ */
