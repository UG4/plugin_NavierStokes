/*
 * filter_impl.h
 *
 *  Created on: 8.04.2013
 *      Author: Christian Wehner
 */

#ifndef __H__UG__NAVIER_STOKES__INCOMPRESSIBLE__FILTER__IMPL__
#define __H__UG__NAVIER_STOKES__INCOMPRESSIBLE__FILTER__IMPL__

#include "common/util/provider.h"

namespace ug{
namespace NavierStokes{

// compute average element size of grid
template<typename TGridFunction>
number ConstantBoxFilter<TGridFunction>::compute_average_element_size(SmartPtr<TGridFunction> u){
	typedef typename TGridFunction::const_element_iterator const_iterator;
	typedef typename TGridFunction::element_type element_type;
	typename TGridFunction::domain_type::position_accessor_type& aaPos
				= u->domain()->position_accessor();
	number averageElemSize = 0;
	std::vector<MathVector<dim> > vCorner;
	size_t nrOfElem = 0;
	const_iterator iter = u->template begin<elem_type>();
	const_iterator iterEnd = u->template end<elem_type>();
	for(; iter != iterEnd; ++iter)
	{
		elem_type* elem = *iter;
		ReferenceObjectID roid = elem->reference_object_id();
		CollectCornerCoordinates(vCorner, *elem, aaPos);
		averageElemSize += ElementSize<dim>(roid, &vCorner[0]);
		nrOfElem++;
	}
	return averageElemSize/(number)nrOfElem;
}

template<int dim,typename TPosAcc>
void getElemCoord(MathVector<dim>& co,Vertex* vrt,TPosAcc posAcc){
	co = posAcc[vrt];
}

template<int dim,typename TPosAcc>
void getElemCoord(MathVector<dim>& co,EdgeBase* edge,TPosAcc posAcc){
	co = posAcc[edge->vertex(0)];
	co+=posAcc[edge->vertex(1)];
	co*=0.5;
}

template<int dim,typename TPosAcc>
void getElemCoord(MathVector<dim>& co,Face* face,TPosAcc posAcc){
	co = posAcc[face->vertex(0)];
	for (size_t i=1;i<face->num_vertices();i++){
		co+=posAcc[face->vertex(i)];
	}
	co /= (number) face->num_vertices();
}

// search neighbor elements in filter area by breadth first search
template <typename TGridFunction>
template <typename TAElem,typename TDofElem>
void ConstantBoxFilter<TGridFunction>::handleNeighbors(std::vector<TAElem> nb,number filterWidth,std::vector<MathVector<dim> > ipCo,
		std::vector<MathVector<dim> > ipVol,std::vector<MathVector<dim> > ipVolVal,
		PeriodicAttachmentAccessor<TDofElem,Attachment<MathVector<dim> > >& aaU,
		PeriodicAttachmentAccessor<TDofElem,Attachment<number > >& aaVol){
	std::vector<TAElem*> nbrCandidates;
	number tol=0.5*filterWidth;
	for (size_t nbIndex=0;nbIndex<nb.size();nbIndex++){
		// check if node is inside
		MathVector<dim> co;
		getElemCoord(co,nb[nbIndex],m_posAcc);
		bool outsideForAll=true;
		for (size_t ip=0;ip<ipCo.size();ip++){
			bool inside = true;
			for (size_t j=0;j<dim;j++){
				if (std::abs(ipCo[j]-co[j])>tol){
					inside = false;
					break;
				}
			}
			if (inside==true){
				outsideForAll=false;
				aaU[ nb[nbIndex] ]+=ipVolVal[ip];
				aaVol[ nb[nbIndex] ]+= ipVol[ip];
			}
		}
		// if element is inside search neighbors and add them
		if (outsideForAll==false){
			CollectNeighbors(nbrCandidates, *m_grid, nb[nbIndex]);
			size_t nbSize = nb.size();
			for (size_t j=0;j<nbrCandidates.size();j++){
				bool newNeighbor=true;
				for (size_t k=0;k<nb.size();k++){
					if (nbrCandidates[j]==nb[k]){
						newNeighbor=false;
						break;
					}
				}
				if (newNeighbor==true){
					nb.push_back(nbrCandidates[j]);
				}
			}
		}
	}
}

template <typename TGridFunction>
template <typename VType>
void ConstantBoxFilter<TGridFunction>::apply_filter(PeriodicAttachmentAccessor<vertex_type,Attachment<VType> >& acUHat,
			   SmartPtr<TGridFunction> u,PeriodicAttachmentAccessor<vertex_type,Attachment<VType> >& aaU){
//	bool useGridFunction = true;
//	if (u==NULL) useGridFunction = false;

	std::vector<DoFIndex> multInd;

	//	get domain of grid function
	domain_type& domain = *m_uInfo->domain().get();
	DimFV1Geometry<dim> geo;

	// volume attachment
	typedef PeriodicAttachmentAccessor<vertex_type,ANumber > aNumber;
	//  volume attachment
	aNumber acVolume;
	ANumber aVolume;
	m_grid->template attach_to<Vertex>(aVolume);
	acVolume.access(*m_grid,aVolume);

	SetAttachmentValues(acVolume , u->template begin<Vertex>(), u->template end<Vertex>(), 0);

	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

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
				coCoord[i] = m_posAcc[vVrt[i]];
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

	/*			for (size_t co=0;co<noc;co++){
					for (int d=0;d<dim;d++){
						u->dof_indices(elem->vertex(co), d, multInd);
						uValue[co][d]=DoFRef(*u,multInd[0]);
					}
				};*/

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
					for (size_t d=0;d<dim;d++)
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
		VertexIterator vertexIter = u->template begin<Vertex>(si);
		VertexIterator vertexIterEnd = u->template end<Vertex>(si);
		for(  ;vertexIter !=vertexIterEnd; vertexIter++)
		{
			Vertex* vert = *vertexIter;
			if (pbm && pbm->is_slave(vert)) continue;
			acUHat[vert]/=(number)acVolume[vert];
		}
	}
};


} // namespace NavierStokes
} // end namespace ug

#endif /* __H__UG__NAVIER_STOKES__INCOMPRESSIBLE__FILTER__IMPL__ */
