/*
 * cr_cuthill_mckee.h
 *
 *  Created on: 5.12.2012
 *      Author: Christian Wehner
 */

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__CR_CUTHILL_MCKEE__
#define __H__UG__LIB_DISC__DOF_MANAGER__CR_CUTHILL_MCKEE__

#include <vector>
#include "lib_disc/function_spaces/approximation_space.h"
#include "common/profiler/profiler.h"

namespace ug{

/**
 * This class is used to provide an ordering for indices. The ordering relation
 * is based on the connectivity-degree, i.e. on the number of connections the
 * index has. The indices with less connections are ordered first.
 */
struct CompareDeg {
///	constructor, passing field with connections for each index
	CompareDeg(const std::vector<size_t>& vInfo) : m_degree(vInfo) {}

///	comparison operator
	bool operator() (size_t i,size_t j)
	{
		return (m_degree[i] < m_degree[j]);
	}

private:
///	storage field of connections of each index
	const std::vector<size_t>& m_degree;
};
	
	
void computeDegree(std::vector<size_t>& degree,const std::vector<std::vector<size_t> >& vvConnections,size_t minpind){
	size_t n=vvConnections.size();
	degree.resize(n);
	for (size_t i=0;i<minpind;i++){
		degree[i]=vvConnections[i].size();
	}
	for (size_t i=minpind;i<n;i++){
		degree[i]=0;
		for (size_t j=0;j<vvConnections[i].size();j++){
			degree[i]+=degree[j];
		}
	}
}
	
// Cuthill-McKee order on Crouzeix-Raviart elements
void ComputeCRCuthillMcKeeOrder(std::vector<size_t>& vNewIndex,
                              std::vector<std::vector<size_t> >& vvConnections,
                              size_t minpind,
                              bool bReverse){
	size_t numInd = vvConnections.size();
	vNewIndex.resize(numInd);
	size_t numP = numInd-minpind;
	std::vector<std::vector<size_t> > pConnections(numP);
	std::vector<size_t> degree(numInd);
	std::vector<size_t>::iterator it;
	// go over all velocity nodes and fill degree vector (number of connections)
	// set up p adjacency graph, two pressure nodes are seen as associated
	// when they have a common associated velocity node
	for (size_t i=0;i<minpind;i=i+2){
		size_t localNumInd=vvConnections[i].size();
		degree[i]=localNumInd;
		degree[i+1]=localNumInd;
		vNewIndex[i]=numInd;
		vNewIndex[i+1]=numInd;
		size_t localNumP=0;
		size_t associatedP[2];
		for (size_t j=0;j<localNumInd;j++){
			if (vvConnections[i][j]>=minpind){
				associatedP[localNumP]=vvConnections[i][j]-minpind;
				localNumP++;
			}
		}
		if (localNumP>1){
			size_t pi = associatedP[0];
			size_t pj = associatedP[1];
			it = find(pConnections[pi].begin(), pConnections[pi].end(), pj);
			if(it == pConnections[pi].end()) pConnections[pi].push_back(pj+minpind);
			//	add the opposite direction (adjInd -> index)
			it = find(pConnections[pj].begin(), pConnections[pj].end(), pi);
			if(it == pConnections[pj].end()) pConnections[pj].push_back(pi+minpind);
		}
	}
	// go over all pressure nodes, fill degree vector (sum of number of connections of all associated velocity nodes)
	// compute minimum degree pressure node
	size_t minpdeg = numInd;
	size_t minpdegind = minpind;
	for (size_t i=minpind;i<numInd;i++){
		degree[i]=0;
		vNewIndex[i]=numInd;
		size_t localNumInd=vvConnections[i].size();
		for (size_t j=0;j<localNumInd;j++){
			if (vvConnections[i][j]!=i) degree[i]+=degree[vvConnections[i][j]];
		}
		if (degree[i]<minpdeg){
			minpdeg=degree[i];
			minpdegind=i;
		}
	}
	// set up list of pressure nodes, first entry is minimum degree index
	std::vector<size_t> plist(1);
	plist[0]=minpdegind;
	size_t count=0;
	//	Sort neighbours by degreee
	CompareDeg myCompDegree(degree);
	// go over pressure nodes and sort
	for (size_t j=0;j<numP;j++){
		size_t i=plist[j];
		// sort adjacent indices by degree (remark: by definition the degree of velocity indices is always lower than
		// the associated pressure index degree)
		std::sort(	vvConnections[i].begin(),
				          	vvConnections[i].end(), myCompDegree);
		// assign new index
		for (size_t k=0;k<vvConnections[i].size();k++){
			if (vNewIndex[vvConnections[i][k]]<numInd) continue;
			vNewIndex[vvConnections[i][k]]=count;
			//debug UG_LOG("index(" << vvConnections[i][k] << ")=" << count << "\n");
			count++;
		}
		size_t pind=i-minpind;
		// sort associated p nodes (using p adjacency list)
		std::sort(	pConnections[pind].begin(),
	          	pConnections[pind].end(), myCompDegree);
		// check if already in p list, if not insert
		for (size_t k=0;k<pConnections[pind].size();k++){
			if (vNewIndex[pConnections[pind][k]]<numInd) continue;
			bool found=false;
			for (size_t iiiiii=j;iiiiii<plist.size();iiiiii++){
				if (plist[iiiiii]==pConnections[pind][k]){
					found=true;
					break;
				}
			}
			if (found==false){
				plist.push_back(pConnections[pind][k]);
			}
		}
	}
}

// compute adjacency graph on Crouzeix-Raviart-elements
template <typename TElem, typename TDD,typename TGridFunction>
void cr_get_connections(std::vector<std::vector<size_t> >& vvConnection,size_t& minpind,TDD& dd,TGridFunction& u){
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;

	///	world dimension
	static const int dim = domain_type::dim;

	typename TDD::template traits<TElem>::const_iterator iter, iterBegin, iterEnd;
	LocalIndices ind;

	size_t num_dd_ind = dd.num_indices();

	//	clear neighbors and resize
	vvConnection.clear();
	vvConnection.resize(num_dd_ind);

	minpind=num_dd_ind;

	for (int si=0;si<u.num_subsets();si++){
		iterBegin = dd.template begin<TElem>(si);
		iterEnd = dd.template end<TElem>(si);
		// 	check if at least one element exist
		if(iterBegin == iterEnd) continue;
		// 	Loop over all elements
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		// 	get Element
			TElem* elem = *iter;
		// 	get global indices
			dd.indices(elem, ind, false);
			size_t num_sides = elem->num_sides();
			size_t num_ind = num_sides * dim + 1;
			std::vector<size_t> globind(num_ind);
		 // velocity component indices
			for (size_t side=0;side<num_sides;side++)
				for (int d=0;d<dim;d++)
					globind[dim*side+d]=ind.index(d,side);
		 // pressure component index
			size_t pind = ind.index(dim,0);
			globind[num_ind-1]=pind;
			if (pind<minpind) minpind=pind;
			//	add couplings to adjacency graph
			std::vector<size_t>::iterator it;
			for (size_t i=0;i<num_ind;i++){
				size_t iIndex = globind[i];
				for (size_t j=0;j<num_ind;j++){
					size_t jIndex = globind[j];
					//	add connection (iIndex->jIndex)
					it = find(vvConnection[iIndex].begin(), vvConnection[iIndex].end(), jIndex);
					if(it == vvConnection[iIndex].end()) vvConnection[iIndex].push_back(jIndex);
					//	add the opposite direction (adjInd -> index)
					it = find(vvConnection[jIndex].begin(), vvConnection[jIndex].end(), iIndex);
					if(it == vvConnection[jIndex].end()) vvConnection[jIndex].push_back(iIndex);
				}
			}
		}
		std::vector<size_t> degree;
		computeDegree(degree,vvConnection,minpind);
		//	Sort neighbours by degreee
		CompareDeg myCompDegree(degree);	
		for (size_t i=0;i<vvConnection.size();i++){
			std::sort(vvConnection[i].begin(),vvConnection[i].end(),myCompDegree);
		};
	}
/*debug	
	for (size_t i=0;i<vvConnection.size();i++){
		UG_LOG(i << " : ");
		for (size_t j=0;j<vvConnection[i].size();j++){
			UG_LOG(vvConnection[i][j] << ",");
		}
		UG_LOG("\n");
	}*/
}

/// orders the all DofDistributions of the ApproximationSpace using Cuthill-McKee
template <typename TDomain,typename TGridFunction>
void OrderCRCuthillMcKee(ApproximationSpace<TDomain>& approxSpace,TGridFunction& u,
                       bool bReverse)
{
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;

	///	world dimension
	static const int dim = domain_type::dim;

	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object elem_type;

	for (int si=0;si<u.num_subsets();++si){
		if (u.num_fct(si)<2){
			if (u.num_fct(si)!=dim+1){
				UG_THROW("Wrong number of components in approximation space. Needed " << dim+1 << ", given " << u.num_fct(si) << ".\n");
			}
		}
		for (int d=0;d<dim;d++){
			if (u.local_finite_element_id(d) != LFEID(LFEID::CROUZEIX_RAVIART, 1)){
				UG_THROW("Component " << d << " in approximation space must be of Crouzeix-Raviart type.");
			}
		}
		if (u.local_finite_element_id(dim) != LFEID(LFEID::PIECEWISE_CONSTANT,0)){
			UG_THROW("Component dim in approximation space must be of piecewise constant type.");
		}
	}

	if(approxSpace.levels_enabled())
		for(size_t lev = 0; lev < approxSpace.num_levels(); ++lev){
			//	get adjacency graph
			std::vector<std::vector<size_t> > vvConnection;
			size_t minpind;
			cr_get_connections<elem_type,LevelDoFDistribution,TGridFunction>(vvConnection,minpind,*approxSpace.level_dof_distribution(lev),u);
			//	get mapping for cuthill-mckee order
			std::vector<size_t> vNewIndex;
			ComputeCRCuthillMcKeeOrder(vNewIndex,vvConnection,minpind,bReverse);
			//	reorder indices
			approxSpace.level_dof_distribution(lev)->permute_indices(vNewIndex);
		}

//	order surface
	if(approxSpace.top_surface_enabled()){
		//	get adjacency graph
		std::vector<std::vector<size_t> > vvConnection;
		size_t minpind;
		cr_get_connections<elem_type,SurfaceDoFDistribution,TGridFunction>(vvConnection,minpind,*approxSpace.surface_dof_distribution(),u);
		//	get mapping for cuthill-mckee order
		std::vector<size_t> vNewIndex;
		ComputeCRCuthillMcKeeOrder(vNewIndex,vvConnection,minpind,bReverse);
		//	reorder indices
		approxSpace.surface_dof_distribution()->permute_indices(vNewIndex);
	}
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__CR_CUTHILL_MCKEE__ */
