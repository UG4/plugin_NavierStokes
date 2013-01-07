/*
 * cr_reorder.h
 *
 *  Created on: 7.1.2013
 *      Author: Christian Wehner
 */

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__CR_REORDER__
#define __H__UG__LIB_DISC__DOF_MANAGER__CR_REORDER__

#include <vector>
#include "lib_disc/function_spaces/approximation_space.h"
#include "common/profiler/profiler.h"
#include "cr_order_cuthill_mckee.h"
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/sloan_ordering.hpp>
#include <boost/graph/king_ordering.hpp>
#include <boost/graph/minimum_degree_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>
#include <boost/graph/graphviz.hpp>
#include "lib_disc/dof_manager/surface_dof_distribution.h"
#include "lib_disc/dof_manager/level_dof_distribution.h"

namespace ug{
	
// extract p graph from connections, two p nodes are defined adjacent if they share a common velocity neighbour
template <class TGraph>
void extractPGraph(TGraph& G,const std::vector<std::vector<size_t> >& vvConnection,size_t pmin,size_t dim,bool directed){
	for (size_t i=0;i<pmin;i=i+dim){
		size_t s=vvConnection[i].size();
		int ind1 = vvConnection[i][s-2]-pmin;
		if (ind1<0) continue;
		int ind2 = vvConnection[i][s-1]-pmin;
		boost::add_edge(ind1,ind2,G);
		if (directed==true) boost::add_edge(ind2,ind1,G);
	}	
}
	
// extract v graph from connections
template <class TGraph>
void extractVGraph(TGraph& G,const std::vector<std::vector<size_t> >& vvConnection,size_t pmin,size_t dim,bool directed){
	for (size_t i=0;i<pmin;i=i+dim){
		size_t ind1 = round((number)i/dim);
		size_t s = vvConnection[i].size();
		for (size_t j=0;j<s-1;j=j+2){
			size_t ind2 = round((number)vvConnection[i][j]/dim);
			if (dim*ind2>=pmin) break;
			if (ind1==ind2) continue;
			boost::add_edge(ind1,ind2,G);
			if (directed==true) boost::add_edge(ind2,ind1,G);
		}
	}
}

	
void computeNewIndicesFromPIndexset(const std::vector<std::vector<size_t> >& vvConnection,std::vector<size_t>& newIndex,size_t pmin,const std::vector<int>& inv_perm,bool bseparate){
	size_t n = vvConnection.size();
	newIndex.resize(n);
	size_t pn = n - pmin;
	size_t undefined = n+1;
	for (size_t i=0;i<n;i++){
		newIndex[i] = undefined;
	}
	// group p indices behind v indices
	if (bseparate==true){
		size_t vindex=0;
		for (size_t i=0;i<pn;i++){
			size_t pind = inv_perm[i]+pmin;
			newIndex[pind]=pmin+i;
			size_t isize = vvConnection[pind].size()-1;
			for (size_t jj=0;jj<isize;jj++){
				size_t j = vvConnection[pind][jj];
				if (newIndex[j]==undefined){
					newIndex[j]=vindex;
					vindex++;
				}
			}
		}
	// mixed ordering, p indices inserted directly behind associated v nodes
	} else {
		size_t index=0;
		for (size_t i=0;i<pn;i++){
			size_t pind = inv_perm[i]+pmin;
			size_t isize = vvConnection[pind].size()-1;
			for (size_t jj=0;jj<isize;jj++){
				size_t j = vvConnection[pind][jj];
				if (newIndex[j]==undefined){
					newIndex[j]=index;
					index++;
				}
			}
			newIndex[pind]=index;
			index++;
		}
	}
}

	
void computeNewIndicesFromVIndexset(const std::vector<std::vector<size_t> >& vvConnection,std::vector<size_t>& newIndex,size_t pmin,size_t dim,const std::vector<int>& inv_perm,bool bseparate,size_t strategy){
	size_t n = vvConnection.size();
	newIndex.resize(n);
	size_t vn = round((number)1.0/dim*pmin);
	size_t undefined = n+1;
	for (size_t i=0;i<n;i++){
		newIndex[i] = undefined;
	}
	if (bseparate==true){
		size_t pindex=0;
		// group p behind using first given v node
		if (strategy==1){
			for (size_t i=0;i<vn;i++){
				size_t vind = dim*inv_perm[i];
				for (size_t d=0;d<dim;d++) newIndex[vind+d]=dim*i+d;
				size_t isize = vvConnection[vind].size();
				for (size_t jj=isize-2;jj<isize;jj++){
					size_t j = vvConnection[vind][jj];
					if (j<pmin) continue;
					if (newIndex[j]==undefined){
						newIndex[j]=pindex+pmin;
						pindex++;
					}
				}
			}
		}
		// group p behind when all associated v nodes are handled
		if (strategy==2){
			for (size_t i=0;i<vn;i++){
				size_t vind = dim*inv_perm[i];
				for (size_t d=0;d<dim;d++) newIndex[vind+d]=dim*i+d;
				size_t isize = vvConnection[vind].size();
				for (size_t jj=isize-2;jj<isize;jj++){
					size_t j = vvConnection[vind][jj];
					if (j<pmin) continue;
					if (newIndex[j]==undefined){
						bool giveindex=true;
						size_t jsize = vvConnection[j].size()-1;
						for (size_t kk=0;kk<jsize;kk=kk+dim){
							size_t k = vvConnection[j][kk];
							if (newIndex[k]==n+1){
								giveindex=false;
								break;
							}
						}
						if (giveindex==true){
							newIndex[j]=pindex+pmin;
							pindex++;
						}
					}
				}
			}
		}
	// separate==false mix v and p indices, p is inserted when all associated v nodes are handled
	} else {
		size_t index=0;
		for (size_t i=0;i<vn;i++){
			size_t vind = dim*inv_perm[i];
			for (size_t d=0;d<dim;d++) newIndex[vind+d]=index+d;
			index=index+dim;
			size_t isize = vvConnection[vind].size();
			for (size_t jj=isize-2;jj<isize;jj++){
				size_t j = vvConnection[vind][jj];
				if (j<pmin) continue;
				if (newIndex[j]==undefined){
					bool giveindex=true;
					size_t jsize = vvConnection[j].size()-1;
					for (size_t kk=0;kk<jsize;kk=kk+dim){
						size_t k = vvConnection[j][kk];
						if (newIndex[k]==undefined){
							giveindex=false;
							break;
						}
					}
					if (giveindex==true){
						newIndex[j]=index;
						index++;
					}
				}
			}
		}
	}
}

/// orders the dof distribution using Cuthill-McKee
template <typename TDD>
void CRCuthillMcKee(TDD& dofDistr,std::vector<std::vector<size_t> > vvConnection,
						 size_t minpind,bool bReverse,bool bseparate,bool orderp,bool orderv)
{
	PROFILE_FUNC();
		
	size_t n = vvConnection.size();
	size_t pnr = n-minpind;
	size_t vnr = round((number)minpind/2.0);
	
	const size_t dim = dofDistr.dim(0);
		
	typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS,
	boost::property<boost::vertex_color_t, boost::default_color_type,
	boost::property<boost::vertex_degree_t,int> > >
	Graph;
	std::vector<int> inv_permP(0);
	std::vector<int> inv_permV(0);
	std::vector<size_t> newIndex(n);

	if (orderp==true){
		Graph pGraph(pnr);
		inv_permP.resize(pnr);
		extractPGraph(pGraph,vvConnection,minpind,dim,false);
		if (bReverse==true)
			cuthill_mckee_ordering(pGraph, inv_permP.rbegin(), boost::get(boost::vertex_color, pGraph),
							   boost::make_degree_map(pGraph));
		else
			cuthill_mckee_ordering(pGraph, inv_permP.begin(), boost::get(boost::vertex_color, pGraph),
								   boost::make_degree_map(pGraph));
		if (orderv==false)
			computeNewIndicesFromPIndexset(vvConnection,newIndex,minpind,inv_permP,bseparate);
	}
	
	if (orderv==true){
		Graph vGraph(vnr);
		inv_permV.resize(vnr);
		extractVGraph(vGraph,vvConnection,minpind,dim,false);
		if (bReverse==true)
			cuthill_mckee_ordering(vGraph, inv_permV.rbegin(), boost::get(boost::vertex_color, vGraph),
							   boost::make_degree_map(vGraph));
		else 
			cuthill_mckee_ordering(vGraph, inv_permV.begin(), boost::get(boost::vertex_color, vGraph),
								   boost::make_degree_map(vGraph));
		size_t strategy = 1;
		if (orderp==false)
			computeNewIndicesFromVIndexset(vvConnection,newIndex,minpind,dim,inv_permV,bseparate,strategy);
	}

	if ((orderp==true)&&(orderv==true)){
		// order from p and v index set
		for (size_t i=0;i<(number)1.0/dim*minpind;i++){
			for (size_t d=0;d<dim;d++) newIndex[dim*inv_permV[i]+d]=dim*i+d;
		}
		for (size_t i=0;i<pnr;i++){
			newIndex[inv_permP[i]+minpind]=i+minpind;
		}
	}
			
	dofDistr.permute_indices(newIndex);
}
	
/// orders the all DofDistributions of the ApproximationSpace using Cuthill-McKee
template <typename TDomain,typename TGridFunction>
void CROrderCuthillMcKee(ApproximationSpace<TDomain>& approxSpace,TGridFunction& u,
                       bool bReverse,bool bseparate,bool orderp,bool orderv)
{
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;
	
	///	world dimension
	static const int dim = domain_type::dim;
	
	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object elem_type;
	
	std::vector<std::vector<size_t> > vvConnection;
	size_t minpind;
	
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
	
	//	order levels
	if(approxSpace.levels_enabled())
		for(size_t lev = 0; lev < approxSpace.num_levels(); ++lev){
			cr_get_connections<elem_type,LevelDoFDistribution,TGridFunction>(vvConnection,minpind,*approxSpace.level_dof_distribution(lev),u);
			CRCuthillMcKee(*approxSpace.level_dof_distribution(lev),vvConnection,minpind,bReverse,bseparate,orderp,orderv);
		}
	
	//	order surface
	if(approxSpace.top_surface_enabled()){
		cr_get_connections<elem_type,SurfaceDoFDistribution,TGridFunction>(vvConnection,minpind,*approxSpace.surface_dof_distribution(),u);	
		CRCuthillMcKee(*approxSpace.surface_dof_distribution(),vvConnection,minpind,bReverse,bseparate,orderp,orderv);
	}
}

/// orders the dof distribution using Sloan	
template <typename TDD>
void CRSloan(TDD& dofDistr,std::vector<std::vector<size_t> > vvConnection,
						 size_t minpind,bool bseparate,bool orderp,bool orderv)
{
	PROFILE_FUNC();
		
	size_t n = vvConnection.size();
	size_t pnr = n-minpind;
	size_t vnr = round((number)minpind/2.0);
	
	const size_t dim = dofDistr.dim(0);
		
	typedef boost::adjacency_list<
		boost::setS, 
		boost::vecS, 
		boost::undirectedS, 
		boost::property<
		boost::vertex_color_t, 
		boost::default_color_type,
		boost::property<
		boost::vertex_degree_t,
		int,
		boost::property<
		boost::vertex_priority_t,
		double > > > > Graph;
		
	typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
	typedef typename boost::property_map<Graph,boost::vertex_degree_t>::type deg_type;
	typedef typename boost::property_map<Graph, boost::vertex_index_t>::type map_type;
		
	std::vector<int> inv_permP(0);
	std::vector<int> inv_permV(0);
	std::vector<size_t> newIndex(n);

	if (orderp==true){
		Graph G(pnr);
		inv_permP.resize(pnr);
		extractPGraph(G,vvConnection,minpind,dim,false);

		//Creating two iterators over the vertices
		boost::graph_traits<Graph>::vertex_iterator ui, ui_end;
	
		//Creating a property_map with the degrees of the degrees of each vertex
		deg_type deg = get(boost::vertex_degree, G);
		for (boost::tie(ui, ui_end) = vertices(G); ui != ui_end; ++ui)
			deg[*ui] = degree(*ui, G);
	
		//Creating a property_map for the indices of a vertex
		map_type index_map = get(boost::vertex_index, G);
	
		//Setting the start node
		Vertex s = vertex(0, G);
		int ecc;   //defining a variable for the pseudoperipheral radius
    
		//Calculating the pseudoeperipheral node and radius
		Vertex e = pseudo_peripheral_pair(G, s, ecc, boost::get(boost::vertex_color, G), boost::get(boost::vertex_degree, G) );

		//Sloan ordering
		sloan_ordering(G, s, e, inv_permP.begin(), get(boost::vertex_color, G), 
				   get(boost::vertex_degree, G), get(boost::vertex_priority, G));
				   
		if (orderv==false)
			computeNewIndicesFromPIndexset(vvConnection,newIndex,minpind,inv_permP,bseparate);
	}
	
	if (orderv==true){
		Graph G(vnr);
		inv_permV.resize(vnr);
		extractVGraph(G,vvConnection,minpind,dim,false);
		
		//Creating two iterators over the vertices
		boost::graph_traits<Graph>::vertex_iterator ui, ui_end;
	
		//Creating a property_map with the degrees of the degrees of each vertex
		deg_type deg = get(boost::vertex_degree, G);
		for (boost::tie(ui, ui_end) = vertices(G); ui != ui_end; ++ui)
			deg[*ui] = degree(*ui, G);
	
		//Creating a property_map for the indices of a vertex
		map_type index_map = get(boost::vertex_index, G);
	
		//Setting the start node
		Vertex s = vertex(0, G);
		int ecc;   //defining a variable for the pseudoperipheral radius
    
		//Calculating the pseudoeperipheral node and radius
		Vertex e = pseudo_peripheral_pair(G, s, ecc, boost::get(boost::vertex_color, G), boost::get(boost::vertex_degree, G) );

		//Sloan ordering
		sloan_ordering(G, s, e, inv_permV.begin(), get(boost::vertex_color, G), 
				   get(boost::vertex_degree, G), get(boost::vertex_priority, G));

		size_t strategy = 1;
		if (orderp==false)
			computeNewIndicesFromVIndexset(vvConnection,newIndex,minpind,dim,inv_permV,bseparate,strategy);
	}

	if ((orderp==true)&&(orderv==true)){
		// order from p and v index set
		for (size_t i=0;i<(number)1.0/dim*minpind;i++){
			for (size_t d=0;d<dim;d++) newIndex[dim*inv_permV[i]+d]=dim*i+d;
		}
		for (size_t i=0;i<pnr;i++){
			newIndex[inv_permP[i]+minpind]=i+minpind;
		}
	}
			
	dofDistr.permute_indices(newIndex);
}
	
/// orders the all DofDistributions of the ApproximationSpace using Sloan	
template <typename TDomain,typename TGridFunction>
void CROrderSloan(ApproximationSpace<TDomain>& approxSpace,TGridFunction& u,
                       bool bseparate,bool orderp,bool orderv)
{
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;
	
	///	world dimension
	static const int dim = domain_type::dim;
	
	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object elem_type;
	
	std::vector<std::vector<size_t> > vvConnection;
	size_t minpind;
	
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
	
	//	order levels
	if(approxSpace.levels_enabled())
		for(size_t lev = 0; lev < approxSpace.num_levels(); ++lev){
			cr_get_connections<elem_type,LevelDoFDistribution,TGridFunction>(vvConnection,minpind,*approxSpace.level_dof_distribution(lev),u);
			CRSloan(*approxSpace.level_dof_distribution(lev),vvConnection,minpind,bseparate,orderp,orderv);
		}
	
	//	order surface
	if(approxSpace.top_surface_enabled()){
		cr_get_connections<elem_type,SurfaceDoFDistribution,TGridFunction>(vvConnection,minpind,*approxSpace.surface_dof_distribution(),u);	
		CRSloan(*approxSpace.surface_dof_distribution(),vvConnection,minpind,bseparate,orderp,orderv);
	}
}

/// orders the dof distribution using King
template <typename TDD>
void CRKing(TDD& dofDistr,std::vector<std::vector<size_t> > vvConnection,
						 size_t minpind,bool bReverse,bool bseparate,bool orderp,bool orderv)
{
	PROFILE_FUNC();
		
	size_t n = vvConnection.size();
	size_t pnr = n-minpind;
	size_t vnr = round((number)minpind/2.0);
	
	const size_t dim = dofDistr.dim(0);
		
	typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS,
		boost::property<boost::vertex_color_t, boost::default_color_type,
		boost::property<boost::vertex_degree_t,int> > >
	Graph;
	
	std::vector<int> inv_permP(0);
	std::vector<int> inv_permV(0);
	std::vector<size_t> newIndex(n);

	if (orderp==true){
		Graph pGraph(pnr);
		inv_permP.resize(pnr);
		extractPGraph(pGraph,vvConnection,minpind,dim,false);
		if (bReverse==true)
			king_ordering(pGraph, inv_permP.rbegin());
		else
			king_ordering(pGraph, inv_permP.begin());
			
		if (orderv==false)
			computeNewIndicesFromPIndexset(vvConnection,newIndex,minpind,inv_permP,bseparate);
	}
	
	if (orderv==true){
		Graph vGraph(vnr);
		inv_permV.resize(vnr);
		extractVGraph(vGraph,vvConnection,minpind,dim,false);
		if (bReverse==true)
			king_ordering(vGraph, inv_permV.rbegin());
		else
			king_ordering(vGraph, inv_permV.begin());
			
		size_t strategy = 1;
		if (orderp==false)
			computeNewIndicesFromVIndexset(vvConnection,newIndex,minpind,dim,inv_permV,bseparate,strategy);
	}

	if ((orderp==true)&&(orderv==true)){
		// order from p and v index set
		for (size_t i=0;i<(number)1.0/dim*minpind;i++){
			for (size_t d=0;d<dim;d++) newIndex[dim*inv_permV[i]+d]=dim*i+d;
		}
		for (size_t i=0;i<pnr;i++){
			newIndex[inv_permP[i]+minpind]=i+minpind;
		}
	}
			
	dofDistr.permute_indices(newIndex);
}
	
/// orders the all DofDistributions of the ApproximationSpace using Cuthill-McKee
template <typename TDomain,typename TGridFunction>
void CROrderKing(ApproximationSpace<TDomain>& approxSpace,TGridFunction& u,
                       bool bReverse,bool bseparate,bool orderp,bool orderv)
{
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;
	
	///	world dimension
	static const int dim = domain_type::dim;
	
	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object elem_type;
	
	std::vector<std::vector<size_t> > vvConnection;
	size_t minpind;
	
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
	
	//	order levels
	if(approxSpace.levels_enabled())
		for(size_t lev = 0; lev < approxSpace.num_levels(); ++lev){
			cr_get_connections<elem_type,LevelDoFDistribution,TGridFunction>(vvConnection,minpind,*approxSpace.level_dof_distribution(lev),u);
			CRKing(*approxSpace.level_dof_distribution(lev),vvConnection,minpind,bReverse,bseparate,orderp,orderv);
		}
	
	//	order surface
	if(approxSpace.top_surface_enabled()){
		cr_get_connections<elem_type,SurfaceDoFDistribution,TGridFunction>(vvConnection,minpind,*approxSpace.surface_dof_distribution(),u);	
		CRKing(*approxSpace.surface_dof_distribution(),vvConnection,minpind,bReverse,bseparate,orderp,orderv);
	}
}

/// orders the dof distribution using minimum degree algorithm	
template <typename TDD>
void CRMinimumDegree(TDD& dofDistr,std::vector<std::vector<size_t> > vvConnection,
						 size_t minpind,bool bseparate,bool orderp,bool orderv)
{
	PROFILE_FUNC();
		
	size_t n = vvConnection.size();
	size_t pnr = n-minpind;
	size_t vnr = round((number)minpind/2.0);
	
	const size_t dim = dofDistr.dim(0);
		
	typedef boost::adjacency_list<boost::setS, boost::vecS, boost::directedS>  Graph;
		
	typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
	typedef typename boost::property_map<Graph,boost::vertex_degree_t>::type deg_type;
	typedef typename boost::property_map<Graph, boost::vertex_index_t>::type map_type;
		
	std::vector<int> inv_permP(0);
	std::vector<int> inv_permV(0);
	std::vector<size_t> newIndex(n);
	
	int delta=0;

	if (orderp==true){
		Graph G(pnr);
		
		inv_permP.resize(pnr,0);
		std::vector<int> perm(pnr,0);
		std::vector<int> degree(pnr,0);
		std::vector<int> supernode_sizes(pnr,1);
		
		extractPGraph(G,vvConnection,minpind,dim,true);
		
		// Creating a property_map for the indices of a vertex
		map_type index_map = get(boost::vertex_index, G);
		
		// boost minimum degree algorithm does not work for 1 edge graph use king instead
		if (pnr<=2)
			king_ordering(G, inv_permP.rbegin());
		else
			minimum_degree_ordering
			(G,
			boost::make_iterator_property_map(&degree[0],index_map, degree[0]),
			perm.begin(),
			inv_permP.begin(),
			boost::make_iterator_property_map(&supernode_sizes[0],index_map, supernode_sizes[0]), 
			delta, index_map);
				   
		if (orderv==false)
			computeNewIndicesFromPIndexset(vvConnection,newIndex,minpind,inv_permP,bseparate);
	}
	
	if (orderv==true){
		Graph G(vnr);
		inv_permV.resize(vnr,0);
		std::vector<int> perm(vnr,0);
		std::vector<int> degree(vnr,0);
		std::vector<int> supernode_sizes(vnr,1);
		
		extractVGraph(G,vvConnection,minpind,dim,true);
		
		// Creating a property_map for the indices of a vertex
		map_type index_map = get(boost::vertex_index, G);
		
		minimum_degree_ordering
		(G,
		boost::make_iterator_property_map(&degree[0],index_map, degree[0]),
		perm.begin(),
		inv_permV.begin(),
		boost::make_iterator_property_map(&supernode_sizes[0],index_map, supernode_sizes[0]), 
		delta, index_map);

		size_t strategy = 1;
		if (orderp==false)
			computeNewIndicesFromVIndexset(vvConnection,newIndex,minpind,dim,inv_permV,bseparate,strategy);
	}

	if ((orderp==true)&&(orderv==true)){
		// order from p and v index set
		for (size_t i=0;i<(number)1.0/dim*minpind;i++){
			for (size_t d=0;d<dim;d++) newIndex[dim*inv_permV[i]+d]=dim*i+d;
		}
		for (size_t i=0;i<pnr;i++){
			newIndex[inv_permP[i]+minpind]=i+minpind;
		}
	}
			
	dofDistr.permute_indices(newIndex);
}
	
/// orders the all DofDistributions of the ApproximationSpace using minimum degree algorithm	
template <typename TDomain,typename TGridFunction>
void CROrderMinimumDegree(ApproximationSpace<TDomain>& approxSpace,TGridFunction& u,
                       bool bseparate,bool orderp,bool orderv)
{
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;
	
	///	world dimension
	static const int dim = domain_type::dim;
	
	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object elem_type;
	
	std::vector<std::vector<size_t> > vvConnection;
	size_t minpind;
	
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
	
	//	order levels
	if(approxSpace.levels_enabled())
		for(size_t lev = 0; lev < approxSpace.num_levels(); ++lev){
			cr_get_connections<elem_type,LevelDoFDistribution,TGridFunction>(vvConnection,minpind,*approxSpace.level_dof_distribution(lev),u);
			CRMinimumDegree(*approxSpace.level_dof_distribution(lev),vvConnection,minpind,bseparate,orderp,orderv);
		}
	
	//	order surface
	if(approxSpace.top_surface_enabled()){
		cr_get_connections<elem_type,SurfaceDoFDistribution,TGridFunction>(vvConnection,minpind,*approxSpace.surface_dof_distribution(),u);	
		CRMinimumDegree(*approxSpace.surface_dof_distribution(),vvConnection,minpind,bseparate,orderp,orderv);
	}
}


} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__CR_REORDER__ */
