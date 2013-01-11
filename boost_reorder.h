/*
 * boost_reorder.h
 *
 *  Created on: 29.12.2012
 *      Author: Christian Wehner
 */

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__REORDER__
#define __H__UG__LIB_DISC__DOF_MANAGER__REORDER__

#include <vector>
#include <list>
#include <queue>
#include "lib_disc/function_spaces/approximation_space.h"
#include "common/profiler/profiler.h"
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
	
template <typename TGraph>
	void get_graph(TGraph graph){}
	
//////////////////////////////////////////////////////////////////////////////////////////
//			Cuthill-McKee	
//////////////////////////////////////////////////////////////////////////////////////////
	
/// orders the dof distribution using Cuthill-McKee
template <typename TDD>
void BoostOrderCuthillMcKee(TDD& dofDistr,
						bool bReverse)
{
	PROFILE_FUNC();
	size_t n = dofDistr.num_indices();
	
	typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS,
	boost::property<boost::vertex_color_t, boost::default_color_type,
	boost::property<boost::vertex_degree_t,int> > >
	Graph;
	
	Graph G(n);
	
	//	get adjacency graph
	dofDistr.get_graph(G);

	typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
	std::vector<Vertex> inv_perm(n);
	std::vector<size_t> newInd(n);
	
	if (bReverse==true)
		cuthill_mckee_ordering(G, inv_perm.rbegin(), boost::get(boost::vertex_color, G),
                           boost::make_degree_map(G));
	else 
		cuthill_mckee_ordering(G, inv_perm.begin(), boost::get(boost::vertex_color, G),
							   boost::make_degree_map(G));
	
	for (size_t i=0;i<n;i++){
		newInd[inv_perm[i]]=i;
	};
	
	dofDistr.permute_indices(newInd);
}
	
/// orders all DofDistributions of the ApproximationSpace using Cuthill-McKee
template <typename TDomain>
void BoostOrderCuthillMcKee(ApproximationSpace<TDomain>& approxSpace,
						   bool bReverse)
{
	//	order levels
	if(approxSpace.levels_enabled())
		for(size_t lev = 0; lev < approxSpace.num_levels(); ++lev)
			BoostOrderCuthillMcKee(*approxSpace.level_dof_distribution(lev), bReverse);
		
	//	order surface
	if(approxSpace.top_surface_enabled())
		BoostOrderCuthillMcKee(*approxSpace.surface_dof_distribution(), bReverse);
}
	
	
//////////////////////////////////////////////////////////////////////////////////////////
//			Sloan	
//////////////////////////////////////////////////////////////////////////////////////////
	
/// orders the dof distribution using Sloan
template <typename TDD>
void BoostOrderSloan(TDD& dofDistr)
{
	PROFILE_FUNC();
	size_t n = dofDistr.num_indices();
		
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
	
	Graph G(n);
		
	//	get adjacency graph
	dofDistr.get_graph(G);
		
	typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
	typedef typename boost::property_map<Graph,boost::vertex_degree_t>::type deg_type;
	typedef typename boost::property_map<Graph, boost::vertex_index_t>::type map_type;
	
	//Creating two iterators over the vertices
	boost::graph_traits<Graph>::vertex_iterator ui, ui_end;
	
	//Creating a property_map with the degrees of the degrees of each vertex
	deg_type deg = get(boost::vertex_degree, G);
	for (boost::tie(ui, ui_end) = vertices(G); ui != ui_end; ++ui)
		deg[*ui] = degree(*ui, G);
	
	
	//Creating a property_map for the indices of a vertex
	map_type index_map = get(boost::vertex_index, G);
	
	std::vector<Vertex> inv_perm(n);
	std::vector<size_t> newInd(n);
	
	//Setting the start node
    Vertex s = vertex(0, G);
    int ecc;   //defining a variable for the pseudoperipheral radius
    
    //Calculating the pseudoeperipheral node and radius
    Vertex e = pseudo_peripheral_pair(G, s, ecc, boost::get(boost::vertex_color, G), boost::get(boost::vertex_degree, G) );

    //Sloan ordering
    sloan_ordering(G, s, e, inv_perm.begin(), get(boost::vertex_color, G), 
				   get(boost::vertex_degree, G), get(boost::vertex_priority, G));
	
	for (size_t i=0;i<n;i++){
		newInd[inv_perm[i]]=i;
	};
		
	dofDistr.permute_indices(newInd);
}
	
/// orders all DofDistributions of the ApproximationSpace using Sloan
template <typename TDomain>
void OrderSloan(ApproximationSpace<TDomain>& approxSpace)
{
	//	order levels
	if(approxSpace.levels_enabled())
		for(size_t lev = 0; lev < approxSpace.num_levels(); ++lev)
			BoostOrderSloan(*approxSpace.level_dof_distribution(lev));
		
	//	order surface
	if(approxSpace.top_surface_enabled())
		BoostOrderSloan(*approxSpace.surface_dof_distribution());
}

//////////////////////////////////////////////////////////////////////////////////////////
//			Minimum degree
//////////////////////////////////////////////////////////////////////////////////////////
//	6 7 8
//	3 4 5
//	0 1 2
template <typename TGraph>
void square3x3(TGraph& G)
{
	boost::add_edge(0,1,G);
	boost::add_edge(0,3,G);
	boost::add_edge(2,1,G);
	boost::add_edge(2,5,G);
	boost::add_edge(6,3,G);
	boost::add_edge(6,7,G);
	boost::add_edge(8,5,G);
	boost::add_edge(8,7,G);
	boost::add_edge(1,0,G);
	boost::add_edge(1,4,G);
	boost::add_edge(1,2,G);
	boost::add_edge(3,0,G);
	boost::add_edge(3,4,G);
	boost::add_edge(3,6,G);
	boost::add_edge(5,2,G);
	boost::add_edge(5,4,G);
	boost::add_edge(5,8,G);
	boost::add_edge(7,4,G);
	boost::add_edge(7,6,G);
	boost::add_edge(7,8,G);
	boost::add_edge(4,1,G);
	boost::add_edge(4,3,G);
	boost::add_edge(4,5,G);
	boost::add_edge(4,7,G);
	typedef typename boost::graph_traits<TGraph>::edge_iterator fuckyou;
	fuckyou fuck,fuck_end;
	boost::tie(fuck, fuck_end) = edges(G);
	if (fuck==fuck_end) UG_LOG("!!!!!what the fuck!!!!!\n");
	
}
	
/// orders the dof distribution using minimum degree algorithm
template <typename TDD>
void BoostOrderMinimumDegree(TDD& dofDistr)
{
	PROFILE_FUNC();
	size_t n = dofDistr.num_indices();
	
	int delta = 0;
		
	typedef boost::adjacency_list<boost::setS, boost::vecS, boost::directedS>  Graph;
		
	Graph G(9);
		
	//	get adjacency graph
	dofDistr.get_graph(G);
//	square3x3(G);
	
/*	typedef typename boost::graph_traits<Graph>::edge_iterator fuckyou;
	fuckyou fuck,fuck_end;
	boost::tie(fuck, fuck_end) = edges(G);
	if (fuck==fuck_end) UG_LOG("what the fuck!!!!!\n");
	
	size_t count=0;
	for (; fuck != fuck_end; ++fuck) {
		UG_LOG(count++ << " (" << source(*fuck,G) << "," << target(*fuck,G) << ")\n");
	}*/
		
	typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
	typedef typename boost::property_map<Graph, boost::vertex_index_t>::type map_type;
		
	// Creating a property_map for the indices of a vertex
	map_type index_map = get(boost::vertex_index, G);
		
	std::vector<int> inv_perm(n,0);
	std::vector<int> perm(n,0);
	std::vector<int> degree(n,0);
	std::vector<int> supernode_sizes(n,1);
	
	std::vector<size_t> newInd(n);
	
	/// minimum degree order
	minimum_degree_ordering
    (G,
     boost::make_iterator_property_map(&degree[0],index_map, degree[0]),
     inv_perm.begin(),
     perm.begin(),
     boost::make_iterator_property_map(&supernode_sizes[0],index_map, supernode_sizes[0]), 
     delta, index_map);
	
	for (size_t i=0;i<n;i++){
		newInd[i]=inv_perm[i];
//		UG_LOG(i << " " << newInd[i] << "\n");
	};
//	UG_LOG("----------------------------------\n");
	
	dofDistr.permute_indices(newInd);
}
	
/// orders all DofDistributions of the ApproximationSpace using Sloan
template <typename TDomain>
void OrderMinimumDegree(ApproximationSpace<TDomain>& approxSpace)
{
	//	order levels
	if(approxSpace.levels_enabled())
		for(size_t lev = 0; lev < approxSpace.num_levels(); ++lev)
			BoostOrderMinimumDegree(*approxSpace.level_dof_distribution(lev));
		
	//	order surface
	if(approxSpace.top_surface_enabled())
		BoostOrderMinimumDegree(*approxSpace.surface_dof_distribution());
}
	
//////////////////////////////////////////////////////////////////////////////////////////
//			King
//////////////////////////////////////////////////////////////////////////////////////////
	
/// orders the dof distribution using King
template <typename TDD>
void BoostOrderKing(TDD& dofDistr,bool reverse)
{
	PROFILE_FUNC();
	size_t n = dofDistr.num_indices();
	
	typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS,
	boost::property<boost::vertex_color_t, boost::default_color_type,
	boost::property<boost::vertex_degree_t,int> > >
	Graph;
	
	Graph G(n);
	
	//	get adjacency graph
	dofDistr.get_graph(G);
	
	typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
	std::vector<Vertex> inv_perm(n);
	std::vector<size_t> newInd(n);
	
	if (reverse==true)
		king_ordering(G, inv_perm.rbegin());
	else 
		king_ordering(G, inv_perm.begin());
	
	for (size_t i=0;i<n;i++){
		newInd[inv_perm[i]]=i;
	};
	
	dofDistr.permute_indices(newInd);
}
	
/// orders all DofDistributions of the ApproximationSpace using King
template <typename TDomain>
void OrderKing(ApproximationSpace<TDomain>& approxSpace,bool reverse)
{
	//	order levels
	if(approxSpace.levels_enabled())
		for(size_t lev = 0; lev < approxSpace.num_levels(); ++lev)
			BoostOrderKing(*approxSpace.level_dof_distribution(lev),reverse);
		
	//	order surface
	if(approxSpace.top_surface_enabled())
		BoostOrderKing(*approxSpace.surface_dof_distribution(),reverse);
}
	
	
} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__REORDER__ */
