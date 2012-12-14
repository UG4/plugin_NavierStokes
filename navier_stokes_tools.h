/*
 * navier_stokes_tools.h
 *
 *  Created on: 13.12.2012
 *      Author: Christian Wehner
 */

#ifndef __H__UG__LIB_DISC__NAVIER_STOKES_TOOLS__
#define __H__UG__LIB_DISC__NAVIER_STOKES_TOOLS__

#include <vector>
#include "lib_disc/function_spaces/approximation_space.h"
#include "common/profiler/profiler.h"

namespace ug{

/// orders the all DofDistributions of the ApproximationSpace using Cuthill-McKee
template <typename TGridFunction>
void vorticity(TGridFunction& vort,TGridFunction& u)
{
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;

	///	grid type
	typedef typename domain_type::grid_type grid_type;

	///	world dimension
	static const int dim = domain_type::dim;

	// get grid
	grid_type& grid = *u.domain()->grid();

	// position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;

	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object elem_type;

    /// side type
	typedef typename elem_type::side side_type;

	//  volume attachment
	typedef typename Grid::AttachmentAccessor<side_type,ANumber > aSideNumber;
	aSideNumber m_acVolume;
	ANumber m_aVolume;

	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// side iterator
	typedef typename TGridFunction::template traits<side_type>::const_iterator SideIterator;

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = u.domain()->position_accessor();

	DimCRFVGeometry<dim> geo;

	grid.template attach_to<side_type>(m_aVolume);
	m_acVolume.access(grid,m_aVolume);

	SetAttachmentValues(m_acVolume,grid.template begin<side_type>(), grid.template end<side_type>(), 0);

	if (vort.local_finite_element_id(0) != LFEID(LFEID::CROUZEIX_RAVIART, 1)){
				UG_THROW("Component " << 0 << " in approximation space of parameter 1 must be of Crouzeix-Raviart type.");
	}

	for (int d=0;d<dim;d++){
		if (u.local_finite_element_id(d) != LFEID(LFEID::CROUZEIX_RAVIART, 1)){
			UG_THROW("Component " << d << " in approximation space of parameter 2 must be of Crouzeix-Raviart type.");
		}
	}

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	//	create Multiindex
	std::vector<MultiIndex<2> > multInd;

	for(int si = 0; si < u.num_subsets(); ++si)
	{
		if (si>0) continue;
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
				// UG_LOG("co_coord(" << i<< "+1,:)=" << coCoord[i] << "\n");
			}
			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), u.domain()->subset_handler().get());

			size_t nofsides = geo.num_scv();

			typename grid_type::template traits<side_type>::secure_container sides;

			UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Only elements of type elem_type are currently supported");

			grid.associated_elements_sorted(sides, static_cast<elem_type*>(elem) );

			static const size_t MaxNumSidesOfElem = 18;

			typedef MathVector<dim> MVD;
			std::vector<MVD> uValue(MaxNumSidesOfElem);

			for (size_t s=0;s < nofsides;s++)
			{
				for (int d=0;d<dim;d++){
					//	get indices of function fct on vertex
					u.multi_indices(sides[s], d, multInd);
					//	read value of index from vector
					uValue[s][d]=DoFRef(u,multInd[0]);
				}
			}

			for (size_t s=0;s < nofsides;s++)
			{
				number localvort = 0;
				const typename DimCRFVGeometry<dim>::SCV& scv = geo.scv(s);

				for(size_t sh = 0 ; sh < nofsides; ++sh){
					localvort += uValue[sh][1] * scv.global_grad(sh)[0] - uValue[sh][0] * scv.global_grad(sh)[1];
				};

				number vol = scv.volume();

				localvort*=vol;

				vort.multi_indices(sides[s], 0, multInd);
				DoFRef(vort,multInd[0])+=localvort;
				m_acVolume[sides[s]] += vol;
			}
		}
		// average vorticity
		SideIterator sideIter = vort.template begin<side_type>(si);
		SideIterator sideIterEnd = vort.template end<side_type>(si);
		number maxvort = 0;
		for(  ;sideIter !=sideIterEnd; sideIter++)
		{
			//	get Elem
			side_type* elem = *sideIter;
			vort.multi_indices(elem, 0, multInd);
			DoFRef(vort,multInd[0])/=m_acVolume[elem];
			if (DoFRef(vort,multInd[0])*DoFRef(vort,multInd[0])>maxvort*maxvort){
				maxvort = DoFRef(vort,multInd[0]);
			}
		}
	}
	grid.template detach_from<side_type>(m_aVolume);
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__NAVIER_STOKES_TOOLS__ */
