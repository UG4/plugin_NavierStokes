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
#include "lib_grid/lg_base.h"
#include "lib_grid/tools/periodic_boundary_manager.h"
#include "common/profiler/profiler.h"

namespace ug{

// compute vorticity
// for velocity field (u,v,w) the vorticity is \partial_x v - \partial_y u
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
	typedef PeriodicAttachmentAccessor<side_type,ANumber > aSideNumber;
	aSideNumber m_acVolume;
	ANumber m_aVolume;

	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// side iterator
	typedef typename TGridFunction::template traits<side_type>::const_iterator SideIterator;

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	position_accessor_type& posAcc = u.domain()->position_accessor();

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

	PeriodicBoundaryManager* pbm = grid.periodic_boundary_manager();

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
			// if periodic slave continue
			if (pbm && pbm->is_slave(elem)) continue;
			vort.multi_indices(elem, 0, multInd);
			DoFRef(vort,multInd[0])/=m_acVolume[elem];
			if (DoFRef(vort,multInd[0])*DoFRef(vort,multInd[0])>maxvort*maxvort){
				maxvort = DoFRef(vort,multInd[0]);
			}
		}
	}
	grid.template detach_from<side_type>(m_aVolume);
}

// array0 = array1
void copyGhiaNumbers(number array0[17],number array1[17]){
	for (size_t i=0;i<17;i++) array0[i]=array1[i];
}

template <typename TGridFunction>
void drivenCavityEvaluation(TGridFunction& u,size_t Re){
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;

	///	grid type
	typedef typename domain_type::grid_type grid_type;

	///	world dimension
	static const int dim = domain_type::dim;

	// position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;

	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object elem_type;

	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	position_accessor_type& posAcc = u.domain()->position_accessor();

	bool foundRe = false;

	// data from Ghia paper, see also cast3m implementation
	number yco[17]={0.0000,0.0547,0.0625,0.0703,0.1016,0.1719,0.2813,0.4531,0.5000,0.6172,0.7344,0.8516,0.9531,0.9609,0.9688,0.9766,1.0000};
	number xco[17]={0.0000,0.0625,0.0703,0.0781,0.0938,0.1563,0.2266,0.2344,0.5000,0.8047,0.8594,0.9063,0.9453,0.9531,0.9609,0.9688,1.0000};
	number xLineReferenceValue[17];
	number yLineReferenceValue[17];
	number xLineReferenceValue100[17]={0.,-0.03717,-0.04192,-0.04775,-0.06434,-0.10150,-0.15662,-0.2109,-0.20581,-0.13641,0.00332,0.23151,0.68717,
						0.73722,0.78871,0.84123,1.};
	number yLineReferenceValue100[17]={0.00000,0.09233,0.10091,0.10890,0.12317,0.16077,0.17507
						,0.17527
						,0.05454,-0.24533,-0.22445,-0.16914,-0.10313,-0.08864
						,-0.07391,-0.05906,0.00000};
	number xLineReferenceValue400[17]={0.,-0.08186,-0.09266,-0.10338,-0.14612,-0.24299,-0.32726,-0.17119,-0.11477,0.02135,0.16256,0.29093,0.55892,
					0.61756,0.68439,0.75837,1.};
				number yLineReferenceValue400[17]={0.00000,0.18360,0.19713,0.20920,0.22965,0.28124,0.30203
						,0.30174
						,0.05186,-0.38598,-0.44993,-0.3827,-0.22847,-0.19254
						,-0.15663,-0.12146,0.00000};
			number xLineReferenceValue1000[17]={0.,-0.18109,-0.20196,-0.22220,-0.29730,-0.38289,-0.27805,-0.10648,-0.06080,0.05702,0.18719,0.33304,0.46604,
					0.51117,0.57492,0.65928,1.};
			number yLineReferenceValue1000[17]={0.00000,0.27485,0.29012,0.30353,0.32627,0.37095,0.33075
					,0.32235
					,0.02526,-0.31966,-0.42665,-0.51550,-0.39188,-0.33714
					,-0.27669,-0.21388,0.00000};
			number xLineReferenceValue3200[17]={0.00000,-0.32407,-0.35344,-0.37827,-0.41933,-0.34323
					,-0.24427,-0.086636
					,-0.04272,0.07156,0.19791,0.34682,0.46101,0.46547,0.48296
					,0.53236,1.00000};
			number yLineReferenceValue3200[17]={0.00000,0.39560,0.40917,0.41906,0.42768,0.37119,0.29030
					,0.28188
					,0.00999,-0.31184,-0.37401,-0.44307,-0.54053,-0.52357
					,-0.47425,-0.39017,0.00000};
		number xLineReferenceValue5000[17]={0.00000,-0.41165,-0.42901,-0.43643,-0.40435,-0.33050
				,-0.22855,-0.07404
				,-0.03039,0.08183,0.20087,0.33556,0.46036,0.45992,0.46120
				,0.48223,1.00000};
		number yLineReferenceValue5000[17]={0.00000,0.42447,0.43329,0.43648,0.42951,0.35368,0.28066
				,0.27280
				,0.00945,-0.30018,-0.36214,-0.41442,-0.52876,-0.55408
				,-0.55069,-0.49774,0.00000};
		number xLineReferenceValue7500[17]={0.00000,-0.43154,-0.43590,-0.43025,-0.38324,-0.32393
				,-0.23176,-0.07503
				,-0.03800,0.08342,0.20591,0.34228,0.47167,0.47323,0.47048
				,0.47244,1.00000};
		number yLineReferenceValue7500[17]={0.00000,0.43979,0.44030,0.43564,0.41824,0.35060,0.28117
				,0.27348
				,0.00824,-0.30448,-0.36213,-0.41050,-0.48590,-0.52347
				,-0.55216,-0.53858,0.00000};
		number xLineReferenceValue10000[17]={0.00000,-0.42735,-0.42537,-0.41657,-0.38000,-0.32709
				,-0.23186,-0.07540
				,-0.03111,0.08344,0.20673,0.34635,0.47804,0.48070,0.47783
				,0.47221,1.00000};
		number yLineReferenceValue10000[17]={0.00000,0.43983,0.43733,0.43124,0.41487,0.35070,0.28003
				,0.27224
				,0.00831,-0.30719,-0.36737,-0.41496,-0.45863,-0.49099
				,-0.52987,-0.54302,0.00000};

	switch (Re){
	case 100:
		foundRe=true;
		copyGhiaNumbers(xLineReferenceValue,xLineReferenceValue100);
		copyGhiaNumbers(yLineReferenceValue,yLineReferenceValue100);
		break;
	case 400:
		foundRe=true;
		copyGhiaNumbers(xLineReferenceValue,xLineReferenceValue400);
		copyGhiaNumbers(yLineReferenceValue,yLineReferenceValue400);
		break;
	case 1000:
		foundRe=true;
		copyGhiaNumbers(xLineReferenceValue,xLineReferenceValue1000);
		copyGhiaNumbers(yLineReferenceValue,yLineReferenceValue1000);
		break;
	case 3200:
		foundRe=true;
		copyGhiaNumbers(xLineReferenceValue,xLineReferenceValue3200);
		copyGhiaNumbers(yLineReferenceValue,yLineReferenceValue3200);
		break;
	case 5000:
		foundRe=true;
		copyGhiaNumbers(xLineReferenceValue,xLineReferenceValue5000);
		copyGhiaNumbers(yLineReferenceValue,yLineReferenceValue5000);
		break;
	case 7500:
		foundRe=true;
		copyGhiaNumbers(xLineReferenceValue,xLineReferenceValue7500);
		copyGhiaNumbers(yLineReferenceValue,yLineReferenceValue7500);
		break;
	case 10000:
		foundRe=true;
		copyGhiaNumbers(xLineReferenceValue,xLineReferenceValue10000);
		copyGhiaNumbers(yLineReferenceValue,yLineReferenceValue10000);
		break;
	}

	if (foundRe==false){
		UG_THROW("Reynolds number " << Re << " not supported.\n");
	}

	MathVector<dim> xLineEvalPos;
	MathVector<dim> yLineEvalPos;
	xLineEvalPos[0]=0.5;
	yLineEvalPos[1]=0.5;

	number xLineValue[17];
	number yLineValue[17];

	static const size_t unhandled = 1000000;

	for (size_t i=0;i<17;i++){
		xLineValue[i]=unhandled;
		yLineValue[i]=unhandled;
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

		static const number bigNumber = 1e+8;

		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			elem_type* elem = *iter;
			//	get vertices and extract corner coordinates
			const size_t numVertices = elem->num_vertices();
			number xmin=bigNumber,xmax=-bigNumber,ymin=bigNumber,ymax=-bigNumber;
			//UG_LOG("-----------------------------------------------------------------\n");
			// find element bounds
			for(size_t i = 0; i < numVertices; ++i){
				vVrt[i] = elem->vertex(i);
				coCoord[i] = posAcc[vVrt[i]];
				//UG_LOG("co(" << i+1 << ",:)=[" << coCoord[i][0] << "," << coCoord[i][1] << "];" << "\n");
				if (coCoord[i][0]<xmin) xmin=coCoord[i][0];
				if (coCoord[i][0]>xmax) xmax=coCoord[i][0];
				if (coCoord[i][1]<ymin) ymin=coCoord[i][1];
				if (coCoord[i][1]>ymax) ymax=coCoord[i][1];
			};
			//UG_LOG("co(" << numVertices+1 << ",:)=[" << coCoord[0][0] << "," << coCoord[0][1] << "];" << "\n");

			//UG_LOG("bound=[" << xmin << " " << xmax << " " << ymin << " " << ymax << "]\n");

			bool checkx=true;
			bool checky=true;

			if (xmin>0.5) checkx=false;
			else if (xmax<0.5) checkx=false;

			if (ymin>0.5) checky=false;
			else if (ymax<0.5) checky=false;

			if ((checkx==false)&&(checky==false)) continue;

			number interpolation;

			if (checkx==true){
				for (size_t i=0;i<17;i++)
				{
					if ((yco[i]>=ymin)&&(yco[i]<=ymax)){
						xLineEvalPos[1]=yco[i];
						//UG_LOG("p=[" << xLineEvalPos[0] << "," << xLineEvalPos[1] << "]\n");
						if ( ContainsPoint(elem, xLineEvalPos, posAcc) )
						{
							//UG_LOG("***************************************************\n");
							//UG_LOG(xLineEvalPos << "\n");
							if (yco[i]==0){
								//UG_LOG("##################################################################\n");
								//UG_LOG("##################################################################\n");
								//UG_LOG("##################################################################\n");
								//UG_LOG("##################################################################\n");
							}
							//	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
							ReferenceObjectID roid = (ReferenceObjectID) elem->reference_object_id();
							//	get the reference mapping for the element using global corners
							DimReferenceMapping<dim, dim>& mapping
								= ReferenceMappingProvider::get<dim, dim>(roid, coCoord);
							MathVector<dim> localPos;

							//	compute global integration points
							mapping.global_to_local(localPos, xLineEvalPos);

							//	local finite element id
							const LFEID m_id = u.local_finite_element_id(0);

							//	get trial space
							const LocalShapeFunctionSet<dim>& rTrialSpace =
								LocalShapeFunctionSetProvider::get<dim>(roid, m_id);

							//	number of dofs on element
							const size_t num_sh = rTrialSpace.num_sh();

							//	get multiindices of element
							std::vector<MultiIndex<2> > ind;  // 	aux. index array
							u.multi_indices(elem, 0, ind);
							// 	compute approximated solution at integration point
							interpolation = 0.0;
							for(size_t sh = 0; sh < num_sh; ++sh)
							{
								//	get value at shape point (e.g. corner for P1 fct)
								const number valSH = DoFRef(u, ind[sh]);

								//	add shape fct at ip * value at shape
								interpolation += valSH * rTrialSpace.shape(sh, localPos);
							}
							// if there is a previous value, then average
							if (xLineValue[i]!=unhandled) xLineValue[i] = 0.5*(xLineValue[i] + interpolation);
							else xLineValue[i]=interpolation;
						}
					}
				}
			}

			if (checky==true){
				for (size_t i=0;i<17;i++)
				{
					if ((xco[i]>=xmin)&&(xco[i]<=xmax)){
						yLineEvalPos[0]=xco[i];
						if ( ContainsPoint(elem, yLineEvalPos, posAcc) )
						{
							//	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
							ReferenceObjectID roid = (ReferenceObjectID) elem->reference_object_id();
							//	get the reference mapping for the element using global corners
							DimReferenceMapping<dim, dim>& mapping
								= ReferenceMappingProvider::get<dim, dim>(roid, coCoord);
							MathVector<dim> localPos;

							//	compute global integration points
							mapping.global_to_local(localPos, yLineEvalPos);

							//	local finite element id
							const LFEID m_id = u.local_finite_element_id(0);

							//	get trial space
							const LocalShapeFunctionSet<dim>& rTrialSpace =
								LocalShapeFunctionSetProvider::get<dim>(roid, m_id);

							//	number of dofs on element
							const size_t num_sh = rTrialSpace.num_sh();

							//	get multiindices of element
							std::vector<MultiIndex<2> > ind;  // 	aux. index array
							u.multi_indices(elem, 1, ind);
							// 	compute approximated solution at integration point
							interpolation = 0.0;
							for(size_t sh = 0; sh < num_sh; ++sh)
							{
								//	get value at shape point (e.g. corner for P1 fct)
								const number valSH = DoFRef(u, ind[sh]);

								//	add shape fct at ip * value at shape
								interpolation += valSH * rTrialSpace.shape(sh, localPos);
							}
							// if there is a previous value, then average
							if (yLineValue[i]!=unhandled) yLineValue[i] = 0.5*(yLineValue[i] + interpolation);
							else yLineValue[i]=interpolation;
						}
					}
				}
			}
		};// for(  ;iter !=iterEnd; ++iter)
	};// for(int si = 0; si < u.num_subsets(); ++si)
	UG_LOG("\nData evaluation for Re=" << Re << ":\n\n");
	UG_LOG("u values on line through x=0.5:" << "\n\n");
	number maxdiff = 0;
	number diffsum = 0;
	for (size_t i=0;i<17;i++){
		number localdiff=abs(xLineReferenceValue[i]-xLineValue[i]);
		UG_LOG("y(" << i+1 << ") = " << yco[i] << "; u(" << i+1 << ") = " << xLineValue[i] << "; u_ghia(" << i+1 << ") = " << xLineReferenceValue[i] << "; udiff(" << i+1 << ") = " << localdiff << ";\n");
		if (localdiff>maxdiff) maxdiff = localdiff;
		diffsum += localdiff;
	}
	UG_LOG("max difference: " << maxdiff << "\naverage difference: " << (number)diffsum/17.0);
	UG_LOG("\n\nv values on line through y=0.5:" << "\n\n");
	maxdiff = 0;
	diffsum = 0;
	for (size_t i=0;i<17;i++){
		number localdiff=abs(yLineReferenceValue[i]-yLineValue[i]);
		UG_LOG("x(" << i+1 << ") = " << xco[i] << "; v(" << i+1 << ") = " << yLineValue[i] << "; v_ghia(" << i+1 << ") = " << yLineReferenceValue[i] << "; vdiff(" << i+1 << ") = " << abs(yLineReferenceValue[i]-yLineValue[i]) << ";\n");
		if (localdiff>maxdiff) maxdiff = localdiff;
		diffsum += localdiff;
	}
	UG_LOG("max difference: " << maxdiff << "\naverage difference: " << (number)diffsum/17.0);
	UG_LOG("\n\n");
}

/*
// Computation of maximum CFL-number over elements
// General CFL number is given as |u|*deltaT/h where u is velocity, deltaT time step length
// and h mesh width.
// Here following approximation for element local CFL-number is used:
// In element compute interpolated velocity v_bary in barycenter
// go over all "edges" between Crouzeix-Raviart nodes c_i and c_j (middle of edge/face)
// approximate velocity between nodes is 1/|c_i-c_j|*|(c_i-c_j)*v_bary|
// step length is |c_i-c_j| Therefore overall local cfl number is
// 1/|c_i-c_j|^2*|(c_i-c_j)*v_bary|
// Function computes maximum over all elements.
 */
template<typename TGridFunction>
void cflNumber(TGridFunction& u,number deltaT){
	/// domain type
	typedef typename TGridFunction::domain_type domain_type;

	///	world dimension
	static const int dim = domain_type::dim;

	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object elem_type;

	/// side type
	typedef typename elem_type::side side_type;

	///	grid type
	typedef typename domain_type::grid_type grid_type;

	// cfl number
	number maxCfl=0;
	//	get domain
	domain_type& domain = *u.domain().get();

	DimCRFVGeometry<dim> geo;

	//	create Multiindex
	std::vector<MultiIndex<2> > multInd;
	
	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();

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

			MathVector<dim> baryV;
			baryV = 0;
			for (size_t s=0;s<nofsides;s++){
				MathVector<dim> localbaryV;
				for (int d=0;d<dim;d++){
					u.multi_indices(sides[s], d, multInd);
					localbaryV[d]=DoFRef(u,multInd[0]);
				}
				localbaryV *= vShape[s];
				baryV += localbaryV;
			}
			for (size_t i=0;i<nofsides;i++){
				const typename DimCRFVGeometry<dim>::SCV& scvi = geo.scv(i);
				MathVector<dim> iCoord = scvi.global_ip();
				for (size_t j=i+1;j<nofsides;j++){
					const typename DimCRFVGeometry<dim>::SCV& scvj = geo.scv(j);
					MathVector<dim> jCoord = scvj.global_ip();
					MathVector<dim> subVec;
					VecSubtract(subVec,iCoord,jCoord);
					number localCfl=deltaT*(number)1.0/VecTwoNormSq(subVec)*abs(subVec*baryV);
					if (localCfl>maxCfl) maxCfl=localCfl;
				}
			}
		}
	}
	UG_LOG("Max CFL number is " << maxCfl << "\n");
}

// Compute 1/|\Omega| \int_{\Omega} \hat{u_i} \hat{u_i} dx , where \Omega is the computational domain.
// Elementwise approximation of integral by interpolation to barycenter.
template<typename TGridFunction>
void kineticEnergy(TGridFunction& u){
	// total kinetic energy
	number totalE=0;
	// total volume
	number totalVol=0;

	/// domain type
	typedef typename TGridFunction::domain_type domain_type;

	///	world dimension
	static const int dim = domain_type::dim;

	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object elem_type;

	/// side type
	typedef typename elem_type::side side_type;

	///	grid type
	typedef typename domain_type::grid_type grid_type;

	//	get domain of grid function
	domain_type& domain = *u.domain().get();
	DimCRFVGeometry<dim> geo;

	 std::vector<MultiIndex<2> > multInd;

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
		ElemIterator iter = u.template begin<elem_type>(si);
		ElemIterator iterEnd = u.template end<elem_type>(si);

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
					u.multi_indices(sides[s], d, multInd);
					localValue[d]=DoFRef(u,multInd[0]);
				}
				localValue *= vShape[s];
				value += localValue;
				elementVolume += scv.volume();
			}
			for (int d=0;d<dim;d++){
				totalE += elementVolume*value[d]*value[d];
			}

			totalVol+=elementVolume;
		}
	}
	// average
	totalE/=(number)totalVol;
	UG_LOG("Total kinetic energy in domain is " << totalE << "\n");
}


} // end namespace ug

#endif /* __H__UG__LIB_DISC__NAVIER_STOKES_TOOLS__ */
