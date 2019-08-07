/*
 * diffusion_interface_handler_local_impl.h
 *
 *  Created on: 15.01.2015
 *      Author: suze
 */

#ifndef INTERFACE_HANDLER_2PF_IMPL_
#define INTERFACE_HANDLER_2PF_IMPL_

#include <math.h>

namespace ug{
namespace NavierStokes{


//hier methods aus moving_particle_impl.h einfuegen

// call constructor of base class
template<int TWorldDim>
InterfaceHandlerLocal2PF<TWorldDim>::InterfaceHandlerLocal2PF(
		SmartPtr<DiffusionInterfaceProvider<dim> > interfaceProvider,
		SmartPtr<CutElementHandlerImmersed<dim> > cutElementHandler,
		number fluidDensity, number fluidKinVisc) :
		InterfaceHandlerLocalDiffusion<TWorldDim>(interfaceProvider, cutElementHandler),
 		m_fluidDensity(fluidDensity),
 		m_fluidKinVisc(fluidKinVisc),
 		m_numFct(0),
 		m_numCo(0),
 		m_shift_DoFIndex_tri(false),
 		m_shift_DoFIndex_quad(false)
{
}
;

template<int TWorldDim>
void InterfaceHandlerLocal2PF<TWorldDim>::
write_solution(const std::vector<double > verticesValues)
{
	this->m_verticesValue.clear();

	for (size_t i = 0; i < verticesValues.size(); ++i)
 		this->m_verticesValue.push_back(verticesValues[i]);


// through error:
	if ( this->m_verticesValue.size() != verticesValues.size() )
	{
		UG_LOG("m_verticesValue.size(): " << this->m_verticesValue.size() << "\n");
		UG_LOG("verticesValues.size(): " << verticesValues.size() << "\n");
		UG_THROW("in InterfaceHandlerLocal2PF::write_solution: wrong size of m_verticesValue!\n");
	}

}

template<int TWorldDim>
void InterfaceHandlerLocal2PF<TWorldDim>::
reset_defect_on_interface(LocalVector& locD, const size_t size)
{
	if ( size > locD.num_all_dof(0) )
    {
        UG_LOG("in 'reset_defect_on_interface()': size = " << size << ", locD.size = " << locD.num_all_dof(0) << "\n");
		UG_THROW("in 'reset_defect_on_interface()': size = " << size << ", locD.size = " << locD.num_all_dof(0) << " => claimed size is NOT equal to size of solution vector!\n");
    }
// loop and set solution in 'solU_tri':
	for (size_t dof = 0; dof < locD.num_all_dof(0); ++dof)
	{
	// if dof_real is index of m_vertex: get solution from class:
		if ( this->lies_onInterface_size(dof, size) )
			locD.value(0, dof) = 0.0;
	}
}

template<int TWorldDim>
void InterfaceHandlerLocal2PF<TWorldDim>::
reset_jacobian_on_interface(LocalMatrix& locJ, const size_t size)
{

	if ( size > locJ.num_all_row_dof(0) )
    {
        UG_LOG("in 'reset_jacobian_on_interface()': size = " << size << ", locJ.size = " << locJ.num_all_row_dof(0) << "\n");
		UG_THROW("in 'reset_jacobian_on_interface()': size = " << size << ", locJ.size = " << locJ.num_all_row_dof(0) << " => claimed size is NOT equal to size of solution vector!\n");
    }
// loop and set solution in 'solU_tri':
	for (size_t dof1 = 0; dof1 < locJ.num_all_row_dof(0); ++dof1)
	{
		if ( this->lies_onInterface_size(dof1, size) )
		{
		// erase all col-values of chosen row dof1:
			for (size_t dof2 = 0; dof2 < locJ.num_all_col_dof(0); ++dof2)
				locJ.value(0, dof1, 0, dof2) = 0.0;
		}

	}

}

template<int TWorldDim>
double InterfaceHandlerLocal2PF<TWorldDim>::
get_jump_value(const MathVector<dim> position)
{
	return 0.0;

	double absValue = position[0]*position[0] + position[1]* position[1];
	double sum = position[0] + position[1];

	double returnValue = log(absValue) - sin(sum);

	return returnValue;
}

template<int TWorldDim>
double InterfaceHandlerLocal2PF<TWorldDim>::
get_jump_value_ex3(const MathVector<dim> position)
{

	//return 0.0;

	if ( this->get_orientation() == 1)
		return 0.0;

	double absValue = position[0]*position[0] + position[1]* position[1];
	double returnValue = exp(-absValue);

	return returnValue;
}

template<int TWorldDim>
double InterfaceHandlerLocal2PF<TWorldDim>::
get_jump_value_const(const MathVector<dim> position)
{
	return 0.0;

	double absValue = position[0]*position[0] + position[1]* position[1];
	double factor = 8*(2*absValue - position[0] - position[1]);

	return factor*exp(-absValue);

	double sum = position[0] + position[1];

	double returnValue = log(absValue) - sin(sum);

	return returnValue;
}

template<int TWorldDim>
double InterfaceHandlerLocal2PF<TWorldDim>::
get_jump_grad_value_ex3(const MathVector<dim> position)
{
	//return 0.0;

	if ( this->get_orientation() == -1)
		return 0.0;

	double absValue = position[0]*position[0] + position[1]* position[1];
	double sum = position[0] + position[1];

	double returnValue = 8*(2*absValue - sum)*exp(-absValue);

	return returnValue;
}

template<int TWorldDim>
double InterfaceHandlerLocal2PF<TWorldDim>::
get_jump_grad_value(const MathVector<dim> position)
{
	if ( this->get_orientation() == 1)
		return 2.0;
	else
		return 0.0;

	double absValue = position[0]*position[0] + position[1]* position[1];

	return -exp(-absValue);

	double sum = position[0] + position[1];
	MathVector<dim> normal;
	normal[0] = position[0];
	normal[1] = position[1];

	double returnValue = 2*(sin(sum) + 2)*(position[0]*normal[0] + position[1]*normal[1])/absValue;
	returnValue -= cos(sum)*(cos(sum)+2)*(normal[0] + normal[1]);

	return returnValue;
}



template<int TWorldDim>
double InterfaceHandlerLocal2PF<TWorldDim>::
get_source_kappa(const MathVector<dim> position)
{
	if ( this->get_orientation() == 1)
	{
		return  16.0*16.0;
	}
	else
	{
		if ( this->get_orientation() != -1)
			UG_THROW("wrong orientation!\n");

		double dist_x = position[0] - 0.1;
		double dist_y = position[1] - 0.2;
	 	double dist = sqrt(dist_x*dist_x+dist_y*dist_y);

	 	return 200*16*dist*dist;

	}

}

template<int TWorldDim>
double InterfaceHandlerLocal2PF<TWorldDim>::
get_source(const MathVector<dim> position)
{
	double absValue = position[0]*position[0] + position[1]* position[1];
	double sum = position[0] + position[1];
//	double dist = sqrt(absValue);

	double returnValue = 0.0;
	if ( this->get_orientation() == 1)
		returnValue = 2*sum*cos(sum)/absValue;
	else
	{
		returnValue = -4*sin(sum)*(cos(sum)+1);

		if ( this->get_orientation() != -1)
			UG_THROW("wrong orientation!\n");
	}
	return 2.0; //returnValue;
}

template<int TWorldDim>
LocalVector InterfaceHandlerLocal2PF<TWorldDim>::
set_source(LocalIndices ind, const size_t size, const bool bElementIsCut)
{
	LocalVector source;
	ind.resize_dof(0, size);
	source.resize(ind);

// loop and set solution in 'solU_tri':
	for (size_t dof = 0; dof < source.num_all_dof(0); ++dof)
	{
/*		if ( dof > 2 )
			source.value(0, dof) = sourceIm[2]; // done during 'add_def_elem_local()'
		else
			source.value(0, dof) = sourceIm[dof]; // done during 'add_def_elem_local()'
*/


	// if dof_real is index of m_vertex: get solution from class:
		if ( this->lies_onInterface_size(dof, size) )
		{
			size_t dof_real = this->real_index_size(dof, size);
 			source.value(0, dof) = get_source_kappa(this->get_VerticesPos(dof_real));
		}
		else
		{
/*
            size_t dof_orig = this->corner_orig(dof);
            source.value(0, dof) = sourceIm[dof_orig]; // done during 'add_def_elem_local()'
*/
		// careful: this->corner() returns corner coords of new corners => use dof-Index, NOT dof_orig-Index?!
 			source.value(0, dof) = get_source_kappa(this->corner(dof));

 		}

	}

	return source;
}


template<int TWorldDim>
LocalVector InterfaceHandlerLocal2PF<TWorldDim>::
set_jump_values(LocalIndices ind, const size_t size)
{
	LocalVector jump;
	ind.resize_dof(0, size);
	jump.resize(ind);

// loop and set solution in 'solU_tri':
	for (size_t dof = 0; dof < jump.num_all_dof(0); ++dof)
	{

	// if dof_real is index of m_vertex: get solution from class:
		if ( this->lies_onInterface_size(dof, size) )
		{
			size_t dof_real = this->real_index_size(dof, size);
 			jump.value(0, dof) = get_jump_value_ex3(this->get_VerticesPos(dof_real));
		}
		else
		{
			jump.value(0, dof) = 0.0;
		}
	}

	return jump;
}


template<int TWorldDim>
LocalVector InterfaceHandlerLocal2PF<TWorldDim>::
set_jump_grad_values(LocalIndices ind, const size_t size)
{
	LocalVector jump_grad;
	ind.resize_dof(0, size);
	jump_grad.resize(ind);

// loop and set solution in 'solU_tri':
	for (size_t dof = 0; dof < jump_grad.num_all_dof(0); ++dof)
	{

	// if dof_real is index of m_vertex: get solution from class:
		if ( this->lies_onInterface_size(dof, size) )
		{
			size_t dof_real = this->real_index_size(dof, size);
 			jump_grad.value(0, dof) = get_jump_grad_value_ex3(this->get_VerticesPos(dof_real));
		}
		else
		{
			jump_grad.value(0, dof) = 0.0;
		}
	}

	return jump_grad;

}


template<int TWorldDim>
void InterfaceHandlerLocal2PF<TWorldDim>::
set_local_sol(LocalVector& solU, const size_t size, const LocalVector& lvec, const int orientation)
{
	if ( size > solU.num_all_dof(0) )
    {
        UG_LOG("in 'set_local_sol()': size = " << size << ", solU.size = " << solU.num_all_dof(0) << "\n");
		UG_THROW("in 'set_local_sol()': size = " << size << ", solU.size = " << solU.num_all_dof(0) << " => claimed size is NOT equal to size of solution vector!\n");
    }
// loop and set solution in 'solU_tri':
 	for(size_t fct=0; fct < solU.num_all_fct(); ++fct)
 		for (size_t dof = 0; dof < solU.num_all_dof(fct); ++dof)
 		{
 			size_t dof_real = this->real_index_size(dof, size);

 		// if dof_real is index of m_vertex: get solution from class:
 			if ( this->lies_onInterface_size(dof, size) )
 			{
 				solU.value(fct, dof) = this->get_sol(dof_real);
 			}
 			else
 			{
 				solU.value(fct, dof) = lvec.value(fct,dof_real);
 			}
 		}

}

/*
template <int TWorldDim>
int InterfaceHandlerLocal2PF<TWorldDim>::
CollectCorners_FlatTop_2d(GridObject* elem)
{
	//////////////////////////////////////////////
	// 1) fill vector with fluid corners:
	//////////////////////////////////////////////

	this->m_vCornerCoords.clear();
	this->m_vInterfaceID.clear();
	this->m_vOriginalCornerID.clear();

// buffer vectors for (cornerCoords, cornerIndex)
	std::vector<std::pair<MathVector<dim>, size_t > > vOutsideCorners;
	std::vector<std::pair<MathVector<dim>, size_t > > vInsideCorners;
	std::vector<std::pair<MathVector<dim>, size_t > > vNearIntCorners;

//	collect all vertices of the element
	std::vector<Vertex*> vVertex;
	CollectVertices(vVertex, *this->m_spMG, elem);

 	bool isFTVertex = false;
 	for(size_t i = 0; i < vVertex.size(); ++i)
	{
	// remember boolian for check, weather there existe at least one vertex, which is FT!
		isFTVertex = this->is_FTVertex(vVertex[i], i);
		if ( isFTVertex )
			break;
	}

	if ( !isFTVertex )
		UG_THROW("Error in 'CollectCorners_FlatTop_2d': no vertex is FTVertex: should be true for at least 1 vertex!\n");

	//	collect all edges of the element
	std::vector<Edge*> vEdges;
	CollectEdgesSorted(vEdges, *this->m_spMG, elem);

	// loop vertices
	//////////////////////////////////////////////
	// REMARK:
	// order is the same as in 'vCornerCoords', therefore we can be sure, that the
	// order of the new 'vCornerIBCoords' will be consistent with the grid standard
	//////////////////////////////////////////////

	bool bNearInterface = false;
	for(size_t i = 0; i < vVertex.size(); ++i)
	{
		Vertex* vrtRoot = vVertex[i];

		//////////////////////////////////////////////
		// case 1:
		// vertex insideDomain
 		if ( !this->is_FTVertex(vrtRoot, i) )
 		{
			if ( this->is_nearInterfaceVertex(vrtRoot, i) )
				UG_THROW("NearInterface BUT !is_FT => neuerdings Fehler!!....\n");

			this->m_vCornerCoords.push_back(this->m_aaPos[vrtRoot]);
			this->m_vOriginalCornerID.push_back(i);

			vInsideCorners.push_back(std::make_pair(this->m_aaPos[vrtRoot], i));
		}
		//////////////////////////////////////////////
  		// case 2:
		// vertex = FT + ON interface
		// 		=> KEINE Berechnung von 'intersectionPoint' notwendig! -> pushen und alten index pushen

		// REMARK: is_nearInterfaceVerx = false per default, if m_vThresholdOnLevel = 0.0
		else if ( this->is_nearInterfaceVertex(vrtRoot, i) )
		{
 			bNearInterface = true;
 			this->m_vCornerCoords.push_back(this->m_aaPos[vrtRoot]);
 			this->m_vOriginalCornerID.push_back(i);
 			this->m_vInterfaceID.push_back(this->m_vCornerCoords.size()-1);  // attention: push AFTER 'm_vCornerCoords.push_back()'!!

			vOutsideCorners.push_back(std::make_pair(this->m_aaPos[vrtRoot], i));
			vNearIntCorners.push_back(std::make_pair(this->m_aaPos[vrtRoot], i));

		}
		//////////////////////////////////////////////
  		// case 3:
		// vertex 'outsideFluid'
		// 		=> NEUE Position berechen+pushen und alten index pushen
		else
		{
 			//////////////////////////////////////////////////////////////////////////////////////////
			// loop alle edges, die interface schneiden und damit einen neuen intersectionPnt
			// beitragen zum damit assoziierten alten index
			for(size_t e = 0; e < vEdges.size(); ++e)
			{
				Edge* edge = vEdges[e];
				std::vector<Vertex*> vVertexEdge;
				CollectVertices(vVertexEdge, *this->m_spMG, edge);
				if ( vVertexEdge.size() != 2 )
					UG_THROW("error in collecting vertices associated to an edge!....EXIT!...\n");

				Vertex* vrt1 = vVertexEdge[0];
				Vertex* vrt2 = vVertexEdge[1];
				size_t vrtInd1 = get_vertex_index(vrt1, elem);
				size_t vrtInd2 = get_vertex_index(vrt2, elem);

				MathVector<dim> intersectionPnt;

	 		///////////////////////////////////////////////////////////////////
 			// lies vrtRoot on a cutted edge?
		 	///////////////////////////////////////////////////////////////////
			// case1: vrtRoot is intersectionPnt with insideCorner = near_interface_corner => remove!
				if ( this->is_nearInterfaceVertex(vrt2, vrtInd2) || this->is_nearInterfaceVertex(vrt1, vrtInd1) )
				{ bNearInterface = true; continue; }
			 // case2: vert2 = outsideParticle && vrt1 = insideParticle:
				else if ( vrtRoot == vrt1 && !this->is_FTVertex(vrt2, vrtInd2) ){
					this->get_intersection_point(intersectionPnt, vrt2, vrt1);
 				}
			// case3: vrt1 = outsideParticle && vrt2 = insideParticle:
				else if ( vrtRoot == vrt2 && !this->is_FTVertex(vrt1, vrtInd1) )
					this->get_intersection_point(intersectionPnt, vrt1, vrt2);
				else
 					continue;

			// check for correct inersectionPnt
				if ( fabs(this->get_LSvalue_byPosition(intersectionPnt)) > 1e-6  )
					UG_THROW("in 'CollectIBCorners2d()': Error in computation of 'intersectionPnt':\n "
							" intersectionPnt = " << intersectionPnt << "\n distance from interace = " << fabs(get_LSvalue_byPosition(intersectionPnt)) << "\n");

	 		///////////////////////////////////////////////////////////////////
	 		// only push_back, if not included yet!
			// 	-> can be ONLY the case, if the intersectionPoint is a node
	 			if ( ! this->isIncluded(this->m_vCornerCoords, intersectionPnt) )
	 			{

	 				this->m_vCornerCoords.push_back(intersectionPnt);
	 				this->m_vOriginalCornerID.push_back(i);
	 				this->m_vInterfaceID.push_back(this->m_vCornerCoords.size()-1);  // attention: push AFTER 'm_vCornerCoords.push_back()'!!

	 				vOutsideCorners.push_back(std::make_pair(intersectionPnt, i));
   	 			}


 			} // end edge-loop

		} // end else-case

 	} // end vrt-loop

////////////////////////////////////////////////////////////////////////////////////////////
// Postprecessing for quadrilaterals ( <=>  vOutsideCorners == 2 )
// (vInsideCorners.size() == 2) && (bNearInterface)	 => ALL nodes insideFluid, BUT one ON surface
//		=> no Quadrilateral, but Triangle!!
////////////////////////////////////////////////////////////////////////////////////////////
	MathVector<dim> normalDir(0.0);
	if ( (this->m_vCornerCoords.size() == 4) && (!bNearInterface) && (dim == 2) )
		this->ResortQuadrilateral(vInsideCorners, vOutsideCorners, normalDir);
	else if ( bNearInterface )
	{
	// Quadrilateral -> Triangle
		if ( vInsideCorners.size() == 1 ) // case 1
		{
			// do nothing, since re-sorting not necessary...???
		}
	// skip whole element, since only FT points are included
		else if ( vInsideCorners.size() == 0 )
			UG_THROW("in 'CollectCorners_FlatTop_2d()': vInsideCorners.size() "
					"= " << vInsideCorners.size() << "not possible!\n");
	}


	return this->m_vCornerCoords.size();

}


// called by geo.update()!!
template <int TWorldDim>
bool InterfaceHandlerLocal2PF<TWorldDim>::
update_elem(GridObject* elem, const MathVector<TWorldDim>* vCornerCoords, int interfaceOrientation)
{
	bool do_update_local = false;
	this->m_vBF.clear();

// computing flat top modus
	this->m_elemModus = compute_element_modus(elem, interfaceOrientation);

	switch(this->m_elemModus)
	{
 		case INSIDE_DOM:	   if ( dim == 2 ) this->set_flat_top_data(elem, vCornerCoords, ROID_TRIANGLE);
 							   if ( dim == 3 ) this->set_flat_top_data(elem, vCornerCoords, ROID_TETRAHEDRON);
							   break;	// usual assembling
		case OUTSIDE_DOM: 	   if ( dim == 2 ) this->set_flat_top_data(elem, vCornerCoords, ROID_TRIANGLE);
							   if ( dim == 3 ) this->set_flat_top_data(elem, vCornerCoords, ROID_TETRAHEDRON);
							   break;	// usual assembling
		case CUT_BY_INTERFACE: this->compute_flat_top_data(elem);
								//if ( m_roid == ROID_PYRAMID )UG_THROW("PYRAMID\n");
							   do_update_local = true;
						  	   break;  // flat top assembling
		default:
			throw(UGError("Error in InterfaceHandlerLocalDiffusion::update(): switch(m_elemModus)!"));
	}

 	return do_update_local;

}
*/

} // namespace NavierStokes
} // end ug namespace

#endif /* INTERFACE_HANDLER_2PF_IMPL_ */
