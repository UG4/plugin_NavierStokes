/*
 * moving_particle_tools.h
 *
 *  Created on: 20.01.2015
 *      Author: suze
 */

#ifndef MOVING_PARTICLE_TOOLS_H_
#define MOVING_PARTICLE_TOOLS_H_

namespace ug{
namespace NavierStokes{


template <typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::
interpolate_point(ConstSmartPtr<DoFDistribution> dd,
		const vector_type& u,
		const MathVector<dim>& evalPos,			// input data
		MathVector<dim+1>& interpolation)		// output data
{
	MathVector<dim> localPos(0.0);

//	get data
	typename TDomain::position_accessor_type aaPos = this->m_spApproxSpace->domain()->position_accessor();

	//	create Function Group
	FunctionGroup fctGrp(dd->function_pattern(), TokenizeString(m_spCutElementHandler->m_fctNames));

	//	get iterators for all elems on subset
	typedef typename DoFDistribution::dim_traits<dim>::grid_base_object grid_base_object;

	typename DoFDistribution::traits<grid_base_object>::const_iterator iter, iterEnd;
	iter = dd->begin<grid_base_object>();
	iterEnd = dd->end<grid_base_object>();


//	loop elements in order to compute 'U_global' and 'omega_global':
    for( ; iter != iterEnd; ++iter)
	{
    //	get element
		grid_base_object* elem = *iter;

		if ( ContainsPoint(elem, evalPos, aaPos) )
		{
        //	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
			ReferenceObjectID roid = (ReferenceObjectID) elem->reference_object_id();

        //	get all corner coordinates
			std::vector<MathVector<dim> > vCorner;
			CollectCornerCoordinates(vCorner, *elem, aaPos, true);

        //	get the reference mapping for the element using global corners
			DimReferenceMapping<dim, dim>& mapping
									= ReferenceMappingProvider::get<dim, dim>(roid, vCorner);


        //	compute global integration points
			mapping.global_to_local(localPos, evalPos);

        //	loop all velocity components
			for(int cmp = 0; cmp < dim+1; ++cmp)
			{
            //	get fct id for compent
				const size_t fct = fctGrp[cmp];

            //	local finite element id
				const LFEID m_id = dd->local_finite_element_id(fct);

            //	get trial space
				const LocalShapeFunctionSet<dim>& rTrialSpace =
							LocalFiniteElementProvider::get<dim>(roid, m_id);

            //	number of dofs on element
				const size_t num_sh = rTrialSpace.num_sh();

            //	get multiindices of element
				std::vector<DoFIndex> ind;  // 	aux. index array
				dd->dof_indices(elem, fct, ind);

            //	check multi indices
				if(ind.size() != num_sh)
					UG_THROW("L2ErrorIntegrand::evaluate: Wrong number of"
							" multi indices.");

            // 	compute approximated solution at integration point
				interpolation[cmp] = 0.0;
				for(size_t sh = 0; sh < num_sh; ++sh)
				{
					//	get value at shape point (e.g. corner for P1 fct)
						const number valSH = DoFRef(u, ind[sh]);
						//UG_LOG("interpolate: valSh = " << valSH << "\n");
					//	add shape fct at ip * value at shape
						interpolation[cmp] += valSH * rTrialSpace.shape(sh, localPos);
						//UG_LOG("interpolate: interpolation[" << cmp << "] = " << interpolation[cmp] << "\n");

				}

			}  // end cmp-loop
		} // end if-ContainsPoint
	} // end element-loop

}


template <typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::
print_pressure_nodal(const vector_type& u, const int topLevel)
{
	FILE *pressure;
	char filename[40];
	sprintf(filename, "pressure_profile_level%d.txt", topLevel);
	pressure = fopen(filename, "w");
	int numIter = 0;
	number dist = 0.0;

	const int levIndex = get_Index(GridLevel(topLevel, GridLevel::LEVEL));
	std::vector<grid_base_object*> ElemListLog = m_spCutElementHandler->m_vvvElemListCut[levIndex][0];
	ConstSmartPtr<DoFDistribution> dd = this->m_spApproxSpace->dof_distribution(GridLevel(topLevel, GridLevel::LEVEL));
	typename TDomain::position_accessor_type aaPos = m_spCutElementHandler->m_aaPos;

 	const MathVector<dim>& center = m_spCutElementHandler->get_center(0);

	typedef typename std::vector<grid_base_object*>::iterator ListIter;

// loop all cut elements in order to print pressure for fluid nodes, which lie near the interface
	for(ListIter listIter = ElemListLog.begin();
		listIter != ElemListLog.end(); ++listIter)
	{
		numIter++;
	//	get element
		grid_base_object* elem = *listIter;

	//	collect all vertices of the element
		std::vector<Vertex*> vVertex;
 		CollectVertices(vVertex, *m_spCutElementHandler->m_spMG, elem);

		MathVector<dim> bufferVector;
		VecSubtract(bufferVector, aaPos[vVertex[0]], aaPos[vVertex[2]]);
		dist = 0.5*VecLength(bufferVector);

	//	loop vertices
		for(size_t v = 0; v < vVertex.size(); ++v)
		{
		//	get vertex
			Vertex* vrt = vVertex[v];

         // print value for vertices in the fluid and near the interface:
 			if ( !m_spCutElementHandler->is_outsideFluid(vrt) )
			{
 				std::vector<DoFIndex>  vInd;
 				if(dd->inner_dof_indices(vrt, dim, vInd) != 1)
					UG_THROW("Only one index expected.");

				MathVector<dim> radialVector;
				radialVector[0] = aaPos[vrt][0]-center[0];
				radialVector[1] = aaPos[vrt][1]-center[1];

				number teta = atan2(aaPos[vrt][1]-center[1], aaPos[vrt][0]-center[0]);
				if ( teta < 0.0 )
					teta += 2*3.1415926;

				fprintf(pressure, "%e \t %e \t %lu # teta, value pressure, KnotenIndex vInd[0] \n", teta, DoFRef(u, vInd[0]), vInd[0][0]);

			}
		}

	}


	fclose(pressure);

}

template <typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::
print_deltaP(const vector_type& u, const int topLevel)
{
	FILE *pressure;
	char filename[40];
	sprintf(filename, "delta_p_level%d.txt", topLevel);
	pressure = fopen(filename, "w");

   	ConstSmartPtr<DoFDistribution> dd = this->m_spApproxSpace->dof_distribution(GridLevel(topLevel, GridLevel::LEVEL));

  	MathVector<dim> printPos1;
  	MathVector<dim> printPos2;
  	printPos1[0] = 0.15;
  	printPos1[1] = 0.2;
  	printPos2[0] = 0.25;
  	printPos2[1] = 0.2;

  	MathVector<dim+1> interpolVal1;
  	MathVector<dim+1> interpolVal2;

  	interpolate_point(dd, u, printPos1, interpolVal1);
  	interpolate_point(dd, u, printPos2, interpolVal2);

  	number deltaP = interpolVal1[dim] - interpolVal2[dim];

   	fprintf(pressure, "%e \t %e \t %e \n# interpolVal1[dim], interpolVal2[dim], deltaP\n", interpolVal1[dim], interpolVal2[dim], deltaP);
	fclose(pressure);

}

template <typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::
print_pressure_teta(const vector_type& u, const int topLevel)
{
	FILE *pressure;
	char filename[40];
	sprintf(filename, "pressure_teta_level%d.txt", topLevel);
	pressure = fopen(filename, "w");

   	ConstSmartPtr<DoFDistribution> dd = this->m_spApproxSpace->dof_distribution(GridLevel(topLevel, GridLevel::LEVEL));
 	const MathVector<dim>& center = m_spCutElementHandler->get_center(0);

	size_t t_max = 500;
 	MathVector<dim> printPos;
	MathVector<dim+1> interpolVal;
	for(size_t t = 0; t < t_max; ++t)
	{
		number teta = t*(2*PI/t_max);
		number co = cos(teta);
		number si = sin(teta);
		printPos[0] = 0.2*co + center[0];
		printPos[1] = 0.2*si + center[1];
		interpolate_point(dd, u, printPos, interpolVal);

		//UG_LOG("after: interpolVal[" << dim<< "] = " << interpolVal[dim] << "\n");
 		fprintf(pressure, "%e \t", teta);
 		fprintf(pressure, "%e \t", interpolVal[dim]);
  		fprintf(pressure, "%e \t %e # teta, value_pressure, printPos[0], printPot[1]\n", printPos[0], printPos[1]);
	}

	fclose(pressure);

}

template <typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::
print_velocity(const vector_type& u, const int topLevel, number time, const char* filename)
{
// Parameter von class 'MovingParticle' extrahieren:
    const int levIndex = get_Index(GridLevel(topLevel, GridLevel::LEVEL));
    const bool isTimedep = is_time_dependent();
    
    size_t numPrt = m_spCutElementHandler->num_particles();
    
    for(size_t p = 0; p < numPrt; ++p)
    {
        
#ifdef UG_PARALLEL
    // use size of member 'CutElementHandler_FlatTop::m_vvvElemListCut' in order to indicate,
    // whether a particle lies on a processor or not
        std::vector<grid_base_object*> ElemList = m_spCutElementHandler->m_vvvElemListCut[levIndex][p];
        UG_LOG("1 MovingParticle::output_velocity: ElemList.size(): " << ElemList.size() << "\n");
        if (ElemList.size() == 0) {
            UG_LOG("2 MovingParticle::output_velocity: ElemList.size(): " << ElemList.size() << " => skip assembling! \n");
            continue;
        }
#endif
        
    // Parameter von class 'CutElementHandler' extrahieren:
    
    // get multiindices for translation and rotation of particle
        std::vector < DoFIndex > transInd = m_spCutElementHandler->get_transInd(levIndex, p);
        std::vector < DoFIndex > rotInd = m_spCutElementHandler->get_rotInd(levIndex, p);
        
        MathVector<dim>  transSol;
        MathVector<dim>  rotSol;
        
        for ( int d = 0; d < dim; ++d )
        {
            transSol[d] = DoFRef(u, transInd[d]);
            rotSol[d]	= DoFRef(u, rotInd[d]);
        }
        
        // auf untergeordnete class 'ParticleProvider' weiterleiten:
        m_spCutElementHandler->print_velocity(transSol, rotSol, p, isTimedep, time, filename);
    }
    
}

} // end namespace NavierStokes
} // end namespace ug



#endif /* MOVING_PARTICLE_TOOLS_H_ */
