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
		const MathVector<dim>& evalPos,				// input data
		MathVector<dim+1>& interpolation)		// output data
{
	MathVector<dim> localPos(0.0);

//	get data
	typename TDomain::position_accessor_type aaPos = this->m_spApproxSpace->domain()->position_accessor();

	//	create Function Group
	FunctionGroup fctGrp(dd->function_pattern(), TokenizeString(m_spParticleHandlerGlobal->m_fctNames));

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

	//UG_LOG("end: interpolation = " << interpolation[dim] << "\n");
 }


template <typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::
print_pressure_nodal(const vector_type& u, const int topLevel)
{
	UG_LOG("START print_pressure_nodal()\n");

	FILE *pressure;
	char filename1[40];
	sprintf(filename1, "pressure_nodal_teta_level%d.txt", topLevel);
	pressure = fopen(filename1, "w");
	FILE *pressure_on_interface;
	char filename2[40];
	sprintf(filename2, "pressure_nodal_teta_level_ON_interface%d.txt", topLevel);
	pressure_on_interface = fopen(filename2, "w");
	int numIter = 0;
	number dist = 0.0;

	const int levIndex = get_Index(GridLevel(topLevel, GridLevel::LEVEL));
	std::vector<grid_base_object*> ElemListLog = m_spParticleHandlerGlobal->m_vvvElemListCut[levIndex][0];
	ConstSmartPtr<DoFDistribution> dd = this->m_spApproxSpace->dof_distribution(GridLevel(topLevel, GridLevel::LEVEL));
	typename TDomain::position_accessor_type aaPos = m_spParticleHandlerGlobal->m_aaPos;

 	const MathVector<dim>& center = m_spParticleHandlerGlobal->get_center(0);

 	UG_LOG("center " << center << "\n");

	typedef typename std::vector<grid_base_object*>::iterator ListIter;

	// loop all elements relevant for prtIndex-th particle
	for(ListIter listIter = ElemListLog.begin();
		listIter != ElemListLog.end(); ++listIter)
	{
		numIter++;
	//	get element
		grid_base_object* elem = *listIter;

	//	collect all vertices of the element
		std::vector<Vertex*> vVertex;
 		CollectVertices(vVertex, *m_spParticleHandlerGlobal->m_spMG, elem);

		MathVector<dim> bufferVector;
		VecSubtract(bufferVector, aaPos[vVertex[0]], aaPos[vVertex[2]]);
		dist = 0.5*VecLength(bufferVector);

	//	loop vertices
		for(size_t v = 0; v < vVertex.size(); ++v)
		{
		//	get vertex
			Vertex* vrt = vVertex[v];

 			if ( !m_spParticleHandlerGlobal->is_outsideFluid(vrt) )
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
 			if ( m_spParticleHandlerGlobal->is_outsideFluid(vrt) )
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

				fprintf(pressure_on_interface, "%e \t %e \t %lu # teta, value pressure, KnotenIndex vInd[0] \n", teta, DoFRef(u, vInd[0]), vInd[0][0]);

			}
		}

	}


	fclose(pressure);
	fclose(pressure_on_interface);

	UG_LOG("print_pressure_nodal() done...\n");
}

template <typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::
print_deltaP(const vector_type& u, const int topLevel)
{
	UG_LOG("START print_pressure()\n");

	FILE *pressure;
	char filename1[40];
	sprintf(filename1, "delta_pressure_level%d.txt", topLevel);
	pressure = fopen(filename1, "w");

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

  	UG_LOG("interpolVal1[" << dim<< "] = " << interpolVal1[dim] << "\n");
  	UG_LOG("interpolVal2[" << dim<< "] = " << interpolVal2[dim] << "\n");
  	UG_LOG("deltaP = " << deltaP << "\n");

   	fprintf(pressure, "%e \t %e \t %e \n# interpolVal1[dim], interpolVal1[dim], deltaP\n", interpolVal1[dim], interpolVal2[dim], deltaP);
	fclose(pressure);

}

template <typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::
print_pressure(const vector_type& u, const int topLevel)
{
	UG_LOG("START print_pressure()\n");

	FILE *pressure;
	char filename1[40];
	sprintf(filename1, "pressure_MovPrtSTdFV_teta_level%d.txt", topLevel);
	pressure = fopen(filename1, "w");

   	ConstSmartPtr<DoFDistribution> dd = this->m_spApproxSpace->dof_distribution(GridLevel(topLevel, GridLevel::LEVEL));
 	const MathVector<dim>& center = m_spParticleHandlerGlobal->get_center(0);

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

	//UG_THROW("Vorsicht bei Enwednugn von 'print_pressure()': früher wurde der Vektor 'interpolVal' an der Stelle dim+1 statt dim ausgelesen!!\n");
}

    /*
template <typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::
compute_error_on_circle(const vector_type& u, const int topLevel, number r)
{
 	FILE *solution;
	char filename1[40];
	sprintf(filename1, "error_MovPrt_radius%e_level%d.txt", r, topLevel);
	solution = fopen(filename1, "w");

 	const MathVector<dim>& center = m_spParticleHandlerGlobal->get_center(0);
 	const number R = m_spParticleHandlerGlobal->get_radius(0);
	ConstSmartPtr<DoFDistribution> dd = this->m_spApproxSpace->dof_distribution(GridLevel(topLevel, GridLevel::LEVEL));

	number error_u = 0.0;
	number error_v = 0.0;
	number error_p = 0.0;

	size_t t_max = 500;
 	MathVector<dim> printPos;
	MathVector<dim+1> interpolVal;
	for(size_t t = 0; t < t_max; ++t)
	{
		number teta = t*(2*PI/t_max);
		number co = cos(teta);
		number si = sin(teta);
		printPos[0] = r*co + center[0];
		printPos[1] = r*si + center[1];
		interpolate_point(dd, u, printPos, interpolVal);

		number sol_u = ((R*R-r*r)*co*co + r*r*log(r/R) + 0.5*(r*r-R*R))/(r*r);
		number sol_v = ((R*R-r*r)*si*co)/(r*r);
		number sol_p = -2*co*r;

		error_u += fabs(interpolVal[0]-sol_u);
		error_v += fabs(interpolVal[1]-sol_v);
		error_p += fabs(interpolVal[2]-sol_p);

		fprintf(solution, "%e \t", teta);
		for ( size_t i = 0; i <= dim; ++i)
			fprintf(solution, "%e \t", interpolVal[i]);

 		fprintf(solution, "%e \t %e # teta, interpolVal[0], interpolVal[1], interpolVal[2], printPos[0], printPot[1]\n", printPos[0], printPos[1]);

	}

	fclose(solution);

	UG_LOG("error_u: " << error_u/500 << "\n");
	UG_LOG("error_v: " << error_v/500 << "\n");
	UG_LOG("error_p: " << error_p/500 << "\n\n");

}
*/
template <typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::
print_velocity(const vector_type& u, const int topLevel, number time, const char* filename)
{
	UG_LOG("MovingParticle::print_velocity(): Start: \n");

   // Parameter von class 'MovingParticle' extrahieren:
    const int levIndex = get_Index(GridLevel(topLevel, GridLevel::LEVEL));
    const bool isTimedep = is_time_dependent();
    
    size_t numPrt = m_spParticleHandlerGlobal->num_particles();
    
//    if ( numPrt > 2 )
//        UG_THROW("Particle:output_velocity: VORSICHT, output nicht implementiert fuer mehr als 2 particle! -> m_bShared[][] ist Problem! ... Exit! \n");
    
    for(size_t p = 0; p < numPrt; ++p)
    {
        
#ifdef UG_PARALLEL
        std::vector<grid_base_object*> ElemList = m_spParticleHandlerGlobal->m_vvvElemListCut[levIndex][p];
        UG_LOG("1 MovingParticle::output_velocity: ElemList.size(): " << ElemList.size() << "\n");
        if (ElemList.size() == 0) {
            UG_LOG("2 MovingParticle::output_velocity: ElemList.size(): " << ElemList.size() << " => skip assembling! \n");
            continue;
        }
#endif
        
        // Parameter von class 'CutElementHandler' extrahieren:
        
        // get multiindices for translation and rotation of particle
        std::vector < DoFIndex > transInd = m_spParticleHandlerGlobal->get_transInd(levIndex, p);
        std::vector < DoFIndex > rotInd = m_spParticleHandlerGlobal->get_rotInd(levIndex, p);
        
        MathVector<dim>  transSol;
        MathVector<dim>  rotSol;
        
        for ( int d = 0; d < dim; ++d )
        {
            transSol[d] = DoFRef(u, transInd[d]);
            rotSol[d]	= DoFRef(u, rotInd[d]);
        }
        
        // auf untergeordnete class 'ParticleProvider' weiterleiten:
        m_spParticleHandlerGlobal->print_velocity(transSol, rotSol, p, isTimedep, time, filename);
    }
    
	UG_LOG("output_velocity....DONE\n");
        
}
    /*

template <typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::
print_velocity_many_particles(const vector_type& u, const int topLevel, number time, const char* filename)
{
	UG_LOG("MovingParticle::print_velocity(): Start: \n");

	size_t numPrt = m_spParticleHandlerGlobal->num_particles();

	const int levIndex = get_Index(GridLevel(topLevel, GridLevel::LEVEL));
	number distance = -10.0;

	// compute distance between 2 particles
	if ( numPrt == 2 )
	{
		const number radius1 = m_spParticleHandlerGlobal->get_radius(0);
		const number radius2 = m_spParticleHandlerGlobal->get_radius(1);

		const number radius = radius1+radius2;

		MathVector<dim>  distVec;
 		VecSubtract(distVec, m_spParticleHandlerGlobal->get_center(1), m_spParticleHandlerGlobal->get_center(0));

		distance = VecDot(distVec, distVec);
		distance = sqrt(distance);
		distance = distance - radius;
	}
	for(size_t p = 0; p < numPrt; ++p)
	{

#ifdef UG_PARALLEL
 		std::vector<grid_base_object*> ElemList = m_spParticleHandlerGlobal->m_vvvElemListCut[levIndex][p];
 		UG_LOG("1 MovingParticle::output_velocity: ElemList.size(): " << ElemList.size() << "\n");
		if (ElemList.size() == 0) {
 			UG_LOG("2 MovingParticle::output_velocity: ElemList.size(): " << ElemList.size() << " => skip assembling! \n");
			continue;
		}
#endif
	// get multiindices for translation and rotation of particle
 		std::vector < DoFIndex > transInd = m_spParticleHandlerGlobal->get_transInd(levIndex, p);
		std::vector < DoFIndex > rotInd = m_spParticleHandlerGlobal->get_rotInd(levIndex, p);

		MathVector<dim>  transSol;
		MathVector<dim>  rotSol;

		MathVector<dim>  center = m_spParticleHandlerGlobal->get_center(p);
		const number radius = m_spParticleHandlerGlobal->get_radius(p);

		for ( int d = 0; d < dim; ++d )
		{
			transSol[d] = DoFRef(u, transInd[d]);
			rotSol[d]	= DoFRef(u, rotInd[d]);
		}

		number f1 = 0.0;
		number gravity = -9.81;
		if ( !is_time_dependent() )
			f1 = transSol[0]*4.0/(radius*radius*gravity);

		std::string name(filename);

		char * cstr = new char [name.size()+1];
		strcpy (cstr, name.c_str());


		if ( !is_time_dependent() )
		{
			FILE* print_velocity = fopen(name.c_str(), "a");


			fprintf(print_velocity,"%e \t %e \t ",radius, f1);
			for ( int d = 0; d < dim; ++d )
				fprintf(print_velocity,"%e \t %e \t ", transSol[d], rotSol[d]);
			for ( int d = 0; d < dim; ++d )
				fprintf(print_velocity,"%e \t ", center[d]);

			fprintf(print_velocity," # radius, f1, transVel[0], rotVel[0], transVel[1], rotVel[1], center_coords, (m_bSharedIP = ");
			fprintf(print_velocity, "0 ), ");

			fprintf(print_velocity," (m_bSharedElem = ");
			fprintf(print_velocity, "0 )  \n ");
	fclose(print_velocity);

		}
		else
		{
			FILE* print_velocity = fopen(name.c_str(), "a");

 			fprintf(print_velocity,"%e \t ",time);

			for ( int d = 0; d < dim; ++d )
				fprintf(print_velocity,"%e \t %e \t ", transSol[d], rotSol[d]);

			for ( int d = 0; d < dim; ++d )
				fprintf(print_velocity,"%e \t ", center[d]);


			fprintf(print_velocity,"%e \t ", distance);

			fprintf(print_velocity," # time, transVel[0], rotVel[0], transVel[1], rotVel[1], center_coords, (m_bSharedIP = ");
            fprintf(print_velocity, "0  ), ");

			fprintf(print_velocity," (m_bSharedElem = ");
		fprintf(print_velocity, "0 )  \n ");
		fclose(print_velocity);

		}

	}

	UG_LOG("output_velocity....DONE\n");

}
*/
} // end namespace NavierStokes
} // end namespace ug



#endif /* MOVING_PARTICLE_TOOLS_H_ */
