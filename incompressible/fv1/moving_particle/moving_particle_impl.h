/*
 * moving_particle_impl.h
 *
 *  Created on: 20.01.2015
 *      Author: suze
 */

#ifndef MOVING_PARTICLE_IMPL_H_
#define MOVING_PARTICLE_IMPL_H_

#include <map>

namespace ug {
namespace NavierStokes {


///////////////////////////////////////////////////////////
// Implementation of the methods class
// 	 		'MovingParticle'
///////////////////////////////////////////////////////////

template<typename TDomain, typename TAlgebra>
MovingParticle<TDomain, TAlgebra>::MovingParticle(
		SmartPtr<IAssemble<TAlgebra> > ass,
		SmartPtr<NavierStokesFV1<TDomain> > spMaster,
		SmartPtr<CutElementHandlerFlatTop<dim> > cutElementHandler,
		number fluidDensity, number fluidKinVisc) :
		m_spParticleHandlerGlobal(cutElementHandler),
		m_spInterfaceHandlerLocal(
				new InterfaceHandlerLocalParticle<dim>(cutElementHandler, fluidDensity, fluidKinVisc) ),
		m_spInterfaceMapper(
				new ParticleMapper<TDomain, TAlgebra>(m_spInterfaceHandlerLocal)),
		m_spInterfaceBndCond(
				new ParticleBndCond<TDomain>(spMaster, m_spInterfaceHandlerLocal))
{
    
	if (cutElementHandler->num_particles() == 0)
		UG_THROW("MovingParticle::Constructor(): no particles initializen in 'cutElementHandler\n");

	// initialize singleton and set local handler
	typedef DimFV1FTGeometry<dim, dim, InterfaceHandlerLocalParticle<dim> > TFVGeom;
	TFVGeom& geo = GeomProvider<TFVGeom>::get(LFEID(LFEID::LAGRANGE, dim, 1), 1);
	geo.set_interface_handler(m_spInterfaceHandlerLocal);

// initialize mapper within domainDisc:
	SmartPtr<AssemblingTuner<TAlgebra> > assAdapt = ass->ass_tuner();
	assAdapt->set_mapping(m_spInterfaceMapper.get());
	assAdapt->enable_modify_solution(true);
    
    
	// => assTuner->modify_LocSol() = mapper->modify_LocSol()
	// see: ass_tuner.h: 114

	//  ToDo: reconstruct modify_LocSol: ???
	//	assAdapt->set_local_modifier();

}

/*
 template <typename TDomain, typename TAlgebra>
 void MovingParticle<TDomain, TAlgebra>::
 clear_solution(SmartPtr<VectorTimeSeries<vector_type> > vSol, const GridLevel& topGrid)
 {
 //approxSpace->dof_distribution(GridLevel(GridLevel::TOP, GridLevel::SURFACE)), time);

 }
 */

template<typename TDomain, typename TAlgebra>
number MovingParticle<TDomain, TAlgebra>::MeanElementDiameter(TDomain& domain,
		int level) {
	typedef typename domain_traits<TDomain::dim>::grid_base_object TElem;

	//typedef typename std::vector<grid_base_object*>::iterator ListIter;
	typedef typename geometry_traits<TElem>::iterator ListIter;

	ListIter iter = domain.grid()->template begin<TElem>(level);
	ListIter iterEnd = domain.grid()->template end<TElem>(level);

	number mean = 0.0;
	size_t numIter = 0;
	for (; iter != iterEnd; ++iter) {
		mean += ElementDiameterSq(*domain.grid(), domain.position_accessor(),
				*iter);
		numIter++;
	}

	mean = mean / numIter;
	UG_LOG("mean = " << std::sqrt(mean) << "\n");

#ifdef UG_PARALLEL
	// share value between all procs
	pcl::ProcessCommunicator com;
	// ToDO: PCL_RO_MIN oder MAX oder doch mean noch berechnen m.H.v. /numProcs???
	//UG_THROW("in MeanElementDiameter(): was macht 'allredue' im Parallelen???\n");
	mean = com.allreduce(mean, PCL_RO_MIN);
#endif

	UG_LOG("nach com.allreduce: mean = " << std::sqrt(mean) << "\n");
	return std::sqrt(mean);
}

template<typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::initialize_threshold(TDomain& domain,
		const int baseLevel, const int topLevel) {
	UG_LOG("----------------- START initialize_threshold() ---------------- \n");

	if (baseLevel < 0)
		UG_THROW(
				"initialize_threshold(): no cast of baselevel from 'int' tp 'size_t' possible! \n");
	if (topLevel < 0)
		UG_THROW(
				"initialize_threshold(): no cast of toplevel from 'int' tp 'size_t' possible! \n");

	typedef typename domain_traits<TDomain::dim>::grid_base_object TElem;

// compute level-dependent value for threshold:
	for (size_t lev = baseLevel; lev <= (size_t) topLevel; ++lev) {
		const number maxLength = MaxElementDiameter(domain, lev);
		const number meanLength = MeanElementDiameter(domain, lev);
		UG_LOG("maxLength = " << maxLength << "\n");
		UG_LOG("meanLength = " << meanLength << "\n");
		UG_LOG("threshold_max = " << maxLength*maxLength << "\n");
		UG_LOG("threshold_mean = " << meanLength*meanLength << "\n");

		m_spParticleHandlerGlobal->set_threshold(lev, meanLength * meanLength);
		//m_spParticleHandlerGlobal->set_threshold(lev, 0.003);

	}

	UG_LOG("----------------- END initialize_threshold() ---------------- \n");

}

template <typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::
get_velocity(MathVector<dim>& transSol, MathVector<dim>& rotSol, const vector_type& u, const int topLevel, number deltaT, const size_t prtIndex)
{
    UG_LOG("MovingParticle::print_velocity(): Start: \n");
        
    size_t numPrt = m_spParticleHandlerGlobal->num_particles();
    
  //  if ( numPrt > 2 )
  //      UG_THROW("Particle:output_velocity: VORSICHT, output nicht implementiert fuer mehr als 2 particle! -> m_bShared[][] ist Problem! ... Exit! \n");
        
    const int levIndex = get_Index(GridLevel(topLevel, GridLevel::LEVEL));

            
#ifdef UG_PARALLEL
    std::vector<grid_base_object*> ElemList = m_spParticleHandlerGlobal->m_vvvElemListCut[levIndex][prtIndex];
    UG_LOG("1 MovingParticle::output_velocity: ElemList.size(): " << ElemList.size() << "\n");
    if (ElemList.size() == 0) {
        UG_LOG("2 MovingParticle::output_velocity: ElemList.size(): " << ElemList.size() << " => skip assembling! \n");
        return;
    }
#endif
    // get multiindices for translation and rotation of particle
    std::vector < DoFIndex > transInd = m_spParticleHandlerGlobal->get_transInd(levIndex, prtIndex);
    std::vector < DoFIndex > rotInd = m_spParticleHandlerGlobal->get_rotInd(levIndex, prtIndex);
    
                
    for ( int d = 0; d < dim; ++d )
    {
        transSol[d] = DoFRef(u, transInd[d]);
        rotSol[d]	= DoFRef(u, rotInd[d]);
        UG_LOG("in get_velocity: transSol: " << transSol[d] << "\t rotSol: " << rotSol[d] << "\n");
    }
            

    
}
    

// transfer direction: DoFRef(u, transInd) ---> m_vvLinearVelocity/m_vvAngularVelocity
template<typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::copy_particle_solution(vector_type& u,
		const int topLevel) {
	const int levIndex = get_Index(GridLevel(topLevel, GridLevel::LEVEL));

	for (size_t p = 0; p < m_spParticleHandlerGlobal->num_particles(); ++p) {

#ifdef UG_PARALLEL
		std::vector<grid_base_object*> ElemList =
				m_spParticleHandlerGlobal->m_vvvElemListCut[levIndex][p];
		UG_LOG(
				"1 MovingParticle::copy_solution: ElemList.size(): " << ElemList.size() << "\n");
		if (ElemList.size() == 0) {
			UG_LOG(
					"2 MovingParticle::copy_solution: ElemList.size(): " << ElemList.size() << " => skip assembling! \n");
			continue;
		}
#endif
#ifdef UG_DEBUG
		UG_LOG(
				"copy_solution(): VORHER transSol(0) = " << m_spParticleHandlerGlobal->get_transSol(p,0) << "\n");
		UG_LOG(
				"copy_solution(): VORHER transSol(1) = " << m_spParticleHandlerGlobal->get_transSol(p, 1) << "\n");
		UG_LOG(
				"copy_solution(): VORHER rotSol(0) = " << m_spParticleHandlerGlobal->get_rotSol(p,0) << "\n");
		UG_LOG(
				"copy_solution(): VORHER rot1Sol(1) = " << m_spParticleHandlerGlobal->get_rotSol(p, 1) << "\n");
#endif
		// 'get_transInd()' returns NOT YET updated indices (updated during 'movingParticle:update()' vie lua)
		std::vector < DoFIndex > transInd =
				m_spParticleHandlerGlobal->get_transInd(levIndex, p);
		std::vector < DoFIndex > rotInd = m_spParticleHandlerGlobal->get_rotInd(
				levIndex, p);

		UG_LOG(" transInd(0) = " << transInd[0] << "\n");
		UG_LOG(" rotInd(0) = " << rotInd[0] << "\n");

		for (int d = 0; d < dim; ++d) {
			number solution = DoFRef(u, transInd[d]);
			m_spParticleHandlerGlobal->set_extraSolTrans(solution, p, 0, d);
			DoFRef(u, transInd[d]) = 0.0;

			solution = DoFRef(u, rotInd[d]);
			m_spParticleHandlerGlobal->set_extraSolRot(solution, p, 0, d);
			DoFRef(u, rotInd[d]) = 0.0;
		}
#ifdef UG_DEBUG
		UG_LOG(
				"copy_solution(): NACHHER transSol(0) = " << m_spParticleHandlerGlobal->get_transSol(p,0) << "\n");
		UG_LOG(
				"copy_solution(): NACHHER transSol(1) = " << m_spParticleHandlerGlobal->get_transSol(p, 1) << "\n");
		UG_LOG(
				"copy_solution(): NACHHER rotSol(0) = " << m_spParticleHandlerGlobal->get_rotSol(p,0) << "\n");
		UG_LOG(
				"copy_solution(): NACHHER rot1Sol(1) = " << m_spParticleHandlerGlobal->get_rotSol(p, 1) << "\n");
#endif
	} // end particle loop

}

template<typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::compute_gradient_local_max(
		vector_type& sol, vector_type& grad,
		SmartPtr<ApproximationSpace<TDomain> > spApproxSpace,
		const int topLevel) {

// get data
	typedef typename domain_traits<dim>::grid_base_object grid_base_object;
	typename DoFDistribution::traits<grid_base_object>::const_iterator iter,
			iterEnd;

	ConstSmartPtr<DoFDistribution> dd = spApproxSpace->dof_distribution(
			GridLevel(topLevel, GridLevel::LEVEL));
	typename TDomain::position_accessor_type aaPos =
			m_spParticleHandlerGlobal->m_aaPos;

	iter = dd->template begin<grid_base_object>();
	iterEnd = dd->template end<grid_base_object>();

	//	loop elements in order to compute the volume and set rhs:
	for (; iter != iterEnd; iter++) {
		//	get element
		grid_base_object* elem = *iter;

		if (m_spParticleHandlerGlobal->m_spOutsideMarker->is_marked(elem))
			continue;

		std::vector < DoFIndex > vInd1;
		std::vector < DoFIndex > vInd2;

		std::vector<Edge*> vEdges;
		CollectEdgesSorted(vEdges, *m_spParticleHandlerGlobal->m_spMG, elem);

		if (m_spParticleHandlerGlobal->m_spCutMarker->is_marked(elem)) {
			for (size_t e = 0; e < vEdges.size(); ++e) {
				Edge* edge = vEdges[e];
				std::vector<Vertex*> vVertexEdge;
				CollectVertices(vVertexEdge, *m_spParticleHandlerGlobal->m_spMG,
						edge);
				if (vVertexEdge.size() != 2)
					UG_THROW(
							"error in collecting vertices associated to an edge!....EXIT!...\n");

				Vertex* vrt1 = vVertexEdge[0];
				Vertex* vrt2 = vVertexEdge[1];

				// no gradient needs to be computed on outside edges:
				if (m_spInterfaceHandlerLocal->is_FTVertex(vrt1)
						&& m_spInterfaceHandlerLocal->is_FTVertex(vrt2))
					continue;

				if (dd->inner_dof_indices(vrt1, dim, vInd1) != 1)
					UG_THROW(
							"MovingParticle::compute_gradient(): Only one index expected.");
				if (dd->inner_dof_indices(vrt2, dim, vInd2) != 1)
					UG_THROW(
							"MovingParticle::compute_gradient(): Only one index expected.");

				// compute local gradient along the edge
				//	UG_LOG("DoFRef(sol,vInd1[0]): " << DoFRef(sol,vInd1[0]) << "\n");
				//	UG_LOG("DoFRef(sol,vInd2[0]): " << DoFRef(sol,vInd2[0]) << "\n");

				number grad_edge = fabs(
						DoFRef(sol, vInd1[0]) - DoFRef(sol, vInd2[0]));
				number edgeLength = VecDistance(aaPos[vrt1], aaPos[vrt2]);

				if (!m_spInterfaceHandlerLocal->is_FTVertex(vrt1)
						&& !m_spInterfaceHandlerLocal->is_FTVertex(vrt2)) { // nothing has to be adapted!
				} else if (m_spInterfaceHandlerLocal->is_FTVertex(vrt1)) {
					if (m_spInterfaceHandlerLocal->is_FTVertex(vrt2))
						UG_THROW("hmmm...wrong in 'compute_gradient()'...\n");

					// compute intersectionPoint:
					MathVector<dim> intersectionPnt;
					if (m_spInterfaceHandlerLocal->is_nearInterfaceVertex(vrt1))
						VecCopy(intersectionPnt, aaPos[vrt1], 0.0);
					else
						m_spInterfaceHandlerLocal->get_intersection_point(
								intersectionPnt, vrt2, vrt1);

					// compute edgeLength:
					edgeLength = VecDistance(aaPos[vrt2], intersectionPnt);
				} else if (m_spInterfaceHandlerLocal->is_FTVertex(vrt2)) {
					if (m_spInterfaceHandlerLocal->is_FTVertex(vrt1))
						UG_THROW("hmmm...wrong in 'compute_gradient()'...\n");

					// compute intersectionPoint:
					MathVector<dim> intersectionPnt;
					if (m_spInterfaceHandlerLocal->is_nearInterfaceVertex(vrt2))
						VecCopy(intersectionPnt, aaPos[vrt2], 0.0);
					else
						m_spInterfaceHandlerLocal->get_intersection_point(
								intersectionPnt, vrt1, vrt2);
					// compute edgeLength:
					edgeLength = VecDistance(aaPos[vrt1], intersectionPnt);
				}

				grad_edge = grad_edge / edgeLength;
				// set local maximum over all gradients along the edge to a given vertex
				DoFRef(grad, vInd1[0]) = std::max(DoFRef(grad, vInd1[0]),
						grad_edge);
				DoFRef(grad, vInd2[0]) = std::max(DoFRef(grad, vInd2[0]),
						grad_edge);

			} // end edge-loop
		} else {
			for (size_t e = 0; e < vEdges.size(); ++e) {
				Edge* edge = vEdges[e];
				std::vector<Vertex*> vVertexEdge;
				CollectVertices(vVertexEdge, *m_spParticleHandlerGlobal->m_spMG,
						edge);
				if (vVertexEdge.size() != 2)
					UG_THROW(
							"error in collecting vertices associated to an edge!....EXIT!...\n");

				Vertex* vrt1 = vVertexEdge[0];
				Vertex* vrt2 = vVertexEdge[1];
				if (dd->inner_dof_indices(vrt1, dim, vInd1) != 1)
					UG_THROW(
							"MovingParticle::compute_gradient(): Only one index expected.");
				if (dd->inner_dof_indices(vrt2, dim, vInd2) != 1)
					UG_THROW(
							"MovingParticle::compute_gradient(): Only one index expected.");

				// compute local gradient along the edge
				number edgeLength = VecDistance(aaPos[vrt1], aaPos[vrt2]);
				number grad_edge = fabs(
						DoFRef(sol, vInd1[0]) - DoFRef(sol, vInd2[0]));
				grad_edge = grad_edge / edgeLength;

				// set local maximum over all gradients along the edge to a given vertex
				DoFRef(grad, vInd1[0]) = std::max(DoFRef(grad, vInd1[0]),
						grad_edge);
				DoFRef(grad, vInd2[0]) = std::max(DoFRef(grad, vInd2[0]),
						grad_edge);

			} // end edge-loop

		}

	} // end elem-loop

}

template<typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::adjust_global_solution(vector_type& u,
		const int topLevel) {
	const int levIndex = get_Index(GridLevel(topLevel, GridLevel::LEVEL));

	for (size_t p = 0; p < m_spParticleHandlerGlobal->num_particles(); ++p) {

#ifdef UG_PARALLEL
		std::vector<grid_base_object*> ElemList =
				m_spParticleHandlerGlobal->m_vvvElemListCut[levIndex][p];
		UG_LOG(
				"1 MovingParticle::copy_and_fill_solution: ElemList.size(): " << ElemList.size() << "\n");
		if (ElemList.size() == 0) {
			UG_LOG(
					"2 MovingParticle::copy_and_fill_solution: ElemList.size(): " << ElemList.size() << " => skip assembling! \n");
			continue;
		}
#endif

		// 'get_transInd()' returns NOT YET updated indices (updated during 'movingParticle:update()' vie lua)
		std::vector < DoFIndex > transInd =
				m_spParticleHandlerGlobal->get_transInd(levIndex, p);
		std::vector < DoFIndex > rotInd = m_spParticleHandlerGlobal->get_rotInd(
				levIndex, p);

		MathVector<dim> transSol = m_spParticleHandlerGlobal->get_transSol(p,
				0);
		MathVector<dim> rotSol = m_spParticleHandlerGlobal->get_rotSol(p, 0);

		UG_LOG(
				"in adjust_global_solution(): transSol(0) = " << m_spParticleHandlerGlobal->get_transSol(p,0) << "\n");
		UG_LOG(
				"in adjust_global_solution(): transSol(1) = " << m_spParticleHandlerGlobal->get_transSol(p, 1) << "\n");
		UG_LOG(
				"in adjust_global_solution(): rotSol(0) = " << m_spParticleHandlerGlobal->get_rotSol(p,0) << "\n");
		UG_LOG(
				"in adjust_global_solution	(): rot1Sol(1) = " << m_spParticleHandlerGlobal->get_rotSol(p, 1) << "\n");

		ConstSmartPtr<DoFDistribution> dd =
				this->m_spApproxSpace->dof_distribution(
						GridLevel(topLevel, GridLevel::LEVEL));
		typename TDomain::position_accessor_type aaPos =
				m_spParticleHandlerGlobal->m_aaPos;

		const MathVector<dim>& center = m_spParticleHandlerGlobal->get_center(
				p);

		typedef typename std::vector<grid_base_object*>::iterator ListIter;

		// loop all elements relevant for prtIndex-th particle
		std::vector<grid_base_object*> ElemListLog =
				m_spParticleHandlerGlobal->m_vvvElemListOutside[levIndex][p];

		for (ListIter listIter = ElemListLog.begin();
				listIter != ElemListLog.end(); ++listIter) {
			//	get element
			grid_base_object* elem = *listIter;

			//	collect all vertices of the element
			std::vector<Vertex*> vVertex;
			CollectVertices(vVertex, *m_spParticleHandlerGlobal->m_spMG, elem);

			//	loop vertices
			for (size_t v = 0; v < vVertex.size(); ++v) {
				//	get vertex
				Vertex* vrt = vVertex[v];

				for (size_t fct = 0; fct < dim; ++fct) {
					//	create multi index
					std::vector < DoFIndex > vInd;
					//	get multi indices
					if (dd->inner_dof_indices(vrt, fct, vInd) != 1)
						UG_THROW("Only one index expected.");

					// set solution: particle velocity
					MathVector<dim> radialCo;
					VecSubtract(radialCo, aaPos[vrt], center);
					MathMatrix<dim, dim> rotationMatCo =
							m_spParticleHandlerGlobal->get_rotationMat(
									radialCo);

					//	set solution
					DoFRef(u, vInd[0]) = transSol[fct];
					for (int d = 0; d < dim; ++d)
						DoFRef(u, vInd[0]) += rotationMatCo[fct][d] * rotSol[d];

					DoFRef(u, vInd[0]) = 0.0;

				}

			} // end vrt-loop
		} // end outsideElem-loop

	} // end prt-loop

}

template<typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::fill_particle_solution(vector_type& u,
		const int topLevel, const number time) {

	bool output = false;

	const char* filename = "freed_solution";
	std::string name(filename);
	char ext[50];
	sprintf(ext, ".txt");
	name.append(ext);
	FILE* outputFile = fopen(name.c_str(), "a");
	if (output) {
		fprintf(outputFile,
				"------------------------------------------------------------- \n");
		fprintf(outputFile, "-- 'fill_particle_solution()' for time = %e -- \n",
				time);
		fprintf(outputFile,
				"------------------------------------------------------------- \n");
	}

	const int levIndex = get_Index(GridLevel(topLevel, GridLevel::LEVEL));

	for (size_t p = 0; p < m_spParticleHandlerGlobal->num_particles(); ++p) {

#ifdef UG_PARALLEL
		std::vector<grid_base_object*> ElemList =
				m_spParticleHandlerGlobal->m_vvvElemListCut[levIndex][p];
		UG_LOG(
				"1 MovingParticle::copy_and_fill_solution: ElemList.size(): " << ElemList.size() << "\n");
		if (ElemList.size() == 0) {
			UG_LOG(
					"2 MovingParticle::copy_and_fill_solution: ElemList.size(): " << ElemList.size() << " => skip assembling! \n");
			continue;
		}
#endif

		// 'get_transInd()' returns NOT YET updated indices (updated during 'movingParticle:update()' vie lua)
		std::vector < DoFIndex > transInd =
				m_spParticleHandlerGlobal->get_transInd(levIndex, p);
		std::vector < DoFIndex > rotInd = m_spParticleHandlerGlobal->get_rotInd(
				levIndex, p);

		MathVector<dim> transSol = m_spParticleHandlerGlobal->get_transSol(p,0);
		MathVector<dim> rotSol = m_spParticleHandlerGlobal->get_rotSol(p, 0);

		UG_LOG(
				"in fill_particle_solution(): transSol(0) = " << m_spParticleHandlerGlobal->get_transSol(p,0) << "\n");
		UG_LOG(
				"in fill_particle_solution(): transSol(1) = " << m_spParticleHandlerGlobal->get_transSol(p, 1) << "\n");
		UG_LOG(
				"in fill_particle_solution(): rotSol(0) = " << m_spParticleHandlerGlobal->get_rotSol(p,0) << "\n");
		UG_LOG(
				"in fill_particle_solution(): rot1Sol(1) = " << m_spParticleHandlerGlobal->get_rotSol(p, 1) << "\n");

		ConstSmartPtr<DoFDistribution> dd =
				this->m_spApproxSpace->dof_distribution(
						GridLevel(topLevel, GridLevel::LEVEL));
		typename TDomain::position_accessor_type aaPos =
				m_spParticleHandlerGlobal->m_aaPos;

		const MathVector<dim>& center = m_spParticleHandlerGlobal->get_center(
				p);

		typedef typename std::vector<grid_base_object*>::iterator ListIter;

		// loop all elements relevant for prtIndex-th particle
		std::vector<grid_base_object*> ElemListLog =
				m_spParticleHandlerGlobal->m_vvvElemListOutside[levIndex][p];

		/*		for(ListIter listIter = ElemListLog.begin();
		 listIter != ElemListLog.end(); ++listIter)
		 {
		 //	get element
		 grid_base_object* elem = *listIter;

		 //	collect all vertices of the element
		 std::vector<Vertex*> vVertex;
		 CollectVertices(vVertex, *m_spParticleHandlerGlobal->m_spMG, elem);

		 //	loop vertices
		 for(size_t v = 0; v < vVertex.size(); ++v)
		 {
		 //	get vertex
		 Vertex* vrt = vVertex[v];


		 for (size_t fct = 0; fct < dim; ++fct)
		 {
		 //	create multi index
		 std::vector<DoFIndex>  vInd;
		 //	get multi indices
		 if(dd->inner_dof_indices(vrt, fct, vInd) != 1)
		 UG_THROW("Only one index expected.");

		 // set solution: particle velocity
		 MathVector<dim> radialCo;
		 VecSubtract(radialCo, aaPos[vrt], center);
		 MathMatrix<dim,dim> rotationMatCo = m_spParticleHandlerGlobal->get_rotationMat(radialCo);

		 //	set solution
		 DoFRef(u,vInd[0]) = transSol[fct];
		 for ( int d = 0; d < dim; ++d )
		 DoFRef(u,vInd[0]) += rotationMatCo[fct][d]*rotSol[d];

		 }


		 } // end vrt-loop
		 } // end outsideElem-loop
		 */
		// loop all elements relevant for prtIndex-th particle
		ElemListLog = m_spParticleHandlerGlobal->m_vvvElemListCut[levIndex][p];

		for (ListIter listIter = ElemListLog.begin();
				listIter != ElemListLog.end(); ++listIter) {
			//	get element
			grid_base_object* elem = *listIter;

			//	collect all vertices of the element
			std::vector<Vertex*> vVertex;
			CollectVertices(vVertex, *m_spParticleHandlerGlobal->m_spMG, elem);

			//	loop vertices
			for (size_t v = 0; v < vVertex.size(); ++v) {
				//	get vertex
				Vertex* vrt = vVertex[v];

				// ---> is_outsideFluid() is a geometrical check w.r.t new particle coords
				// ---> m_spOutsideMarker was marked w.r.t old particle coords
				//		===> !is_inside && is_outside = is_freed :-)
				if (m_spParticleHandlerGlobal->is_outsideFluid_prtIndex(p, vrt)
						&& m_spParticleHandlerGlobal->m_spOutsideMarker->is_marked(
								vrt)) {
					if (1) {
						fprintf(outputFile, " freed node: %e \t %e \t",
								aaPos[vrt][0], aaPos[vrt][1]);
						if (dim == 3)
							fprintf(outputFile, " node: %e", aaPos[vrt][2]);
						fprintf(outputFile, " time: %e \n\n", time);
					}

					for (size_t fct = 0; fct < dim; ++fct) {
						//	create multi index
						std::vector < DoFIndex > vInd;
						//	get multi indices
						if (dd->inner_dof_indices(vrt, fct, vInd) != 1)
							UG_THROW("Only one index expected.");

						// set solution: particle velocity
						MathVector<dim> radialCo;
						VecSubtract(radialCo, aaPos[vrt], center);
						MathMatrix<dim, dim> rotationMatCo =
								m_spParticleHandlerGlobal->get_rotationMat(
										radialCo);

						if (output)
                        {
							fprintf(outputFile, "fct = %lu:\n", fct);
							fprintf(outputFile," transSol: %e, \t rotSol: %e\n", transSol[fct], rotSol[fct]);
							fprintf(outputFile, " vorher: %e, \t", DoFRef(u, vInd[0]));
						}
						//	set solution
						DoFRef(u, vInd[0]) = transSol[fct];

                        if ( output)
                        {
                            UG_LOG(" fill_sol_: transInd(0) = " << transInd[0] << "\n");
                            UG_LOG(" fill_sol_: rotInd(0) = " << rotInd[0] << "\n");
                            UG_LOG(" fill_sol_: transSol[" << fct << "]: " << transSol[fct] << "\n");
                            UG_LOG(" fill_sol_: rotSol[" << fct << "]: " << rotSol[fct] << "\n");
                        }
                        
						for (int d = 0; d < dim; ++d)
							DoFRef(u, vInd[0]) += rotationMatCo[fct][d]
									* rotSol[d];

						if (output)
							fprintf(outputFile, " nachher: %e \n\n",
									DoFRef(u, vInd[0]));

					}
					// wenn node ganz im Inneren vom Partikel liegt (!FlatTopVrtMarker), dann hat sich das Partikel zu schnell
					// weiterbewegt -> CFL-Schranke verletzt!
					if (!m_spParticleHandlerGlobal->m_spFlatTopVrtMarker->is_marked(
							vrt)
							&& m_spParticleHandlerGlobal->m_spOutsideMarker->is_marked(
									vrt))
						UG_THROW(
								"in 'fill_particle_solution()': CFL-Schranke verletzt!\n");
				}

			} // end vrt-loop
		} // end cutElem-loop

	} // end prt-loop

	fclose(outputFile);

}

template<typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::clear_solution(vector_type& u,
		const int topLevel) {
//	const int levIndex = get_Index(GridLevel(topLevel, GridLevel::LEVEL));

}

// call via lua BEFORE applying 'movingParticle:update()', which changes the global indices of trandSol/rotSol
// necessary for re-writing solution to NEW indices of vector 'u' during 'update_solution'
// 			=> handling instance/buffering data := m_spParticleHandlerGlobal->(m_vSolTrans/m_vSolRot)

// transfer direction:  m_vSolTrans ---> DoFRef(u, transInd)
template<typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::update_particle_solution(vector_type& u, const int topLevel)
{
	const int levIndex = get_Index(GridLevel(topLevel, GridLevel::LEVEL));

	for (size_t p = 0; p < m_spParticleHandlerGlobal->num_particles(); ++p) {

#ifdef UG_PARALLEL
		std::vector<grid_base_object*> ElemList =
				m_spParticleHandlerGlobal->m_vvvElemListCut[levIndex][p];
		UG_LOG(
				"1 MovingParticle::update_solution:  ElemList.size(): " << ElemList.size() << "\n");
		if (ElemList.size() == 0) {
			UG_LOG(
					"2 MovingParticle::update_solution: ElemList.size(): " << ElemList.size() << " => skip assembling! \n");
			continue;
		}
#endif

		// 'get_transInd()' returns allready updated indices (updated during 'movingParticle:update()' vie lua)
		std::vector < DoFIndex > transInd = m_spParticleHandlerGlobal->get_transInd(levIndex, p);
		std::vector < DoFIndex > rotInd   = m_spParticleHandlerGlobal->get_rotInd(levIndex, p);

		UG_LOG("in 'update_particle_solution()': transInd = " << transInd[0] << "\n");
		UG_LOG("in 'update_particle_solution()': rotInd = " << rotInd[0] << "\n");

		for (int d = 0; d < dim; ++d) {
#ifdef UG_DEBUG
			UG_LOG(
					"VORHER: in 'update_solution()': DoFRef(u, transInd[d]) = " << DoFRef(u, transInd[d]) << "\n");
			UG_LOG(
					"VORHER: in 'update_solution()': DoFRef(u, rotInd[d]) = " << DoFRef(u, rotInd[d]) << "\n");
#endif
			DoFRef(u, transInd[d]) = m_spParticleHandlerGlobal->get_transSol(p,0)[d];
			DoFRef(u, rotInd[d])   = m_spParticleHandlerGlobal->get_rotSol(p, 0)[d];
#ifdef UG_DEBUG
			UG_LOG(
					"NACHHER: in 'update_solution()': DoFRef(u, transInd[d]) = " << DoFRef(u, transInd[d]) << "\n");
			UG_LOG(
					"NACHHER: in 'update_solution()': DoFRef(u, rotInd[d]) = " << DoFRef(u, rotInd[d]) << "\n");
#endif
		}

	} // end particle loop

}

    
    template<typename TDomain, typename TAlgebra>
    bool MovingParticle<TDomain, TAlgebra>::mark_collision_area(SmartPtr<AdaptiveRegularRefiner_MultiGrid> refiner, int level)
    {
        
        //typedef typename TDomain::position_accessor_type	position_accessor_type;
        typedef typename GeomObjBaseTypeByDim<TDomain::dim>::base_obj_type elem_t;
        //typedef typename ElementStorage<elem_t>::SectionContainer::iterator BaseElemIterator;
        //typedef typename elem_t::side side_t;
        //typedef typename std::vector<elem_t>::iterator ElemIterator;
        typedef typename std::map<elem_t*,int>::iterator ElemMapIterator;
        
        typedef typename geometry_traits<elem_t>::iterator ElemIterator;
        
        bool elem_is_cut_by_2 = false;
        
        std::map<elem_t*,int> num_cuts;
        
        //	access the grid and the position attachment
        MultiGrid* mg = (refiner->multi_grid());
        //position_accessor_type& aaPos = dom.position_accessor();
        
        for (ElemIterator elemIter = mg->begin<elem_t>(); elemIter != mg->end<elem_t>(); ++elemIter) {
            elem_t* elem = *elemIter;
            for (size_t p = 0; p < m_spParticleHandlerGlobal->num_particles(); ++p) {
                //std::vector<elem_t> ElemListLog = m_vvvElemListCut[level][p];
                //for (ElemIterator elem = ElemListLog.begin();
                //elem != ElemListLog.end(); ++elem) {
                bool cut = false;
                for (size_t v = 0; v < elem->num_vertices(); ++v) {
                    const MathVector<dim>& vrtPos = m_spParticleHandlerGlobal->m_aaPos[elem->vertex(v)];
                    const MathVector<dim>& center = m_spParticleHandlerGlobal->get_center(p);
                    const number radius = m_spParticleHandlerGlobal->get_radius(p);
                    
                    if (VecDistance(vrtPos, center) <= radius) {
                        cut = true;
                    }
                }
                if (cut) {
                    num_cuts[elem]++;
                }
            }
        }
        for (ElemMapIterator cut_elem = num_cuts.begin(); cut_elem != num_cuts.end(); ++cut_elem) {
            if (cut_elem->second >= 2 && !mg->has_children(cut_elem->first)) {
                refiner->mark(cut_elem->first);
                elem_is_cut_by_2 = true;
            }
        }
        return elem_is_cut_by_2;
    }
    
    template<typename TDomain, typename TAlgebra>
    double MovingParticle<TDomain, TAlgebra>::estimate_repulsive_force_parameters(vector_type& u, int topLevel, double maxElemDiameter, double deltaT){
        // Get and synchronize velocities
        const int levIndex = get_Index(GridLevel(topLevel, GridLevel::LEVEL));
        
        std::vector<MathVector<dim> > transVel;
        int num_particles = m_spParticleHandlerGlobal->num_particles();
        transVel.resize(num_particles);
        
        bool verbose = false;
        MathVector<dim> transVelP;

#ifdef UG_PARALLEL
        pcl::ProcessCommunicator com;
        
        for (size_t p = 0; p < num_particles; ++p) {
            
            std::vector<grid_base_object*> ElemList =
            m_spParticleHandlerGlobal->m_vvvElemListCut[levIndex][p];
            if (ElemList.size() == 0) {
                for (int d = 0; d < dim; ++d) {
                    transVelP[d] = 0.0;
                }
            } else {
                std::vector < DoFIndex > transInd =
                m_spParticleHandlerGlobal->get_transInd(levIndex, p);
                
                for (int d = 0; d < dim; ++d) {
                    transVelP[d] = DoFRef(u, transInd[d]);
                }
                
            }
            if (verbose) {
                UG_LOG(pcl::ProcRank()<< "sends (" << transVelP[0] << "," << transVelP[1] << ")\n");
            }
            com.allreduce(&transVelP[0],&transVel[p][0], dim, MPI_DOUBLE, PCL_RO_SUM);
            if (verbose) {
                if (dim == 1) {
                    UG_LOG("transVel[p] = ("<< transVel[p][0] <<")\n");
                } else{
                    if (dim == 2) {
                        UG_LOG("transVel[p] = ("<< transVel[p][0] << "," << transVel[p][1] <<")\n");
                    } else {
                        if (dim == 3) {
                            UG_LOG("transVel[p] = ("<< transVel[p][0] << "," << transVel[p][1] << "," << transVel[p][2] <<")\n");
                        }
                    }
                }
            }
        }
#else
        for (size_t p = 0; p < num_particles; ++p) {
            std::vector < DoFIndex > transInd =
            m_spParticleHandlerGlobal->get_transInd(levIndex, p);
            
            for (int d = 0; d < dim; ++d) {
                transVelP[d] = DoFRef(u, transInd[d]);
            }
        }
#endif
        
        double eps = 0.0;
        //const int levIndex = get_Index(GridLevel(13957 /*topLevel*/, GridLevel::LEVEL));
        
        for (size_t p = 0; p < num_particles; ++p){
            // Get values for particle p
            MathVector<dim> center_p = m_spParticleHandlerGlobal->get_center(p);
            VecScaleAdd(center_p, 1.0, center_p, deltaT, transVel[p]);
            number radius_p = m_spParticleHandlerGlobal->get_radius(p);
            //double Mass_p = m_spInterfaceMapper->Mass(13957,p);
            for (size_t q = 0; q < num_particles; ++q){
                if (p == q) {
                    continue;
                }
                if (verbose)
                    UG_LOG("p = "<<p<<", q = "<<q<<"\n");
                // Get values for particle q
                MathVector<dim> center_q = m_spParticleHandlerGlobal->get_center(q);
                VecScaleAdd(center_q, 1.0, center_q, deltaT, transVel[q]);
                number radius_q = m_spParticleHandlerGlobal->get_radius(q);
                //double Mass_q = m_spInterfaceMapper->Mass(13957,q);
                // Calculate force
                if (verbose) {
                    if (dim == 1) {
                        UG_LOG("center_p = ("<< center_p[0] <<")\n");
                    } else{
                        if (dim == 2) {
                            UG_LOG("center_p = ("<< center_p[0] << "," << center_p[1] <<")\n");
                        } else {
                            if (dim == 3) {
                                UG_LOG("center_p = ("<< center_p[0] << "," << center_p[1] << "," << center_p[2] <<")\n");
                            }
                        }
                    }
                    if (dim == 1) {
                        UG_LOG("center_q = ("<< center_q[0] <<")\n");
                    } else{
                        if (dim == 2) {
                            UG_LOG("center_q = ("<< center_q[0] << "," << center_q[1] <<")\n");
                        } else {
                            if (dim == 3) {
                                UG_LOG("center_q = ("<< center_q[0] << "," << center_q[1] << "," << center_q[2] <<")\n");
                            }
                        }
                    }
                }
                number centerDist = VecDistance(center_q, center_p);
                if (verbose)
                    UG_LOG("centerDist = "<< centerDist <<"\n");
                number s_ij = centerDist - radius_p - radius_q;
                if (verbose)
                    UG_LOG("centerDist - radius_p - radius_q = "<< s_ij <<"\n");
                number dist = maxElemDiameter - s_ij;
                if (verbose)
                    UG_LOG("maxElemDiameter - s_ij  = " << dist <<"\n");
                eps = std::max(eps,dist);
            }
        }
        return eps;
    }


} // end namespace NavierStokes
} // end namespace ug

#endif /* MOVING_PARTICLE_IMPL_H_ */
