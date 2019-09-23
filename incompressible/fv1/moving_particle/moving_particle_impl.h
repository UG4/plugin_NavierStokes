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
		SmartPtr<NavierStokesFV1_cutElem<TDomain> > spMaster,
		SmartPtr<CutElementHandler_FlatTop<dim> > cutElementHandler,
		number fluidDensity, number fluidKinVisc)
:       m_spCutElementHandler(cutElementHandler),
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
    
// needs to be enabled, in order to call 'spAssTuner->modify_LocalData()' during element disc assembling
	assAdapt->enable_modify_solution(true);

}

    
template<typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::
init(vector_type& u, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, const int baseLevel,
     const int topLevel)
{
 // get data
    m_spApproxSpace = spApproxSpace;
    ConstSmartPtr<DoFDistribution> dd = spApproxSpace->dof_distribution(GridLevel(topLevel, GridLevel::LEVEL));
    
// call of 'CutElementHadlerFT::update_interface_data(): for all given levels
    m_spCutElementHandler->template init<TDomain>(dd, baseLevel, topLevel);
    
//  transfer particle velocity:
//    FROM  m_vvLinearVelocity/m_vvAngularVelocity  (initialized during 'ParticleProvider::add()'
//    TO    DoFRef(u, transInd)/DoFRef(u, rotInd)   (for the later solution process)
     write_particle_velocity(u, topLevel);
  
}
 
    
template<typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::
update(vector_type& u, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, const int baseLevel,
       const int topLevel, const number time, const number deltaT)
{
        int topLev = spApproxSpace->num_levels()-1;
        if ( topLev != topLevel )
            UG_THROW("NavierStokes::update: parameter 'topLevel' = " << topLevel << " != "
                     << topLev << "current top leven! \n");
        
    // A. copy particle velocity to data: DoFRef(u, transInd/rotInd) ---> m_vvLinearVelocity/m_vvAngularVelocity
    // reverse direction done during 'NavierStokes::write_particle_velocity()' --> see F.
        store_particle_velocity(u, topLevel);
        
#ifdef UG_PARALLEL
    // B. synchronize particle data when everything in ParticleProvider is up to date.
        const int levIndex = get_Index(GridLevel(topLevel, GridLevel::LEVEL));
        m_spCutElementHandler->synchronize_particles(levIndex);
#endif
    
    // C. calculate new particle coordinates
        m_spCutElementHandler->update_prtCoords(topLevel, deltaT);
    
    // D. write solution to 'FREED' nodes, i.e. nodes which were inside particle in the last
    //  time step and are outside the particle (i.e. freed) in the current time step
    //      --> needs to be called AFTER 'update_prtCoords(), since being inside particle refers
    //          already to new coordinates
    //      --> needs to be called BEFORE 'm_spCutElementHandler->template init(), since being
    //          inside particle refers to OLD BoolMarker
        fill_freed_nodes(u, topLevel, time);
        
    // E. update marker and global indices => new location for (u, transInd), (u, rotInd) !
        ConstSmartPtr<DoFDistribution> dd = spApproxSpace->dof_distribution(GridLevel(topLevel, GridLevel::LEVEL));
        m_spCutElementHandler->template init<TDomain>(dd, baseLevel, topLevel);
        
    // F. write solution into new location of global indices 'transInd' and 'rotInd':
    //  transfer direction:  m_vvLinearVelocity/m_vvAngularVelocity ---> DoFRef(u, transInd)/DoFRef(u, rotInd)
        write_particle_velocity(u, topLevel);
        
    // reset volume:
        m_spInterfaceMapper->reset_volume();
}
    
template<typename TDomain, typename TAlgebra>
int MovingParticle<TDomain, TAlgebra>::
get_Index(const GridLevel& gridLevel)
{
    ConstSmartPtr<DoFDistribution> dd = m_spApproxSpace->dof_distribution(gridLevel);
    
    const int levIndex = m_spCutElementHandler->get_Index(gridLevel, dd);
    
    return levIndex;
}
    
template<typename TDomain, typename TAlgebra>
number MovingParticle<TDomain, TAlgebra>::
MeanElementDiameter(TDomain& domain, int level)
{
    
	typedef typename domain_traits<TDomain::dim>::grid_base_object TElem;
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
void MovingParticle<TDomain, TAlgebra>::
initialize_threshold(TDomain& domain, const int baseLevel, const int topLevel)
{
	if (baseLevel < 0)
		UG_THROW("initialize_threshold(): no cast of baselevel from 'int' tp 'size_t' possible! \n");
	if (topLevel < 0)
		UG_THROW("initialize_threshold(): no cast of toplevel from 'int' tp 'size_t' possible! \n");

	typedef typename domain_traits<TDomain::dim>::grid_base_object TElem;

// compute level-dependent value for threshold:
	for (size_t lev = baseLevel; lev <= (size_t) topLevel; ++lev) {
		const number maxLength = MaxElementDiameter(domain, lev);
		const number meanLength = MeanElementDiameter(domain, lev);
		UG_LOG("maxLength = " << maxLength << "\n");
		UG_LOG("meanLength = " << meanLength << "\n");
		UG_LOG("threshold_max = " << maxLength*maxLength << "\n");
		UG_LOG("threshold_mean = " << meanLength*meanLength << "\n");

		set_threshold(lev, meanLength * meanLength);

	}

	UG_LOG("----------------- END initialize_threshold() ---------------- \n");

}
    

// transfer direction: DoFRef(u, transInd) ---> m_vvLinearVelocity/m_vvAngularVelocity
template<typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::
store_particle_velocity(vector_type& u, const int topLevel)
{
    bool output = false;
	const int levIndex = get_Index(GridLevel(topLevel, GridLevel::LEVEL));

	for (size_t p = 0; p < m_spCutElementHandler->num_particles(); ++p) {

#ifdef UG_PARALLEL
    // use size of member 'CutElementHandler_FlatTop::m_vvvElemListCut' in order
    // to indicate, whether a particle lies on a processor or not
		std::vector<grid_base_object*> ElemList = m_spCutElementHandler->m_vvvElemListCut[levIndex][p];
		UG_LOG("1 MovingParticle::store_particle_velocity: ElemList.size(): " << ElemList.size() << "\n");
		if (ElemList.size() == 0) {
			UG_LOG("2 MovingParticle::store_particle_velocity: ElemList.size(): "
                   << ElemList.size() << " => skip assembling! \n");
			continue;
		}
#endif
        if ( output ){
		UG_LOG("copy_solution(): VORHER transSol(0) = " << m_spCutElementHandler->get_transSol(p,0) << "\n");
		UG_LOG("copy_solution(): VORHER transSol(1) = " << m_spCutElementHandler->get_transSol(p, 1) << "\n");
		UG_LOG("copy_solution(): VORHER rotSol(0) = " << m_spCutElementHandler->get_rotSol(p,0) << "\n");
		UG_LOG("copy_solution(): VORHER rot1Sol(1) = " << m_spCutElementHandler->get_rotSol(p, 1) << "\n");
        }
		
    // 'get_transInd()' returns NOT YET updated indices (updated during 'movingParticle:update()' vie lua)
		std::vector < DoFIndex > transInd =
				m_spCutElementHandler->get_transInd(levIndex, p);
		std::vector < DoFIndex > rotInd   =
                m_spCutElementHandler->get_rotInd(levIndex, p);

		UG_LOG(" transInd(0) = " << transInd[0] << "\n");
		UG_LOG(" rotInd(0) = " << rotInd[0] << "\n");

    // finally: write solution from m_vvLinearVelocity/m_vvAngularVelocity to DoFRef and set
    //          solution in DoFRef to zero
		for (int d = 0; d < dim; ++d)
        {
			number solution = DoFRef(u, transInd[d]);
			m_spCutElementHandler->set_extraSolTrans(solution, p, 0, d);
			DoFRef(u, transInd[d]) = 0.0;

			solution = DoFRef(u, rotInd[d]);
			m_spCutElementHandler->set_extraSolRot(solution, p, 0, d);
			DoFRef(u, rotInd[d]) = 0.0;
		}
        
        if ( output ){
		UG_LOG("copy_solution(): NACHHER transSol(0) = " << m_spCutElementHandler->get_transSol(p,0) << "\n");
		UG_LOG("copy_solution(): NACHHER transSol(1) = " << m_spCutElementHandler->get_transSol(p, 1) << "\n");
		UG_LOG("copy_solution(): NACHHER rotSol(0) = " << m_spCutElementHandler->get_rotSol(p,0) << "\n");
		UG_LOG("copy_solution(): NACHHER rot1Sol(1) = " << m_spCutElementHandler->get_rotSol(p, 1) << "\n");
        }
        
	} // end particle loop

}

    
    
template <typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::
get_velocity(MathVector<dim>& transSol, MathVector<dim>& rotSol, const vector_type& u,
                 const int topLevel, number deltaT, const size_t prtIndex)
{
    // get data
        const int levIndex = get_Index(GridLevel(topLevel, GridLevel::LEVEL));
        
#ifdef UG_PARALLEL
        // use size of member 'CutElementHandler_FlatTop::m_vvvElemListCut' in order
        // to indicate, whether a particle lies on a processor or not
        std::vector<grid_base_object*> ElemList = m_spCutElementHandler->m_vvvElemListCut[levIndex][prtIndex];
        UG_LOG("1 MovingParticle::get_velocity: ElemList.size(): " << ElemList.size() << "\n");
        if (ElemList.size() == 0) {
            UG_LOG("2 MovingParticle::get_velocity: ElemList.size(): " << ElemList.size() <<
                   " => skip assembling! \n");
            return;
        }
#endif
    // get multiindices for translation and rotation of particle
        std::vector < DoFIndex > transInd = m_spCutElementHandler->get_transInd(levIndex, prtIndex);
        std::vector < DoFIndex > rotInd = m_spCutElementHandler->get_rotInd(levIndex, prtIndex);
        
        
        for ( int d = 0; d < dim; ++d )
        {
            transSol[d] = DoFRef(u, transInd[d]);
            rotSol[d]	= DoFRef(u, rotInd[d]);
        }
        
}
    
template<typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::
compute_gradient_local_max(vector_type& sol, vector_type& grad,
		SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, const int topLevel)
{

// get data
	typedef typename domain_traits<dim>::grid_base_object grid_base_object;
	typename DoFDistribution::traits<grid_base_object>::const_iterator iter, iterEnd;

	ConstSmartPtr<DoFDistribution> dd = spApproxSpace->dof_distribution(
			GridLevel(topLevel, GridLevel::LEVEL));
	typename TDomain::position_accessor_type aaPos =
			m_spCutElementHandler->m_aaPos;

	iter = dd->template begin<grid_base_object>();
	iterEnd = dd->template end<grid_base_object>();

//	loop elements in order to compute the volume and set rhs:
	for (; iter != iterEnd; iter++) {
    //	get element
		grid_base_object* elem = *iter;

		if (m_spCutElementHandler->m_spOutsideMarker->is_marked(elem))
			continue;

		std::vector < DoFIndex > vInd1;
		std::vector < DoFIndex > vInd2;

		std::vector<Edge*> vEdges;
		CollectEdgesSorted(vEdges, *m_spCutElementHandler->m_spMG, elem);

		if (m_spCutElementHandler->m_spCutMarker->is_marked(elem)) {
			for (size_t e = 0; e < vEdges.size(); ++e) {
				Edge* edge = vEdges[e];
				std::vector<Vertex*> vVertexEdge;
				CollectVertices(vVertexEdge, *m_spCutElementHandler->m_spMG,
						edge);
				if (vVertexEdge.size() != 2)
					UG_THROW("error in collecting vertices associated to an edge!...EXIT!...\n");

				Vertex* vrt1 = vVertexEdge[0];
				Vertex* vrt2 = vVertexEdge[1];

            // no gradient needs to be computed on outside edges:
				if (m_spInterfaceHandlerLocal->is_onInterfaceVertex(vrt1)
						&& m_spInterfaceHandlerLocal->is_onInterfaceVertex(vrt2))
					continue;

				if (dd->inner_dof_indices(vrt1, dim, vInd1) != 1)
					UG_THROW("MovingParticle::compute_gradient(): Only one index expected.");
				if (dd->inner_dof_indices(vrt2, dim, vInd2) != 1)
					UG_THROW("MovingParticle::compute_gradient(): Only one index expected.");

            // compute local gradient along the edge
				number grad_edge = fabs(DoFRef(sol, vInd1[0]) - DoFRef(sol, vInd2[0]));
				number edgeLength = VecDistance(aaPos[vrt1], aaPos[vrt2]);

				if (!m_spInterfaceHandlerLocal->is_onInterfaceVertex(vrt1)
						&& !m_spInterfaceHandlerLocal->is_onInterfaceVertex(vrt2)) { // nothing has to be adapted!
				} else if (m_spInterfaceHandlerLocal->is_onInterfaceVertex(vrt1)) {
					if (m_spInterfaceHandlerLocal->is_onInterfaceVertex(vrt2))
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
				} else if (m_spInterfaceHandlerLocal->is_onInterfaceVertex(vrt2)) {
					if (m_spInterfaceHandlerLocal->is_onInterfaceVertex(vrt1))
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
				DoFRef(grad, vInd1[0]) = std::max(DoFRef(grad, vInd1[0]), grad_edge);
				DoFRef(grad, vInd2[0]) = std::max(DoFRef(grad, vInd2[0]), grad_edge);

			} // end edge-loop
		} else {
			for (size_t e = 0; e < vEdges.size(); ++e) {
				Edge* edge = vEdges[e];
				std::vector<Vertex*> vVertexEdge;
				CollectVertices(vVertexEdge, *m_spCutElementHandler->m_spMG,
						edge);
				if (vVertexEdge.size() != 2)
					UG_THROW(
							"error in collecting vertices associated to an edge!....EXIT!...\n");

				Vertex* vrt1 = vVertexEdge[0];
				Vertex* vrt2 = vVertexEdge[1];
				if (dd->inner_dof_indices(vrt1, dim, vInd1) != 1)
					UG_THROW("MovingParticle::compute_gradient(): Only one index expected.");
				if (dd->inner_dof_indices(vrt2, dim, vInd2) != 1)
					UG_THROW("MovingParticle::compute_gradient(): Only one index expected.");

				// compute local gradient along the edge
				number edgeLength = VecDistance(aaPos[vrt1], aaPos[vrt2]);
				number grad_edge = fabs(
						DoFRef(sol, vInd1[0]) - DoFRef(sol, vInd2[0]));
				grad_edge = grad_edge / edgeLength;

				// set local maximum over all gradients along the edge to a given vertex
				DoFRef(grad, vInd1[0]) = std::max(DoFRef(grad, vInd1[0]), grad_edge);
				DoFRef(grad, vInd2[0]) = std::max(DoFRef(grad, vInd2[0]), grad_edge);

			} // end edge-loop

		}

	} // end elem-loop

}

template<typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::
adjust_global_solution(vector_type& u, const int topLevel)
{
// get data
    const int levIndex = get_Index(GridLevel(topLevel, GridLevel::LEVEL));

// loop all particles
	for (size_t p = 0; p < m_spCutElementHandler->num_particles(); ++p)
    {

#ifdef UG_PARALLEL
    // use size of member 'CutElementHandler_FlatTop::m_vvvElemListCut' in order
    // to indicate, whether a particle lies on a processor or not
		std::vector<grid_base_object*> ElemList = m_spCutElementHandler->m_vvvElemListCut[levIndex][p];
		UG_LOG("1 MovingParticle::adjust_global_solution: ElemList.size(): " << ElemList.size() << "\n");
		if (ElemList.size() == 0) {
			UG_LOG("2 MovingParticle::adjust_global_solution: ElemList.size(): "
                   << ElemList.size() << " => skip assembling! \n");
			continue;
		}
#endif

		MathVector<dim> transSol = m_spCutElementHandler->get_transSol(p,0);
		MathVector<dim> rotSol = m_spCutElementHandler->get_rotSol(p, 0);

		ConstSmartPtr<DoFDistribution> dd = this->m_spApproxSpace->dof_distribution(GridLevel(topLevel, GridLevel::LEVEL));
		typename TDomain::position_accessor_type aaPos = m_spCutElementHandler->m_aaPos;
        typedef typename std::vector<grid_base_object*>::iterator ListIter;

		const MathVector<dim>& center = m_spCutElementHandler->get_center(p);


    // loop all elements relevant for p-th particle
		std::vector<grid_base_object*> ElemListLog = m_spCutElementHandler->m_vvvElemListOutside[levIndex][p];

		for (ListIter listIter = ElemListLog.begin(); listIter != ElemListLog.end(); ++listIter)
        {
        //	get element
			grid_base_object* elem = *listIter;

        //	collect all vertices of the element
			std::vector<Vertex*> vVertex;
			CollectVertices(vVertex, *m_spCutElementHandler->m_spMG, elem);

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
							m_spCutElementHandler->get_rotationMat(
									radialCo);

                //	set solution: u = U + w x r
					DoFRef(u, vInd[0]) = transSol[fct];
					for (int d = 0; d < dim; ++d)
						DoFRef(u, vInd[0]) += rotationMatCo[fct][d] * rotSol[d];

			//		DoFRef(u, vInd[0]) = 0.0;
				}
			} // end vrt-loop
		} // end outsideElem-loop

	} // end prt-loop

}

template<typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::
fill_freed_nodes(vector_type& u, const int topLevel, const number time)
{
// get data
	const int levIndex = get_Index(GridLevel(topLevel, GridLevel::LEVEL));

// loop all particles
	for (size_t p = 0; p < m_spCutElementHandler->num_particles(); ++p)
    {

#ifdef UG_PARALLEL
    // use size of member 'CutElementHandler_FlatTop::m_vvvElemListCut' in order
    // to indicate, whether a particle lies on a processor or not
		std::vector<grid_base_object*> ElemList =
				m_spCutElementHandler->m_vvvElemListCut[levIndex][p];
		UG_LOG("1 MovingParticle::fill_freed_nodes: ElemList.size(): " << ElemList.size() << "\n");
		if (ElemList.size() == 0) {
			UG_LOG("2 MovingParticle::fill_freed_nodes: ElemList.size(): "
                   << ElemList.size() << " => skip assembling! \n");
			continue;
		}
#endif

		ConstSmartPtr<DoFDistribution> dd = this->m_spApproxSpace->dof_distribution(GridLevel(topLevel, GridLevel::LEVEL));
		typename TDomain::position_accessor_type aaPos = m_spCutElementHandler->m_aaPos;
        typedef typename std::vector<grid_base_object*>::iterator ListIter;

		const MathVector<dim>& center = m_spCutElementHandler->get_center(p);
        
        MathVector<dim> transSol = m_spCutElementHandler->get_transSol(p,0);
        MathVector<dim> rotSol = m_spCutElementHandler->get_rotSol(p, 0);


    // loop all cut elements relevant for p-th particle
        std::vector<grid_base_object*> ElemListLog = m_spCutElementHandler->m_vvvElemListCut[levIndex][p];

		for (ListIter listIter = ElemListLog.begin();
				listIter != ElemListLog.end(); ++listIter)
        {
        //	get element
			grid_base_object* elem = *listIter;

        //	collect all vertices of the element
			std::vector<Vertex*> vVertex;
			CollectVertices(vVertex, *m_spCutElementHandler->m_spMG, elem);

        //	loop vertices
			for (size_t v = 0; v < vVertex.size(); ++v) {
            //	get vertex
				Vertex* vrt = vVertex[v];

            // ---> is_insideParticle_with_given_index() is a geometrical check w.r.t NEW particle coords
            //       (update_prtCoords() needs to be called BEFORE calling this method!)
            // ---> m_spOutsideMarker was marked w.r.t OLD particle coords
            //       (update_interface_data() has to be called AFTER calling this method!)
            //		===> !is_inside && is_outside = is_freed :-)
				if ( m_spCutElementHandler->is_insideParticle_with_given_index(p, vrt)
                  && m_spCutElementHandler->m_spOutsideMarker->is_marked(vrt) )
                {
					for (size_t fct = 0; fct < dim; ++fct) {
                    //	create multi index
						std::vector < DoFIndex > vInd;
                    //	get multi indices
						if (dd->inner_dof_indices(vrt, fct, vInd) != 1)
							UG_THROW("Only one index expected.");

                    // set solution: particle velocity
						MathVector<dim> radialCo;
						VecSubtract(radialCo, aaPos[vrt], center);
						MathMatrix<dim, dim> rotationMatCo = m_spCutElementHandler->get_rotationMat(radialCo);

                    //	set solution
						DoFRef(u, vInd[0]) = transSol[fct];
						for (int d = 0; d < dim; ++d)
							DoFRef(u, vInd[0]) += rotationMatCo[fct][d] * rotSol[d];


					}
                // if vertex lies outside particle and NOT near interface, i.e. !FlatTopVrtMarker = true,
                // then the particle moved to fast in the last timestep
                //      ==> CFL criterion not satisfied
					if (!m_spCutElementHandler->m_spInterfaceVrtMarker->is_marked(vrt)
                      && m_spCutElementHandler->m_spOutsideMarker->is_marked(vrt))
						UG_THROW("in 'fill_freed_nodes()': CFL criterion not satisfied!\n");
				}

			} // end vrt-loop
		} // end cutElem-loop

	} // end prt-loop

}


// call via lua BEFORE applying 'movingParticle:update()', which changes the global indices of trandSol/rotSol
// necessary for re-writing solution to NEW indices of vector 'u' during 'update_solution'
// 			=> handling instance/buffering data := m_spCutElementHandler->(m_vSolTrans/m_vSolRot)

// transfer direction:  m_vSolTrans ---> DoFRef(u, transInd)
template<typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::
write_particle_velocity(vector_type& u, const int topLevel)
{
	const int levIndex = get_Index(GridLevel(topLevel, GridLevel::LEVEL));

	for (size_t p = 0; p < m_spCutElementHandler->num_particles(); ++p)
    {
#ifdef UG_PARALLEL
    // use size of member 'CutElementHandler_FlatTop::m_vvvElemListCut' in order to indicate,
    // whether a particle lies on a processor or not
		std::vector<grid_base_object*> ElemList =
				m_spCutElementHandler->m_vvvElemListCut[levIndex][p];
		UG_LOG("1 MovingParticle::write_particle_velocity:  ElemList.size(): " << ElemList.size() << "\n");
		if (ElemList.size() == 0) {
			UG_LOG("2 MovingParticle::write_particle_velocity: ElemList.size(): "
                   << ElemList.size() << " => skip assembling! \n");
			continue;
		}
#endif

    // 'get_transInd()' returns allready updated indices (updated during 'movingParticle:update()' via lua)
		std::vector < DoFIndex > transInd = m_spCutElementHandler->get_transInd(levIndex, p);
		std::vector < DoFIndex > rotInd   = m_spCutElementHandler->get_rotInd(levIndex, p);

//		UG_LOG("in 'write_particle_velocity()': transInd = " << transInd[0] << "\n");
//		UG_LOG("in 'write_particle_velocity()': rotInd = " << rotInd[0] << "\n");

		for (int d = 0; d < dim; ++d)
        {
        // write solution to the GridFunction 'u' into the vertices, where their DoFs were allocated
			DoFRef(u, transInd[d]) = m_spCutElementHandler->get_transSol(p,0)[d];
			DoFRef(u, rotInd[d])   = m_spCutElementHandler->get_rotSol(p, 0)[d];
            
		}

	} // end particle loop

}

    
template<typename TDomain, typename TAlgebra>
bool MovingParticle<TDomain, TAlgebra>::
mark_collision_area(SmartPtr<AdaptiveRegularRefiner_MultiGrid> refiner, int level)
{
        
    typedef typename GeomObjBaseTypeByDim<TDomain::dim>::base_obj_type elem_t;
    typedef typename std::map<elem_t*,int>::iterator ElemMapIterator;
    typedef typename geometry_traits<elem_t>::iterator ElemIterator;
        
    bool elem_is_cut_by_2 = false;
        
    std::map<elem_t*,int> num_cuts;
        
//	access the grid and the position attachment
    MultiGrid* mg = (refiner->multi_grid());
    
    for (ElemIterator elemIter = mg->begin<elem_t>(); elemIter != mg->end<elem_t>(); ++elemIter) {
            elem_t* elem = *elemIter;
            for (size_t p = 0; p < m_spCutElementHandler->num_particles(); ++p) {
                //std::vector<elem_t> ElemListLog = m_vvvElemListCut[level][p];
                //for (ElemIterator elem = ElemListLog.begin();
                //elem != ElemListLog.end(); ++elem) {
                bool cut = false;
                for (size_t v = 0; v < elem->num_vertices(); ++v) {
                    const MathVector<dim>& vrtPos = m_spCutElementHandler->m_aaPos[elem->vertex(v)];
                    const MathVector<dim>& center = m_spCutElementHandler->get_center(p);
                    const number radius = m_spCutElementHandler->get_radius(p);
                    
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
double MovingParticle<TDomain, TAlgebra>::
estimate_repulsive_force_parameters(vector_type& u, int topLevel, double maxElemDiameter, double deltaT)
{
    
    // Get and synchronize velocities
        const int levIndex = get_Index(GridLevel(topLevel, GridLevel::LEVEL));
        
        std::vector<MathVector<dim> > transVel;
        size_t num_particles = m_spCutElementHandler->num_particles();
        transVel.resize(num_particles);
        
        bool verbose = false;
        MathVector<dim> transVelP;

#ifdef UG_PARALLEL
        pcl::ProcessCommunicator com;
        
        for (size_t p = 0; p < num_particles; ++p) {
            
            std::vector<grid_base_object*> ElemList =
            m_spCutElementHandler->m_vvvElemListCut[levIndex][p];
            if (ElemList.size() == 0) {
                for (int d = 0; d < dim; ++d) {
                    transVelP[d] = 0.0;
                }
            } else {
                std::vector < DoFIndex > transInd =
                m_spCutElementHandler->get_transInd(levIndex, p);
                
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
            m_spCutElementHandler->get_transInd(levIndex, p);
            
            for (int d = 0; d < dim; ++d) {
                transVelP[d] = DoFRef(u, transInd[d]);
            }
        }
#endif
        
        double eps = 0.0;
        //const int levIndex = get_Index(GridLevel(13957 /*topLevel*/, GridLevel::LEVEL));
        
        for (size_t p = 0; p < num_particles; ++p){
            // Get values for particle p
            MathVector<dim> center_p = m_spCutElementHandler->get_center(p);
            VecScaleAdd(center_p, 1.0, center_p, deltaT, transVel[p]);
            number radius_p = m_spCutElementHandler->get_radius(p);
            //double Mass_p = m_spInterfaceMapper->Mass(13957,p);
            for (size_t q = 0; q < num_particles; ++q){
                if (p == q) {
                    continue;
                }
                if (verbose)
                    UG_LOG("p = "<<p<<", q = "<<q<<"\n");
                // Get values for particle q
                MathVector<dim> center_q = m_spCutElementHandler->get_center(q);
                VecScaleAdd(center_q, 1.0, center_q, deltaT, transVel[q]);
                number radius_q = m_spCutElementHandler->get_radius(q);
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

#ifdef UG_PARALLEL
template<typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::
pre_balancing_update(vector_type& u, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace,
                          const int baseLevel, const int topLevel, const number time, number deltaT)
{
        int topLev = spApproxSpace->num_levels()-1;
        if ( topLev != topLevel )
            UG_THROW("NavierStokes::update: parameter 'topLevel' = " << topLevel << " != "
                     << topLev << "current top leven! \n");
        
        const char* filename = "update_times";
        std::string name(filename);
        char ext[50];
        sprintf(ext, "_%d.txt", pcl::ProcRank());
        //sprintf(ext, ".txt");
        name.append(ext);
        FILE* outputFile = fopen(name.c_str(), "a");
        
        std::clock_t begin = std::clock();
    // A. copy solution to data: DoFRef(u, transInd) ---> m_vSolTrans
    // reverse direction done during 'NavierStokes::write_particle_velocity()'
        store_particle_velocity(u, topLevel);
        std::clock_t end = std::clock();
        fprintf(outputFile,"store_particle_velocity %f", double(end - begin) / CLOCKS_PER_SEC);
        
        begin = std::clock();
        
    // C. calculate new particle coordinates
        m_spCutElementHandler->update_prtCoords(topLevel, deltaT);
        end = std::clock();
        fprintf(outputFile,"update_prtCoords %f", double(end - begin) / CLOCKS_PER_SEC);
        
        begin = std::clock();
    
    // B. synchronize particle data when everything in ParticleProvider is up to date.
        const int levIndex = get_Index(GridLevel(topLevel, GridLevel::LEVEL));
        m_spCutElementHandler->synchronize_particles(levIndex);
        
        end = std::clock();
        fprintf(outputFile,"synchronize_particles %f", double(end - begin) / CLOCKS_PER_SEC);
        
        begin = std::clock();
        
    // D. fill particle nodes with their real solution
    // 		--> FIRST store_particle_velocity() necessary, since eventually during
    //			fill extraDoF-nodes will be overwriten
        fill_freed_nodes(u, topLevel, time);
        end = std::clock();
        fprintf(outputFile,"fill_freed_nodes %f", double(end - begin) / CLOCKS_PER_SEC);
        
        fclose(outputFile);
        
}
    
    
template<typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::
post_balancing_update(vector_type& u, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, const int baseLevel, const int topLevel, const number time, number deltaT)
{
    // E. update data => new node for (u, transInd) !
        ConstSmartPtr<DoFDistribution> dd = spApproxSpace->dof_distribution(GridLevel(topLevel, GridLevel::LEVEL));
        m_spCutElementHandler->template init<TDomain>(dd, baseLevel, topLevel);
        
    // F. update solution from data: m_vSolTransDoFRef(u, transInd)
    //  transfer direction:  m_vSolTrans ---> DoFRef(u, transInd)
        write_particle_velocity(u, topLevel);
        
    // reset volume:
        m_spInterfaceMapper->reset_volume();
}
#endif  /* UG_PARALLEL */


    
} // end namespace NavierStokes
} // end namespace ug

#endif /* MOVING_PARTICLE_IMPL_H_ */
