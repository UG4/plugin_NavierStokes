/*
 * moving_particle.h
 *
 *  Created on: 20.01.2015
 *      Author: suze
 */

#ifndef MOVING_PARTICLE_H_
#define MOVING_PARTICLE_H_

#ifdef UG_PARALLEL
    #include "lib_grid/parallelization/load_balancer_util.h"
    #include "../../../../Parmetis/src/unificator_interface.h"
    #include "lib_grid/grid/neighborhood.h"
    #include "lib_grid/grid/neighborhood_util.h"              // for GetConnectedNeighbor
    #include "lib_grid/grid/grid_util.h"
    #include "lib_grid/grid/grid.h"
    #include "lib_grid/grid/grid_base_object_traits.h"        // for geometry_traits
#endif

#include "../../incompressible_navier_stokes_base.h"

#include "lib_disc/spatial_disc/immersed_util/interface_handler/interface_handler_flat_top_cut/interface_handler_particle.h"
#include "loc_to_glob_mapper_particle.h"
#include "immersed_bnd_cond_particle.h"

// for method 'mark_collision_area()':
#include "lib_grid/refinement/adaptive_regular_mg_refiner.h"

#include <memory>
#include <ctime>


namespace ug{
namespace NavierStokes{



template <	typename TDomain, typename TAlgebra>
class MovingParticle
		: public IImmersedInterface<TDomain, TAlgebra>
{
	public:
	///	world Dimension
		static const int dim = TDomain::dim;

	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

		typedef typename domain_traits<dim>::grid_base_object grid_base_object;

 		MovingParticle(SmartPtr<IAssemble<TAlgebra> > ass,
                       SmartPtr<NavierStokesFV1_cutElem<TDomain> > spMaster,
                       SmartPtr<CutElementHandler_FlatTop<dim> > cutElementHandler,
                       number fluidDensity, number fluidKinVisc);

 	// destructor
		virtual ~MovingParticle(){};

    //////////////////////////////////////////////////////////////////////////////////////////
    // main methods
    //////////////////////////////////////////////////////////////////////////////////////////
    
    //  general initialisation of set up data;
    //  most important: call of 'update_interface_data()' during init() of CutElementHandler
    //          --> marks the cut elements and interface vertices as
    //              INSIDE, OUTSIDE, CUT_BY_INTERFACE(element)/ON_INTERFACE(vertx)
    //          --> sets up the element lists CutElementHandlerParticle::m_vvvElemList,
    //                  ::m_vvvElemListCut, ::m_vvvElemListOutside
	/// called via .lua:
        void init(vector_type& u, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, const int baseLevel,
              const int topLevel);
    
    // mainly updates the coordinates of the particles for the next time step and updates the marker accordingly
    //      --> for that the solution in the DoFs for the particle motion needs to be
    //          (1) copied from DoFRef(u, ) into external storage:     via call of 'store_particle_velocity()'
    //          (2) written back from storage to DoFRef(u, ):          via call of 'write_particle_velocity()'
        void update(vector_type& u, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, const int baseLevel,
                const int topLevel, const number time, const number deltaT);
    
#ifdef UG_PARALLEL
        void pre_balancing_update(vector_type& u, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace,
                              const int baseLevel, const int topLevel, const number time, number deltaT);
    
        void post_balancing_update(vector_type& u, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace,
                               const int baseLevel, const int topLevel, const number time, number deltaT);
#endif
    
    //////////////////////////////////////////////////////////////////////////////////////////
    // helper methods for init() and update()
    //////////////////////////////////////////////////////////////////////////////////////////
    
    // transfer direction:
    //  DoFRef(u, transInd) ---> ParticleProvider::m_vvLinearVelocity/m_vvAngularVelocity
        void store_particle_velocity(vector_type& u, const int topLevel);

    // transfer direction:
    //  ParticleProvider::m_vvLinearVelocity/m_vvAngularVelocity ---> DoFRef(u, transInd)
        void write_particle_velocity(vector_type& u, const int topLevel);
    
    //////////////////////////////////////////////////////////////////////////////////////////
    // write solution to 'FREED' nodes, i.e. nodes which were inside particle in the last
    //  time step and are outside the particle (i.e. freed) in the current time step
    //		--> call method via .lua BEFORE 'solTimeSeries:push_discard_oldest(oldestSol, time)':
    //			(if not: in case that outside nodes are inside the domain AFTER updating the
    //           particle coordinates, NO solution will be defined here!)
        void fill_freed_nodes(vector_type& u, const int topLevel, const number time);
    
    //////////////////////////////////////////////////////////////////////////////////////////
    // helper methods for init() and update()
    //////////////////////////////////////////////////////////////////////////////////////////
    
    // currently not used
        void compute_gradient_local_max(vector_type& sol, vector_type& grad,
                                        SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, const int topLevel);
   
    /// writes the physical velocity values to the nodes lying in the particle domain
    //   (call via lua for visualisation purposes)
		void adjust_global_solution(vector_type& u, const int topLevel);


    //////////////////////////////////////////////////////////////////////////////////////////
    /// methods for setting up the forces for collision (by Jonas Simon, see Master Thesis)
    //////////////////////////////////////////////////////////////////////////////////////////
    
        void set_repulsive_force(bool repForce,  number forceValue)
            { m_spInterfaceMapper->set_repulsive_force(repForce, forceValue); }
        void set_glowinski_repulsive_force(bool maxRepForce, number rho, number eps)
            { m_spInterfaceMapper->set_glowinski_repulsive_force(maxRepForce, rho, eps); }
        void set_minimum_correction_force(bool EquiRepForce, number repulsiveDistance)
            { m_spInterfaceMapper->set_minimum_correction_force(EquiRepForce, repulsiveDistance); }
    
        bool mark_collision_area(SmartPtr<AdaptiveRegularRefiner_MultiGrid> refiner, int level);
        double estimate_repulsive_force_parameters(vector_type& u, int topLevel, double MaxElemDiameter, double deltaT);
        void set_forceLog(bool val)  { m_spInterfaceMapper->set_forceLog(val); }
        void set_mpi_routine(int val){ m_spCutElementHandler->set_mpi_routine(val); }
    
    
    //////////////////////////////////////////////////////////////////////////////////////////
    // some getter and setter methods
    //////////////////////////////////////////////////////////////////////////////////////////

        SmartPtr<ParticleBndCond<TDomain> > get_BndCond() { return m_spInterfaceBndCond; }

    /// checks if grid data is updated and returns the 'levIndex' of the 'gridLevel' via access to 'CutElementHandlerBase::m_Map'
        int get_Index(const GridLevel& gridLevel);
    
    /// returns true, if computation is time dependent
        bool is_time_dependent() { return m_spInterfaceMapper->is_time_dependent();}
    

    //////////////////////////////////////////////////////////////////////////////////////////
    /// lua-methods for set up:
    //////////////////////////////////////////////////////////////////////////////////////////
    
    // the 'threshold' defines the bandwidth around the immersed interface, in which a node
    //  counts as 'OUTSIDE' or 'ON_INTERFACE' during call of 'CutElementHandler::is_outside()'
    //  and 'CutElementHandler::is_nearInterface()
        void initialize_threshold(TDomain& domain, const int baseLevel, const int topLevel);
        void set_threshold(size_t level, const number threshold)
            { m_spCutElementHandler->set_threshold(level, threshold); }
        number MeanElementDiameter(TDomain& domain, int level);
    
    // If StdFV-assembling is ON, NO new 'vCornerCoords' will be computed on the cut elements.
    // The original nodes ACCROSS the interface (ON the euclidian mesh) will be chosen for the
    // computation of the solution of the interface
    //  ==> the standard shape functions w.r.t. the non-conforming mesh will be used
    //       (similar as in common 'ficticious domain' methods)
    //  ==> the shape functions will NOT be 1 ON the interface and the gradient will NOT
    //          point normal to the interface
        void set_StdFV_assembling(bool bValue)              { m_spInterfaceHandlerLocal->set_StdFV_assembling(bValue);}
        bool StdFV_assembling()                             { return m_spInterfaceHandlerLocal->StdFV_assembling(); }
    
    // setting and getting flag for printing of cut-element data into file:
        void set_print_cutElemData(bool bValue) { m_spInterfaceHandlerLocal->set_print_cutElemData(bValue); }

        void set_element_diameter(double val)               { m_spInterfaceMapper->set_element_diameter(val); }
        void set_volume_comp_mode(bool bVolumeCompMode)     { m_spInterfaceMapper->set_volume_comp_mode(bVolumeCompMode); }
        void set_gravity(bool gravity, number gravityConst) { m_spInterfaceMapper->set_gravity(gravity, gravityConst); }
        void set_time_step(number dt)                       { m_spInterfaceMapper->set_time_step(dt); }

    // writes the velocity of the prtIndex-th particle into 'transSol' and 'rotSol' for further computatino or output
        void get_velocity(MathVector<dim>& transSol, MathVector<dim>& rotSol, const vector_type& u,
                      const int topLevel, number deltaT, const size_t prtIndex);
    

    //////////////////////////////////////////////////////////////////////////////////////////
    // lua-methods for output
    //////////////////////////////////////////////////////////////////////////////////////////

    /// writing data to file; called via .lua
        void print_velocity(const vector_type& u, const int topLevel, number time, const char* filename);
    
    // writes the pressure on the front and back of the cylinder and the delta of it into a file
        void print_deltaP(const vector_type& u, const int topLevel);
    
    // writes the pressure value along the circular interface in grid nodes into file
        void print_pressure_nodal(const vector_type& u, const int topLevel);
    
    // writes the pressure value along the circular interface for equidistant angle teta into file
        void print_pressure_teta(const vector_type& u, const int topLevel);

    
    // interpolates the value at a non-grid point (used for the computation of deltaP in 'print_deltaP()'
        void interpolate_point(ConstSmartPtr<DoFDistribution> dd, const vector_type& u,
                           const MathVector<dim>& evalPos, MathVector<dim+1>& interpolation);
    
    // returns the number of DoFs on the original grid
        size_t get_numCutElements(const int gridlevel, const size_t prtIndex)
            {
                ConstSmartPtr<DoFDistribution> dd = m_spApproxSpace->dof_distribution(GridLevel(gridlevel, GridLevel::LEVEL));
                const int levIndex = m_spCutElementHandler->get_Index(gridlevel, dd);
                
                return m_spCutElementHandler->get_numCutElements(levIndex, prtIndex);
            }

    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    // class member
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

	private:
    // new member
		SmartPtr<CutElementHandler_FlatTop<dim> > m_spCutElementHandler;

	// member from base class
        SmartPtr<InterfaceHandlerLocalParticle<dim> > m_spInterfaceHandlerLocal;    // contains class member
                                                                                    //  'CutElementHandler_FlatTop'
		SmartPtr<ParticleMapper<TDomain, TAlgebra> > m_spInterfaceMapper;           // contains class member
                                                                                    //  'IInterfaceHandlerLocal'
		SmartPtr<ParticleBndCond<TDomain> > m_spInterfaceBndCond;                   // contains class member
                                                                                    //  'IInterfaceHandlerLocal'

	///	current ApproxSpace
		SmartPtr<ApproximationSpace<TDomain> > m_spApproxSpace;


};  // end class MovingParticle

    
    
    
#ifdef UG_PARALLEL
template<typename TDomain>
class ParticleUnificator : public parmetis::IUnificator<typename GeomObjBaseTypeByDim<TDomain::dim>::base_obj_type>
{

    public:
        typedef typename GeomObjBaseTypeByDim<TDomain::dim>::base_obj_type elem_t;
        typedef typename elem_t::side side_t;
        typedef Attachment<int> AElemIndex;
        typedef Attachment<std::vector<int> > AElemIndices;
        
        ///	world Dimension
        static const int dim = TDomain::dim;
        
        ///	Type of position coordinates
        typedef MathVector<dim> position_type;
        
        ///	Type of Position Attachment
        typedef Attachment<position_type> position_attachment_type;
        
        ///	Type of Accessor to the Position Data Attachment
        typedef Grid::VertexAttachmentAccessor<position_attachment_type>
        position_accessor_type;
        
        /// constructor
        ParticleUnificator(SmartPtr<TDomain> spDom)
        : m_spMG(spDom->grid())
        {
            //	get position attachment
            m_aPos = GetDefaultPositionAttachment<position_attachment_type>();
            
            // 	let position accessor access Vertex Coordinates
            if(!m_spMG->has_attachment<Vertex>(m_aPos))
                m_spMG->attach_to<Vertex>(m_aPos);
            m_aaPos.access(*m_spMG, m_aPos);
        }
        
        void update_particles(SmartPtr<ParticleProvider<dim> > provider )
        {
            m_vParticleCoord.clear();
            
            size_t numPrt = provider->num_particles();
            for ( size_t i = 0; i < numPrt; ++i )
                m_vParticleCoord.push_back(std::make_pair(provider->get_radius(i), provider->get_center(i)));
            
            UG_LOG("in update_particles:\n")
            for ( size_t i = 0; i < m_vParticleCoord.size(); ++i )
                UG_LOG("m_vParticleCoord: " << m_vParticleCoord[i].second << "\n");
            
            UG_LOG("in update_particles end....\n\n\n");
            
        }
        
        virtual void unify(
                           MultiGrid* mg,
                           int lvl,
                           int localOffset,
                           const Grid::AttachmentAccessor<elem_t, AElemIndex>& aaElemInd, // local indices!
                           const Grid::AttachmentAccessor<side_t, AElemIndices>& aaSideElemInd, // global indices!
                           std::vector<std::pair<int, int> >& unificationPairs) const // global indices!
        {
            const DistributedGridManager& dgm = *mg->distributed_grid_manager();
            
            std::vector<elem_t*> vNeighbors;
            
            //UG_LOG("unify() on lvl " << lvl << "\n");
            
            typename geometry_traits<elem_t>::iterator iterEnd, iter;
            iter = mg->begin<elem_t>(lvl);
            iterEnd = mg->end<elem_t>(lvl);
            
            // iterate all elements
            for(; iter != iterEnd; ++iter) {
                elem_t* elem = *iter;
                
                int locInd = aaElemInd[elem];
                
                //UG_LOG("element is " << locInd << "..\n");
                
                // skip ghosts
                if (locInd == -1) {
                    continue;
                }
                
                // only process particle elements.
                if (is_PartOfParticle(elem)) {
                    //UG_LOG("local index is " << locInd << " and global index is " << locInd + localOffset << "\n");
                    
                    //std::vector<side_t*> sides;
                    typename Grid::traits<side_t>::secure_container sides;
                    //CollectAssociated(sides, *mg, elem, true);
                    mg->associated_elements(sides, elem);
                    //UG_LOG("sides has length " << sides.size() << "\n");
                    // iterate over current element's sides
                    for (size_t i = 0; i < sides.size(); ++i) {
                        //UG_LOG("Check side between elements " << aaSideElemInd[sides[i]][0] << " and " << aaSideElemInd[sides[i]][1] << ".\n");
                        // check if neighbor on the other side is stored on another process
                        // if not use aaElemInd with local indices
                        if (!dgm.is_in_horizontal_interface(sides[i]))	// QUESTION: Check elem_t or side_t here?
                        {
                            // get neighbors of element
                            //CollectNeighbors(vNeighbors, elem, *mg);
                            
                            // iterate over neighbor elements
                            /*for (elem_t* neighbor : vNeighbors) {
                             
                             // check if neighbor element is inside particle too
                             if (is_PartOfParticle(neighbor)) {
                             UG_LOG("Add " << aaElemInd[elem]+localOffset << "(" << localOffset << ") and " << aaElemInd[neighbor]+localOffset << " to unificationPairs\n");
                             unificationPairs.push_back(std::make_pair(aaElemInd[elem] + localOffset, aaElemInd[neighbor] + localOffset));
                             }
                             }*/
                            
                            // get back the elements to current side and add them to unificationPairs
                            typename Grid::traits<elem_t>::secure_container elemList;
                            mg->associated_elements(elemList, sides[i]);
                            UG_COND_THROW(elemList.size() > 2, "More than two " << elemList[0]->reference_object_id()
                                          << "s associated with " << sides[i]->reference_object_id() << ".");
                            
                            if (elemList.size() == 2)
                            {
                                int locInd0 = aaElemInd[elemList[0]];
                                int locInd1 = aaElemInd[elemList[1]];
                                
                                // exclude ghosts
                                if (locInd0 == -1 || locInd1 == -1)
                                    continue;
                                
                                // add to unification pairs
                                //UG_LOG("(1) Add " << locInd0+localOffset << "(" << localOffset << ") and " << locInd1+localOffset << " to unificationPairs\n");
                                unificationPairs.push_back(std::make_pair(locInd0+localOffset, locInd1+localOffset));
                            } else {
                                //UG_LOG("elemList.size() = " << elemList.size() << " != 2 shouldn't appear.\n");
                            }
                        }
                        // if other processes have to be addressed for neighbors, use aaSideElemInd with global indices instead
                        else{
                            //UG_LOG("Elements on different processes\n");
                            const std::vector<int>& inds = aaSideElemInd[sides[i]];
                            
                            // if element side is part of global boundary aaSideElemInd[elem] holds only the element itself.
                            if (inds.size() == 1)
                                continue;
                            
                            UG_ASSERT(inds.size() == 2, "Indices attachment vector has more or less than 2 entries.");
                            //UG_LOG("(2) Add " << inds[0] << " and " << inds[1] << " to unificationPairs\n");
                            unificationPairs.push_back(std::make_pair(inds[0], inds[1]));
                        }
                    }
                }
            }
            //UG_LOG("unify() finished.\n");
        }
        
        /**
         * Function that estimates if any particle is close to any process border.
         */
        bool rebalance(MultiGrid* mg, int lvl){
            typename geometry_traits<elem_t>::iterator iterEnd, iter;
            iter = mg->begin<elem_t>(lvl);
            iterEnd = mg->end<elem_t>(lvl);
            
            int repart = 0;
            
            // iterate over all elements
            for(; iter != iterEnd; ++iter) {
                elem_t* elem = *iter;
                typename Grid::traits<side_t>::secure_container sides;
                mg->associated_elements(sides, elem);
                
                // iterate over current element's sides
                for (size_t i = 0; i < sides.size(); ++i) {
                    // get back the elements to current side and add them to unificationPairs
                    typename Grid::traits<elem_t>::secure_container elemList;
                    mg->associated_elements(elemList, sides[i]);
                    // if border element
                    //        			UG_LOG("elemList.size(): " << elemList.size())
                    if (elemList.size() < 2) {
                        if (is_PartOfParticle(elemList[0])) {
                            repart = 1;
                        }
                    }
                }
            }
            
            pcl::ProcessCommunicator com;
            int repart_all = 0;
            com.allreduce(&repart, &repart_all, 1, MPI_INT, PCL_RO_SUM);
            
            return repart_all;
        }
        
    private:
        bool is_PartOfParticle(elem_t* elem) const
        {
            bool part = false;
            for (size_t p = 0; p < m_vParticleCoord.size(); ++p) {
                // special case if particle is smaller than one element of the multigrid's base level -> ist das überhaupt relevant? kann ein Baselevelelement partitioniert werden?
                if (ContainsPoint(elem, m_vParticleCoord[p].second, m_aaPos)) {
                    //UG_LOG("Particle center is inside this element.\n");
                    return true;
                }
                // check if at least one vertex is inside particle
                for(size_t i = 0; i < elem->num_vertices(); ++i) {
                    if (is_insideParticle(elem->vertex(i),p)) {
                        return true;
                    }
                }
                //				// exclude triangle-circle overlap
                //				for(size_t i = 0; i < elem->num_vertices(); ++i) {
                //					for(size_t j = i; j < elem->num_vertices(); ++j) {
                //						if (intersection(m_vParticleCoord[p].second,m_vParticleCoord[p].first,m_aaPos[elem->vertex(i)],m_aaPos[elem->vertex(j)])) {
                //							return true;
                //						}
                //					}
                //				}
                
                MathVector<dim> elemCenter = m_aaPos[elem->vertex(0)];
                for(size_t i = 1; i < elem->num_vertices(); ++i) {
                    //					const MathVector<dim>& vrtPos = m_aaPos[elem->vertex(i)];
                    //					for (size_t j = 0; j < dim; ++j) {
                    //						elemCenter[j] += vrtPos[j];
                    //					}
                    elemCenter += m_aaPos[elem->vertex(i)];
                }
                elemCenter /= elem->num_vertices();
                //UG_LOG("Element Center is (" << elemCenter[0] << ", " <<elemCenter[1] << ").\n");
                number maxDistance = 0.0;
                for(size_t i = 1; i < elem->num_vertices(); ++i) {
                    number distance = VecDistance(m_aaPos[elem->vertex(i)],elemCenter);
                    maxDistance = (maxDistance > distance) ? maxDistance : distance;
                }
                //UG_LOG("maxDistance is " << maxDistance << ".\n");
                const number radius = m_vParticleCoord[p].first;
                const MathVector<dim>& center = m_vParticleCoord[p].second;
                number centerDistance = VecDistance(elemCenter, center);
                //UG_LOG("centerDistance is " << centerDistance << ".\n");
                part = (part or (radius + maxDistance > centerDistance));
            }
            return part;
            //return false;
        }
        
        bool is_insideParticle(Vertex* vrt, size_t prtIndex) const
        {
            // get distance to first radius as initial value
            const MathVector<dim>& vrtPos = m_aaPos[vrt];
            const MathVector<dim>& center = m_vParticleCoord[prtIndex].second;
            const number radius = m_vParticleCoord[prtIndex].first;
            
            if (VecDistance(vrtPos, center) <= radius) {
                return true;
            }
            return false; //set_nearInterface(vrt, prtIndex);
        }
        
        bool is_outside(Vertex* vrt)
        {
            //	loop over all centers and pick the index with minimal distance
            for (size_t p = 0; p < m_vParticleCoord.size(); ++p)
                if ( !is_insideParticle(vrt, p) )
                    return true;
            return false;
        }
        
        bool conn_is_cut_by_interface(side_t* conn)
        {
            // init data
            bool insideFluid = false;
            bool outsideFluid = false;
            
            //	loop vertices
            for(size_t i = 0; i < conn->num_vertices(); ++i)
            {
                //	check if inside a particle
                if (is_outside(conn->vertex(i))) outsideFluid = true;
                else 				 insideFluid  = true;
                
                //if (insideFluid && is_nearInterface(vrt))
                //	UG_THROW("in 'get_elem_modus()': case 'set_nearInterface(vrt) = true' not possible!\n");
                
            } // vertex loop
            
            // return elem type
            if (insideFluid && outsideFluid) {
                return true;
                //UG_LOG("edge " << conn->grid_data_index() << " should have weight " << m_weight);
            }
            return false;
            
        }
        
        bool intersection(MathVector<dim> center, double radius, MathVector<dim> vertex1, MathVector<dim> vertex2) const {
            number alpha;
            
            MathVector<dim> lineDir;
            MathVector<dim> rayDir;
            // lineDir = vertex1 - vertex2;
            VecSubtract(lineDir, vertex1, vertex2);
            // rayDir = vertex2 - center;
            VecSubtract(rayDir, vertex2, center);
            
            const number a = VecDot(lineDir, lineDir);
            const number b = 2.0 * VecDot(lineDir, rayDir);
            const number c = VecDot(vertex2, vertex2) + VecDot(center, center)
            - 2 * VecDot(vertex2, center) - radius * radius;
            
            const number discriminant = b * b - 4 * a * c;
            
            // check that 'vrtPosOut' and 'vrtPosIn' really lie on different sides of the circle:
            //if (discriminant < -1e-8)
            //	UG_LOG("discriminant is " << discriminant << " < -1e-8\n");
            //UG_THROW("Value of discriminant = " << discriminant << "\n");
            
            // discriminant = 0!
            //const number alpha1 = (-b - sqrt(discriminant)) / (2.0 * a);
            //const number alpha2 = (-b + sqrt(discriminant)) / (2.0 * a);
            
            //			if (alpha1 <= alpha2)
            //				alpha = alpha1;
            //			else
            //				alpha = alpha2;
            
            if (discriminant < 0)
                return false;
            else
                return true;
        }
        
        SmartPtr<MultiGrid> m_spMG;
        position_attachment_type m_aPos;	///<Position Attachment
        position_accessor_type m_aaPos;		///<Accessor
        
        std::vector<std::pair<number, MathVector<dim> > > m_vParticleCoord;

}; // end class ParticleUnificator
    
#endif

} // end namespace MovingParticle
} // end namespace ug


#include "moving_particle_impl.h"
#include "moving_particle_tools.h"



#endif /* MOVING_PARTICLE_H_ */
