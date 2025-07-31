/*
 * moving_particle.h
 *
 *  Created on: 20.01.2015
 *      Author: suze
 */

#ifndef PARTICLE_MAPPER_H_
#define PARTICLE_MAPPER_H_


#ifdef UG_PARALLEL
 	#include "lib_grid/parallelization/load_balancer_util.h"
#endif

#include "lib_disc/spatial_disc/immersed_util/immersed_interface_base.h"

namespace ug{
namespace NavierStokes{



template <typename TDomain, typename TAlgebra>
class ParticleMapper : public IInterfaceMapper<TAlgebra>
{
	public:
	///	World dimension
		static const int dim = TDomain::dim;

	/// abbreviation for pressure
		static const size_t _P_ = dim;

	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	///	Type of geometric base object
		typedef typename domain_traits<dim>::grid_base_object grid_base_object;

	/// used boundary face type
		typedef typename DimFV1FTGeometry<dim, dim, InterfaceHandlerLocalParticle<dim> >::BF interfaceBF;

	public:
	///	call base class constructor
 		ParticleMapper(SmartPtr<InterfaceHandlerLocalParticle<dim> > localHandler)
			 : m_spInterfaceHandlerLocal(localHandler),
			   m_spCutElementHandler(m_spInterfaceHandlerLocal->get_cutElementHandler()),
			   m_gravityConst(0.0),
               m_bGravity(true),
               m_bRepulsiveForce(false),
               m_bGlowRepulsiveForce(false),
               m_bMinimumCorrectionForce(false),
               m_bForceLog(false),
               m_rho(0.0),
               m_epsilon(0.0),
               m_repulsiveForce(0.0),
               m_repulsiveDistance(0.0),
               m_dt(0.0),
               m_bTimeDep(false),
               m_meanDiameter(0.0),
               m_volume(0.0),
               m_bVolumeCompExact(true),
               m_bUsualAss(false)
		{
 		// init boolian arrays for all registered particles
  			size_t numPrt = num_particles();
   			m_bFlagGravity.resize(numPrt, false);
  			m_bFlagInertial.resize(numPrt, false);
		};

 	///	destructor
 		virtual ~ParticleMapper() {};

	///////////////////////////////////////////////////////////////////////////////
 	///
 	///  base class methods not needed for this class
 	///
 	///////////////////////////////////////////////////////////////////////////////

        void adjust_mat_global(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd){};

    ///////////////////////////////////////////////////////////////////////////////
  	///
    ///  A. 	send local entries to global vector
    ///
    ///////////////////////////////////////////////////////////////////////////////

	/// base method
		virtual void add_local_vec_to_global(vector_type& vec, const LocalVector& lvec, ConstSmartPtr<DoFDistribution> dd);
    
 	/// methods called by base method for cut element case
        virtual void add_local_vec_to_global_interface(vector_type& vec, const LocalVector& lvec,
                                                       ConstSmartPtr<DoFDistribution> dd);
    
    /// methods called by base method for the case, that TWO elements are cut by interface
    ///         --> tag: CUT_BY_2_INTERFACE
    /// REMARK: not finally tested!!
 		virtual void add_local_vec_to_global_interface_for2(vector_type& vec, const LocalVector& lvec,
                                                       ConstSmartPtr<DoFDistribution> dd);
    
    /// central method to write/map the entries for fluid-particle coupling
        virtual void add_local_vec_to_global_FT(vector_type& vec, const LocalVector& lvec, std::vector<DoFIndex> transInd,
                                            std::vector<DoFIndex> rotInd);
    
    /// sends local entries to particle DoFs in case of two particles, cutting an element
    ///         --> tag: CUT_BY_2_INTERFACE
    /// REMARK: not finally tested!!
        virtual void add_local_vec_to_global_FT_for2_StdFV(vector_type& vec, const LocalVector& lvec,
                                                    std::vector<DoFIndex> transInd1, std::vector<DoFIndex> rotInd1,
                                                    std::vector<DoFIndex> transInd2, std::vector<DoFIndex> rotInd2);
        virtual void add_local_vec_to_global_FT_for2(vector_type& vec, const LocalVector& lvec,
                                                     std::vector<DoFIndex> transInd1, std::vector<DoFIndex> rotInd1,
                                                     std::vector<DoFIndex> transInd2, std::vector<DoFIndex> rotInd2);
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    ///  B. assemble components of the rhs for the particle
    ///
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // for time dependent case
 		void add_mass_part_def(vector_type& vec, std::vector<DoFIndex> transInd,
                    std::vector<DoFIndex> rotInd, const int levIndex, const int prtIndex);
    
    // called during 'add_local_vec_to_global_interface()': calls according components for assemgling the rhs
    // (see methods below)
 		void add_rhs(vector_type& vec, std::vector<DoFIndex> transInd, std::vector<DoFIndex> rotInd,
                    const int levIndex, const int prtIndex);
 	
    // additional methods called by 'add_rhs()':
        void set_gravitational_rhs(vector_type& vec, std::vector<DoFIndex> transInd,
                        std::vector<DoFIndex> rotInd, const int levIndex, const int prtIndex);
        void add_repulsive_force_rhs(vector_type& vec, std::vector<DoFIndex> transInd,
                        std::vector<DoFIndex> rotInd, const int levIndex, const int prtIndex);
        void add_glowinski_repulsive_force_rhs(vector_type& vec, std::vector<DoFIndex> transInd,
                        std::vector<DoFIndex> rotInd, const int levIndex, const int prtIndex);
        void add_minimum_correction_force_rhs(vector_type& vec, std::vector<DoFIndex> transInd,
                        std::vector<DoFIndex> rotInd, const int levIndex, const int prtIndex);

    
   	///////////////////////////////////////////////////////////////////////////////
 	/// helper function for 'add_local_vec_to_global_FT_for2()':
    ///         --> tag: CUT_BY_2_INTERFACE
    //  REMARK: not finally tested!
    
    /// special mapping of local DoFs due to element, cut by two particles
 		size_t map_for2(size_t dof);
 
 	/// called during 'add_local_vec_to_global_FT_for2()':
 		void assemble_QuadriCorners(vector_type& vec, const LocalVector& lvec,
 				std::vector<DoFIndex> transInd1, std::vector<DoFIndex> rotInd1,
 				std::vector<DoFIndex> transInd2, std::vector<DoFIndex> rotInd2);

    
 	///////////////////////////////////////////////////////////////////////////////
    ///
 	///  C. local matrix to global
    ///
    ///////////////////////////////////////////////////////////////////////////////

    
    /// base method
		virtual void add_local_mat_to_global(matrix_type& mat, const LocalMatrix& lmat,
				ConstSmartPtr<DoFDistribution> dd);
    
    /// methods called by base method for cut element case
 		virtual void add_local_mat_to_global_interface(matrix_type& mat, const LocalMatrix& lmat,
                                                       ConstSmartPtr<DoFDistribution> dd);
    
    /// methods called by base method for the case, that TWO elements are cut by interface
    ///         --> tag: CUT_BY_2_INTERFACE
    /// REMARK: not finally tested!!
 		virtual void add_local_mat_to_global_interface_for2(matrix_type& mat, const LocalMatrix& lmat,
                                                       ConstSmartPtr<DoFDistribution> dd);
    
    /// central method to write/map the entries for fluid-particle coupling --> tag: CUT_BY_INTERFACE
 		virtual void add_local_mat_to_global_FT(matrix_type& mat, const LocalMatrix& lmat,
                std::vector<DoFIndex> transInd, std::vector<DoFIndex> rotInd);
    
    /// sends local entries to particle DoFs in case of two particles, cutting an element
    ///         --> tag: CUT_BY_2_INTERFACE
    /// REMARK: not finally tested!!
 		virtual void add_local_mat_to_global_FT_for2_StdFV(matrix_type& mat, const LocalMatrix& lmat,
 				std::vector<DoFIndex> transInd1, std::vector<DoFIndex> rotInd1,
 				std::vector<DoFIndex> transInd2, std::vector<DoFIndex> rotInd2);
 		virtual void add_local_mat_to_global_FT_for2(matrix_type& mat, const LocalMatrix& lmat,
 				std::vector<DoFIndex> transInd1, std::vector<DoFIndex> rotInd1,
 				std::vector<DoFIndex> transInd2, std::vector<DoFIndex> rotInd2);

    
    // for time dependent case
        void add_mass_part_jac(matrix_type& mat, std::vector<DoFIndex> transInd,
                           std::vector<DoFIndex> rotInd, const int levIndex, const int prtIndex);
    
    /// sets dirichlet rows for DoFs lying outside the fluid
    /// 	->  IFF element does NOT include transInd/rotInd-node! (which also lie outside fluid,
    ///         i.e. inside the particle domain
        void set_identity_mat_constraint(matrix_type& mat, const LocalMatrix& lmat,
                                             ConstSmartPtr<DoFDistribution> dd);

 	///////////////////////////////////////////////////////////////////////////////
 	/// REMARK:
 	///	During DomainDiscretization::assemble_jacobian:
 	/// calling
 	/// 		--->  m_spAssTuner->modify_LocalData(pModifyMemory, vSol, dd);
 	///////////////////////////////////////////////////////////////////////////////

	///	resizes local solution vector for the assemlby on cut elements, which
    ///     potentially have more nodes than the original element
 		virtual void modify_LocalData(LocalMatrix& locJ, LocalVector& locU,
                                      ConstSmartPtr<DoFDistribution> dd);
 		virtual void modify_LocalData(LocalVectorTimeSeries& uT, LocalMatrix& locJ,
                                      LocalVector& locU, ConstSmartPtr<DoFDistribution> dd);
    
    ///	(1) resizes local solution vector for the assemlby on cut elements, which
    ///     potentially have more nodes than the original element
    /// (2) writes the particle velocities in the designated, local particle indices
 		virtual void modify_LocalData(LocalVector& locD, LocalVector& tmpLocD, LocalVector& locU,
                                      ConstSmartPtr<DoFDistribution> dd);
 		virtual void modify_LocalData(LocalVectorTimeSeries& uT, LocalVector& locD, LocalVector& tmpLocD,
                                      LocalVector& locU, ConstSmartPtr<DoFDistribution> dd, size_t t);

 	/// method called during 'modify_LocalData()' to set local sol for CUT_BY_2_INTERFACE
 		void set_QuadriSol(LocalVector& locU, LocalVector& locD);

 	/// called during 'modify_LocalData(locU)':
    /// copies old local vector to the new local algebra due to cut element
 		void map_local_data(LocalVector& d);

 	/// called during modify_LocalData() for resizing local data on cut elements:
 		void resize_local_indices(LocalVector& locU)
            { m_spInterfaceHandlerLocal->resize_local_indices(locU); }
 		void resize_local_indices(LocalVector& locU, size_t numCo)
            { m_spInterfaceHandlerLocal->resize_local_indices(locU, numCo); }

    //  writes/copies the linear and angular velocity of the particle to data storage
    // in 'ParticleProvider:m_vvLinearVelocity' (called during 'modify_GlobalSol')
  		void set_extraSol(vector_type& vec, std::vector<DoFIndex> transInd, std::vector<DoFIndex> rotInd,
                          const int timeIndex, const int prtIndex);

	///////////////////////////////////////////////////////////////////////////////
	/// REMARK:
	///	During DomainDiscretization::assemble_jacobian:
	/// calling
	/// 		--->  m_spAssTuner->modify_GlobalSol(pModifyMemory, vSol, dd);
	///////////////////////////////////////////////////////////////////////////////

    /// the method 'modify_GlobalSol()' does not 'modify' the solution, but: the
    /// computed solution at the trans und rot DoFs gets written/stored into data of
    /// class 'ParticleProvider::m_vvLinearVelocity/m_vvAngularVelocity',
    /// since they will be overwritten during the local computation of defect
    /// (for each call of domainDisc, i.e. newton step)

		virtual void modify_GlobalSol(vector_type& uMod, const vector_type& u, ConstSmartPtr<DoFDistribution> dd);

		virtual void modify_GlobalSol(SmartPtr<VectorTimeSeries<vector_type> > vSolMod,
				ConstSmartPtr<VectorTimeSeries<vector_type> > vSol, ConstSmartPtr<DoFDistribution> dd);

	///////////////////////////////////////////////////////////////////////////////
	/// setter methods:
	///////////////////////////////////////////////////////////////////////////////

		void set_gravity(bool gravity, number gravityConst)
            { m_bGravity = gravity; m_gravityConst = gravityConst;}
        void set_repulsive_force(bool repulsive, number forceValue)
            {m_bRepulsiveForce = repulsive; m_repulsiveForce = forceValue;}
        void set_glowinski_repulsive_force(bool repulsive, number rho, number epsilon)
            {m_bGlowRepulsiveForce = repulsive; m_rho = rho; m_epsilon = epsilon;}
        void set_minimum_correction_force(bool repulsive, std::vector<MathVector<dim> > equiDist)
            {m_bMinimumCorrectionForce = repulsive; m_repulsiveDistance = equiDist;}

        void set_bUsualAss(bool UsualAss)               { m_bUsualAss = UsualAss; }
		void set_volume_comp_mode(bool bVolumeCompMode) { m_bVolumeCompExact = bVolumeCompMode;}
        void set_time_step(number dt)                   { m_dt = dt; set_time_dependent(true);}
		void set_time_dependent(bool bTimeDep)          { m_bTimeDep = bTimeDep; }
        void set_element_diameter(double diameter)      { m_meanDiameter = diameter;}
    
    // flag for output
        void set_forceLog(bool val)                     {m_bForceLog = val;}

    
    ///////////////////////////////////////////////////////////////////////////////
    /// getter methods:
    ///////////////////////////////////////////////////////////////////////////////
    
    // gets the current index of the particle, which cuts the current element;
    // stored in 'InterfaceHandlerLocalParticle'
        int get_prtIndex()                          { return this->m_spInterfaceHandlerLocal->get_prtIndex(); }
    
    // IFF two particles cut an element: gets the index of the particle, to which
    //  the node with local index 'dof' belongs
        int getPrtIndex(size_t dof)                 { return m_spInterfaceHandlerLocal->getPrtIndex(dof);}

        size_t num_particles() const                { return m_spCutElementHandler->num_particles();}
        bool UsualAss()                             { return m_bUsualAss; }
		bool is_time_dependent()                    { return m_bTimeDep;}
        bool gravitation_force()                    { return m_bGravity; }
    
	    number get_time_step()
            { if ( !is_time_dependent() )
                UG_THROW("Call for time step, BUT: not timedependent computation!\n");
                                                        return m_dt; }
    
    /// checks if grid data is updated and returns the 'levIndex' of the 'gridLevel'
    /// via access to 'CutElementHandlerBase::m_Map'
        int get_Index(const GridLevel& gridLevel, ConstSmartPtr<DoFDistribution> dd)
            { return m_spCutElementHandler->get_Index(gridLevel, dd); }

    /////////////////////////////////////////////////////////////////////////////////////
    // the mapping modus inherits the coupling categorie 
    //  --> pair combinations (from,to) of: velDoF in fluid, velDoF in particle,
    //                          pressureDoF in fluid, pressureDof on particle boundary
    //
    // => for the case 'pressureDof on particle boundary': is treated as usual DoF
    // => ONLY for the case 'velDoF in particle' mapping to the global index needs to be
    //       changed (for the row or column) and in case of rotational DoF also different
    //       values (multiplied by rotationMat) need to be choosen
    /////////////////////////////////////////////////////////////////////////////////////

        int getMappingModus(size_t fct1, size_t dof1, size_t fct2, size_t dof2);
        int getMappingModus_for2(size_t fct1, size_t dof1, size_t fct2, size_t dof2);
    

	 // access to 'LocalMatrix' data within 'ParticleBndCond' class
     // (used during 'add_local_mat_to_global_FT()' in case the boundary
     // conditions were computed by the class 'ParticleBndCond')
		number get_rotJ_ind(size_t fct1, size_t dof1, size_t fct2, size_t dof2)
			{ return m_spInterfaceHandlerLocal->get_rotJ_ind(fct1, dof1, fct2, dof2); }
		number get_rotJ_rot(size_t fct1, size_t dof1, size_t fct2, size_t dof2)
			{ return m_spInterfaceHandlerLocal->get_rotJ_rot(fct1, dof1, fct2, dof2); }

		number get_rotD(size_t fct, size_t dof)
			{ return m_spInterfaceHandlerLocal->get_rotD(fct, dof); }

    
    ///////////////////////////////////////////////////////////////////////////////
    /// Forwarding computational methods to the interface provider:
    ///////////////////////////////////////////////////////////////////////////////

    // analytical computation of the volume of a dics/sphere
  		number Volume(int levIndex, size_t prtIndex)
            { return m_spCutElementHandler->Volume(levIndex, prtIndex); }
 		number Mass(const int levIndex, const int prtIndex)
            { return m_spCutElementHandler->Mass(levIndex, prtIndex, m_spInterfaceHandlerLocal->get_density_fluid()); }
 		number Mass(const int levIndex, const int prtIndex, const number volume)
            { return m_spCutElementHandler->Mass(levIndex, prtIndex, volume, m_spInterfaceHandlerLocal->get_density_fluid()); }
        number MomOfInertia(const int levIndex, const int prtIndex)
            { return m_spCutElementHandler->MomOfInertia(levIndex, prtIndex, m_spInterfaceHandlerLocal->get_density_fluid()); }
 		number MomOfInertia(const int levIndex, const int prtIndex, const number volume)
            { return m_spCutElementHandler->MomOfInertia(levIndex, prtIndex, volume, m_spInterfaceHandlerLocal->get_density_fluid()); }

    // if m_bVolumeCompExact = false: the volume will computed based on the volume of the parts
    // covered by in the particle, instead of using the analytical formular for a disc/sphere
        number compute_volume(int levIndex, size_t prtIndex);
        void reset_volume(){m_volume = 0.0;}
    
    ///////////////////////////////////////////////////////////////////////////////
    /// class member
    ///////////////////////////////////////////////////////////////////////////////

	private:
	// member from base class
 		SmartPtr<InterfaceHandlerLocalParticle<dim> > m_spInterfaceHandlerLocal;
 	// new member
 		SmartPtr<CutElementHandler_FlatTop<dim> > m_spCutElementHandler;


     // gravityConst for call during 'set_gravitational_rhs()'
     	number m_gravityConst;

    // boolian to add gravity force to global defect during 'add_local_vec_to_global_interface()'
     	bool m_bGravity;                            // default = false;

    /// handles the call of'set_gravity()' exactly ONCE during 'modify_LocalData()'
        std::vector<bool> m_bFlagGravity; 			// default = false
    
    /// handles the call of'set_gravity()' exactly ONCE during 'modify_LocalData()'
        std::vector<bool> m_bFlagInertial; 			// default = false

    
    // boolian to add repulsive force in 'add_rhs()'
        bool m_bRepulsiveForce;
        bool m_bGlowRepulsiveForce;
        bool m_bMinimumCorrectionForce;
        bool m_bForceLog;
        number m_rho;
        number m_epsilon;
        number m_repulsiveForce;
        std::vector<MathVector<dim> > m_repulsiveDistance;

    // used within 'set_gravitational_rhs()' for computation of rhs
        number m_dt;                                // default = 0.0;
        bool m_bTimeDep;
    
    // used within 'add_repulsive_force' for computation of rhs
        double m_meanDiameter;
    
    // used for the iterative computation of the volume, if not computet analytically
    	number m_volume;
    	bool m_bVolumeCompExact;  // default = true; if false, the volume will computed based
                                  // on the volume of the parts covered by in the particle,
                                  // instead of using the analytical formular for a disc/sphere
    
     	bool m_bUsualAss;         // default = false;


};


} // end namespace NavierStokes
} // end namespace ug


#include "loc_to_glob_mapper_particle_impl.h"



#endif /* PARTICLE_MAPPER_H_ */
