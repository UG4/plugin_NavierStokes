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
			   m_spParticleHandlerGlobal(m_spInterfaceHandlerLocal->get_cutElementHandler()),
			   m_gravityConst(0.0), m_bGravity(true),
               m_dt(0.0), m_bTimeDep(false),
               m_volume(0.0), m_bVolumeCompExact(true),
               m_bUsualAss(false), m_meanDiameter(0.0),
               m_bRepulsiveForce(false),
               m_bGlowRepulsiveForce(false),
               m_bMinimumCorrectionForce(false),
               m_repulsiveDistance(0.0), m_rho(0.0),
               m_repulsiveForce(0.0), m_epsilon(0.0),
               m_bForceLog(false)

		{
 		// init boolian arrays for all registered particles
  			size_t numPrt = num_particles();
   			m_bFlagGravity.resize(numPrt, false);
  			m_bFlagInertial.resize(numPrt, false);
//  		particleValues.data.resize(numPrt);
		};

 	///	destructor
 		~ParticleMapper() {};

	///////////////////////////////////////////////////////////////////////////////
 	///
 	///  base class methods and helper methods called by them
 	///
 	///////////////////////////////////////////////////////////////////////////////

 	///////////////////////////////////////////////////////////////////////////////
  	///  A. local vector to global

	///	send local entries to global vector
		virtual void add_local_vec_to_global(vector_type& vec, const LocalVector& lvec, ConstSmartPtr<DoFDistribution> dd);
 		virtual void add_local_vec_to_global_interface(vector_type& vec, const LocalVector& lvec, ConstSmartPtr<DoFDistribution> dd);
 		virtual void add_local_vec_to_global_interface_for2(vector_type& vec, const LocalVector& lvec, ConstSmartPtr<DoFDistribution> dd);

		void adjust_mat_global(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd){};

 		void add_mass_part_def(vector_type& vec, std::vector<DoFIndex> transInd, std::vector<DoFIndex> rotInd, const int levIndex, const int prtIndex);
 		void add_rhs(vector_type& vec, std::vector<DoFIndex> transInd, std::vector<DoFIndex> rotInd, const int levIndex, const int prtIndex);
 		void set_gravitational_rhs(vector_type& vec, std::vector<DoFIndex> transInd, std::vector<DoFIndex> rotInd, const int levIndex, const int prtIndex);
        void add_repulsive_force_rhs(vector_type& vec, std::vector<DoFIndex> transInd, std::vector<DoFIndex> rotInd, const int levIndex, const int prtIndex);
        void add_glowinski_repulsive_force_rhs(vector_type& vec, std::vector<DoFIndex> transInd, std::vector<DoFIndex> rotInd, const int levIndex, const int prtIndex);
        void add_minimum_correction_force_rhs(vector_type& vec, std::vector<DoFIndex> transInd, std::vector<DoFIndex> rotInd, const int levIndex, const int prtIndex);

		void map_all_local_vec_to_transDoF(vector_type& vec, const LocalVector& lvec, ConstSmartPtr<DoFDistribution> dd);
 		virtual void add_local_vec_to_global_FT(vector_type& vec, const LocalVector& lvec, std::vector<DoFIndex> transInd, std::vector<DoFIndex> rotInd);
 		virtual void add_local_vec_to_global_FT_for2_StdFV(vector_type& vec, const LocalVector& lvec,
 				std::vector<DoFIndex> transInd1, std::vector<DoFIndex> rotInd1,
 				std::vector<DoFIndex> transInd2, std::vector<DoFIndex> rotInd2);
 		virtual void add_local_vec_to_global_FT_for2(vector_type& vec, const LocalVector& lvec,
 				std::vector<DoFIndex> transInd1, std::vector<DoFIndex> rotInd1,
 				std::vector<DoFIndex> transInd2, std::vector<DoFIndex> rotInd2);
 	/// helper function for 'add_local_vec_to_global_FT_for2()':
 		size_t map_for2(size_t dof);

  	/// called during 'add_local_vec_to_global_FT_for2()':
 		void assemble_fluid_nodes(vector_type& vec, const LocalVector& lvec);

 	/// called during 'add_local_vec_to_global_FT_for2()':
 		void assemble_QuadriCorners(vector_type& vec, const LocalVector& lvec,
 				std::vector<DoFIndex> transInd1, std::vector<DoFIndex> rotInd1,
 				std::vector<DoFIndex> transInd2, std::vector<DoFIndex> rotInd2);

 	///////////////////////////////////////////////////////////////////////////////
 	///  B. local matrix to global

	///	send local entries to global rhs
		virtual void add_local_mat_to_global(matrix_type& mat, const LocalMatrix& lmat,
				ConstSmartPtr<DoFDistribution> dd);
 		virtual void add_local_mat_to_global_interface(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd);
 		virtual void add_local_mat_to_global_interface_for2(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd);

 		void add_mass_part_jac(matrix_type& mat, std::vector<DoFIndex> transInd, std::vector<DoFIndex> rotInd, const int levIndex, const int prtIndex);
  		void map_all_local_mat_to_transDoF(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd);
 		virtual void add_local_mat_to_global_FT(matrix_type& mat, const LocalMatrix& lmat,
						std::vector<DoFIndex> transInd, std::vector<DoFIndex> rotInd);
 		virtual void add_local_mat_to_global_FT_for2_StdFV(matrix_type& mat, const LocalMatrix& lmat,
 				std::vector<DoFIndex> transInd1, std::vector<DoFIndex> rotInd1,
 				std::vector<DoFIndex> transInd2, std::vector<DoFIndex> rotInd2);
 		virtual void add_local_mat_to_global_FT_for2(matrix_type& mat, const LocalMatrix& lmat,
 				std::vector<DoFIndex> transInd1, std::vector<DoFIndex> rotInd1,
 				std::vector<DoFIndex> transInd2, std::vector<DoFIndex> rotInd2);

 		/// instead of calling base class method 'set_identity_mat()':
 	 	/// 	->  IFF element does NOT include transInd/rotInd-node!
 	 		void set_identity_mat_constraint(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd);

 	///////////////////////////////////////////////////////////////////////////////
 	/// REMARK:
 	///	During DomainDiscretization::assemble_jacobian:
 	/// calling
 	/// 		--->  m_spAssTuner->modify_LocalData(pModifyMemory, vSol, dd);
 	///////////////////////////////////////////////////////////////////////////////

	///	modifies local solution vector for adapted defect computation
 		virtual void modify_LocalData(LocalMatrix& locJ, LocalVector& locU, ConstSmartPtr<DoFDistribution> dd);
 		virtual void modify_LocalData(LocalVectorTimeSeries& uT, LocalMatrix& locJ, LocalVector& locU, ConstSmartPtr<DoFDistribution> dd);

 		virtual void modify_LocalData(LocalVector& locD, LocalVector& tmpLocD, LocalVector& locU, ConstSmartPtr<DoFDistribution> dd);
 		virtual void modify_LocalData(LocalVectorTimeSeries& uT, LocalVector& locD, LocalVector& tmpLocD, LocalVector& locU,
 				ConstSmartPtr<DoFDistribution> dd, size_t t);

 	/// method called during 'modify_LocalData()' to set local sol for CUT_BY_2_INTERFACE
 		void set_QuadriSol(LocalVector& locU, LocalVector& locD);

 	/// called during 'modify_LocalData(locU)'
 		void map_local_data(LocalVector& d);

 	/// called during modify_LocalData() for resizing:
 		void resize_local_indices(LocalVector& locU)
		{ m_spInterfaceHandlerLocal->resize_local_indices(locU); }
 		void resize_local_indices(LocalVector& locU, size_t numCo)
		{ m_spInterfaceHandlerLocal->resize_local_indices(locU, numCo); }

	/// calls m_spParticleHandlerGlobal->set_extraSolTrans/Rot
	/// called during 'modify_GlobalSol'
  		void set_extraSol(vector_type& vec, std::vector<DoFIndex> transInd, std::vector<DoFIndex> rotInd, const int timeIndex, const int prtIndex);

	///////////////////////////////////////////////////////////////////////////////
	/// REMARK:
	///	During DomainDiscretization::assemble_jacobian:
	/// calling
	/// 		--->  m_spAssTuner->modify_GlobalSol(pModifyMemory, vSol, dd);
	/// instead of calling
	///			---->  m_vConstraint[i]->modify_solution(pModifyMemory, vSol, dd);
	///////////////////////////////////////////////////////////////////////////////

	///	modifies local solution vector for adapted defect computation
		virtual void modify_GlobalSol(vector_type& uMod, const vector_type& u, ConstSmartPtr<DoFDistribution> dd);

		virtual void modify_GlobalSol(SmartPtr<VectorTimeSeries<vector_type> > vSolMod,
				ConstSmartPtr<VectorTimeSeries<vector_type> > vSol, ConstSmartPtr<DoFDistribution> dd);

	///////////////////////////////////////////////////////////////////////////////
	/// Further helper methods:
	///////////////////////////////////////////////////////////////////////////////

		int getPrtIndex(size_t dof){ return m_spInterfaceHandlerLocal->getPrtIndex(dof);}
		int getMappingModus(size_t fct1, size_t dof1, size_t fct2, size_t dof2);
		int getMappingModus_for2(size_t fct1, size_t dof1, size_t fct2, size_t dof2);
		int get_prtIndex() { return this->m_spInterfaceHandlerLocal->get_prtIndex(); }
		size_t num_particles() const { return m_spParticleHandlerGlobal->num_particles();}
		bool UsualAss() { return m_bUsualAss; }

		void set_bUsualAss(bool UsualAss) { m_bUsualAss = UsualAss; }

		void set_gravity(bool gravity, number gravityConst) { m_bGravity = gravity; m_gravityConst = gravityConst;}
        void set_repulsive_force(bool repulsive, number forceValue) {m_bRepulsiveForce = repulsive; m_repulsiveForce = forceValue;}
        void set_glowinski_repulsive_force(bool repulsive, number rho, number epsilon) {m_bGlowRepulsiveForce = repulsive; m_rho = rho; m_epsilon = epsilon;}
        void set_minimum_correction_force(bool repulsive, number equiDist) {m_bMinimumCorrectionForce = repulsive; m_repulsiveDistance = equiDist;}

		void set_volume_comp_mode(bool bVolumeCompMode) { m_bVolumeCompExact = bVolumeCompMode;}

		bool gravitation_force() { return m_bGravity; }

        void set_time_step(number dt) { m_dt = dt; set_time_dependent(true);}

		void set_time_dependent(bool bTimeDep) { m_bTimeDep = bTimeDep; }
		bool is_time_dependent() { return m_bTimeDep;}

        void set_element_diameter(double diameter){m_meanDiameter = diameter;UG_LOG("Mean element diameter is " << m_meanDiameter);}

	    number get_time_step()
	    { if ( !is_time_dependent() ) UG_THROW("Call for time step, BUT: not timedependent computation!\n");
          return m_dt; }

	 // access methods for '..._FT()':
		number get_rotJ_ind(size_t fct1, size_t dof1, size_t fct2, size_t dof2)
			{ return m_spInterfaceHandlerLocal->get_rotJ_ind(fct1, dof1, fct2, dof2); }
		number get_rotJ_rot(size_t fct1, size_t dof1, size_t fct2, size_t dof2)
			{ return m_spInterfaceHandlerLocal->get_rotJ_rot(fct1, dof1, fct2, dof2); }

		number get_rotD(size_t fct, size_t dof)
			{ return m_spInterfaceHandlerLocal->get_rotD(fct, dof); }

	    int get_Index(const GridLevel& gridLevel, ConstSmartPtr<DoFDistribution> dd)
	    { return m_spParticleHandlerGlobal->get_Index(gridLevel, dd); }

  		void reset_volume(){m_volume = 0.0;}

  		number Volume(int levIndex, size_t prtIndex)
            { return m_spParticleHandlerGlobal->Volume(levIndex, prtIndex); }
 		number Mass(const int levIndex, const int prtIndex)
            { return m_spParticleHandlerGlobal->Mass(levIndex, prtIndex, m_spInterfaceHandlerLocal->get_density_fluid()); }
 		number Mass(const int levIndex, const int prtIndex, const number volume)
            { return m_spParticleHandlerGlobal->Mass(levIndex, prtIndex, volume, m_spInterfaceHandlerLocal->get_density_fluid()); }
        number MomOfInertia(const int levIndex, const int prtIndex)
            { return m_spParticleHandlerGlobal->MomOfInertia(levIndex, prtIndex, m_spInterfaceHandlerLocal->get_density_fluid()); }
 		number MomOfInertia(const int levIndex, const int prtIndex, const number volume)
            { return m_spParticleHandlerGlobal->MomOfInertia(levIndex, prtIndex, volume, m_spInterfaceHandlerLocal->get_density_fluid()); }

 		number compute_volume(int levIndex, size_t prtIndex);

        void set_forceLog(bool val) {
            m_bForceLog = val;
        }

	private:
	// member from base class
 		SmartPtr<InterfaceHandlerLocalParticle<dim> > m_spInterfaceHandlerLocal;
 	// new member
 		SmartPtr<CutElementHandlerFlatTop<dim> > m_spParticleHandlerGlobal;


     // gravityConst for call during 'set_gravitational_rhs()'
     	number m_gravityConst;

    // boolian to add gravity force to global defect during 'add_local_vec_to_global_interface()'
     	bool m_bGravity;	// default = false;

    // boolian to add repulsive force in 'add_rhs()'
        bool m_bRepulsiveForce;
        bool m_bGlowRepulsiveForce;
        bool m_bMinimumCorrectionForce;
        number m_rho;
        number m_epsilon;
        number m_repulsiveForce;
        number m_repulsiveDistance;

    // used within 'set_gravitational_rhs()' for computation of rhs
        number m_dt;	// default = 0.0;
        bool m_bTimeDep;
    
    // used within 'add_repulsive_force' for computation of rhs
        double m_meanDiameter;
    
    
    /// handles the call of'set_gravity()' exactly ONCE during 'modify_LocalData()'
     	std::vector<bool> m_bFlagGravity; 			// default = false

    /// handles the call of'set_gravity()' exactly ONCE during 'modify_LocalData()'
    	std::vector<bool> m_bFlagInertial; 			// default = false

    	number m_volume;
    	bool m_bVolumeCompExact;
     	bool m_bUsualAss;	// default = false;

        bool m_bForceLog;

};


} // end namespace NavierStokes
} // end namespace ug


#include "loc_to_glob_mapper_particle_impl.h"
#include "loc_to_glob_mapper_particle_tools.h"



#endif /* PARTICLE_MAPPER_H_ */
