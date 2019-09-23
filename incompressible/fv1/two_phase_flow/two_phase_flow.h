/*
 * diffusion_interface.h
 *
 *  Created on: 24.08.2017
 *      Author: suze
 */

#ifndef PF2_INTERFACE_H_
#define PF2_INTERFACE_H_


#ifdef UG_PARALLEL
 	#include "lib_grid/parallelization/load_balancer_util.h"
#endif

#include "../navier_stokes_fv1.h"
#include "../../../navier_stokes_base.h"
#include "interface_handler_2pf.h"
#include "loc_to_glob_mapper_2pf.h"
#include "lib_disc/spatial_disc/immersed_util/immersed_interface_base.h"

namespace ug{
namespace NavierStokes{



template <	typename TDomain, typename TAlgebra>
class MovingInterface2PF
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

 		MovingInterface2PF(
 					   SmartPtr<IAssemble<TAlgebra> > ass,
 					   SmartPtr<NavierStokesFV1<TDomain> > spMaster,
 					   SmartPtr<DiffusionInterfaceProvider<dim> > interfaceProvider,
 					   SmartPtr<CutElementHandler_TwoSided<dim> > cutElementHandler,
 					   number fluidDensity1, number fluidDensity2);

		void set_StdFV_assembling(bool bValue) { m_spInterfaceHandlerLocal->set_StdFV_assembling(bValue);}
     	bool StdFV_assembling() { return m_spInterfaceHandlerLocal->StdFV_assembling(); }


 	// destructor
		~MovingInterface2PF(){};

	/// called via .lua:
		void initialize_threshold(TDomain& domain, const int baseLevel, const int topLevel);
	    void set_threshold(size_t level, const number threshold)
	    { m_spCutElementHandler->set_threshold(level, threshold); }

	//////////////////////////////////////////////////////////////////////////////////
	/// Info - 'initialize_interface()':
	///
	/// computes vertices on intersection of cut element edges and interface:
	/// for 2d: instead of computing intersections: count number of cut elements!
	/// 	-> #cutElements == #m_vertices
	/// 	-> called during init()
	//////////////////////////////////////////////////////////////////////////////////

	    const size_t initialize_interface(vector_type& u, ConstSmartPtr<DoFDistribution> dd);


		number MeanElementDiameter(TDomain& domain, int level);

		//////////////////////////////////////////////////////////////////////////////////
		// ---> .lua:  update(u, deltaT, u:grid_level()): update global indices (transInd...)
		// 				=> A. copy_solution(topLev)
		// 				   B. update(baseLev-topLev)
		// 				   C. update_solution(topLev)
 		//////////////////////////////////////////////////////////////////////////////////
		// write solution to nodes outside fluid with particle velocities
		//		--> call method vie .lua BEFORE 'solTimeSeries:push_discard_oldest(oldestSol, time)':
		//			=> in case that outside nodes are inside AFTER update_prtCoords: NO solution defined here!


		/// call of the method via lua to set the real velocity values within the particle domain
		void adjust_global_solution(vector_type& u, const int topLevel);
		void fill_particle_solution(vector_type& u, const int topLevel, const number time);
		void update(vector_type& u, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, const int baseLevel, const int topLevel, const number time)
		{
			int topLev = spApproxSpace->num_levels()-1;
			if ( topLev != topLevel )
				UG_THROW("MovingParticle::update: parameter 'topLevel' = " << topLevel << " != "
								 << topLev << "current top leven! \n");

		// fill particle nodes with their real solution
 			fill_particle_solution(u, topLevel, time);

		// update data: bool_marker
			ConstSmartPtr<DoFDistribution> dd = spApproxSpace->dof_distribution(GridLevel(topLevel, GridLevel::LEVEL));
			m_spCutElementHandler->template init<TDomain>(dd, baseLevel, topLevel);

  		}

		
		void set_analytic_solution(vector_type& u, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, SmartPtr<MultiGrid> mg, const int topLevel);
		void adjust_for_error(vector_type& u, vector_type& uCopy, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, SmartPtr<MultiGrid> mg, const int topLevel);

		double compute_solution_value(const MathVector<dim>& vrtPos);

		number get_L2Error()
		{ return m_spInterfaceHandlerLocal->get_L2Error(); }

 		void init(vector_type& u, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, const int baseLevel, const int topLevel, bool bScaleDoFs)
		{
   			m_spApproxSpace = spApproxSpace;

			ConstSmartPtr<DoFDistribution> dd = spApproxSpace->dof_distribution(GridLevel(topLevel, GridLevel::LEVEL));


			size_t numDoFs = u.size();
			const size_t num_interfaceDoFs = initialize_interface(u, dd);
			const size_t num_newDoFs = 4*num_interfaceDoFs;

			m_spInterfaceMapper->set_numDoFs(numDoFs);
			m_spInterfaceMapper->set_numNewDoFs(num_interfaceDoFs);

			UG_LOG("________________ numDoFs = " << numDoFs << "\n");
			UG_LOG("________________ num_interfaceDoFs = " << num_interfaceDoFs << "\n");
			UG_LOG("________________ num_newDoFs = " << num_newDoFs << "\n");

		// values for new DoFs are set to 0.0 by the 'resize()'-method (see vector.h):
			if ( bScaleDoFs )
				u.resize(numDoFs + 2*num_newDoFs);
			else
				u.resize(numDoFs + num_newDoFs);

			UG_LOG("AGAIN: in init(): numALLDoFs = " <<  u.size() << "\n");

			m_spInterfaceMapper->set_bScaleDoFs(bScaleDoFs);
			m_spInterfaceHandlerLocal->set_bScaleDoFs(bScaleDoFs);

			m_spInterfaceHandlerLocal->L2Error_init();

		// not necessary anymore: only local evaluations within diffusion problem!
			//m_spCutElementHandler->template init_marker<TDomain>(dd, baseLevel, topLevel);
		}

	   /// checks if grid data is updated and returns 'levIndex'-pair for 'gridLevel' in 'm_Map'
	    int get_Index(const GridLevel& gridLevel)
	    {
	    	ConstSmartPtr<DoFDistribution> dd = m_spApproxSpace->dof_distribution(gridLevel);

 	 		const int levIndex = m_spCutElementHandler->get_Index(gridLevel, dd);

	    	return levIndex;
	    }
    
	    //ToDo: method needed?
	    void update_interface( const int topLevel, number deltaT);

		bool is_time_dependent() { return m_spInterfaceHandlerLocal->is_time_dependent();}

 		void interpolate_point(ConstSmartPtr<DoFDistribution> dd, const vector_type& u,
							   const MathVector<dim>& evalPos, MathVector<dim+1>& interpolation);
 
		void print_deltaP(const vector_type& u, const int topLevel);
		void print_pressure(const vector_type& u, const int topLevel);
		void print_pressure_nodal(const vector_type& u, const int topLevel);

	/// writing data to file; called via .lua
		void print_velocity(const vector_type& u, const int topLevel, number time, const char* filename);
 
	private:
	///	current ApproxSpace
		SmartPtr<ApproximationSpace<TDomain> > m_spApproxSpace;

		SmartPtr<DiffusionInterfaceProvider<dim> > m_spInterfaceProvider;
		SmartPtr<CutElementHandler_TwoSided<dim> > m_spCutElementHandler;
		SmartPtr<InterfaceHandlerLocal2PF<dim> > m_spInterfaceHandlerLocal;

		SmartPtr<InterfaceMapper2PF<TDomain, TAlgebra> > m_spInterfaceMapper;

};

} // end namespace NavierStokes
} // end namespace ug


#include "two_phase_flow_impl.h"

#endif /* PF2_INTERFACE_H_ */




