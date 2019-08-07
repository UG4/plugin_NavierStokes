/*
 * diffusion_interface_handler_local.h
 *
 *  Created on: 19.01.2015
 *      Author: suze
 */

#ifndef INTERFACE_HANDLER_2PF_H_
#define INTERFACE_HANDLER_2PF_H_

#include "lib_disc/spatial_disc/immersed_util/interface_handler/interface_handler_two_sided_cut/interface_handler_diffusion.h"

namespace ug{
namespace NavierStokes{
    

template <int TWorldDim>
class InterfaceHandlerLocal2PF : public InterfaceHandlerLocalDiffusion<TWorldDim>
{

	public:
	///	World dimension
		static const int dim = TWorldDim;

	/// used boundary face type
		typedef typename DimFV1CutGeometry<dim, dim, InterfaceHandlerLocalDiffusion<dim> >::BF interfaceBF;

        typedef typename domain_traits<dim>::grid_base_object grid_base_object;

  		InterfaceHandlerLocal2PF(SmartPtr<DiffusionInterfaceProvider<dim> > interfaceProvider,
 									   SmartPtr<CutElementHandlerImmersed<dim> > cutElementHandler,
 									   number fluidDensity, number fluidKinVisc);

    	virtual ~InterfaceHandlerLocal2PF()	{}

	///////////////////////////////////////////////////////////////
    /// redefine base calss methods:
	///////////////////////////////////////////////////////////////
    //	bool update_elem(GridObject* elem, const MathVector<TWorldDim>* vCornerCoords, int interfaceOrientation);

    //	int CollectCorners_FlatTop_2d(GridObject* elem);

    //    size_t get_vertex_index(Vertex* vrt, GridObject* elem);


   ///////////////////////////////////////////////////////////////
   /// new methods:
   ///////////////////////////////////////////////////////////////

		size_t num_particles() const { return this->m_spInterfaceProvider->num_particles();}

	/// get solution values
		MathVector<dim> get_solution(size_t prtIndex, size_t timeSeriesInd)
			{ return this->m_spInterfaceProvider->get_solution(prtIndex, timeSeriesInd); }


	    number get_density() { return this->m_spInterfaceProvider->get_density(0); }
	    number get_density_fluid() { return m_fluidDensity; }
	    number get_kinVisc_fluid() { return m_fluidKinVisc; }

		const LocalIndices& get_local_indices() const { return m_ind; }

		void set_jacobian_tri(const LocalMatrix locJ) { m_locJ_tri = locJ; }
		void set_jacobian_quad(const LocalMatrix locJ){ m_locJ_quad = locJ; }

		void set_defect_tri(const LocalVector locD) { m_locD_tri = locD; }
		void set_defect_quad(const LocalVector locD){ m_locD_quad = locD; }

		void reset_defect_on_interface(LocalVector& locD, const size_t size);

		void reset_jacobian_on_interface(LocalMatrix& locJ, const size_t size);


		void set_solution_tri(const LocalVector locU) { m_locU_tri = locU; }
		void set_solution_quad(const LocalVector locU){ m_locU_quad = locU; }

		void set_DoF_tag_tri(const bool bFactor2_for_DoFIndex)
		{ m_shift_DoFIndex_tri = bFactor2_for_DoFIndex; return; }
		void set_DoF_tag_quad(const bool bFactor2_for_DoFIndex)
		{ m_shift_DoFIndex_quad = bFactor2_for_DoFIndex; return; }

		const size_t get_index_shift_tri() const { return m_shift_DoFIndex_tri; }
		const size_t get_index_shift_quad() const { return m_shift_DoFIndex_quad; }

		LocalMatrix& get_local_jacobian_tri()  { return m_locJ_tri; }
		LocalMatrix& get_local_jacobian_quad() { return m_locJ_quad; }

		LocalVector& get_local_defect_tri()  { return m_locD_tri; }
		LocalVector& get_local_defect_quad() { return m_locD_quad; }

		LocalVector& get_local_solution_tri()  { return m_locU_tri; }
		LocalVector& get_local_solution_quad() { return m_locU_quad; }

		void set_local_sol(LocalVector& solU, const size_t size, const LocalVector& lvec, const int orientation);

		LocalVector set_jump_values(LocalIndices ind, const size_t size);
		LocalVector set_jump_grad_values(LocalIndices ind, const size_t size);
		LocalVector set_source(const std::vector<double> sourceIm, LocalIndices ind, const size_t size, const bool bElementIsCut)
		{LocalVector dummy; return dummy;}
		LocalVector set_source(LocalIndices ind, const size_t size, const bool bElementIsCut);

		double get_jump_value_const(const MathVector<dim> position);
		double get_jump_value(const MathVector<dim> position);
		double get_jump_value_ex3(const MathVector<dim> position);
		double get_jump_grad_value(const MathVector<dim> position);
		double get_jump_grad_value_ex3(const MathVector<dim> position);

		double get_source(const MathVector<dim> position);
		double get_source_kappa(const MathVector<dim> position);

	// writes solution of global vector vec into this->m_verticesValue-array:
		void write_solution(const std::vector<double > verticesValues);


 		void resize_local_data(LocalVector locU)
		{
 			LocalIndices ind = locU.get_indices();

            for (size_t fct = 0; fct < (dim+1); ++fct)
            {
            // resize for cut triangle
                ind.resize_dof(fct, 3);
                m_locU_tri.resize(ind);
                m_locD_tri.resize(ind);
                m_locJ_tri.resize(ind);
            }
            for (size_t fct = 0; fct < (dim+1); ++fct)
            {
            // resize for cut quadrilateral
                ind.resize_dof(fct, 4);
                m_locU_quad.resize(ind);
                m_locD_quad.resize(ind);
                m_locJ_quad.resize(ind);
            }
            
            UG_LOG("check LocalIndices ind:" << ind << "\n");

			return;
		}

	private:

    /// fuid parameter imported by Constructor()
     	number m_fluidDensity;
     	number m_fluidKinVisc;

	/// size of local algebra for flat top element: 'm_numFct' x 'm_numCo'
		size_t m_numFct;
	/// number of corners of flat top element
		size_t m_numCo;

	/// new local algebra for resized flat top element
		LocalIndices m_ind;

	// local data for assembling:
		LocalMatrix m_locJ_tri;
		LocalMatrix m_locJ_quad;
		LocalVector m_locD_tri;
		LocalVector m_locD_quad;
		LocalVector m_locU_tri;
		LocalVector m_locU_quad;

	// scale factor for access to DoFIndex on triangle or quadri as cut element
	// --> for call during 'add_local_def/jac_to_global_interface()':
		bool m_shift_DoFIndex_tri;
		bool m_shift_DoFIndex_quad;

};

}// namespace NavierStokes
} // end ug namespace

#include "interface_handler_2pf_impl.h"


#endif /* INTERFACE_HANDLER_2PF_H_ */
