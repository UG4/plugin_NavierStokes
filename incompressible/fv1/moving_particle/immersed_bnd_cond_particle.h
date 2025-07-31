/*
 * particle_bnd_cond.h
 *
 *  Created on: 30.01.2015
 *      Author: suze
 */

#ifndef PARTICLE_BND_COND_H_
#define PARTICLE_BND_COND_H_

#ifdef UG_PARALLEL
 	#include "lib_grid/parallelization/load_balancer_util.h"
#endif

#include "../../incompressible_navier_stokes_base.h"


#include "lib_disc/spatial_disc/immersed_util/immersed_interface_base.h"

namespace ug{
namespace NavierStokes{

template <typename TDomain>
class ParticleBndCond : public IInterfaceBndCond<TDomain>
{
	public:
	///	World dimension
		static const int dim = TDomain::dim;

	/// abbreviation for pressure
		static const size_t _P_ = dim;

	/// used boundary face type
		typedef typename DimFV1FTGeometry<dim, dim, InterfaceHandlerLocalParticle<dim> >::BF interfaceBF;

	///	call base class constructor

        ParticleBndCond(SmartPtr<NavierStokesFV1_cutElem<TDomain> > spMaster,
					SmartPtr<InterfaceHandlerLocalParticle<dim> > localHandler);

    
    /// destructor
        virtual ~ParticleBndCond(){};
    
    ///////////////////////////////////////////////////////////////////////////////
    ///
    ///  base class methods
    ///
    ///////////////////////////////////////////////////////////////////////////////

	public:
	///	type of trial space for each function used
		void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);
 
	///	prepares the element loop
		template <typename TElem, typename TFVGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si){
		}

	///	prepares the element for evaluation
		template <typename TElem, typename TFVGeom>
		void prep_elem(const LocalVector& u, GridObject* elem,
				const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[]);

	///	finishes the element loop
		template <typename TElem, typename TFVGeom>
		void fsh_elem_loop(){ }

	///	adds the stiffness part to the local jacobian
		template <typename TElem, typename TFVGeom>
		void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem,
				const MathVector<dim> vCornerCoords[]);

    /// helper method called by 'add_jac_A_elem()'
		void add_jac_A_elem_Quadri_for2(LocalMatrix& J, const LocalVector locU);

	///	adds the stiffness part to the local defect
		template <typename TElem, typename TFVGeom>
		void add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem,
				const MathVector<dim> vCornerCoords[]);
    
    /// helper method called by 'add_def_A_elem()'
		void add_def_A_elem_Quadri_for2(LocalVector& locD, const LocalVector locU);
		void add_quadri_to_defect(LocalVector& d, const LocalVector& quadriD);

		template <typename TElem, typename TFVGeom>
		void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem,
				const MathVector<dim> vCornerCoords[]);

		template <typename TElem, typename TFVGeom>
		void add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem,
				const MathVector<dim> vCornerCoords[]);

		template <typename TElem, typename TFVGeom>
		void add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[]);


    // computes contributions to local stiffness matrix due to the stresses on the particle boundary
		void diffusive_flux_Jac(const size_t ip, const interfaceBF& bf, LocalMatrix& J, const LocalVector& u);
		void diffusive_flux_Jac_for2(const size_t ip, const interfaceBF& bf, LocalMatrix& J, const LocalVector& u);
		void diffusive_flux_Jac_rot(const size_t ip, const interfaceBF& bf, LocalMatrix& J, const LocalVector& u, number importDensity);

    // computes contributions to local defect due to the stresses on the particle boundary
		void diffusive_flux_defect(const size_t ip, const interfaceBF& bf, LocalVector& d, const LocalVector& u);
		void diffusive_flux_defect_rot(const size_t ip, const interfaceBF& bf, LocalVector& d, const LocalVector& u);

    // reset entries for the momentum equations for interface-corner:
    // set entries in local data to zero BEFORE adding the computed contributions from the boundary condition
		void remove_equations(LocalVector& d, std::vector<size_t> vFctID, std::vector<size_t> vDofID);
		void remove_equations(LocalMatrix& J, std::vector<size_t> vFctID, std::vector<size_t> vDofID);

    // gets the current index of the particle, which cuts the current element;
    // stored in 'InterfaceHandlerLocalParticle'
		int get_prtIndex() { return m_spInterfaceHandlerLocal->get_prtIndex(); }

	    number get_density(int prtIndex) { return m_spInterfaceHandlerLocal->get_density(prtIndex); }
	    number get_density_fluid()       { return m_spInterfaceHandlerLocal->get_density_fluid(); }
	    number get_kinVisc_fluid()       { return m_spInterfaceHandlerLocal->get_kinVisc_fluid(); }

    // method called by 'ParticleMapper' class for acces to boundary conditions, assembled by this class
	    void copy_local_couplings_jac()
            { m_spInterfaceHandlerLocal->set_local_couplings_jac(rotJ_ind, rotJ_rot); }
	    void copy_local_couplings_def()
            { m_spInterfaceHandlerLocal->set_local_couplings_def(rotD); }

	    void write_QuadriSol(const LocalVector origU);

    /// remaps entries for Quadri+Tri-combination of CUT_BY_2_INTERFACE element
	    size_t remap_for2(size_t dof);

    
    ///////////////////////////////////////////////////////////////////////////////
    /// class member
    ///////////////////////////////////////////////////////////////////////////////

	protected:
	// master element disc
		SmartPtr<NavierStokesFV1_cutElem<TDomain> > m_spMaster;

	// member from base class
		SmartPtr<InterfaceHandlerLocalParticle<dim> > m_spInterfaceHandlerLocal;

	// local data for assembling:
		LocalMatrix rotJ_ind;
		LocalMatrix rotJ_rot;
		LocalVector rotD;


	protected:
		void register_all(bool bHang);
		template<typename TElem, typename TFVGeom>
		void register_func();

};



} // end namespace NavierStokes
} // end namespace ug

#include "immersed_bnd_cond_particle_impl.h"



#endif /* PARTICLE_BND_COND_H_ */
