/*
 * navier_stokes_fe.h
 *
 *  Created on: 20.09.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__PLUGINS__NAVIER_STOKES__NAVIER_STOKES_FE__
#define __H__UG__PLUGINS__NAVIER_STOKES__NAVIER_STOKES_FE__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"

#include "../navier_stokes_base.h"

namespace ug{
namespace NavierStokes{

/// \ingroup lib_disc_elem_disc
/// @{

/// Finite Element Discretization for the incompressible Navier-Stokes Equation
template<	typename TDomain>
class NavierStokesFE
	: public NavierStokesBase<TDomain>
{
	protected:
	///	Base class type
		typedef NavierStokesBase<TDomain> base_type;

	///	own type
		typedef NavierStokesFE<TDomain> this_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	public:
	///	Constructor (setting default values)
	/// \{
		NavierStokesFE(const char* functions, const char* subsets);
		NavierStokesFE(const std::vector<std::string>& vFct, const std::vector<std::string>& vSubset);
	/// \}

	///	sets the kinematic viscosity
		void set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > user);

	///	returns kinematic viscosity
		SmartPtr<CplUserData<number, dim> > kinematic_viscosity() {return m_imKinViscosity.user_data ();}

	///	sets the density
		void set_density(SmartPtr<CplUserData<number, dim> > user);

	///	returns density
		SmartPtr<CplUserData<number, dim> > density() {return m_imDensity.user_data ();}

	///	sets the source function
		void set_source(SmartPtr<CplUserData<MathVector<dim>, dim> > user);

	///	sets the stabilization parameter
		void set_stabilization(number alpha) {m_stabParam = alpha;}

	public:
	///	type of trial space for each function used
		bool request_finite_element_id(const std::vector<LFEID>& vLfeID);

	///	switches between non-regular and regular grids
		virtual bool request_non_regular_grid(bool bNonRegular);

	///	returns string identifying disc type
		virtual std::string disc_type() const {return "fe";};

	protected:
	///	stabilization parameter
		number m_stabParam;

	///	quadrature order
		int m_quadOrder;

	///	current shape function set
		LFEID m_vLFEID;
		LFEID m_pLFEID;

	///	current element
		GeometricObject* m_pElem;

		void init();

	private:
	///	Data import for source
		DataImport<MathVector<dim>, dim> m_imSource;

	///	Data import for kinematic viscosity
		DataImport<number, dim> m_imKinViscosity;

	///	Data import for kinematic viscosity
		DataImport<number, dim> m_imDensity;

	/// position access
		const MathVector<dim>* m_vCornerCoords;

	/// abbreviation for pressure
		static const size_t _P_ = dim;

		using base_type::m_bPecletBlend;
		using base_type::m_bExactJacobian;
		using base_type::m_bStokes;
		using base_type::m_bLaplace;

	public:
		template<typename TElem, typename VGeom, typename PGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);

		template<typename TElem, typename VGeom, typename PGeom>
		void prep_elem(TElem* elem, const LocalVector& u);

	///	finishes the loop over all elements
		template<typename TElem, typename VGeom, typename PGeom>
		void fsh_elem_loop();

	///	assembles the local stiffness matrix using a finite volume scheme
		template<typename TElem, typename VGeom, typename PGeom>
		void add_jac_A_elem(LocalMatrix& J, const LocalVector& u);

	///	assembles the local mass matrix using a finite volume scheme
		template<typename TElem, typename VGeom, typename PGeom>
		void add_jac_M_elem(LocalMatrix& J, const LocalVector& u);

	///	assembles the stiffness part of the local defect
		template<typename TElem, typename VGeom, typename PGeom>
		void add_def_A_elem(LocalVector& d, const LocalVector& u);

	///	assembles the mass part of the local defect
		template<typename TElem, typename VGeom, typename PGeom>
		void add_def_M_elem(LocalVector& d, const LocalVector& u);

	///	assembles the local right hand side
		template<typename TElem, typename VGeom, typename PGeom>
		void add_rhs_elem(LocalVector& d);

	// 	FVHO Assemblings
		void register_all_funcs(const LFEID& vLfeID, const LFEID& pLfeID, const int quadOrder);
		template<typename TElem, typename VGeom, typename PGeom> void register_func();

};

/// @}

} // namespace NavierStokes
} // end namespace ug

#endif /*__H__UG__PLUGINS__NAVIER_STOKES__NAVIER_STOKES_FE__*/
