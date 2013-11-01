/*
 * navier_stokes_fv.h
 *
 *  Created on: 20.09.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FV__NAVIER_STOKES_FV__
#define __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FV__NAVIER_STOKES_FV__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"

#include "../incompressible_navier_stokes_base.h"

namespace ug{
namespace NavierStokes{

/// \ingroup lib_disc_elem_disc
/// @{

template<	typename TDomain>
class NavierStokesFV
	: public IncompressibleNavierStokesBase<TDomain>
{
	private:
	///	Base class type
		typedef IncompressibleNavierStokesBase<TDomain> base_type;

	///	own type
		typedef NavierStokesFV<TDomain> this_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	public:
	///	Constructor (setting default values)
	/// \{
		NavierStokesFV(const char* functions, const char* subsets);
		NavierStokesFV(const std::vector<std::string>& vFct, const std::vector<std::string>& vSubset);
	/// \}

	///	sets the kinematic viscosity
		void set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > user);

	///	returns kinematic viscosity
		SmartPtr<CplUserData<number, dim> > kinematic_viscosity() {return m_imKinViscosity.user_data ();}

	///	sets the density
		void set_density(SmartPtr<CplUserData<number, dim> > user);

	///	returns density
		SmartPtr<CplUserData<number, dim> > density() {return m_imDensitySCVF.user_data ();}

	///	sets the source function
		void set_source(SmartPtr<CplUserData<MathVector<dim>, dim> > user);

	public:
	///	type of trial space for each function used
		void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

	///	returns string identifying disc type
		virtual std::string disc_type() const {return "fv";};

	protected:
	///	current shape function set
		LFEID m_vLFEID;
		LFEID m_pLFEID;

	///	quadrature order
		int m_quadOrder;

		void init();

	private:
	///	Data import for source
		DataImport<MathVector<dim>, dim> m_imSource;

	///	Data import for kinematic viscosity
		DataImport<number, dim> m_imKinViscosity;

	///	Data import for density
		DataImport<number, dim> m_imDensitySCVF;
		DataImport<number, dim> m_imDensitySCVFp;
		DataImport<number, dim> m_imDensitySCV;

	/// abbreviation for pressure
		static const size_t _P_ = dim;

		using base_type::m_bPecletBlend;
		using base_type::m_bFullNewtonFactor;
		using base_type::m_bStokes;
		using base_type::m_bLaplace;

	public:
		template<typename TElem, typename VGeom, typename PGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);

		template<typename TElem, typename VGeom, typename PGeom>
		void prep_elem(const LocalVector& u, GeometricObject* elem, const MathVector<dim> vCornerCoords[]);

	///	finishes the loop over all elements
		template<typename TElem, typename VGeom, typename PGeom>
		void fsh_elem_loop();

	///	assembles the local stiffness matrix using a finite volume scheme
		template<typename TElem, typename VGeom, typename PGeom>
		void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GeometricObject* elem, const MathVector<dim> vCornerCoords[]);

	///	assembles the local mass matrix using a finite volume scheme
		template<typename TElem, typename VGeom, typename PGeom>
		void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GeometricObject* elem, const MathVector<dim> vCornerCoords[]);

	///	assembles the stiffness part of the local defect
		template<typename TElem, typename VGeom, typename PGeom>
		void add_def_A_elem(LocalVector& d, const LocalVector& u, GeometricObject* elem, const MathVector<dim> vCornerCoords[]);

	///	assembles the mass part of the local defect
		template<typename TElem, typename VGeom, typename PGeom>
		void add_def_M_elem(LocalVector& d, const LocalVector& u, GeometricObject* elem, const MathVector<dim> vCornerCoords[]);

	///	assembles the local right hand side
		template<typename TElem, typename VGeom, typename PGeom>
		void add_rhs_elem(LocalVector& d, GeometricObject* elem, const MathVector<dim> vCornerCoords[]);

	protected:
	///	register util
		void register_all_funcs(const LFEID& vLfeID, const LFEID& pLfeID);
		template<typename TElem, typename VGeom, typename PGeom>
		void register_func();

		std::vector<std::vector<number> > m_vvPShape;
		std::vector<std::vector<number> > m_vvVShape;
};

/// @}

} // namespace NavierStokes
} // end namespace ug

#endif /*__H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FV__NAVIER_STOKES_FV__*/
