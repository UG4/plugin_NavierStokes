/*
 * navier_stokes_fvcr.h
 *
 *  Created on: 20.09.2010
 *      Author: cwehner
 */

#ifndef __H__UG__PLUGINS__NAVIER_STOKES__NAVIER_STOKES_FVCR__
#define __H__UG__PLUGINS__NAVIER_STOKES__NAVIER_STOKES_FVCR__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"

#include "../navier_stokes_base.h"
#include "../upwind_interface.h"

namespace ug{
namespace NavierStokes{

/// \ingroup lib_disc_elem_disc
/// @{

template<	typename TDomain>
class NavierStokesFVCR
	: public NavierStokesBase<TDomain>
{
	private:
	///	Base class type
		typedef NavierStokesBase<TDomain> base_type;

	///	own type
		typedef NavierStokesFVCR<TDomain> this_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	public:
	///	Constructor (setting default values)
	/// \{
		NavierStokesFVCR(const char* functions, const char* subsets);
		NavierStokesFVCR(const std::vector<std::string>& vFct, const std::vector<std::string>& vSubset);
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

		void set_defect_upwind(bool defectUpwind) { m_bDefectUpwind = defectUpwind;}
		bool get_defect_upwind() {return m_bDefectUpwind; }

	///	sets an upwinding for the convective term of momentum equation
		void set_upwind(SmartPtr<INavierStokesUpwind<dim> > spUpwind)
			{m_spConvUpwind = spUpwind;}

	///	sets the upwind based on a string identifier
		void set_upwind(const std::string& name)
			{m_spConvUpwind = CreateNavierStokesUpwind<dim>(name);}

	public:
	///	type of trial space for each function used
		void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

	///	returns string identifying disc type
		virtual std::string disc_type() const {return "fvcr";};

	protected:
		void init();

	public:
	/// mixed upwind for Crouzeix-Raviart elements
		template <typename TFVGeom>
		inline number peclet_blend(MathVector<dim>& UpwindVel, const TFVGeom& geo, size_t ip,
		                           const MathVector<dim>& StdVel, number kinVisco);

	protected:
	/// flag if using upwind in defect computation
		bool m_bDefectUpwind;

	///	Data import for source
		DataImport<MathVector<dim>, dim> m_imSource;

	///	Data import for kinematic viscosity
		DataImport<number, dim> m_imKinViscosity;

	///	Data import for kinematic viscosity
		DataImport<number, dim> m_imDensitySCVF;
		DataImport<number, dim> m_imDensitySCV;

	///	Upwinding for velocity in convective term of momentum equation
		SmartPtr<INavierStokesUpwind<dim> > m_spConvUpwind;

	/// position access
		const MathVector<dim>* m_vCornerCoords;

	/// abbreviation for pressure
		static const size_t _P_ = dim;

		using base_type::m_bPecletBlend;
		using base_type::m_bExactJacobian;
		using base_type::m_bStokes;
		using base_type::m_bLaplace;

	public:
		template <typename TElem, typename TFVGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);

		template <typename TElem, typename TFVGeom>
		void prep_elem(TElem* elem, LocalVector& u);

		template <typename TElem, typename TFVGeom>
		void fsh_elem_loop();
		
		template <typename TElem, typename TFVGeom>
		void add_jac_A_elem(LocalMatrix& J, const LocalVector& u);
		
		template <typename TElem, typename TFVGeom>
		void add_def_M_elem(LocalVector& d, const LocalVector& u);

		template <typename TElem, typename TFVGeom>
		void add_def_A_elem(LocalVector& d, const LocalVector& u);

		template <typename TElem, typename TFVGeom>
		void add_jac_M_elem(LocalMatrix& J, const LocalVector& u);

		template <typename TElem, typename TFVGeom>
		void add_rhs_elem(LocalVector& d);

	private:
		template <typename TElem, typename TFVGeom>
		void register_func();
		void register_all_funcs(bool bHang);
};

/// @}

} // namespace NavierStokes
} // end namespace ug

#endif /*__H__UG__PLUGINS__NAVIER_STOKES__NAVIER_STOKES_FVCR__*/
