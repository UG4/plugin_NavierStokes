/*
 * no_normal_stress_outflow_fv.h
 *
 *  Created on: 27.03.2012
 *  D. Logashenko, A. Vogel
 */

#ifndef __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FV__BND__NO_NORMAL_STRESS_OUTFLOW_FV_
#define __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FV__BND__NO_NORMAL_STRESS_OUTFLOW_FV_

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"

#include "../../bnd/no_normal_stress_outflow_base.h"

namespace ug{
namespace NavierStokes{


/// \ingroup lib_disc_elem_disc
/// @{

/// The zero-stress (neutral) outflow boundary condition for the incompressible NS equation
/**
 * This class implements the so-called neutral boundary condition for outflow
 * boundaries. Note that this class can be used only with the stabilized
 * vertex-centered discretization of the Navier-Stokes equations.
 *
 * This boundary condition imposes two equations on the unknown functions
 * at the outflow boundary \f$ F \f$:
 * <ul>
 * <li> \f$ \int_F p ds = 0 \f$, and
 * <li> \f$ \sigma \vec{n} \cdot \vec{n} = 0 \f$ on \f$ F \f$.
 * </ul>
 * where \f$ \vec{n} \f$ is the outer normal at \f$ F \f$ and
 * \f$ \sigma = \mu (\nabla \vec{u} + (\nabla \vec{u})^T) \f$ the stress tensor.
 */
template<	typename TDomain>
class NavierStokesNoNormalStressOutflowFV
	: public NavierStokesNoNormalStressOutflowBase<TDomain>
{
	private:
	///	Base class type
		typedef NavierStokesNoNormalStressOutflowBase<TDomain> base_type;

	///	own type
		typedef NavierStokesNoNormalStressOutflowFV<TDomain> this_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	public:
	///	Constructor (setting default values)
		NavierStokesNoNormalStressOutflowFV(SmartPtr< IncompressibleNavierStokesBase<TDomain> > spMaster);

	protected:
	///	sets the kinematic viscosity
		virtual void set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > data)
			{m_imKinViscosity.set_data(data);}

	///	sets the density
		virtual void set_density(SmartPtr<CplUserData<number, dim> > data)
			{m_imDensity.set_data(data); m_imDensityP.set_data(data);}

	public:
	///	type of trial space for each function used
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

	public:
	///	prepares the element loop
		template<typename TElem, typename VGeom, typename PGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);

	///	prepares the element for evaluation
		template<typename TElem, typename VGeom, typename PGeom>
		void prep_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	finishes the element loop
		template<typename TElem, typename VGeom, typename PGeom>
		void fsh_elem_loop();

	///	adds the stiffness part to the local jacobian
		template<typename TElem, typename VGeom, typename PGeom>
		void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	adds the stiffness part to the local defect
		template<typename TElem, typename VGeom, typename PGeom>
		void add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	public:
	///	dummy implementations
	///	\{
		template<typename TElem, typename VGeom, typename PGeom>
		void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]){}
		template<typename TElem, typename VGeom, typename PGeom>
		void add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]){}
		template<typename TElem, typename VGeom, typename PGeom>
		void add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[]){}
	/// \}

	private:
	/// adds the diffusive part of the local Jacobian of the momentum equation
		template <typename BF>
		inline void diffusive_flux_Jac
		(
			const size_t i,
			const size_t ip,
			const BF& bf,
			LocalMatrix& J,
			const LocalVector& u
		);
	/// adds the diffusive part of the local defect of the momentum equation
		template <typename BF>
		inline void diffusive_flux_defect
		(
			const size_t i,
			const size_t ip,
			const BF& bf,
			LocalVector& d,
			const LocalVector& u
		);
	/// adds the convective part of the local Jacobian of the momentum equation
		template <typename BF>
		inline void convective_flux_Jac
		(
			const size_t i,
			const size_t ip,
			const BF& bf,
			LocalMatrix& J,
			const LocalVector& u
		);
	/// adds the convective part of the local defect of the momentum equation
		template <typename BF>
		inline void convective_flux_defect
		(
			const size_t i,
			const size_t ip,
			const BF& bf,
			LocalVector& d,
			const LocalVector& u
		);

	protected:
	/// abbreviation for pressure
		static const size_t _P_ = dim;

		using base_type::m_spMaster;
		using base_type::m_vBndSubSetIndex;

	///	Data import for kinematic viscosity
		DataImport<number, dim> m_imKinViscosity;

	/// Data import for density
		DataImport<number, dim> m_imDensity;
		DataImport<number, dim> m_imDensityP;

	/// Boundary integration points of the viscosity and the density
		std::vector<MathVector<dim> > m_vLocIPv;
		std::vector<MathVector<dim> > m_vGloIPv;
		std::vector<MathVector<dim> > m_vLocIPp;
		std::vector<MathVector<dim> > m_vGloIPp;

	protected:
	///	current shape function set
		LFEID m_vLFEID;
		LFEID m_pLFEID;

	///	quadrature order
		int m_quadOrder;

	protected:
	///	register util
		void register_all_funcs(const LFEID& vLfeID, const LFEID& pLfeID);
		template<typename TElem, typename VGeom, typename PGeom>
		void register_func();

		std::vector<std::vector<number> > m_vvVShape;
};

/// @}

} // namespace NavierStokes
} // end namespace ug

#endif /*__H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FV__BND__NO_NORMAL_STRESS_OUTFLOW_FV_*/
