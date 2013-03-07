/*
 * no_normal_stress_outflow.h
 *
 *  Created on: 27.03.2012
 *  D. Logashenko, A. Vogel
 */

#ifndef __H__UG__PLUGINS__NAVIER_STOKES__NO_NORMAL_STRESS_OUTFLOW_FV1_
#define __H__UG__PLUGINS__NAVIER_STOKES__NO_NORMAL_STRESS_OUTFLOW_FV1_

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
class NavierStokesNoNormalStressOutflowFV1
	: public NavierStokesNoNormalStressOutflowBase<TDomain>
{
	private:
	///	Base class type
		typedef NavierStokesNoNormalStressOutflowBase<TDomain> base_type;

	///	own type
		typedef NavierStokesNoNormalStressOutflowFV1<TDomain> this_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	public:
	///	Constructor (setting default values)
		NavierStokesNoNormalStressOutflowFV1(SmartPtr< NavierStokesBase<TDomain> > spMaster);

	public:
	///	type of trial space for each function used
		virtual bool request_finite_element_id(const std::vector<LFEID>& vLfeID);

	///	switches between non-regular and regular grids
		virtual bool request_non_regular_grid(bool bNonRegular);

	public:
	///	prepares the element loop
		template <typename TElem, typename TFVGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);

	///	prepares the element for evaluation
		template <typename TElem, typename TFVGeom>
		void prep_elem(TElem* elem, const LocalVector& u);

	///	finishes the element loop
		template <typename TElem, typename TFVGeom>
		void fsh_elem_loop();

	///	adds the stiffness part to the local jacobian
		template <typename TElem, typename TFVGeom>
		void add_jac_A_elem(LocalMatrix& J, const LocalVector& u);

	///	adds the stiffness part to the local defect
		template <typename TElem, typename TFVGeom>
		void add_def_A_elem(LocalVector& d, const LocalVector& u);

	public:
	///	dummy implementations
	///	\{
		template <typename TElem, typename TFVGeom>
		void add_jac_M_elem(LocalMatrix& J, const LocalVector& u){}
		template <typename TElem, typename TFVGeom>
		void add_def_M_elem(LocalVector& d, const LocalVector& u){}
		template <typename TElem, typename TFVGeom>
		void add_rhs_elem(LocalVector& d){}
	/// \}

	private:
	/// adds the diffusive part of the local Jacobian of the momentum equation
		template <typename BF>
		inline void diffusive_flux_Jac
		(
			const size_t ip,
			const BF& bf,
			LocalMatrix& J,
			const LocalVector& u
		);
	/// adds the diffusive part of the local defect of the momentum equation
		template <typename BF>
		inline void diffusive_flux_defect
		(
			const size_t ip,
			const BF& bf,
			LocalVector& d,
			const LocalVector& u
		);
	/// adds the convective part of the local Jacobian of the momentum equation
		template <typename BF>
		inline void convective_flux_Jac
		(
			const size_t ip,
			const BF& bf,
			LocalMatrix& J,
			const LocalVector& u
		);
	/// adds the convective part of the local defect of the momentum equation
		template <typename BF>
		inline void convective_flux_defect
		(
			const size_t ip,
			const BF& bf,
			LocalVector& d,
			const LocalVector& u
		);
	
	protected:
	/// abbreviation for pressure
		static const size_t _P_ = dim;

		using base_type::m_imDensity;
		using base_type::m_imKinViscosity;
		using base_type::m_spMaster;

		using base_type::m_vLocIP;
		using base_type::m_vGloIP;
		using base_type::m_vBndSubSetIndex;

	protected:
		void register_all_funcs(bool bHang);
		template<typename TElem, typename TFVGeom>
		void register_func();

};

/// @}

} // namespace NavierStokes
} // end namespace ug

#endif /*__H__UG__PLUGINS__NAVIER_STOKES__NO_NORMAL_STRESS_OUTFLOW_FV1_*/
