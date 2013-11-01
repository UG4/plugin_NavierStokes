/*
 * compressible_navier_stokes_base.h
 *
 *  Created on: 31.10.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UG__PLUGINS__NAVIER_STOKES__COMPRESSIBLE__COMPRESSIBLE_NAVIER_STOKES_BASE__
#define __H__UG__PLUGINS__NAVIER_STOKES__COMPRESSIBLE__COMPRESSIBLE_NAVIER_STOKES_BASE__

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

/// Finite Volume Element Discretization for the compressible Navier-Stokes Equation
/**
 * This class implements the IElemDisc interface to provide element local
 * assemblings for the compressible Navier-Stokes equation.
 *
 * The unknowns of the equation are:
 * <ul>
 * <li> \f$ \vec{u} \f$ velocity
 * <li> \f$ p \f$ pressure
 * <li> \f$ \rho \f$ density.
 * </ul>
 *
 * The equation takes the form
 * \f{align*}
 * 	\frac{\partial \rho \vec{u}}{\partial t}
 * 	- \nabla \left( \rho \nu (\nabla \vec{u} + (\nabla \vec{u})^T)
 * 		- \frac{2}{3} Id * \nabla \vec{u} \right)
 * 	+ \nabla \cdot \left( \rho \vec{u} \vec{u}^T \right)
 *  + \nabla p
 *  &= \vec{f}_m\\
 *
 *  \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \vec{u}) &= f_c \\
 *
 *  \frac{1}{\gamma - 1}\frac{\partial \rho}{\partial t}
 *  + \frac{1}{2} \frac{\partial \rho \vec{u}^2}{\partial t}
 *  + \frac{\gamma}{\gamma - 1} \nabla \cdot (\rho \vec{u}^T)
 *  + \frac{1}{2} \nabla \cdot (\rho \vec{u}^2 \vec{u}^T)
 *  - \nabla \left( \vec{u}^T \cdot \tau \right) &= f_e \\
 *
 * \f}
 *
 * with
 * <ul>
 * <li> \f$ \tau \f$ is the shear stress tensor \f$ \tau = \rho \nu (\nabla \vec{u} + (\nabla \vec{u})^T)
 * 		- \frac{2}{3} Id * \nabla \vec{u} \f$
 * <li>	\f$ \nu \f$ is the kinematic viscosity (temporarily constant, implement)
 * <li>	\f$ \vec{f} \equiv f(\vec{x},t) \f$ is a Source Term (not implemented yet)
 * </ul>
 *
 * The first equation models the conservation of momentum and is therefore
 * referred to as the <b>Momentum equation</b>. The second equation models the
 * conservation of mass and is known as the <b>Mass continuity equation</b> or
 * simply <b>Continuity equation</b>.
 *
 * In component-wise notation, the equation reads
 *
 * \f{align*}
 * \frac{\partial \rho u_{d_1}}{\partial t}
 * - \sum_{d_2 = 1}^{\text{dim}} \frac{\partial}{\partial x_{d_2}}
 * 		\left( \rho \nu \left( \frac{\partial u_{d_1}}{\partial x_{d_2}}
 * 			+ \frac{\partial u_{d_2}}{\partial x_{d_1}} \right)\right)
 * + \sum_{d_2 = 1}^{\text{dim}}  \frac{\partial}{\partial x_{d_2}}
 * 		\left( \rho u_{d_2} u_{d_1} \right)
 * + \frac{\partial p}{\partial x_{d_1}}
 * &= f_{d_1} \qquad \forall \quad d_1 = 1, \dots, \text{dim}\\
 *
 * \sum_{d = 1}^{\text{dim}}  \frac{\partial \rho u_d}{\partial x_{d}} &= 0
 * \f}
 *
 * Since the density is assumed to be constant, it set to \f$\rho \equiv 1 \f$.
 *
 * The finite volume element discretization uses a dual mesh consisting of
 * Control Volumes \f$\mathcal{B}_h\f$, that cover the domain. For each
 * control volume \f$B \in \mathcal{B}_h\f$ the following equation is solved
 *
 * \f{align*}
 * 	\frac{\partial}{\partial t} \int_{B} \begin{pmatrix} \vec{u} \\ 0 \end{pmatrix} dV
 *  + \int_{\partial B}
 *  	\begin{pmatrix}
 *  		- \rho \nu \left(\nabla \vec{u} + (\nabla \vec{u})^T \right) \vec{n}
 *  		+ \vec{u} \vec{u}^T \vec{n} + p \vec{n} \\
 * 		 	\vec{u} \vec{n}
 * 		\end{pmatrix} dS
 * = \int_B \begin{pmatrix} \vec{f} \\ 0 \end{pmatrix} dV
 * \f}
 *
 * where \f$ \vec{n}\f$ is the outer normal on the boundary of the control volume
 * and \f$\int_B \cdot \; dV \f$ and \f$\int_{\partial B} \cdot \; dS \f$ are the
 * integration over the control volume and the integration over the boundary of
 * the control volume. By \$f B_{sh} \$f is denoted the control volume associated
 * to a shape function (i.e. vertices, since we use P1 Lagrange ansatz functions),
 * T usually denotes the Element.
 *
 * In order to number the local unknowns, we use the following notation.
 * <ul>
 * <li> \f$ sh = 0, \dots, n_{sh} - 1\f$ loop the \f$n_{sh}\f$ shape functions. For
 * 		this implementation these are the corners of the element
 * <li> \f$ d = 0, \dots, \text{dim} - 1 \f$ are the velocity component of each
 * 		corner
 * <li> _P_ := dim indicates the pressure component
 * </ul>
 * The access to local unknowns of local vectors and matrices is now given by
 * \f{align*}
 * 	&\mathbf{d}(d_1, sh), \mathbf{d}(\text{\_P\_}, sh),\\
 *  &\mathcal{J}(d_1, sh_1, d_2, sh_2), \mathcal{J}(d_1, sh_1, \text{\_P\_}, sh_2),
 * \f}
 *
 * \tparam	TDomain		Domain
 * \tparam	TAlgebra	Algebra
 */
template<	typename TDomain>
class CompressibleNavierStokesBase
	: public NavierStokesBase<TDomain>
{
	protected:
	///	Base class type
		typedef NavierStokesBase<TDomain> base_type;

	///	own type
		typedef CompressibleNavierStokesBase<TDomain> this_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	public:
	///	Constructor (setting default values)
	/// \{
		CompressibleNavierStokesBase(const char* functions, const char* subsets);
		CompressibleNavierStokesBase(const std::vector<std::string>& vFct, const std::vector<std::string>& vSubset);
	/// \}

	///	sets the kinematic viscosity
	/**
	 * This method sets the kinematic viscosity value.
	 *
	 * \param[in]	data		kinematic Viscosity
	 */
	///	\{
		virtual void set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > user) = 0;

	///	returns kinematic viscosity
		virtual SmartPtr<CplUserData<number, dim> > kinematic_viscosity() = 0;

	///	sets the adiabatic index (also known as 'heat capacity ratio' or 'ratio of specific heats')
		virtual void set_adiabatic_index(SmartPtr<CplUserData<number, dim> > user) = 0;
		void set_adiabatic_index(number val);
#ifdef UG_FOR_LUA
		void set_adiabatic_index(const char* fctName);
#endif

	///	returns adiabatic index
		virtual SmartPtr<CplUserData<number, dim> > adiabatic_index() = 0;

	///	sets the source function
	/**
	 * This method sets the source value. A zero value is assumed as default.
	 * \param[in]	data		source data
	 */
	///	\{
		virtual void set_source(SmartPtr<CplUserData<MathVector<dim>, dim> > user) = 0;

	///	sets if Mach-number blending is used in momentum equation
		void set_mach_number_blend(bool machNrBlend) {m_bMachNrBlend = machNrBlend;}

  	public:
	///	returns if local time series is needed
		virtual bool requests_local_time_series() {return true;}

	///	returns string identifying disc type
		virtual std::string disc_type() const = 0;

	protected:
	///	flag if using Mach-number Blending
		bool m_bMachNrBlend;

	///	factor for exact jacobian, (1 for exact jacobian, 0 for fix point)
		using base_type::m_bFullNewtonFactor;
};

/// @}

} // namespace NavierStokes
} // end namespace ug

#endif /*__H__UG__PLUGINS__NAVIER_STOKES__COMPRESSIBLE__COMPRESSIBLE_NAVIER_STOKES_BASE__*/


