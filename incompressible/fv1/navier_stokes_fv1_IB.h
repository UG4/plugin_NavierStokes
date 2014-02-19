/*
 * navier_stokes.h
 *
 *  Created on: 20.09.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__PLUGINS__PARTICLE_LADEN_FLOW__NAVIER_STOKES_FVIB__
#define __H__UG__PLUGINS__PARTICLE_LADEN_FLOW__NAVIER_STOKES_FVIB__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"

#include "../incompressible_navier_stokes_base.h"
#include "../../upwind_interface.h"
#include "stabilization.h"
#include "navier_stokes_fv1.h"

#include "lib_disc/spatial_disc/disc_util/fv1ib_geom.h"


namespace ug{
namespace NavierStokes{

/// \ingroup lib_disc_elem_disc
/// @{

/// Finite Volume Element Discretization for the incompressible Navier-Stokes Equation
/**
 * This class implements the IElemDisc interface to provide element local
 * assemblings for the incompressible Navier-Stokes equation.
 *
 * The unknowns of the equation are:
 * <ul>
 * <li> \f$ \vec{u} \f$ velocity
 * <li> \f$ p \f$ pressure.
 * </ul>
 *
 * The equation takes the form
 * \f{align*}
 * 	\frac{\partial \rho \vec{u}}{\partial t}
 * 	- \nabla \left( \rho \nu (\nabla \vec{u} + (\nabla \vec{u})^T) \right)
 * 	+ \nabla \cdot \left( \rho \vec{u} \vec{u}^T \right)
 *  + \nabla p
 *  &= \vec{f}\\
 *
 *  \nabla \cdot (\rho \vec{u}) &= 0
 * \f}
 *
 * with
 * <ul>
 * <li>	\f$ \rho \f$ is the constant density
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
class NavierStokesFV1IB
 	 : public NavierStokesFV1<TDomain>

{
	protected:
	///	Base class type
		typedef NavierStokesFV1<TDomain> base_type;

	///	own type
		typedef NavierStokesFV1IB<TDomain> this_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	public:
	///	Constructor (setting default values)
	/// \{
		NavierStokesFV1IB(const char* functions, const char* subsets);
 		NavierStokesFV1IB(const std::vector<std::string>& vFct, const std::vector<std::string>& vSubset) : NavierStokesFV1<TDomain>(vFct, vSubset){};
	/// \}

		template <typename TElem>
 		void adapt_FVGeometry(FV1IBGeometry<TElem, dim> geo, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]){};

		template <typename TElem, typename TFVGeom>
		void prep_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	/// copy of the original 'add_jac_A_elem()' implementation with adaptions due to the inner boundary:
	/// ImpEq-P1: No pressure term for the momentum equation for scvf being part of the IB
	/// ImpEq-P2: Pressure value for computation of the pressure term in ImpEq in ip's of the remaining
	///		amount of scvf's (= 1scvf (2d); = 3scvf (3d)) depend ONLY on
	/// 	the "participating" corners of the ip, i.e. pressure = average of the pressure in corners
	/// ContEq-V1: Pure mass flux for the continuity equation for scvf being part of the IB
	/// ContEq-V2: Velocity value for computation of the mass flux in ContEq in ip's of the remaining
	///		amount of scvf's (= 1scvf (2d); = 3scvf (3d)) depend ONLY on
	/// 	the "participating" corners of the ip, i.e. pressure = average of the pressure in corners
		template <typename TElem, typename TFVGeom>
		void add_jac_A_elem_IB(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		template <typename TElem, typename TFVGeom>
		void add_def_A_elem_IB(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		template <typename TElem, typename TFVGeom>
		void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{
			bool elemIsCut_byIB = false;

			if ( elemIsCut_byIB )
				add_jac_A_elem_IB<TElem, TFVGeom>(J, u, elem, vCornerCoords);
			else
				this->NavierStokesFV1<TDomain>::template add_jac_A_elem<TElem, TFVGeom>(J, u, elem, vCornerCoords);
				//this->NavierStokesFV1<TDomain>::template add_jac_A_elem<TElem, FV1Geometry<TElem, dim> >(J, u, elem, vCornerCoords);

		}


		template <typename TElem, typename TFVGeom>
		void add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{
			UG_LOG("hier Aufruf: add_def_A_elem()...\n");
			UG_LOG("jeppa!...\n");

			bool elemIsCut_byIB = false;

			if ( elemIsCut_byIB )
				add_def_A_elem_IB<TElem, TFVGeom>(d, u, elem, vCornerCoords);
			else
				this->NavierStokesFV1<TDomain>::template add_def_A_elem<TElem,TFVGeom>(d, u, elem, vCornerCoords);
				//this->NavierStokesFV1<TDomain>::template add_def_A_elem<TElem, FV1Geometry<TElem, dim> >(d, u, elem, vCornerCoords);

 		}

		void init();

	///	type of trial space for each function used
		void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

	private:
		void register_substitution_funcs(bool bHang);
		template <typename TElem, typename TFVGeom>
		void register_substitutes();



};


/// @}

} // namespace NavierStokes
} // end namespace ug

#include "navier_stokes_fv1_IB_impl.h"


#endif /*__H__UG__PLUGINS__PARTICLE_LADEN_FLOW__NAVIER_STOKES_FVIB__*/
