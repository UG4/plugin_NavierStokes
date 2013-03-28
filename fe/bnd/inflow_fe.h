/*
 * inflow_fe.h
 *
 *  Created on: 01.09.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__NAVIER_STOKES__BND__INFLOW_FE__
#define __H__UG__NAVIER_STOKES__BND__INFLOW_FE__

#include "../../bnd/inflow_base.h"
#include "../navier_stokes_fe.h"

#include "lib_disc/spatial_disc/elem_disc/neumann_boundary/neumann_boundary_base.h"
#include "lib_disc/spatial_disc/constraints/dirichlet_boundary/lagrange_dirichlet_boundary.h"

namespace ug{
namespace NavierStokes{

template <	typename TDomain, typename TAlgebra>
class NavierStokesInflowFE
	: public NavierStokesInflowBase<TDomain, TAlgebra>
{
	private:
		static const int dim = TDomain::dim;

	public:
	///	returns the number of element discs
		virtual size_t num_elem_disc() const {return 0;}

	///	returns the element disc
		virtual SmartPtr<IElemDisc<TDomain> > elem_disc(size_t i) {return m_spNeumannDisc;}

	///	returns the number of constraints
		virtual size_t num_constraint() const {return 1;}

	///	returns an element disc
		virtual SmartPtr<IDomainConstraint<TDomain, TAlgebra> > constraint(size_t i) {return m_spDirichletConstraint;}

	public:
	///	Constructor
		NavierStokesInflowFE(SmartPtr< NavierStokesFE<TDomain> > spMaster);

	///	sets the velocity to a given value
		void add(SmartPtr<UserData<MathVector<dim>, dim> > user, const char* subsetsBND);

	protected:
	///	neumann disc for pressure equation
		SmartPtr<NeumannBoundaryBase<TDomain> > m_spNeumannDisc;

	///	dirichlet disc for velocity components
		SmartPtr<DirichletBoundary<TDomain,TAlgebra> > m_spDirichletConstraint;

	/// The master discretization:
		SmartPtr< NavierStokesFE<TDomain> > m_spMaster;
};

} // namespace NavierStokes
} // end namespace ug

#include "inflow_fe_impl.h"

#endif /* __H__UG__NAVIER_STOKES__BND__INFLOW_FE__ */
