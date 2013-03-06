/*
 * inflow.h
 *
 *  Created on: 01.09.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__NAVIER_STOKES__BND__INFLOW__
#define __H__UG__NAVIER_STOKES__BND__INFLOW__

#include "lib_disc/spatial_disc/disc_item.h"

#include "lib_disc/spatial_disc/elem_disc/neumann_boundary/neumann_boundary_base.h"
#include "lib_disc/spatial_disc/constraints/dirichlet_boundary/lagrange_dirichlet_boundary.h"

namespace ug{
namespace NavierStokes{

template <	typename TDomain, typename TAlgebra>
class NavierStokesInflow
	: public IDiscretizationItem<TDomain, TAlgebra>
{
	private:
		const static int dim = TDomain::dim;

	public:
	///	returns the number of element discs
		virtual size_t num_elem_disc() const {
			if (m_spMaster->disc_scheme() != "staggered") return 1;
			else return 0;
		}

	///	returns the element disc
		virtual SmartPtr<IDomainElemDisc<TDomain> > elem_disc(size_t i) {return m_spNeumannDisc;}

	///	returns the number of constraints
		virtual size_t num_constraint() const {return 1;}

	///	returns an element disc
		virtual SmartPtr<IDomainConstraint<TDomain, TAlgebra> > constraint(size_t i) {return m_spDirichletConstraint;}

	public:
	///	Constructor
		NavierStokesInflow(SmartPtr< NavierStokes<TDomain> > spMaster);

	///	sets the velocity to a given value
	///	\{
		void add(SmartPtr<UserData<MathVector<dim>, dim> > user, const char* subsetsBND);
		void add(const std::vector<number>& vVel, const char* subsetsBND);
#ifdef UG_FOR_LUA
		void add(const char* luaFctName, const char* subsetsBND);
#endif
	///	\}

	protected:
	///	neumann disc for pressure equation
		SmartPtr<NeumannBoundaryBase<TDomain> > m_spNeumannDisc;

	///	dirichlet disc for velocity components
		SmartPtr<DirichletBoundary<TDomain,TAlgebra> > m_spDirichletConstraint;

	///	name of velocity+pressure components
		std::vector<std::string> m_vFctName;

	/// The master discretization:
		SmartPtr< NavierStokes<TDomain> > m_spMaster;
};

} // namespace NavierStokes
} // end namespace ug

#include "inflow_impl.h"

#endif /* __H__UG__NAVIER_STOKES__BND__INFLOW__ */
