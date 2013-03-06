/*
 * wall.h
 *
 *  Created on: 01.09.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__NAVIER_STOKES__BND__WALL__
#define __H__UG__NAVIER_STOKES__BND__WALL__

#include "lib_disc/spatial_disc/disc_item.h"

#include "lib_disc/spatial_disc/constraints/dirichlet_boundary/lagrange_dirichlet_boundary.h"

namespace ug{
namespace NavierStokes{

template <	typename TDomain, typename TAlgebra>
class NavierStokesWall
	: public IDiscretizationItem<TDomain,TAlgebra>
{
	public:
	///	returns the number of element discs
		virtual size_t num_elem_disc() const {return 0;}

	///	returns the element disc
		virtual SmartPtr<IDomainElemDisc<TDomain> > elem_disc(size_t i) {return NULL;}

	///	returns the number of constraints
		virtual size_t num_constraint() const {return 1;}

	///	returns an element disc
		virtual SmartPtr<IDomainConstraint<TDomain, TAlgebra> > constraint(size_t i) {return m_spDirichletConstraint;}

	public:
	///	Constructor
		NavierStokesWall(SmartPtr< NavierStokes<TDomain> > spMaster);

	///	sets the velocity to a given value
		void add(const char* subsetsBND);

	protected:
	///	dirichlet disc for velocity components
		SmartPtr<DirichletBoundary<TDomain,TAlgebra> > m_spDirichletConstraint;

	///	name of velocity+pressure components
		std::vector<std::string> m_vFctName;

	/// The master discretization:
		SmartPtr< NavierStokes<TDomain> > m_spMaster;
};

} // namespace NavierStokes
} // end namespace ug

#include "wall_impl.h"

#endif /* __H__UG__NAVIER_STOKES__BND__WALL__ */
