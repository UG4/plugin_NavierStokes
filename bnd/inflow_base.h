/*
 * inflow_base.h
 *
 *  Created on: 01.09.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__NAVIER_STOKES__BND__INFLOW_BASE__
#define __H__UG__NAVIER_STOKES__BND__INFLOW_BASE__

#include "lib_disc/spatial_disc/disc_item.h"

namespace ug{
namespace NavierStokes{

template <	typename TDomain, typename TAlgebra>
class NavierStokesInflowBase
	: public IDiscretizationItem<TDomain, TAlgebra>
{
	private:
		static const int dim = TDomain::dim;

	public:
	///	returns the number of element discs
		virtual size_t num_elem_disc() const = 0;

	///	returns the element disc
		virtual SmartPtr<IElemDisc<TDomain> > elem_disc(size_t i) = 0;

	///	returns the number of constraints
		virtual size_t num_constraint() const = 0;

	///	returns an element disc
		virtual SmartPtr<IDomainConstraint<TDomain, TAlgebra> > constraint(size_t i) = 0;

	public:
	///	sets the velocity to a given value
	///	\{
		virtual void add(SmartPtr<UserData<MathVector<dim>, dim> > user, const char* subsetsBND) = 0;
		void add(const std::vector<number>& vVel, const char* subsetsBND);
#ifdef UG_FOR_LUA
		void add(const char* luaFctName, const char* subsetsBND);
#endif
	///	\}
};

} // namespace NavierStokes
} // end namespace ug

#include "inflow_base_impl.h"

#endif /* __H__UG__NAVIER_STOKES__BND__INFLOW_BASE__ */
