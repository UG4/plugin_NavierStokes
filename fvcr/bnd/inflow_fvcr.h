/*
 * inflow_fvcr.h
 *
 *  Created on: 01.09.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__NAVIER_STOKES__BND__INFLOW_FVCR__
#define __H__UG__NAVIER_STOKES__BND__INFLOW_FVCR__

#include "../../bnd/inflow_base.h"
#include "../navier_stokes_fvcr.h"

#include "lib_disc/spatial_disc/constraints/dirichlet_boundary/lagrange_dirichlet_boundary.h"

namespace ug{
namespace NavierStokes{

template <	typename TDomain, typename TAlgebra>
class NavierStokesInflowFVCR
	: public NavierStokesInflowBase<TDomain, TAlgebra>
{
	private:
		static const int dim = TDomain::dim;

	public:
	///	returns the number of element discs
		virtual size_t num_elem_disc() const {return 0;}

	///	returns the element disc
		virtual SmartPtr<IElemDisc<TDomain> > elem_disc(size_t i) {return NULL;}

	///	returns the number of constraints
		virtual size_t num_constraint() const {return 1;}

	///	returns an element disc
		virtual SmartPtr<IDomainConstraint<TDomain, TAlgebra> > constraint(size_t i) {return m_spDirichletConstraint;}

	public:
	///	Constructor
		NavierStokesInflowFVCR(SmartPtr< NavierStokesFVCR<TDomain> > spMaster);

	///	sets the velocity to a given value
		void add(SmartPtr<CplUserData<MathVector<dim>, dim> > user, const char* subsetsBND);

	protected:
	///	dirichlet disc for velocity components
		SmartPtr<DirichletBoundary<TDomain,TAlgebra> > m_spDirichletConstraint;

	/// The master discretization:
		SmartPtr< NavierStokesFVCR<TDomain> > m_spMaster;
};

} // namespace NavierStokes
} // end namespace ug

#include "inflow_fvcr_impl.h"

#endif /* __H__UG__NAVIER_STOKES__BND__INFLOW_FVCR__ */
