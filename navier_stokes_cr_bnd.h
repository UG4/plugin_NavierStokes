/*
 * navier_stokes_cr_bnd.h
 *
 *  Created on: 23.07.2012
 *      Author: Christian Wehner
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__NAVIER_STOKES_CR_BND__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__NAVIER_STOKES_CR_BND__

#include "lib_disc/spatial_disc/disc_item.h"

#include "lib_disc/spatial_disc/elem_disc/neumann_boundary/neumann_boundary.h"
#include "lib_disc/spatial_disc/constraints/dirichlet_boundary/lagrange_dirichlet_boundary.h"

namespace ug{
namespace NavierStokes{

template <	typename TDomain,
			typename TAlgebra>
class CRNavierStokesInflow
	: public IDiscretizationItem<TDomain, TAlgebra>
{
	private:
		const static int dim = TDomain::dim;

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
		CRNavierStokesInflow(const char* functions, const char* subsets);

	///	sets the velocity to a given value
	///	\{
		void add(SmartPtr<UserData<MathVector<dim>, dim> > user, const char* subsetsBND);
		void add(number vel_x, const char* subsetsBND);
		void add(number vel_x, number vel_y, const char* subsetsBND);
		void add(number vel_x, number vel_y, number vel_z, const char* subsetsBND);
#ifdef UG_FOR_LUA
		void add(const char* fctName, const char* subsetsBND);
#endif
	///	\}

	protected:
	///	sets the symbolic names for velocity
		void set_functions(const char* functions);

	///	dirichlet disc for velocity components
		SmartPtr<DirichletBoundary<TDomain,TAlgebra> > m_spDirichletConstraint;

	///	name of velocity components
		std::string m_velNames;

};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <	typename TDomain,
			typename TAlgebra>
class CRNavierStokesWall
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

	///	virtual destructor
		~CRNavierStokesWall() {}

	public:
	///	Constructor
		CRNavierStokesWall(const char* functions);

	///	sets the velocity to a given value
		void add(const char* subsetsBND);

	protected:
	///	sets the symbolic names for velocity
		void set_functions(const char* functions);

	///	dirichlet disc for velocity components
		SmartPtr<DirichletBoundary<TDomain,TAlgebra> > m_spDirichletConstraint;

	///	name of velocity components
		std::vector<std::string> m_velNames;
};

} // namespace NavierStokes
} // end namespace ug

#include "navier_stokes_cr_bnd_impl.h"

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__NAVIER_STOKES_CR_BND__ */
