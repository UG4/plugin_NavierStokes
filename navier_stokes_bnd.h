/*
 * navier_stokes_bnd.h
 *
 *  Created on: 01.09.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__NAVIER_STOKES_BND__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__NAVIER_STOKES_BND__

#include "lib_disc/spatial_disc/disc_item.h"

#include "lib_disc/spatial_disc/elem_disc/neumann_boundary/neumann_boundary.h"
#include "lib_disc/spatial_disc/constraints/dirichlet_boundary/lagrange_dirichlet_boundary.h"

namespace ug{
namespace NavierStokes{

template <	typename TDomain,
			typename TAlgebra>
class NavierStokesInflow
	: public IDiscretizationItem<TDomain, TAlgebra>
{
	private:
		const static int dim = TDomain::dim;

	public:
	///	returns the number of element discs
		virtual size_t num_elem_disc() const {return 1;}

	///	returns the element disc
		virtual SmartPtr<IDomainElemDisc<TDomain> > elem_disc(size_t i) {return m_spNeumannDisc;}

	///	returns the number of constraints
		virtual size_t num_constraint() const {return 1;}

	///	returns an element disc
		virtual SmartPtr<IDomainConstraint<TDomain, TAlgebra> > constraint(size_t i) {return m_spDirichletConstraint;}

	public:
	///	Constructor
		NavierStokesInflow(const char* functions, const char* subsets);

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
	///	sets the symbolic names for velocity and pressure
		void set_functions(const char* functions);

	///	neumann disc for pressure equation
		SmartPtr<NeumannBoundary<TDomain> > m_spNeumannDisc;

	///	dirichlet disc for velocity components
		SmartPtr<DirichletBoundary<TDomain,TAlgebra> > m_spDirichletConstraint;

	///	name of velocity components
		std::string m_velNames;

	///	name of pressure components
		std::string m_pressureName;
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <	typename TDomain,
			typename TAlgebra>
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

	///	virtual destructor
		~NavierStokesWall() {}

	public:
	///	Constructor
		NavierStokesWall(const char* functions);

	///	sets the velocity to a given value
		void add(const char* subsetsBND);

	protected:
	///	sets the symbolic names for velocity and pressure
		void set_functions(const char* functions);

	///	dirichlet disc for velocity components
		SmartPtr<DirichletBoundary<TDomain,TAlgebra> > m_spDirichletConstraint;

	///	name of velocity components
		std::vector<std::string> m_velNames;
};

} // namespace NavierStokes
} // end namespace ug

#include "navier_stokes_bnd_impl.h"

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__NAVIER_STOKES_BND__ */
