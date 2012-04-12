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
		NavierStokesInflow(const char* functions, const char* subsets)
			: m_spNeumannDisc(new NeumannBoundary<TDomain>(subsets)),
			  m_spDirichletConstraint(new DirichletBoundary<TDomain,TAlgebra>)
		{
			set_functions(functions);
		}

	///	sets the velocity to a given value
		bool add(SmartPtr<IPData<MathVector<dim>, dim> > user, const char* subsetsBND)
		{
			if(m_velNames.empty() || m_pressureName.empty())
			{
				UG_LOG("ERROR in 'NavierStokesInflow::add': Symbolic names for"
						" velocity and pressure not set. Please set them first.\n");
				return false;
			}
			m_spNeumannDisc->add(user, m_pressureName.c_str(), subsetsBND);
			m_spDirichletConstraint->add(user, m_velNames.c_str(), subsetsBND);

			return true;
		}

	protected:
	///	sets the symbolic names for velocity and pressure
		void set_functions(const char* functions)
		{
			std::string strFunctions(functions);
			std::vector<std::string> tokens;

			TokenizeString(strFunctions, tokens, ',');

			if(tokens.size() != TDomain::dim + 1)
			{
				UG_THROW("ERROR in 'NavierStokesInflow::set_functions': This Boundary "
						"Condition works on exactly dim+1 (velocity+pressure) "
						"components, but "<<tokens.size()<<"components given.\n");
			}

			m_velNames.clear();
			for(int i=0;i<TDomain::dim; ++i)
			{
				if(i>0) m_velNames.append(",");
				m_velNames.append(tokens[i]);
			}

			m_pressureName = tokens[TDomain::dim];
		}

	///	neumann disc for pressure equation
		SmartPtr<NeumannBoundary<TDomain> > m_spNeumannDisc;

	///	dirichlet disc for velocity components
		SmartPtr<DirichletBoundary<TDomain,TAlgebra> > m_spDirichletConstraint;

	///	name of velocity components
		std::string m_velNames;

	///	name of pressure components
		std::string m_pressureName;
};

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
		NavierStokesWall(const char* functions)
			: m_spDirichletConstraint(new DirichletBoundary<TDomain,TAlgebra>)
		{
			set_functions(functions);
		}

	///	sets the velocity to a given value
		void add(const char* subsetsBND)
		{
			if(m_velNames.empty())
			{
				UG_THROW("ERROR in 'NavierStokesWall::add': Symbolic names for"
						" velocity and pressure not set. Please set them first.\n");
			}

			for(int i = 0; i < TDomain::dim; ++i)
			{
				m_spDirichletConstraint->add(0.0, m_velNames[i].c_str(), subsetsBND);
			}
		}

	protected:
	///	sets the symbolic names for velocity and pressure
		void set_functions(const char* functions)
		{
			std::string strFunctions(functions);

			m_velNames.clear();
			TokenizeString(strFunctions, m_velNames, ',');

			if(m_velNames.size() != TDomain::dim + 1)
			{
				UG_THROW("ERROR in 'NavierStokesWall::set_functions': This Boundary "
						"Condition works on exactly dim+1 (velocity+pressure) "
						"components, but "<<m_velNames.size()<<"components given.\n");
			}
		}

	///	dirichlet disc for velocity components
		SmartPtr<DirichletBoundary<TDomain,TAlgebra> > m_spDirichletConstraint;

	///	name of velocity components
		std::vector<std::string> m_velNames;
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__NAVIER_STOKES_BND__ */
