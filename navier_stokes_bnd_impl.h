/*
 * navier_stokes_bnd_impl.h
 *
 *  Created on: 12.04.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__NAVIER_STOKES_BND_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__NAVIER_STOKES_BND_IMPL__

#include "navier_stokes_bnd.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	NavierStokesInflow
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
NavierStokesInflow<TDomain,TAlgebra>::
NavierStokesInflow(const char* functions, const char* subsets)
	: m_spNeumannDisc(new NeumannBoundary<TDomain>(subsets)),
	  m_spDirichletConstraint(new DirichletBoundary<TDomain,TAlgebra>)
{
	set_functions(functions);
}


template <typename TDomain, typename TAlgebra>
void NavierStokesInflow<TDomain,TAlgebra>::
set_functions(const char* functions)
{
	std::string strFunctions(functions);
	std::vector<std::string> tokens;

	TokenizeString(strFunctions, tokens, ',');

	if(tokens.size() != TDomain::dim + 1)
		UG_THROW_FATAL("NavierStokesInflow::set_functions: This Boundary "
				"Condition works on exactly dim+1 (velocity+pressure) "
				"components, but "<<tokens.size()<<"components given.");

	m_velNames.clear();
	for(int i=0;i<TDomain::dim; ++i)
	{
		if(i>0) m_velNames.append(",");
		m_velNames.append(tokens[i]);
	}

	m_pressureName = tokens[TDomain::dim];
}

template <typename TDomain, typename TAlgebra>
void NavierStokesInflow<TDomain,TAlgebra>::
add(SmartPtr<IPData<MathVector<dim>, dim> > user, const char* subsetsBND)
{
	if(m_velNames.empty() || m_pressureName.empty())
		UG_THROW_FATAL("NavierStokesInflow::add: Symbolic names for"
				" velocity and pressure not set. Please set them first.");

	m_spNeumannDisc->add(user, m_pressureName.c_str(), subsetsBND);
	m_spDirichletConstraint->add(user, m_velNames.c_str(), subsetsBND);
}

template <typename TDomain, typename TAlgebra>
void NavierStokesInflow<TDomain,TAlgebra>::
add(number vel_x, const char* subsetsBND)
{
	if(dim != 1)
		UG_THROW_FATAL("NavierStokesInflow: Setting velocity vector of dimension 1"
						" to a Discretization for world dim " << dim);

	SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
	f->set_entry(0, vel_x);
	this->add(f, subsetsBND);
}

template <typename TDomain, typename TAlgebra>
void NavierStokesInflow<TDomain,TAlgebra>::
add(number vel_x, number vel_y, const char* subsetsBND)
{
	if(dim != 2)
		UG_THROW_FATAL("NavierStokesInflow: Setting velocity vector of dimension 2"
						" to a Discretization for world dim " << dim);

	SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
	f->set_entry(0, vel_x);
	f->set_entry(1, vel_y);
	add(f, subsetsBND);
}

template <typename TDomain, typename TAlgebra>
void NavierStokesInflow<TDomain,TAlgebra>::
add(number vel_x, number vel_y, number vel_z, const char* subsetsBND)
{
	if(dim != 3)
		UG_THROW_FATAL("NavierStokesInflow: Setting velocity vector of dimension 3"
						" to a Discretization for world dim " << dim);

	SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
	f->set_entry(0, vel_x);
	f->set_entry(1, vel_y);
	f->set_entry(2, vel_z);
	add(f, subsetsBND);
}


#ifdef UG_FOR_LUA
template <typename TDomain, typename TAlgebra>
void NavierStokesInflow<TDomain,TAlgebra>::
add(const char* fctName, const char* subsetsBND)
{
	add(CreateSmartPtr(new LuaUserData<MathVector<dim>, dim>(fctName)), subsetsBND);
}
#endif


////////////////////////////////////////////////////////////////////////////////
//	NavierStokesWall
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
NavierStokesWall<TDomain,TAlgebra>::
NavierStokesWall(const char* functions)
	: m_spDirichletConstraint(new DirichletBoundary<TDomain,TAlgebra>)
{
	set_functions(functions);
}

template <typename TDomain, typename TAlgebra>
void NavierStokesWall<TDomain,TAlgebra>::
add(const char* subsetsBND)
{
	if(m_velNames.empty())
		UG_THROW_FATAL("NavierStokesWall::add: Symbolic names for"
				" velocity and pressure not set. Please set them first.");

	for(int i = 0; i < TDomain::dim; ++i)
	{
		m_spDirichletConstraint->add(0.0, m_velNames[i].c_str(), subsetsBND);
	}
}

template <typename TDomain, typename TAlgebra>
void NavierStokesWall<TDomain,TAlgebra>::
set_functions(const char* functions)
{
	std::string strFunctions(functions);

	m_velNames.clear();
	TokenizeString(strFunctions, m_velNames, ',');

	if(m_velNames.size() != TDomain::dim + 1)
		UG_THROW_FATAL("NavierStokesWall::set_functions: This Boundary "
				"Condition works on exactly dim+1 (velocity+pressure) "
				"components, but "<<m_velNames.size()<<"components given.");
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__NAVIER_STOKES_BND_IMPL__ */
