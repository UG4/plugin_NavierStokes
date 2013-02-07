/*
 * navier_stokes_common.h
 *
 *  Created on: 04.07.2012
 *      Author: Christian Wehner
 */

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

#include "navier_stokes.h"
#include <vector>
#include <string>

namespace ug{
namespace NavierStokes{

////////////////////////////////////////////////////////////////////////////////
//	Constructor - set default values
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
NavierStokes<TDomain>::NavierStokes(const char* functions,
                                    const char* subsets,
                                    const char* discType)
: IDomainElemDisc<TDomain>(functions, subsets),
  m_bStokes(false),
  m_bLaplace(false),
  m_bDefectUpwind(true),
  m_spStab(NULL),
  m_spConvStab(NULL),
  m_spConvUpwind(NULL)
{
	init();

	if(discType != NULL) set_disc_scheme(discType);
	else set_disc_scheme("fv1");
};

template<typename TDomain>
NavierStokes<TDomain>::NavierStokes(const std::vector<std::string>& vFct,
                                    const std::vector<std::string>& vSubset,
                                    const char* discType)
: IDomainElemDisc<TDomain>(vFct, vSubset),
  m_bStokes(false),
  m_bLaplace(false),
  m_bDefectUpwind(true),
  m_spStab(NULL),
  m_spConvStab(NULL),
  m_spConvUpwind(NULL)
{
	init();

	if(discType != NULL) set_disc_scheme(discType);
	else set_disc_scheme("fv1");
};

template<typename TDomain>
void NavierStokes<TDomain>::init()
{
//	check number of functions
	if(this->num_fct() != dim+1)
		UG_THROW("Wrong number of functions: The ElemDisc 'NavierStokes'"
					   " needs exactly "<<dim+1<<" symbolic function.");

//	set default options
	set_stokes(false);
	set_peclet_blend(false);
	set_exact_jacobian(false);
	set_defect_upwind(true);

//	register imports
	this->register_import(m_imSource);
	this->register_import(m_imKinViscosity);
	this->register_import(m_imDensitySCVF);
	this->register_import(m_imDensitySCVFp);
	this->register_import(m_imDensitySCV);

	m_imSource.set_rhs_part();
	m_imDensitySCV.set_mass_part();

//	default value for density
	set_density(1.0);

//	set defaults
	m_quadOrder = -1;
	m_discScheme = "fv1";
}

template<typename TDomain>
bool NavierStokes<TDomain>::request_non_regular_grid(bool bNonRegular)
{
//	switch, which assemble functions to use.
	if(bNonRegular)
	{
		UG_LOG("ERROR in 'NavierStokes::request_non_regular_grid':"
				" Non-regular grid not implemented.\n");
		return false;
	}

//	this disc supports regular grids
	return true;
}

template<typename TDomain>
void NavierStokes<TDomain>::
set_disc_scheme(const char* c_scheme)
{
	//	convert to string
	std::string scheme = c_scheme;

	//	check
	if(	scheme != std::string("fe") &&
		scheme != std::string("fv1") &&
		scheme != std::string("fv") &&
		scheme != std::string("staggered"))
	{
		UG_THROW("NavierStokes: Supported schemes: fe, fv1, fv, staggered.");
	}

	//	remember
	m_discScheme = scheme;

	//	update assemble functions
	set_ass_funcs();
}

template<typename TDomain>
bool NavierStokes<TDomain>::
request_finite_element_id(const std::vector<LFEID>& vLfeID)
{
//	check number
	if(vLfeID.size() != dim+1)
	{
		UG_LOG("NavierStokes:"
				" Wrong number of functions given. Need exactly "<<dim+1<<"\n");
		return false;
	}

	for(int d = 1; d < dim; ++d)
		if(vLfeID[0] != vLfeID[d])
		{
			UG_LOG("NavierStokes: trial spaces for velocity expected to be"
					" identical for all velocity components.\n");
			return false;
		}

//	staggered
	if(m_discScheme == "staggered"){
		for(int d = 0; d < dim; ++d)
			if(vLfeID[d].type() != LFEID::CROUZEIX_RAVIART)
			{
				UG_LOG("NavierStokes: 'staggered' expects Crouzeix-Raviart trial"
						" space for velocity.\n");
				return false;
			}
		if(vLfeID[dim].type() != LFEID::PIECEWISE_CONSTANT)
		{
			UG_LOG("NavierStokes: 'staggered' expects piecewise constant trial"
					" space for pressure.\n");
			return false;
		}
	}
//	fv1
	else if(m_discScheme == "fv1")
	{
		for(int d = 0; d <= dim; ++d)
			if(vLfeID[d].type() != LFEID::LAGRANGE || vLfeID[d].order() != 1)
			{
				UG_LOG("NavierStokes: 'fv1' expects Lagrange P1 trial space "
						"for velocity and pressure.\n");
				return false;
			}
	}
//	fv
	else if(m_discScheme == "fv")
	{
		for(int d = 0; d <= dim; ++d)
			if(vLfeID[d].type() != LFEID::LAGRANGE)
			{
				UG_LOG("NavierStokes: 'fv' expects Lagrange trial space "
						"for velocity and pressure.\n");
				return false;
			}
		for(int d = 0; d < dim; ++d)
			if(vLfeID[d].order() != vLfeID[dim].order() + 1)
			{
				UG_LOG("NavierStokes: 'fv' expects Lagrange trial space "
						"P_k for velocity and P_{k-1} pressure.\n");
				return false;
			}
	}
//	fe
	else if(m_discScheme == "fe")
	{
		// all admissible
	}
//	wrong scheme
	else {
		UG_LOG("NavierStokes: disc scheme '"<<m_discScheme<<"' not reconized.\n");
		return false;
	}

//	remember lfeID;
	m_vLFEID = vLfeID[0];
	m_pLFEID = vLfeID[dim];

//	set order
	m_vorder = vLfeID[0].order();
	m_porder = vLfeID[dim].order();

	//	update assemble functions
	set_ass_funcs();

	//	is supported
	return true;
}

template<typename TDomain>
void NavierStokes<TDomain>::
set_ass_funcs()
{
//	set default quadrature order if not set by user

	m_quadOrder = 2* m_vorder + 1;

	//	switch, which assemble functions to use; both supported.
	if(m_discScheme == "fv1") register_all_fv1_funcs(false);
	else if(m_discScheme == "fe") register_all_fe_funcs(m_vorder);
	else if(m_discScheme == "fv") register_all_fvho_funcs(m_vorder);
	else if(m_discScheme == "staggered") register_all_cr_funcs(false);
	else {
		UG_THROW("NavierStokes: disc scheme '"<<m_discScheme<<"' not reconized.\n");
	}

}

/////////// kinematic Viscosity

template<typename TDomain>
void NavierStokes<TDomain>::
set_kinematic_viscosity(SmartPtr<UserData<number, dim> > data)
{
	m_imKinViscosity.set_data(data);
}

template<typename TDomain>
void NavierStokes<TDomain>::
set_kinematic_viscosity(number val)
{
	set_kinematic_viscosity(CreateSmartPtr(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void NavierStokes<TDomain>::
set_kinematic_viscosity(const char* fctName)
{
	set_kinematic_viscosity(LuaUserDataFactory<number, dim>::create(fctName));
}
#endif

/////////// density

template<typename TDomain>
void NavierStokes<TDomain>::
set_density(SmartPtr<UserData<number, dim> > data)
{
	m_imDensitySCVF.set_data(data);
	m_imDensitySCVFp.set_data(data);
	m_imDensitySCV.set_data(data);
}

template<typename TDomain>
void NavierStokes<TDomain>::
set_density(number val)
{
	set_density(CreateSmartPtr(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void NavierStokes<TDomain>::
set_density(const char* fctName)
{
	set_density(LuaUserDataFactory<number, dim>::create(fctName));
}
#endif

/////////// Source

template<typename TDomain>
void NavierStokes<TDomain>::
set_source(SmartPtr<UserData<MathVector<dim>, dim> > data)
{
	m_imSource.set_data(data);
}

template<typename TDomain>
void NavierStokes<TDomain>::
set_source(const std::vector<number>& vSource)
{
	SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>(vSource));
	set_source(f);
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void NavierStokes<TDomain>::
set_source(const char* fctName)
{
	set_source(LuaUserDataFactory<MathVector<dim>, dim>::create(fctName));
}
#endif

} // namespace NavierStokes
} // end namespace ug

