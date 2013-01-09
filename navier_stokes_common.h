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
  m_spStab(NULL),
  m_spConvStab(NULL),
  m_spConvUpwind(NULL)
{
	init();

	if(discType != NULL) set_disc_scheme(discType);
	else set_disc_scheme("stab");
}

template<typename TDomain>
NavierStokes<TDomain>::NavierStokes(const std::vector<std::string>& vFct,
                                    const std::vector<std::string>& vSubset,
                                    const char* discType)
: IDomainElemDisc<TDomain>(vFct, vSubset),
  m_bStokes(false),
  m_bLaplace(false),
  m_spStab(NULL),
  m_spConvStab(NULL),
  m_spConvUpwind(NULL)
{
	init();

	if(discType != NULL) set_disc_scheme(discType);
	else set_disc_scheme("stab");
}

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

//	register imports
	this->register_import(m_imSource);
	this->register_import(m_imKinViscosity);
	this->register_import(m_imDensitySCVF);
	this->register_import(m_imDensitySCV);

	m_imSource.set_rhs_part();
	m_imDensitySCV.set_mass_part();

//	default value for density
	set_density(1.0);
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
	if(	scheme != std::string("stab") &&
		scheme != std::string("staggered"))
	{
		UG_THROW("NavierStokes: Only 'stab' and 'staggered' supported.");
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
	return true;// todo check
	if(m_discScheme != "staggered"){
	//	check number
		if(vLfeID.size() != dim+1) return false;
		//	check that Lagrange 1st order
		for(size_t i = 0; i < vLfeID.size(); ++i)
			if(vLfeID[i] != LFEID(LFEID::LAGRANGE, 1)) return false;
			return true;
	}
	m_lfeID = vLfeID[0];

	//	update assemble functions
	set_ass_funcs();

	//	is supported
	return true;
}

template<typename TDomain>
void NavierStokes<TDomain>::
set_ass_funcs()
{
	//	switch, which assemble functions to use; both supported.
	if(m_discScheme == "stab") register_all_fv1_funcs(false);
	else if(m_discScheme == "staggered") register_all_cr_funcs(false);
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
set_source(number f_x)
{
	UG_THROW("NavierStokes: Setting source vector of dimension 1"
					" to a Discretization for world dim " << dim);
}

template<>
void NavierStokes<Domain1d>::
set_source(number f_x)
{
	SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
	f->set_entry(0, f_x);
	set_source(f);
}

template<typename TDomain>
void NavierStokes<TDomain>::
set_source(number f_x, number f_y)
{
	UG_THROW("NavierStokes: Setting source vector of dimension 2"
					" to a Discretization for world dim " << dim);
}

template<>
void NavierStokes<Domain2d>::
set_source(number f_x, number f_y)
{
	SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
	f->set_entry(0, f_x);
	f->set_entry(1, f_y);
	set_source(f);
}

template<typename TDomain>
void NavierStokes<TDomain>::
set_source(number f_x, number f_y, number f_z)
{
	UG_THROW("NavierStokes: Setting source vector of dimension 3"
					" to a Discretization for world dim " << dim);
}

template<>
void NavierStokes<Domain3d>::
set_source(number f_x, number f_y, number f_z)
{
	SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
	f->set_entry(0, f_x);
	f->set_entry(1, f_y);
	f->set_entry(2, f_z);
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

