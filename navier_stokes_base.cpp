/*
 * navier_stokes_base.cpp
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

#include "navier_stokes_base.h"
#include <vector>
#include <string>

namespace ug{
namespace NavierStokes{

////////////////////////////////////////////////////////////////////////////////
//	Constructor - set default values
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
NavierStokesBase<TDomain>::NavierStokesBase(const char* functions,
                                    const char* subsets)
: IDomainElemDisc<TDomain>(functions, subsets),
  m_bPecletBlend(false),
  m_bExactJacobian(false),
  m_bStokes(false),
  m_bLaplace(false)
{
};

template<typename TDomain>
NavierStokesBase<TDomain>::NavierStokesBase(const std::vector<std::string>& vFct,
                                    const std::vector<std::string>& vSubset)
: IDomainElemDisc<TDomain>(vFct, vSubset),
  m_bPecletBlend(false),
  m_bExactJacobian(false),
  m_bStokes(false),
  m_bLaplace(false)
{
};


/////////// kinematic Viscosity

template<typename TDomain>
void NavierStokesBase<TDomain>::
set_kinematic_viscosity(number val)
{
	set_kinematic_viscosity(CreateSmartPtr(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void NavierStokesBase<TDomain>::
set_kinematic_viscosity(const char* fctName)
{
	set_kinematic_viscosity(LuaUserDataFactory<number, dim>::create(fctName));
}
#endif

/////////// density

template<typename TDomain>
void NavierStokesBase<TDomain>::
set_density(number val)
{
	set_density(CreateSmartPtr(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void NavierStokesBase<TDomain>::
set_density(const char* fctName)
{
	set_density(LuaUserDataFactory<number, dim>::create(fctName));
}
#endif

/////////// Source

template<typename TDomain>
void NavierStokesBase<TDomain>::
set_source(const std::vector<number>& vSource)
{
	SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>(vSource));
	set_source(f);
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void NavierStokesBase<TDomain>::
set_source(const char* fctName)
{
	set_source(LuaUserDataFactory<MathVector<dim>, dim>::create(fctName));
}
#endif

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class NavierStokesBase<Domain1d>;
#endif
#ifdef UG_DIM_2
template class NavierStokesBase<Domain2d>;
#endif
#ifdef UG_DIM_3
template class NavierStokesBase<Domain3d>;
#endif

} // namespace NavierStokes
} // end namespace ug

