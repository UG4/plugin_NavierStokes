/*
 * incompressible_navier_stokes_base.cpp
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

#include "incompressible_navier_stokes_base.h"
#include <vector>
#include <string>

namespace ug{
namespace NavierStokes{

////////////////////////////////////////////////////////////////////////////////
//	Constructor - set default values
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
IncompressibleNavierStokesBase<TDomain>::IncompressibleNavierStokesBase(const char* functions,
                                    const char* subsets)
: NavierStokesBase<TDomain>(functions, subsets),
  m_bPecletBlend(false),
  m_bStokes(false),
  m_bLaplace(false),
  m_bBingham(false)
{
};

template<typename TDomain>
IncompressibleNavierStokesBase<TDomain>::IncompressibleNavierStokesBase(const std::vector<std::string>& vFct,
                                    const std::vector<std::string>& vSubset)
: NavierStokesBase<TDomain>(vFct, vSubset),
  m_bPecletBlend(false),
  m_bStokes(false),
  m_bLaplace(false),
  m_bBingham(false)
{
};


/////////// density

template<typename TDomain>
void IncompressibleNavierStokesBase<TDomain>::
set_density(number val)
{
	set_density(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void IncompressibleNavierStokesBase<TDomain>::
set_density(const char* fctName)
{
	set_density(LuaUserDataFactory<number, dim>::create(fctName));
}
#endif

/////////// bingham

template<typename TDomain>
void IncompressibleNavierStokesBase<TDomain>::
set_bingham_viscosity(number val)
{
  set_bingham_viscosity(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void IncompressibleNavierStokesBase<TDomain>::
set_bingham_viscosity(const char* fctName)
{
  set_bingham_viscosity(LuaUserDataFactory<number, dim>::create(fctName));
}
#endif


template<typename TDomain>
void IncompressibleNavierStokesBase<TDomain>::
set_yield_stress(number val)
{
  set_yield_stress(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void IncompressibleNavierStokesBase<TDomain>::
set_yield_stress(const char* fctName)
{
  set_yield_stress(LuaUserDataFactory<number, dim>::create(fctName));
}
#endif


////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class IncompressibleNavierStokesBase<Domain1d>;
#endif
#ifdef UG_DIM_2
template class IncompressibleNavierStokesBase<Domain2d>;
#endif
#ifdef UG_DIM_3
template class IncompressibleNavierStokesBase<Domain3d>;
#endif

} // namespace NavierStokes
} // end namespace ug



