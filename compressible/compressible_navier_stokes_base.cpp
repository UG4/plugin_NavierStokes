/*
 * compressible_navier_stokes_base.cpp
 *
 *  Created on: 31.10.2013
 *      Author: raphaelprohl
 */

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

#include "compressible_navier_stokes_base.h"
#include <vector>
#include <string>

namespace ug{
namespace NavierStokes{

////////////////////////////////////////////////////////////////////////////////
//	Constructor - set default values
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
CompressibleNavierStokesBase<TDomain>::CompressibleNavierStokesBase(const char* functions,
                                    const char* subsets)
: NavierStokesBase<TDomain>(functions, subsets),
  m_bMachNrBlend(false)
{
};

template<typename TDomain>
CompressibleNavierStokesBase<TDomain>::CompressibleNavierStokesBase(const std::vector<std::string>& vFct,
                                    const std::vector<std::string>& vSubset)
: NavierStokesBase<TDomain>(vFct, vSubset),
  m_bMachNrBlend(false)
{
};

/////////// adiabatic Index

template<typename TDomain>
void CompressibleNavierStokesBase<TDomain>::
set_adiabatic_index(number val)
{
	set_adiabatic_index(CreateSmartPtr(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void CompressibleNavierStokesBase<TDomain>::
set_adiabatic_index(const char* fctName)
{
	set_adiabatic_index(LuaUserDataFactory<number, dim>::create(fctName));
}
#endif


////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class CompressibleNavierStokesBase<Domain1d>;
#endif
#ifdef UG_DIM_2
template class CompressibleNavierStokesBase<Domain2d>;
#endif
#ifdef UG_DIM_3
template class CompressibleNavierStokesBase<Domain3d>;
#endif

} // namespace NavierStokes
} // end namespace ug




