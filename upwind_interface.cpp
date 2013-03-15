/*
 * upwind.cpp
 *
 *  Created on: 10.03.2011
 *      Author: andreasvogel
 */

// for minimum
#include <limits>
#include <algorithm>
#include <locale>

#include "upwind.h"

namespace ug {
namespace NavierStokes{

template <int dim>
SmartPtr<INavierStokesUpwind<dim> > CreateNavierStokesUpwind(const std::string& name)
{
	std::string n = TrimString(name);
	std::transform(n.begin(), n.end(), n.begin(), ::tolower);

	if(n == "no") return SmartPtr<INavierStokesUpwind<dim> >(new NavierStokesNoUpwind<dim>());
	if(n == "full") return SmartPtr<INavierStokesUpwind<dim> >(new NavierStokesFullUpwind<dim>());
	if(n == "weighted") return SmartPtr<INavierStokesUpwind<dim> >(new NavierStokesWeightedUpwind<dim>());
	if(n == "skewed") return SmartPtr<INavierStokesUpwind<dim> >(new NavierStokesSkewedUpwind<dim>());
	if(n == "linearprofileskewed" ||
	   n == "lps") return SmartPtr<INavierStokesUpwind<dim> >(new NavierStokesLinearProfileSkewedUpwind<dim>());
	if(n == "positive" ||
	   n == "pos") return SmartPtr<INavierStokesUpwind<dim> >(new NavierStokesPositiveUpwind<dim>());
	if(n == "regular" ||
	   n == "reg") return SmartPtr<INavierStokesUpwind<dim> >(new NavierStokesRegularUpwind<dim>());

	UG_THROW("NavierStokes: upwind type '"<<name<<"' not found. Options are: "
	         "no, full, weighted, skewed, LinearProfileSkewed (lps), positive (pos), regular (reg)");
}

////////////////////////////////////////////////////////////////////////////////
//	explicit instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_2
template SmartPtr<INavierStokesUpwind<2> >CreateNavierStokesUpwind<2>(const std::string& name);
#endif
#ifdef UG_DIM_3
template SmartPtr<INavierStokesUpwind<3> >CreateNavierStokesUpwind<3>(const std::string& name);
#endif

} // namespace NavierStokes
} // end namespace ug
