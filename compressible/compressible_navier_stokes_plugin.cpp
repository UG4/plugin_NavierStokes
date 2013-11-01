/*
 * compressible_navier_stokes_plugin.cpp
 *
 *  Created on: 01.11.2013
 *      Author: raphaelprohl
 */


#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"

#include "../register_navier_stokes.h"

#include "compressible_navier_stokes_base.h"

#include "fv1/register_fv1.h"


using namespace std;
using namespace ug::bridge;

namespace ug{
namespace NavierStokes{

/**
 * Class exporting the functionality of the plugin. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct FunctionalityComp
{

/**
 * Function called for the registration of Domain and Algebra dependent parts
 * of the plugin. All Functions and Classes depending on both Domain and Algebra
 * are to be placed here when registering. The method is called for all
 * available Domain and Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain, typename TAlgebra>
static void DomainAlgebra(Registry& reg, string grp)
{}

/**
 * Function called for the registration of Algebra dependent parts.
 * All Functions and Classes depending on Algebra
 * are to be placed here when registering. The method is called for all
 * available Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TAlgebra>
static void Algebra(Registry& reg, string grp)
{}

/**
 * Function called for the registration of Domain dependent parts
 * of the plugin. All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain>
static void Domain(Registry& reg, string grp)
{
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();

//	Compressible Navier-Stokes Base
	{
		typedef CompressibleNavierStokesBase<TDomain> T;
		typedef NavierStokesBase<TDomain> TBase;
		string name = string("CompressibleNavierStokesBase").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("set_mach_number_blend", &T::set_mach_number_blend, "", "Set mach number blending");
		reg.add_class_to_group(name, "CompressibleNavierStokesBase", tag);
	}
}

/**
 * Function called for the registration of Dimension dependent parts
 * of the plugin. All Functions and Classes depending on the Dimension
 * are to be placed here when registering. The method is called for all
 * available Dimension types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <int dim>
static void Dimension(Registry& reg, string grp)
{}

}; // end Functionality
} // end namespace NavierStokes

/**
 * This function is called when the plugin is loaded.
 */
extern "C" void
InitUGPlugin_CompressibleNavierStokes(Registry* reg, string grp)
{
	grp.append("SpatialDisc/CompressibleNavierStokes/");
	typedef NavierStokes::FunctionalityComp Functionality;

	try{
		//RegisterDimension2d3dDependent<Functionality>(*reg,grp);
		RegisterDomain2d3dDependent<Functionality>(*reg,grp);
//		RegisterAlgebraDependent<Functionality>(*reg,grp);
		//RegisterDomain2d3dAlgebraDependent<Functionality>(*reg,grp);

		Init___NavierStokes(reg, grp);
		Init___CompressibleNavierStokes___FV1(reg, grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}


}// namespace ug
