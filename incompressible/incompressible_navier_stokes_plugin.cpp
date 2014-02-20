/**
 * File for registration of navier stokes routines.
 *
 * created by Sebastian Reiter
 * s.b.reiter@googlemail.com
 * 14.09.2011 (m,d,y)
 */

#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"

#include "../register_navier_stokes.h"

#include "lib_disc/function_spaces/grid_function.h"
#include "lib_algebra/operator/debug_writer.h"

#include "incompressible_navier_stokes_base.h"
#include "navier_stokes_tools.h"
#include "filter.h"

#include "bnd/inflow_base.h"
#include "bnd/wall.h"
#include "bnd/no_normal_stress_outflow_base.h"

#include "fv1/register_fv1.h"
#include "fv/register_fv.h"
#include "fvcr/register_fvcr.h"
#include "fe/register_fe.h"

#include "transforming_smoother.h"


using namespace std;
using namespace ug::bridge;

namespace ug{
namespace NavierStokes{

/**
 * Class exporting the functionality of the plugin. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct FunctionalityIncomp
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
{
	static const int dim = TDomain::dim;
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

	typedef ApproximationSpace<TDomain> approximation_space_type;
	typedef GridFunction<TDomain, TAlgebra> function_type;

	//	NavierStokesInflowBase
	{
		typedef NavierStokesInflowBase<TDomain, TAlgebra> T;
		typedef IDiscretizationItem<TDomain, TAlgebra> TBase;
		string name = string("NavierStokesInflowBase").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
#ifdef UG_FOR_LUA
			.add_method("add", static_cast<void (T::*)(const char*, const char*)>(&T::add), "", "Velocity")
#endif
			.add_method("add", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >, const char*)>(&T::add), "", "Velocity, Subset")
			.add_method("add", static_cast<void (T::*)(const std::vector<number>&, const char*)>(&T::add), "", "Velocity, Subset");
		reg.add_class_to_group(name, "NavierStokesInflowBase", tag);
	}

//	NavierStokesWall
	{
		typedef NavierStokesWall<TDomain, TAlgebra> T;
		typedef IDiscretizationItem<TDomain, TAlgebra> TBase;
		string name = string("NavierStokesWall").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr< IncompressibleNavierStokesBase<TDomain> >)>("MasterElemDisc")
			.add_method("add", &T::add)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesWall", tag);
	}

//  ConstantTophatFilter
	{
		string name = string("ConstantBoxFilter").append(suffix);
		typedef ConstantBoxFilter<function_type> T;
		reg.add_class_<T>(name,grp)
			.template add_constructor<void (*)(number)>("filter width")
			.template add_constructor<void (*)(SmartPtr<function_type>,number)>("domain,element size scale")
			.add_method("apply",static_cast<void (T::*)(SmartPtr<function_type>)>(&T::apply), "", "apply filter");
		reg.add_class_to_group(name, "ConstantBoxFilter", tag);
	}
	
	// CFL number computation
	{
		reg.add_function("cflNumber",static_cast<void (*)(function_type&,number)>(&cflNumber),grp);
	}

	// vorticity computation
	{
		reg.add_function("vorticity", static_cast<void (*)(function_type&,function_type&)>(&vorticity), grp);
	}

	// kinetic energy
	{
		reg.add_function("kineticEnergy", static_cast<number (*)(function_type&)>(&kineticEnergy), grp);
	}

	// interpolate cr to lagrange 1
	{
		reg.add_function("CRToLagrange", static_cast<void (*)(function_type&,function_type&)>(&interpolateCRToLagrange), grp);
	}

	// driven cavity data evaluation
	{
		reg.add_function("DrivenCavityLinesEval", static_cast<void (*)(SmartPtr<function_type>, vector<string>,size_t)>(&DrivenCavityLinesEval<function_type>), grp);
	}

	// filter data
	{
		reg.add_function("filter", static_cast<void (*)(SmartPtr<function_type>,const std::string&)>(&filter), grp);
	}

	// filter data
	{
		reg.add_function("DragLift", &DragLift<function_type>, grp);
	}


//	AssembledTransformingSmoother
/*	{
		typedef AssembledTransformingSmoother<TDomain, TAlgebra> T;
		typedef ILinearIterator<typename TAlgebra::vector_type> TBase;
		typedef DebugWritingObject<TAlgebra> TBase2;
		string name = string("AssembledTransformingSmoother").append(suffix);
		reg.add_class_<T, TBase, TBase2>(name, grp)
			.template add_constructor<void (*)(SmartPtr<IAssemble<TAlgebra> >, SmartPtr<IAssemble<TAlgebra> >, SmartPtr<ILinearIterator<typename TAlgebra::vector_type> >)>("RightTrafoAss, TrafoSystemAss, Smoother")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "AssembledTransformingSmoother", tag);
	}*/
}

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
{
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();
}

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
	const int dim = TDomain::dim;
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();

//	Incompressible Navier-Stokes Base
	{
		typedef IncompressibleNavierStokesBase<TDomain> T;
		typedef NavierStokesBase<TDomain> TBase;
		string name = string("IncompressibleNavierStokesBase").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("set_density", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_density), "", "Density")
			.add_method("set_density", static_cast<void (T::*)(number)>(&T::set_density), "", "Density")
#ifdef UG_FOR_LUA
			.add_method("set_density", static_cast<void (T::*)(const char*)>(&T::set_density), "", "Density")
#endif
			.add_method("set_peclet_blend", &T::set_peclet_blend)
			.add_method("set_grad_div", static_cast<void (T::*)(number)>(&T::set_grad_div), "", "GradDivFactor")
			.add_method("set_laplace", &T::set_laplace)
			.add_method("set_stokes", &T::set_stokes);
		reg.add_class_to_group(name, "IncompressibleNavierStokesBase", tag);
	}

	//	NavierStokesNoNormalStressOutflowBase
	{
		typedef NavierStokesNoNormalStressOutflowBase<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("NavierStokesNoNormalStressOutflowBase").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("add", &T::add, "", "Subset(s)");
		reg.add_class_to_group(name, "NavierStokesNoNormalStressOutflowBase", tag);
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


static void Common(Registry& reg, string grp)
{
	// write numbers into file
	{
		reg.add_function("writeNumbers", static_cast<void (*)(std::string filename,const size_t step,const number t,const number data)>(&writeNumbers), grp);
	}
	// write numbers into file
	{
		reg.add_function("clearFile", static_cast<void (*)(std::string filename)>(&clearFile), grp);
	}
}

}; // end Functionality
} // end namespace NavierStokes


/**
 * This function is called when the plugin is loaded.
 */
extern "C" void
InitUGPlugin_IncompressibleNavierStokes(Registry* reg, string grp)
{
	grp.append("SpatialDisc/NavierStokes/");
	typedef NavierStokes::FunctionalityIncomp Functionality;

	try{
		//RegisterDimension2d3dDependent<Functionality>(*reg,grp);
		RegisterDomain2d3dDependent<Functionality>(*reg,grp);
//		RegisterAlgebraDependent<Functionality>(*reg,grp);
		RegisterDomain2d3dAlgebraDependent<Functionality>(*reg,grp);
		RegisterCommon<Functionality>(*reg,grp);

		Init___NavierStokes(reg, grp);
		Init___NavierStokes___FV1(reg, grp);
		Init___NavierStokes___FVCR(reg, grp);
		Init___NavierStokes___FV(reg, grp);
		Init___NavierStokes___FE(reg, grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}



}// namespace ug
