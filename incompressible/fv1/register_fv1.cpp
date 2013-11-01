
#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"

#include "navier_stokes_fv1.h"
#include "bnd/inflow_fv1.h"
#include "bnd/no_normal_stress_outflow_fv1.h"
#include "bnd/symmetric_boundary_fv1.h"
#include "stabilization.h"

#include "turbulent_viscosity_fv1.h"

#include "lib_disc/function_spaces/grid_function.h"

using namespace std;
using namespace ug::bridge;

namespace ug{
namespace NavierStokes{

/**
 * Class exporting the functionality of the plugin. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct FunctionalityFV1
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
	//static const int dim = TDomain::dim;
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

	//	NavierStokesInflow FV1
	{
		typedef NavierStokesInflowFV1<TDomain, TAlgebra> T;
		typedef NavierStokesInflowBase<TDomain, TAlgebra> TBase;
		string name = string("NavierStokesInflowFV1").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr< NavierStokesFV1<TDomain> >)>("MasterElemDisc")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesInflowFV1", tag);
	}
	
	typedef ug::GridFunction<TDomain, TAlgebra> TFct;
	static const int dim = TDomain::dim;
	
	// Turbulent viscosity data
	// Smagorinsky model
	{
		string name = string("FV1SmagorinskyTurbViscData").append(suffix);
		typedef FV1SmagorinskyTurbViscData<TFct> T;
		typedef CplUserData<number, dim> TBase;
		typedef INewtonUpdate TBase2;
		reg.add_class_<T, TBase,TBase2>(name, grp)
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >,SmartPtr<TFct>,number)>("Approximation space, grid function, model parameter")
			.add_method("set_model_parameter", &T::set_model_parameter)
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_kinematic_viscosity), "", "KinematicViscosity")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(number)>(&T::set_kinematic_viscosity), "", "KinematicViscosity")
		#ifdef UG_FOR_LUA
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(const char*)>(&T::set_kinematic_viscosity), "", "KinematicViscosity")
		#endif
			.add_method("set_turbulence_zero_bnd", &T::setTurbulenceZeroBoundaries)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "FV1SmagorinskyTurbViscData", tag);
	}

	// Dynamic model
	{
		string name = string("FV1DynamicTurbViscData").append(suffix);
		typedef FV1DynamicTurbViscData<TFct> T;
		typedef CplUserData<number, dim> TBase;
		typedef INewtonUpdate TBase2;
		reg.add_class_<T, TBase,TBase2>(name, grp)
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >,SmartPtr<TFct>)>("Approximation space, grid function")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_kinematic_viscosity), "", "KinematicViscosity")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(number)>(&T::set_kinematic_viscosity), "", "KinematicViscosity")
		#ifdef UG_FOR_LUA
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(const char*)>(&T::set_kinematic_viscosity), "", "KinematicViscosity")
		#endif
			.add_method("set_turbulence_zero_bnd", &T::setTurbulenceZeroBoundaries)
			.add_method("set_time_filter", &T::set_time_filter)
			.add_method("set_time_filter_eps", &T::set_time_filter_eps)
			.add_method("set_space_filter", &T::set_time_filter)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "FV1DynamicTurbViscData", tag);
	}

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
	static const int dim = TDomain::dim;
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();

	//	Navier-Stokes FV1
	{
		typedef NavierStokesFV1<TDomain> T;
		typedef IncompressibleNavierStokesBase<TDomain> TBase;
		string name = string("NavierStokesFV1").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Functions#Subset(s)")
			.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Functions#Subset(s)")
			.add_method("set_stabilization",  static_cast<void (T::*)(SmartPtr<INavierStokesFV1Stabilization<dim> >)>(&T::set_stabilization))
			.add_method("set_stabilization",  static_cast<void (T::*)(const std::string&)>(&T::set_stabilization))
			.add_method("set_stabilization",  static_cast<void (T::*)(const std::string&, const std::string&)>(&T::set_stabilization))
			.add_method("set_upwind",  static_cast<void (T::*)(SmartPtr<INavierStokesFV1Stabilization<dim> >)>(&T::set_upwind))
			.add_method("set_upwind",  static_cast<void (T::*)(SmartPtr<INavierStokesUpwind<dim> >)>(&T::set_upwind))
			.add_method("set_upwind",  static_cast<void (T::*)(const std::string&)>(&T::set_upwind))
			.add_method("set_pac_upwind", &T::set_pac_upwind, "", "Set pac upwind")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesFV1", tag);
	}

	//	NavierStokesNoNormalStressOutflow FV1
	{
		typedef NavierStokesNoNormalStressOutflowFV1<TDomain> T;
		typedef NavierStokesNoNormalStressOutflowBase<TDomain> TBase;
		string name = string("NavierStokesNoNormalStressOutflowFV1").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr< IncompressibleNavierStokesBase<TDomain> >)>("MasterDisc")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesNoNormalStressOutflowFV1", tag);
	}

	//	NavierStokesSymBCFV1
	{
		typedef NavierStokesSymBCFV1<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("NavierStokesSymBCFV1").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr< IncompressibleNavierStokesBase<TDomain> >)>("MasterDisc")
			.add_method("add", &T::add, "", "Subset(s)")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesSymBCFV1", tag);
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
{
	string suffix = GetDimensionSuffix<dim>();
	string tag = GetDimensionTag<dim>();

//	INavierStokesFV1Stabilization
	{
		typedef INavierStokesFV1Stabilization<dim> T;
		string name = string("INavierStokesFV1Stabilization").append(suffix);
		reg.add_class_<T>(name, grp)
			.add_method("set_upwind", &T::set_upwind)
			.add_method("set_diffusion_length", &T::set_diffusion_length);
		reg.add_class_to_group(name, "INavierStokesFV1Stabilization", tag);
	}

//	NavierStokesFIELDSStabilization
	{
		typedef NavierStokesFIELDSStabilization<dim> T;
		typedef INavierStokesFV1Stabilization<dim> TBase;
		string name = string("NavierStokesFIELDSStabilization").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesFIELDSStabilization", tag);
	}

//	NavierStokesFLOWStabilization
	{
		typedef NavierStokesFLOWStabilization<dim> T;
		typedef INavierStokesFV1Stabilization<dim> TBase;
		string name = string("NavierStokesFLOWStabilization").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesFLOWStabilization", tag);
	}

}

}; // end Functionality
} // end namespace NavierStokes


/**
 * This function is called when the plugin is loaded.
 */
void Init___NavierStokes___FV1(Registry* reg, string grp)
{
	grp.append("SpatialDisc/IncompressibleNavierStokes/");
	typedef NavierStokes::FunctionalityFV1 Functionality;

	try{
		RegisterDimension2d3dDependent<Functionality>(*reg,grp);
		RegisterDomain2d3dDependent<Functionality>(*reg,grp);
//		RegisterAlgebraDependent<Functionality>(*reg,grp);
		RegisterDomain2d3dAlgebraDependent<Functionality>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// namespace ug
