
#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"

#include "navier_stokes_fvcr.h"

#include "pressure_separation.h"
#include "turbulent_viscosity_fvcr.h"

#include "bnd/inflow_fvcr.h"
#include "bnd/no_normal_stress_outflow_fvcr.h"

#include "bnd/symmetric_boundary_fvcr.h"
#include "cr_reorder.h"
#include "cr_ilut.h"

#include "lib_disc/function_spaces/grid_function.h"

using namespace std;
using namespace ug::bridge;

namespace ug{
namespace NavierStokes{

/**
 * Class exporting the functionality of the plugin. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct FunctionalityFVCR
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

	typedef ApproximationSpace<TDomain> approximation_space_type;
	typedef GridFunction<TDomain, TAlgebra> function_type;

	//	NavierStokesInflow FVCR
	{
		typedef NavierStokesInflowFVCR<TDomain, TAlgebra> T;
		typedef NavierStokesInflowBase<TDomain, TAlgebra> TBase;
		string name = string("NavierStokesInflowFVCR").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr< NavierStokesFVCR<TDomain> >)>("MasterElemDisc")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesInflowFVCR", tag);
	}

	//	Order CR-Cuthill-McKee
	{
		reg.add_function("OrderCRCuthillMcKee", static_cast<void (*)(approximation_space_type&,function_type&, bool)>(&OrderCRCuthillMcKee), grp);
	}
	
	//	Order CR-Cuthill-McKee
	{
		reg.add_function("CROrderCuthillMcKee", static_cast<void (*)(approximation_space_type&,function_type&, bool,bool,bool,bool)>(&CROrderCuthillMcKee), grp);
	}
	
	//	Order CR-Sloan
	{
		reg.add_function("CROrderSloan", static_cast<void (*)(approximation_space_type&,function_type&, bool,bool,bool)>(&CROrderSloan), grp);
	}
	
	//	Order CR-King
	{
		reg.add_function("CROrderKing", static_cast<void (*)(approximation_space_type&,function_type&, bool,bool,bool,bool)>(&CROrderKing), grp);
	}
	
	//	Order CR-Minimum Degree
	{
		reg.add_function("CROrderMinimumDegree", static_cast<void (*)(approximation_space_type&,function_type&, bool,bool,bool)>(&CROrderMinimumDegree), grp);
	}

	typedef ug::GridFunction<TDomain, TAlgebra> TFct;
	static const int dim = TDomain::dim;

	// Turbulent viscosity data
	// Smagorinsky model
	{
		string name = string("CRSmagorinskyTurbViscData").append(suffix);
		typedef CRSmagorinskyTurbViscData<TFct> T;
		typedef UserData<number, dim> TBase;
		typedef INewtonUpdate TBase2;
		reg.add_class_<T, TBase,TBase2>(name, grp)
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >,SmartPtr<TFct>,number)>("Approximation space, grid function, model parameter")
			.add_method("set_model_parameter", &T::set_model_parameter)
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(SmartPtr<UserData<number, dim> >)>(&T::set_kinematic_viscosity), "", "KinematicViscosity")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(number)>(&T::set_kinematic_viscosity), "", "KinematicViscosity")
		#ifdef UG_FOR_LUA
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(const char*)>(&T::set_kinematic_viscosity), "", "KinematicViscosity")
		#endif
			.add_method("set_turbulence_zero_bnd", &T::setTurbulenceZeroBoundaries)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "CRSmagorinskyTurbViscData", tag);
	}

	// Dynamic model
	{
		string name = string("CRDynamicTurbViscData").append(suffix);
		typedef CRDynamicTurbViscData<TFct> T;
		typedef UserData<number, dim> TBase;
		typedef INewtonUpdate TBase2;
		reg.add_class_<T, TBase,TBase2>(name, grp)
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >,SmartPtr<TFct>)>("Approximation space, grid function")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(SmartPtr<UserData<number, dim> >)>(&T::set_kinematic_viscosity), "", "KinematicViscosity")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(number)>(&T::set_kinematic_viscosity), "", "KinematicViscosity")
		#ifdef UG_FOR_LUA
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(const char*)>(&T::set_kinematic_viscosity), "", "KinematicViscosity")
		#endif
			.add_method("set_turbulence_zero_bnd", &T::setTurbulenceZeroBoundaries)
			.add_method("set_time_filter", &T::set_time_filter)
			.add_method("set_time_filter_eps", &T::set_time_filter_eps)
			.add_method("set_space_filter", &T::set_time_filter)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "CRDynamicTurbViscData", tag);
	}

	// SeparatedPressureSource
	{
		string name = string("SeparatedPressureSource").append(suffix);
		typedef SeparatedPressureSource<TFct> T;
		typedef UserData<MathVector<dim>, dim> TBase;
		typedef INewtonUpdate TBase2;
		reg.add_class_<T, TBase,TBase2>(name, grp)
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >,SmartPtr<TFct>)>("Approximation space, grid function")
				.add_method("set_source", static_cast<void (T::*)(SmartPtr<UserData<MathVector<dim>, dim> >)>(&T::set_source), "", "Source")
				.add_method("set_source", static_cast<void (T::*)(number)>(&T::set_source), "", "F_x")
				.add_method("set_source", static_cast<void (T::*)(number,number)>(&T::set_source), "", "F_x, F_y")
				.add_method("set_source", static_cast<void (T::*)(number,number,number)>(&T::set_source), "", "F_x, F_y, F_z")
			#ifdef UG_FOR_LUA
				.add_method("set_source", static_cast<void (T::*)(const char*)>(&T::set_source), "", "Source Vector")
			#endif
				.add_method("update", &T::update)
		.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "SeparatedPressureSource", tag);
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
	
	//	CR ILU Threshold
	{
		typedef CRILUTPreconditioner<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("CRILUT").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Incomplete LU Decomposition with threshold")
		.add_constructor()
		.add_method("set_threshold",static_cast<void (T::*)(number,number,number,number)>(&T::set_threshold),
					"", "threshold", "sets threshold of incomplete LU factorisation")
		.add_method("set_threshold",static_cast<void (T::*)(number)>(&T::set_threshold),
					"", "threshold", "sets threshold of incomplete LU factorisation")
		.add_method("set_info", &T::set_info,
					"", "info", "sets storage information output")
		.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "CRILUT", tag);
	}
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

	//	Navier-Stokes FVCR
	{
		typedef NavierStokesFVCR<TDomain> T;
		typedef NavierStokesBase<TDomain> TBase;
		string name = string("NavierStokesFVCR").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Functions#Subset(s)")
			.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Functions#Subset(s)")
			.add_method("set_upwind",  static_cast<void (T::*)(SmartPtr<INavierStokesUpwind<dim> >)>(&T::set_upwind))
			.add_method("set_upwind",  static_cast<void (T::*)(const std::string&)>(&T::set_upwind))
			.add_method("set_defect_upwind", &T::set_defect_upwind)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesFVCR", tag);
	}

	//	NavierStokesNoNormalStressOutflow FVCR
	{
		typedef NavierStokesNoNormalStressOutflowFVCR<TDomain> T;
		typedef NavierStokesNoNormalStressOutflowBase<TDomain> TBase;
		string name = string("NavierStokesNoNormalStressOutflowFVCR").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr< NavierStokesBase<TDomain> >)>("MasterDisc")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesNoNormalStressOutflowFVCR", tag);
	}

	//	CRNavierStokesSymBC
	{
		typedef CRNavierStokesSymBC<TDomain> T;
		typedef IDomainElemDisc<TDomain> TBase;
		string name = string("CRNavierStokesSymBC").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr< NavierStokesBase<TDomain> >)>("MasterDisc")
			.add_method("add", &T::add, "", "Subset(s)")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "CRNavierStokesSymBC", tag);
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
}

}; // end Functionality
} // end namespace NavierStokes


/**
 * This function is called when the plugin is loaded.
 */
void Init___NavierStokes___FVCR(Registry* reg, string grp)
{
	grp.append("SpatialDisc/NavierStokes/");
	typedef NavierStokes::FunctionalityFVCR Functionality;

	try{
//		RegisterDimensionDependent<Functionality>(*reg,grp);
		RegisterDomainDependent<Functionality>(*reg,grp);
		RegisterAlgebraDependent<Functionality>(*reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// namespace ug
