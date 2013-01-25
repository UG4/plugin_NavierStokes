/**
 * File for registration of navier stokes routines.
 *
 * created by Sebastian Reiter
 * s.b.reiter@googlemail.com
 * 14.09.2011 (m,d,y)
 */

#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"

#include "navier_stokes.h"
#include "upwind.h"
#include "stabilization.h"
#include "navier_stokes_tools.h"

#include "bnd/inflow.h"
#include "bnd/wall.h"
#include "bnd/no_normal_stress_outflow.h"
#include "bnd/symmetric_boundary.h"

#include "turbulent_viscosity_data.h"
#include "cr_reorder.h"
#include "cr_ilut.h"
#include "lib_disc/operator/non_linear_operator/newton_solver/newton.h"


using namespace std;
using namespace ug::bridge;

namespace ug{
namespace NavierStokes{

/**
 * Class exporting the functionality of the plugin. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
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
	typedef GridFunction<TDomain, SurfaceDoFDistribution, TAlgebra> function_type;

//	NavierStokesInflow
	{
		typedef NavierStokesInflow<TDomain, TAlgebra> T;
		typedef IDiscretizationItem<TDomain, TAlgebra> TBase;
		string name = string("NavierStokesInflow").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr< NavierStokes<TDomain> >)>("MasterElemDisc")

			.add_method("add", static_cast<void (T::*)(SmartPtr<UserData<MathVector<dim>, dim> >, const char*)>(&T::add), "", "Velocity, Subset")
			.add_method("add", static_cast<void (T::*)(const std::vector<number>&, const char*)>(&T::add), "", "Velocity, Subset")
#ifdef UG_FOR_LUA
			.add_method("add", static_cast<void (T::*)(const char*, const char*)>(&T::add), "", "Velocity")
#endif

			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesInflow", tag);
	}

//	NavierStokesWall
	{
		typedef NavierStokesWall<TDomain, TAlgebra> T;
		typedef IDiscretizationItem<TDomain, TAlgebra> TBase;
		string name = string("NavierStokesWall").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr< NavierStokes<TDomain> >)>("MasterElemDisc")
			.add_method("add", &T::add)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesWall", tag);
	}
	
	typedef ug::GridFunction<TDomain, SurfaceDoFDistribution, TAlgebra> TFct;
	
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
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "CRDynamicTurbViscData", tag);
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

	// vorticity computation
	{
		reg.add_function("vorticity", static_cast<void (*)(function_type&,function_type&)>(&vorticity), grp);
	}

	// driven cavity data evaluation
	{
		reg.add_function("dcevaluation", static_cast<void (*)(function_type&,size_t)>(&drivenCavityEvaluation), grp);
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

//	Navier-Stokes
	{
		typedef NavierStokes<TDomain> T;
		typedef IDomainElemDisc<TDomain> TBase;
		string name = string("NavierStokes").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Functions#Subset(s)")
			.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Functions#Subset(s)")
			.template add_constructor<void (*)(const char*,const char*, const char*)>("Functions#Subset(s)#DiscScheme")
			.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&, const char*)>("Functions#Subset(s)#DiscScheme")

			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(SmartPtr<UserData<number, dim> >)>(&T::set_kinematic_viscosity), "", "KinematicViscosity")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(number)>(&T::set_kinematic_viscosity), "", "KinematicViscosity")
#ifdef UG_FOR_LUA
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(const char*)>(&T::set_kinematic_viscosity), "", "KinematicViscosity")
#endif

			.add_method("set_density", static_cast<void (T::*)(SmartPtr<UserData<number, dim> >)>(&T::set_density), "", "Density")
			.add_method("set_density", static_cast<void (T::*)(number)>(&T::set_density), "", "Density")
#ifdef UG_FOR_LUA
			.add_method("set_density", static_cast<void (T::*)(const char*)>(&T::set_density), "", "Density")
#endif

			.add_method("set_source", static_cast<void (T::*)(SmartPtr<UserData<MathVector<dim>, dim> >)>(&T::set_source), "", "Source")
			.add_method("set_source", static_cast<void (T::*)(const std::vector<number>&)>(&T::set_source), "", "Source")
#ifdef UG_FOR_LUA
			.add_method("set_source", static_cast<void (T::*)(const char*)>(&T::set_source), "", "Source")
#endif

			.add_method("set_stabilization", &T::set_stabilization)
			.add_method("set_conv_upwind",  static_cast<void (T::*)(SmartPtr<INavierStokesFV1Stabilization<dim> >)>(&T::set_conv_upwind))
			.add_method("set_conv_upwind",  static_cast<void (T::*)(SmartPtr<INavierStokesUpwind<dim> >)>(&T::set_conv_upwind))
			.add_method("set_peclet_blend", &T::set_peclet_blend)
			.add_method("set_exact_jacobian", &T::set_exact_jacobian)
			.add_method("set_laplace", &T::set_laplace)
			.add_method("set_stokes", &T::set_stokes)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokes", tag);
	}

//	NavierStokesNoNormalStressOutflow
	{
		typedef NavierStokesNoNormalStressOutflow<TDomain> T;
		typedef IDomainElemDisc<TDomain> TBase;
		string name = string("NavierStokesNoNormalStressOutflow").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr< NavierStokes<TDomain> >)>("MasterDisc")
			.add_method("add", &T::add, "", "Subset(s)")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesNoNormalStressOutflow", tag);
	}

//	CRNavierStokesSymBC
	{
		typedef CRNavierStokesSymBC<TDomain> T;
		typedef IDomainElemDisc<TDomain> TBase;
		string name = string("CRNavierStokesSymBC").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr< NavierStokes<TDomain> >)>("MasterDisc")
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

/////////////////////////////////////////////////////////////////////////////
// Upwind
/////////////////////////////////////////////////////////////////////////////
	
//  FV1 collocated

//	INavierStokesUpwind
	{
		typedef INavierStokesUpwind<dim> T;
		string name = string("INavierStokesUpwind").append(suffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "INavierStokesUpwind", tag);
	}

//	NavierStokesNoUpwind
	{
		typedef NavierStokesNoUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		string name = string("NavierStokesNoUpwind").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesNoUpwind", tag);
	}

//	NavierStokesFullUpwind
	{
		typedef NavierStokesFullUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		string name = string("NavierStokesFullUpwind").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesFullUpwind", tag);
	}

//	NavierStokesSkewedUpwind
	{
		typedef NavierStokesSkewedUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		string name = string("NavierStokesSkewedUpwind").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesSkewedUpwind", tag);
	}

//	NavierStokesLinearProfileSkewedUpwind
	{
		typedef NavierStokesLinearProfileSkewedUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		string name = string("NavierStokesLinearProfileSkewedUpwind").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesLinearProfileSkewedUpwind", tag);
	}

//	NavierStokesPositiveUpwind
	{
		typedef NavierStokesPositiveUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		string name = string("NavierStokesPositiveUpwind").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesPositiveUpwind", tag);
	}

//	NavierStokesRegularUpwind
	{
		typedef NavierStokesRegularUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		string name = string("NavierStokesRegularUpwind").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesRegularUpwind", tag);
	}

//	NavierStokesWeightedUpwind
	{
		typedef NavierStokesWeightedUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		string name = string("NavierStokesWeightedUpwind").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
		.template add_constructor<void (*)(number)>("Weight factor")
		.add_method("set_weight", &T::set_weight, "", "")
		.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesWeightedUpwind", tag);
	}


/////////////////////////////////////////////////////////////////////////////
// Stabilization
/////////////////////////////////////////////////////////////////////////////


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
extern "C" void
InitUGPlugin_NavierStokes(Registry* reg, string grp)
{
	grp.append("/SpatialDisc/NavierStokes");
	typedef NavierStokes::Functionality Functionality;

	try{
		RegisterDimensionDependent<Functionality>(*reg,grp);
		RegisterDomainDependent<Functionality>(*reg,grp);
		RegisterAlgebraDependent<Functionality>(*reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// namespace ug
