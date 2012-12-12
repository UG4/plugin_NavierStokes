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
#include "navier_stokes_bnd.h"
#include "navier_stokes_cr_bnd.h"
#include "no_normal_stress_outflow.h"
#include "turbulent_viscosity_data.h"
#include "cr_order_cuthill_mckee.h"
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
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")

			.add_method("add", static_cast<void (T::*)(SmartPtr<UserData<MathVector<dim>, dim> >, const char*)>(&T::add), "", "Velocity, Subset")
			.add_method("add", static_cast<void (T::*)(number, const char*)>(&T::add), "", "Vel_x, Subset")
			.add_method("add", static_cast<void (T::*)(number,number, const char*)>(&T::add), "", "Vel_x, Vel_y, Subset")
			.add_method("add", static_cast<void (T::*)(number,number,number, const char*)>(&T::add), "", "Vel_x, Vel_y, Vel_z, Subset")
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
			.template add_constructor<void (*)(const char*)>("Function(s)")
			.add_method("add", &T::add)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesWall", tag);
	}
//	CRNavierStokesInflow
	{
		typedef CRNavierStokesInflow<TDomain, TAlgebra> T;
		typedef IDiscretizationItem<TDomain, TAlgebra> TBase;
		string name = string("CRNavierStokesInflow").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")

			.add_method("add", static_cast<void (T::*)(SmartPtr<UserData<MathVector<dim>, dim> >, const char*)>(&T::add), "", "Velocity, Subset")
			.add_method("add", static_cast<void (T::*)(number, const char*)>(&T::add), "", "Vel_x, Subset")
			.add_method("add", static_cast<void (T::*)(number,number, const char*)>(&T::add), "", "Vel_x, Vel_y, Subset")
			.add_method("add", static_cast<void (T::*)(number,number,number, const char*)>(&T::add), "", "Vel_x, Vel_y, Vel_z, Subset")
#ifdef UG_FOR_LUA
			.add_method("add", static_cast<void (T::*)(const char*, const char*)>(&T::add), "", "Velocity")
#endif

			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "CRNavierStokesInflow", tag);
	}

//	CRNavierStokesWall
	{
		typedef CRNavierStokesWall<TDomain, TAlgebra> T;
		typedef IDiscretizationItem<TDomain, TAlgebra> TBase;
		string name = string("CRNavierStokesWall").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(const char*)>("Function(s)")
			.add_method("add", &T::add)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "CRNavierStokesWall", tag);
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
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "CRDynamicTurbViscData", tag);
	}
	
	//	Order CR-Cuthill-McKee
	{
		reg.add_function("OrderCRCuthillMcKee", static_cast<void (*)(approximation_space_type&,function_type&, bool)>(&OrderCRCuthillMcKee), grp);
	}
	
	//	CR ILU Threshold
	{
		typedef CRILUTPreconditioner<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("CRILUT").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Incomplete LU Decomposition with threshold")
		.add_constructor()
		.add_method("set_threshold", &T::set_threshold,
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
			.add_method("set_source", static_cast<void (T::*)(number)>(&T::set_source), "", "F_x")
			.add_method("set_source", static_cast<void (T::*)(number,number)>(&T::set_source), "", "F_x, F_y")
			.add_method("set_source", static_cast<void (T::*)(number,number,number)>(&T::set_source), "", "F_x, F_y, F_z")
#ifdef UG_FOR_LUA
			.add_method("set_source", static_cast<void (T::*)(const char*)>(&T::set_source), "", "Velocity Field")
#endif

			.add_method("set_stabilization", &T::set_stabilization)
			.add_method("set_conv_upwind",  static_cast<void (T::*)(SmartPtr<INavierStokesFV1Stabilization<dim> >)>(&T::set_conv_upwind))
			.add_method("set_conv_upwind",  static_cast<void (T::*)(SmartPtr<INavierStokesFV1Upwind<dim> >)>(&T::set_conv_upwind))
			.add_method("set_conv_upwind",  static_cast<void (T::*)(SmartPtr<INavierStokesCRUpwind<dim> >)>(&T::set_conv_upwind))
			.add_method("set_peclet_blend", &T::set_peclet_blend)
			.add_method("set_exact_jacobian", &T::set_exact_jacobian)
			.add_method("set_laplace", &T::set_laplace)
			.add_method("set_stokes", &T::set_stokes)
			.add_method("set_disc_scheme", &T::set_disc_scheme)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokes", tag);
	}

//	FVNavierStokesNoNormalStressOutflow
	{
		typedef FVNavierStokesNoNormalStressOutflow<TDomain> T;
		typedef IDomainElemDisc<TDomain> TBase;
		string name = string("FVNavierStokesNoNormalStressOutflow").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr< NavierStokes<TDomain> >)>("MasterDisc")
			.add_method("add", &T::add, "", "Subset(s)")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "FVNavierStokesNoNormalStressOutflow", tag);
	}

//	CRNavierStokesNoNormalStressOutflow
	{
		typedef CRNavierStokesNoNormalStressOutflow<TDomain> T;
		typedef IDomainElemDisc<TDomain> TBase;
		string name = string("CRNavierStokesNoNormalStressOutflow").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr< NavierStokes<TDomain> >)>("MasterDisc")
			.add_method("add", &T::add, "", "Subset(s)")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "CRNavierStokesNoNormalStressOutflow", tag);
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

//	INavierStokesFV1Upwind
	{
		typedef INavierStokesFV1Upwind<dim> T;
		string name = string("INavierStokesFV1Upwind").append(suffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "INavierStokesFV1Upwind", tag);
	}

//	NavierStokesNoUpwind
	{
		typedef NavierStokesNoUpwind<dim> T;
		typedef INavierStokesFV1Upwind<dim> TBase;
		string name = string("NavierStokesNoUpwind").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesNoUpwind", tag);
	}

//	NavierStokesFullUpwind
	{
		typedef NavierStokesFullUpwind<dim> T;
		typedef INavierStokesFV1Upwind<dim> TBase;
		string name = string("NavierStokesFullUpwind").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesFullUpwind", tag);
	}

//	NavierStokesSkewedUpwind
	{
		typedef NavierStokesSkewedUpwind<dim> T;
		typedef INavierStokesFV1Upwind<dim> TBase;
		string name = string("NavierStokesSkewedUpwind").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesSkewedUpwind", tag);
	}

//	NavierStokesLinearProfileSkewedUpwind
	{
		typedef NavierStokesLinearProfileSkewedUpwind<dim> T;
		typedef INavierStokesFV1Upwind<dim> TBase;
		string name = string("NavierStokesLinearProfileSkewedUpwind").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesLinearProfileSkewedUpwind", tag);
	}

//	NavierStokesPositiveUpwind
	{
		typedef NavierStokesPositiveUpwind<dim> T;
		typedef INavierStokesFV1Upwind<dim> TBase;
		string name = string("NavierStokesPositiveUpwind").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesPositiveUpwind", tag);
	}

//	NavierStokesRegularUpwind
	{
		typedef NavierStokesRegularUpwind<dim> T;
		typedef INavierStokesFV1Upwind<dim> TBase;
		string name = string("NavierStokesRegularUpwind").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesRegularUpwind", tag);
	}
	
// CR staggered	
	
	//	INavierStokesCRUpwind
	{
		typedef INavierStokesCRUpwind<dim> T;
		string name = string("INavierStokesCRUpwind").append(suffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "INavierStokesCRUpwind", tag);
	}
	
	//	NavierStokesCRNoUpwind
	{
		typedef NavierStokesCRNoUpwind<dim> T;
		typedef INavierStokesCRUpwind<dim> TBase;
		string name = string("NavierStokesCRNoUpwind").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
		.add_constructor()
		.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesCRNoUpwind", tag);
	}
	
	//	NavierStokesCRFullUpwind
	{
		typedef NavierStokesCRFullUpwind<dim> T;
		typedef INavierStokesCRUpwind<dim> TBase;
		string name = string("NavierStokesCRFullUpwind").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
		.add_constructor()
		.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesCRFullUpwind", tag);
	}
	
	//	NavierStokesCRWeightedUpwind
	{
		typedef NavierStokesCRWeightedUpwind<dim> T;
		typedef INavierStokesCRUpwind<dim> TBase;
		string name = string("NavierStokesCRWeightedUpwind").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
		.template add_constructor<void (*)(number)>("Weight factor")
		.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesCRWeightedUpwind", tag);
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
		RegisterDomainAlgebraDependent<Functionality>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// namespace ug
