// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 14.09.2011 (m,d,y)

// extern headers
#include <iostream>
#include <sstream>
#include <string>

#include "bridge/bridge.h"

#include "navier_stokes.h"
#include "upwind.h"
#include "stabilization.h"
#include "navier_stokes_bnd.h"
#include "no_normal_stress_outflow.h"

#include "lib_algebra/cpu_algebra_types.h"

using namespace std;

namespace ug{

using namespace ug::bridge;

template <typename TDomain, typename TAlgebra>
static void Register__Algebra_Domain(bridge::Registry& reg, string parentGroup)
{
//	typedef
	static const int dim = TDomain::dim;

//	group string
	stringstream grpSS; grpSS << parentGroup << "/" << dim << "d";
	string grp = grpSS.str();

//	suffix and tag
	string dimAlgSuffix = bridge::GetDomainSuffix<TDomain>();
	dimAlgSuffix.append(GetAlgebraSuffix<TAlgebra>());

	string dimAlgTag = GetDomainTag<TDomain>();
	dimAlgTag.append(GetAlgebraTag<TAlgebra>());

//	NavierStokesInflow
	{
		typedef NavierStokesInflow<TDomain, TAlgebra> T;
		typedef IDiscretizationItem<TDomain, TAlgebra> TBase;
		string name = string("NavierStokesInflow").append(dimAlgSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")

			.add_method("add", static_cast<void (T::*)(SmartPtr<IPData<MathVector<dim>, dim> >, const char*)>(&T::add), "", "Velocity, Subset")
			.add_method("add", static_cast<void (T::*)(number, const char*)>(&T::add), "", "Vel_x, Subset")
			.add_method("add", static_cast<void (T::*)(number,number, const char*)>(&T::add), "", "Vel_x, Vel_y, Subset")
			.add_method("add", static_cast<void (T::*)(number,number,number, const char*)>(&T::add), "", "Vel_x, Vel_y, Vel_z, Subset")
#ifdef UG_FOR_LUA
			.add_method("add", static_cast<void (T::*)(const char*, const char*)>(&T::add), "", "Velocity")
#endif

			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesInflow", dimAlgTag);
	}

//	NavierStokesWall
	{
		typedef NavierStokesWall<TDomain, TAlgebra> T;
		typedef IDiscretizationItem<TDomain, TAlgebra> TBase;
		string name = string("NavierStokesWall").append(dimAlgSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(const char*)>("Function(s)")
			.add_method("add", &T::add)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesWall", dimAlgTag);
	}

}

template <typename TDomain>
static void Register__Domain(bridge::Registry& reg, string grp)
{

//	dimension of domain
	static const int dim = TDomain::dim;

//	suffix and tag
	string dimSuffix = GetDomainSuffix<dim>();
	string dimTag = GetDomainTag<dim>();

//	Navier-Stokes
	{
		typedef NavierStokes<TDomain> T;
		typedef IDomainElemDisc<TDomain> TBase;
		string name = string("NavierStokes").append(dimSuffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Functions#Subset(s)")

			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(SmartPtr<IPData<number, dim> >)>(&T::set_kinematic_viscosity), "", "KinematicViscosity")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(number)>(&T::set_kinematic_viscosity), "", "KinematicViscosity")
#ifdef UG_FOR_LUA
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(const char*)>(&T::set_kinematic_viscosity), "", "KinematicViscosity")
#endif

			.add_method("set_source", static_cast<void (T::*)(SmartPtr<IPData<MathVector<dim>, dim> >)>(&T::set_source), "", "Source")
			.add_method("set_source", static_cast<void (T::*)(number)>(&T::set_source), "", "F_x")
			.add_method("set_source", static_cast<void (T::*)(number,number)>(&T::set_source), "", "F_x, F_y")
			.add_method("set_source", static_cast<void (T::*)(number,number,number)>(&T::set_source), "", "F_x, F_y, F_z")
#ifdef UG_FOR_LUA
			.add_method("set_source", static_cast<void (T::*)(const char*)>(&T::set_source), "", "Velocity Field")
#endif

			.add_method("set_stabilization", &T::set_stabilization)
			.add_method("set_conv_upwind",  static_cast<void (T::*)(SmartPtr<INavierStokesStabilization<dim> >)>(&T::set_conv_upwind))
			.add_method("set_conv_upwind",  static_cast<void (T::*)(SmartPtr<INavierStokesUpwind<dim> >)>(&T::set_conv_upwind))
			.add_method("set_peclet_blend", &T::set_peclet_blend)
			.add_method("set_exact_jacobian", &T::set_exact_jacobian)
			.add_method("set_laplace", &T::set_laplace)
			.add_method("set_stokes", &T::set_stokes)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokes", dimTag);
	}

//	FVNavierStokesNoNormalStressOutflow
	{
		typedef FVNavierStokesNoNormalStressOutflow<TDomain> T;
		typedef IDomainElemDisc<TDomain> TBase;
		string name = string("FVNavierStokesNoNormalStressOutflow").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr< NavierStokes<TDomain> >)>("MasterDisc")
			.add_method("add", &T::add, "", "Subset(s)")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "FVNavierStokesNoNormalStressOutflow", dimTag);
	}

/////////////////////////////////////////////////////////////////////////////
// Upwind
/////////////////////////////////////////////////////////////////////////////


//	INavierStokesUpwind
	{
		typedef INavierStokesUpwind<dim> T;
		string name = string("INavierStokesUpwind").append(dimSuffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "INavierStokesUpwind", dimTag);
	}

//	NavierStokesNoUpwind
	{
		typedef NavierStokesNoUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		string name = string("NavierStokesNoUpwind").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesNoUpwind", dimTag);
	}

//	NavierStokesFullUpwind
	{
		typedef NavierStokesFullUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		string name = string("NavierStokesFullUpwind").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesFullUpwind", dimTag);
	}

//	NavierStokesSkewedUpwind
	{
		typedef NavierStokesSkewedUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		string name = string("NavierStokesSkewedUpwind").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesSkewedUpwind", dimTag);
	}

//	NavierStokesLinearProfileSkewedUpwind
	{
		typedef NavierStokesLinearProfileSkewedUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		string name = string("NavierStokesLinearProfileSkewedUpwind").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesLinearProfileSkewedUpwind", dimTag);
	}

//	NavierStokesPositiveUpwind
	{
		typedef NavierStokesPositiveUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		string name = string("NavierStokesPositiveUpwind").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesPositiveUpwind", dimTag);
	}

//	NavierStokesRegularUpwind
	{
		typedef NavierStokesRegularUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		string name = string("NavierStokesRegularUpwind").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesRegularUpwind", dimTag);
	}

/////////////////////////////////////////////////////////////////////////////
// Stabilization
/////////////////////////////////////////////////////////////////////////////


//	INavierStokesStabilization
	{
		typedef INavierStokesStabilization<dim> T;
		string name = string("INavierStokesStabilization").append(dimSuffix);
		reg.add_class_<T>(name, grp)
			.add_method("set_upwind", &T::set_upwind)
			.add_method("set_diffusion_length", &T::set_diffusion_length);
		reg.add_class_to_group(name, "INavierStokesStabilization", dimTag);
	}

//	NavierStokesFIELDSStabilization
	{
		typedef NavierStokesFIELDSStabilization<dim> T;
		typedef INavierStokesStabilization<dim> TBase;
		string name = string("NavierStokesFIELDSStabilization").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesFIELDSStabilization", dimTag);
	}

//	NavierStokesFLOWStabilization
	{
		typedef NavierStokesFLOWStabilization<dim> T;
		typedef INavierStokesStabilization<dim> TBase;
		string name = string("NavierStokesFLOWStabilization").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NavierStokesFLOWStabilization", dimTag);
	}

}



template <typename TAlgebra>
static bool Register__Algebra(bridge::Registry& reg, string parentGroup)
{
//	get group string
	string grp = parentGroup; grp.append("/Discretization");

	try
	{

#ifdef UG_DIM_1
//	Domain dependent part 1D
	{
		typedef Domain<1, MultiGrid, MGSubsetHandler> domain_type;
		Register__Algebra_Domain<domain_type, TAlgebra>(reg, grp);
	}
#endif

#ifdef UG_DIM_2
//	Domain dependent part 2D
	{
		typedef Domain<2, MultiGrid, MGSubsetHandler> domain_type;
		Register__Algebra_Domain<domain_type, TAlgebra>(reg, grp);
	}
#endif

#ifdef UG_DIM_3
//	Domain dependent part 3D
	{
		typedef Domain<3, MultiGrid, MGSubsetHandler> domain_type;
		Register__Algebra_Domain<domain_type, TAlgebra>(reg, grp);
	}
#endif

	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in Register__Algebra_DoFDistribution: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}

	return true;
}

extern "C" void
InitUGPlugin_NavierStokesPlugin(ug::bridge::Registry* reg, std::string parentGroup)
{
	std::string grp(parentGroup); grp.append("NavierStokes/");

	bool bReturn = true;
#ifdef UG_CPU_1
	bReturn &= Register__Algebra<CPUAlgebra>(*reg, grp);
#endif
#ifdef UG_CPU_2
	bReturn &= Register__Algebra<CPUBlockAlgebra<2> >(*reg, grp);
#endif
#ifdef UG_CPU_3
	bReturn &= Register__Algebra<CPUBlockAlgebra<3> >(*reg, grp);
#endif
#ifdef UG_CPU_4
	bReturn &= Register__Algebra<CPUBlockAlgebra<4> >(*reg, grp);
#endif
#ifdef UG_CPU_VAR
	bReturn &= Register__Algebra<CPUVariableBlockAlgebra >(*reg, grp);
#endif

	try
	{
#ifdef UG_DIM_1
			Register__Domain<Domain1d>(*reg, grp);
#endif
#ifdef UG_DIM_2
			Register__Domain<Domain2d>(*reg, grp);
#endif
#ifdef UG_DIM_3
			Register__Domain<Domain3d>(*reg, grp);
#endif
	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterLibDisc_ElemDisc: "
				"Registration failed (using name " << ex.name << ").\n");
		return;
	}

	return;
}


}// end of namespace
