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

#include "lib_algebra/cpu_algebra_types.h"

#include "lib_disc/dof_manager/conform/conform.h"
#include "lib_disc/dof_manager/p1conform/p1conform.h"

using namespace std;

namespace ug{

using namespace ug::bridge;

template <typename TDomain, typename TAlgebra, typename TDoFDistribution>
static void RegisterLibDiscDomain__Algebra_DoFDistribution_Domain(bridge::Registry& reg, string parentGroup)
{
//	typedef
	static const int dim = TDomain::dim;
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef ApproximationSpace<TDomain, TDoFDistribution, TAlgebra> approximation_space_type;

#ifdef UG_PARALLEL
		typedef ParallelGridFunction<GridFunction<TDomain, TDoFDistribution, TAlgebra> > function_type;
#else
		typedef GridFunction<TDomain, TDoFDistribution, TAlgebra> function_type;
#endif

//	group string
	stringstream grpSS; grpSS << parentGroup << "/" << dim << "d";
	string grp = grpSS.str();

//	suffix and tag
	string dimAlgDDSuffix = bridge::GetDomainSuffix<TDomain>();
	dimAlgDDSuffix.append(GetAlgebraSuffix<TAlgebra>());
	dimAlgDDSuffix.append(GetDoFDistributionSuffix<TDoFDistribution>());

	string dimAlgDDTag = GetDomainTag<TDomain>();
	dimAlgDDTag.append(GetAlgebraTag<TAlgebra>());
	dimAlgDDTag.append(GetDoFDistributionTag<TDoFDistribution>());

//	NavierStokesInflow
	{
		typedef boost::function<void (MathVector<dim>& value, const MathVector<dim>& x, number time)> VectorFunctor;
		typedef NavierStokesInflow<TDomain, TDoFDistribution, TAlgebra> T;
		typedef IDiscretizationItem<TDomain, TDoFDistribution, TAlgebra> TBase;
		string name = string("NavierStokesInflow").append(dimAlgDDSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.add_method("add", static_cast<bool (T::*)(VectorFunctor&, const char*)>(&T::add));
		reg.add_class_to_group(name, "NavierStokesInflow", dimAlgDDTag);
	}

//	NavierStokesInflow
	{
		typedef NavierStokesWall<TDomain, TDoFDistribution, TAlgebra> T;
		typedef IDiscretizationItem<TDomain, TDoFDistribution, TAlgebra> TBase;
		string name = string("NavierStokesWall").append(dimAlgDDSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(const char*)>("Function(s)")
			.add_method("add", &T::add);
		reg.add_class_to_group(name, "NavierStokesWall", dimAlgDDTag);
	}

}

template <typename TDomain>
static void RegisterIElemDiscs(bridge::Registry& reg, string grp)
{

//	dimension of domain
	static const int dim = TDomain::dim;

//	suffix and tag
	string dimSuffix = GetDomainSuffix<dim>();
	string dimTag = GetDomainTag<dim>();

//	Navier-Stokes
	{
		typedef FVNavierStokesElemDisc<TDomain> T;
		typedef IDomainElemDisc<TDomain> TBase;
		string name = string("FV1NavierStokes").append(dimSuffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Functions#Subset(s)")
			.add_method("set_kinematic_viscosity", &T::set_kinematic_viscosity)
			.add_method("set_stabilization", &T::set_stabilization)
			.add_method("set_conv_upwind",  static_cast<void (T::*)(INavierStokesStabilization<dim>&)>(&T::set_conv_upwind))
			.add_method("set_conv_upwind",  static_cast<void (T::*)(INavierStokesUpwind<dim>&)>(&T::set_conv_upwind))
			.add_method("set_peclet_blend", &T::set_peclet_blend)
			.add_method("set_exact_jacobian", &T::set_exact_jacobian);
		reg.add_class_to_group(name, "FV1NavierStokes", dimTag);
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
			.add_constructor();
		reg.add_class_to_group(name, "NavierStokesNoUpwind", dimTag);
	}

//	NavierStokesFullUpwind
	{
		typedef NavierStokesFullUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		string name = string("NavierStokesFullUpwind").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor();
		reg.add_class_to_group(name, "NavierStokesFullUpwind", dimTag);
	}

//	NavierStokesSkewedUpwind
	{
		typedef NavierStokesSkewedUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		string name = string("NavierStokesSkewedUpwind").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor();
		reg.add_class_to_group(name, "NavierStokesSkewedUpwind", dimTag);
	}

//	NavierStokesLinearProfileSkewedUpwind
	{
		typedef NavierStokesLinearProfileSkewedUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		string name = string("NavierStokesLinearProfileSkewedUpwind").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor();
		reg.add_class_to_group(name, "NavierStokesLinearProfileSkewedUpwind", dimTag);
	}

//	NavierStokesPositiveUpwind
	{
		typedef NavierStokesPositiveUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		string name = string("NavierStokesPositiveUpwind").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor();
		reg.add_class_to_group(name, "NavierStokesPositiveUpwind", dimTag);
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
			.add_constructor();
		reg.add_class_to_group(name, "NavierStokesFIELDSStabilization", dimTag);
	}

}


template <typename TAlgebra, typename TDoFDistribution>
static bool RegisterLibDiscDomain__Algebra_DoFDistribution(bridge::Registry& reg, string parentGroup)
{
//	get group string
	string grp = parentGroup; grp.append("/Discretization");

	try
	{

#ifdef UG_DIM_1
//	Domain dependent part 1D
	{
		typedef Domain<1, MultiGrid, MGSubsetHandler> domain_type;
		RegisterLibDiscDomain__Algebra_DoFDistribution_Domain<domain_type, TAlgebra, TDoFDistribution>(reg, grp);
	}
#endif

#ifdef UG_DIM_2
//	Domain dependent part 2D
	{
		typedef Domain<2, MultiGrid, MGSubsetHandler> domain_type;
		RegisterLibDiscDomain__Algebra_DoFDistribution_Domain<domain_type, TAlgebra, TDoFDistribution>(reg, grp);
	}
#endif

#ifdef UG_DIM_3
//	Domain dependent part 3D
	{
		typedef Domain<3, MultiGrid, MGSubsetHandler> domain_type;
		RegisterLibDiscDomain__Algebra_DoFDistribution_Domain<domain_type, TAlgebra, TDoFDistribution>(reg, grp);
	}
#endif

	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterLibDiscDomain__Algebra_DoFDistribution: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}

	return true;
}

template <typename TAlgebra>
static bool RegisterLibDiscDomain__Algebra(bridge::Registry& reg, string parentGroup)
{
	bool bReturn = true;
#ifdef DOF_P1
	bReturn &= RegisterLibDiscDomain__Algebra_DoFDistribution<TAlgebra, P1DoFDistribution>(reg, parentGroup);
#endif
#ifdef DOF_GEN
	bReturn &= RegisterLibDiscDomain__Algebra_DoFDistribution<TAlgebra, DoFDistribution >(reg, parentGroup);
#endif

	return bReturn;
}


extern "C" void InitUGPlugin(ug::bridge::Registry* reg, std::string parentGroup)
{
	std::string grp(parentGroup); grp.append("NavierStokes/");

	bool bReturn = true;
#ifdef UG_CPU_1
	bReturn &= RegisterLibDiscDomain__Algebra<CPUAlgebra>(*reg, grp);
#endif
#ifdef UG_CPU_2
	bReturn &= RegisterLibDiscDomain__Algebra<CPUBlockAlgebra<2> >(*reg, grp);
#endif
#ifdef UG_CPU_3
	bReturn &= RegisterLibDiscDomain__Algebra<CPUBlockAlgebra<3> >(*reg, grp);
#endif
#ifdef UG_CPU_4
	bReturn &= RegisterLibDiscDomain__Algebra<CPUBlockAlgebra<4> >(*reg, grp);
#endif
#ifdef UG_CPU_VAR
	bReturn &= RegisterLibDiscDomain__Algebra<CPUVariableBlockAlgebra >(*reg, grp);
#endif

	try
	{
#ifdef UG_DIM_1
	//	Domain dependend part 1D
			RegisterIElemDiscs<Domain1d>(*reg, grp);
#endif

#ifdef UG_DIM_2
	//	Domain dependend part 2D
			RegisterIElemDiscs<Domain2d>(*reg, grp);
#endif

#ifdef UG_DIM_3
	//	Domain dependend part 3D
			RegisterIElemDiscs<Domain3d>(*reg, grp);
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
