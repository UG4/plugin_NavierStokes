/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 *
 * This file is part of UG4.
 *
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 *
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 *
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 *
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifdef UG_PARALLEL
    #include "../../Parmetis/src/unificator_interface.h"
#endif

#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"

#include "common/log.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/dof_manager/dof_distribution.h"
#include "fv1/navier_stokes_fv1.h"
#include "incompressible_navier_stokes_base.h"
#include "fv1/bnd/inflow_fv1.h"
#include "bnd/inflow_base.h"


#include "lib_disc/spatial_disc/local_to_global/local_to_global_mapper.h"
#include "fv1/immersed_interface_base/gmg_transfer/particle_transfer.h"
#include "fv1/moving_particle/interface_handler_particle.h"
#include "fv1/moving_particle/moving_particle.h"
#include "fv1/diffusion_interface/interface_handler_diffusion.h"
#include "fv1/diffusion_interface/diffusion_interface.h"
//#include "fv1/two_phase_flow/two_phase_flow.h"


using namespace std;
using namespace ug::bridge;

namespace ug{
namespace ParticleLadenFlow{

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
 	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

	typedef ApproximationSpace<TDomain> approximation_space_type;
	typedef GridFunction<TDomain, TAlgebra> function_type;

//	IMovingInterfaceBase
    {
        typedef MovingInterfaceBase::IMovingInterface<TDomain, TAlgebra> T;
        string name = string("IMovingInterfaceBase").append(suffix);
        reg.add_class_<T>(name, grp);
        reg.add_class_to_group(name, "IMovingInterfaceBase", tag);
        
    }
/*
    //	MovingInterface2PF
    {
        typedef MovingInterface2PF::MovingInterface2PF<TDomain, TAlgebra> T;
        typedef MovingInterfaceBase::IMovingInterface<TDomain, TAlgebra> TBase;
        string name = string("MovingInterface2PF").append(suffix);
        reg.add_class_<T, TBase>(name, grp)
        .template add_constructor<void (*)(SmartPtr<IAssemble<TAlgebra> > ass,
                                           SmartPtr<NavierStokes::NavierStokesFV1<TDomain> > spMaster,															SmartPtr<MovingInterfaceBase::DiffusionInterfaceProvider<TDomain::dim> > interfaceProvider,
                                           SmartPtr<MovingInterfaceBase::CutElementHandlerImmersed<TDomain::dim> > cutElementHandler,
                                           number fluidDensity1, number fluidDensity2)>("domain disc, global handler")
        .add_method("init", &T::init)
        .add_method("get_integral", &T::get_integral)
        .add_method("adjust_for_error", &T::adjust_for_error)
        .add_method("initialize_threshold", &T::initialize_threshold)
        .add_method("set_threshold", &T::set_threshold, "", "Set Threshold")
        .add_method("set_analytic_solution", &T::set_analytic_solution, "", "Set Threshold")
        .set_construct_as_smart_pointer(true);
        reg.add_class_to_group(name, "MovingInterface2PF", tag);
    }
    
    //	MovingInterfaceDiffusionFE
    {
        typedef MovingInterfaceDiffusionFE::MovingInterfaceDiffusionFE<TDomain, TAlgebra> T;
        typedef MovingInterfaceBase::IMovingInterface<TDomain, TAlgebra> TBase;
        string name = string("MovingInterfaceDiffusionFE").append(suffix);
        reg.add_class_<T, TBase>(name, grp)
        .template add_constructor<void (*)(SmartPtr<IAssemble<TAlgebra> > ass,
                                           SmartPtr<ConvectionDiffusionPlugin::ConvectionDiffusionFE<TDomain> > spMaster,
                                           SmartPtr<MovingInterfaceBase::DiffusionInterfaceProvider<TDomain::dim> > interfaceProvider,
                                           SmartPtr<MovingInterfaceBase::CutElementHandlerImmersed<TDomain::dim> > cutElementHandler)>("domain disc, global handler")
        .add_method("init", &T::init)
        .set_construct_as_smart_pointer(true);
        reg.add_class_to_group(name, "MovingInterfaceDiffusionFE", tag);
    }
    */
    //	MovingInterfaceDiffusion
    {
        
        typedef MovingInterfaceDiffusion::MovingInterfaceDiffusion<TDomain, TAlgebra> T;
        typedef MovingInterfaceBase::IMovingInterface<TDomain, TAlgebra> TBase;
        string name = string("MovingInterfaceDiffusion").append(suffix);
        reg.add_class_<T, TBase>(name, grp)
        .template add_constructor<void (*)(SmartPtr<IAssemble<TAlgebra> > ass,
                                           SmartPtr<ConvectionDiffusionPlugin::ConvectionDiffusionFV1<TDomain> > spMaster,
                                           SmartPtr<MovingInterfaceBase::DiffusionInterfaceProvider<TDomain::dim> > interfaceProvider,
                                           SmartPtr<MovingInterfaceBase::CutElementHandlerImmersed<TDomain::dim> > cutElementHandler)>("domain disc, global handler")
        .add_method("init", &T::init)
        .add_method("get_integral", &T::get_integral)
        .add_method("get_numDoFs", &T::get_numDoFs)
        .add_method("set_Nitsche", &T::set_Nitsche)
        .add_method("adjust_for_error", &T::adjust_for_error)
        .add_method("initialize_threshold", &T::initialize_threshold)
        .add_method("set_threshold", &T::set_threshold, "", "Set Threshold")
        .add_method("set_analytic_solution", &T::set_analytic_solution, "", "Set Threshold")
        .set_construct_as_smart_pointer(true);
        reg.add_class_to_group(name, "MovingInterfaceDiffusion", tag);
    }
    
//	MovingParticle
    {
        typedef MovingParticle::MovingParticle<TDomain, TAlgebra> T;
        typedef MovingInterfaceBase::IMovingInterface<TDomain, TAlgebra> TBase;
        string name = string("MovingParticle").append(suffix);
        reg.add_class_<T, TBase>(name, grp)
        .template add_constructor<void (*)(SmartPtr<IAssemble<TAlgebra> > ass,
                                           SmartPtr<NavierStokes::NavierStokesFV1<TDomain> > spMaster,
                                           SmartPtr<MovingInterfaceBase::CutElementHandlerFlatTop<TDomain::dim> > cutElementHandler,
                                           number fluidDensity, number fluidKinVisc)>("domain disc, global handler")
        .add_method("init", &T::init)
     //   .add_method("compute_functional_combined", &T::compute_functional_combined)
        .add_method("compute_functional", &T::compute_functional)
        .add_method("compute_functional_all", &T::compute_functional_all)
        .add_method("gradient_descent", &T::gradient_descent)
        .add_method("project_directions", &T::project_directions)
        .add_method("rescale_directions", &T::rescale_directions)
        .add_method("print_velocity", &T::print_velocity)
     //   .add_method("print_velocity_many_particles", &T::print_velocity_many_particles)
        .add_method("print_pressure", &T::print_pressure)
        .add_method("print_deltaP", &T::print_deltaP)
        .add_method("print_pressure_nodal", &T::print_pressure_nodal)
   //     .add_method("compute_error_on_circle", &T::compute_error_on_circle)
        .add_method("update", &T::update)
        .add_method("get_velocity", &T::get_velocity)
        .add_method("adjust_global_solution", &T::adjust_global_solution)
        .add_method("compute_gradient_local_max", &T::compute_gradient_local_max)
        .add_method("set_gravity", &T::set_gravity)
        .add_method("set_time_step", &T::set_time_step)
        .add_method("set_volume_comp_mode", &T::set_volume_comp_mode)
        .add_method("set_StdFV_assembling", &T::set_StdFV_assembling)
        .add_method("initialize_threshold", &T::initialize_threshold)
        .add_method("get_BndCond", &T::get_BndCond)
//      .add_method("get_particles", &T::get_particles)
    // methods for parallel many-particle simulations:
#ifdef UG_PARALLEL
        .add_method("pre_balancing_update", &T::pre_balancing_update)
        .add_method("post_balancing_update", &T::post_balancing_update)
#endif
        .add_method("set_repulsive_force", &T::set_repulsive_force)
        .add_method("set_glowinski_repulsive_force", &T::set_glowinski_repulsive_force)
        .add_method("set_minimum_correction_force", &T::set_minimum_correction_force)
        .add_method("set_element_diameter", &T::set_element_diameter)
        .add_method("MeanElementDiameter", &T::MeanElementDiameter)
        .add_method("set_forceLog", &T::set_forceLog)
        .add_method("set_mpi_routine", &T::set_mpi_routine)
        .add_method("estimate_repulsive_force_parameters", &T::estimate_repulsive_force_parameters)
        .set_construct_as_smart_pointer(true);
        reg.add_class_to_group(name, "MovingParticle", tag);
    }
    
    //	ParticleTransfer
    {
        typedef ParticleTransfer<TDomain, TAlgebra> T;
        typedef ITransferOperator<TDomain, TAlgebra> TBase;
        string name = string("ParticleTransfer").append(suffix);
        reg.add_class_<T, TBase>(name, grp)
        .template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> > approxSpace, SmartPtr<MovingInterfaceBase::CutElementHandlerFlatTop<TDomain::dim> > cutElementHandler)>("approxSpace, globalHandler")
        .add_method("set_debug", &T::set_debug, "", "")
        .add_method("set_use_transposed", &T::set_use_transposed, "", "")
        .set_construct_as_smart_pointer(true);
        reg.add_class_to_group(name, "ParticleTransfer", tag);
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
 	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();

    // ParticleBndCond
    {
        typedef MovingParticle::ParticleBndCond<TDomain> T;
        typedef IElemDisc<TDomain> TBase;
        string name = string("ParticleBndCond").append(suffix);
        reg.add_class_<T, TBase>(name, grp),
        reg.add_class_to_group(name, "ParticleBndCond", tag);
    }
#ifdef UG_PARALLEL
    // ParticleUnificator
    {
        typedef MovingParticle::ParticleUnificator<TDomain> T;
        typedef typename GeomObjBaseTypeByDim<TDomain::dim>::base_obj_type TBaseElem;
        typedef parmetis::IUnificator<TBaseElem> TBase;
        string name = string("ParticleUnificator").append(suffix);
        reg.add_class_<T, TBase>(name, grp)
        .template add_constructor<void (*)(SmartPtr<TDomain>)>("")
        .add_method("update_particles", &T::update_particles)
        .add_method("rebalance", &T::rebalance)
        .set_construct_as_smart_pointer(true);
        //            .add_method("get_weight", &T::get_weight);
        //           .add_method("reweight", &T::reweight);
        reg.add_class_to_group(name, "ParticleUnificator", tag);
    }
#endif
    
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

    // CutElementHandlerFlatTop
    {
        typedef MovingInterfaceBase::CutElementHandlerFlatTop<dim> T;
        string name = string("CutElementHandlerFlatTop").append(suffix);
        reg.add_class_<T>(name, grp)
        .template add_constructor<void (*)(SmartPtr<MultiGrid> mg, const char*, SmartPtr<MovingParticle::ParticleProviderSphere<dim> >)>("multigrid, fct names")
        .template add_constructor<void (*)(SmartPtr<MultiGrid> mg, const char*, SmartPtr<MovingParticle::ParticleProviderEllipse<dim> >)>("multigrid, fct names")
   //     .add_method("update_prtCoords", &T::update_prtCoords)
        .set_construct_as_smart_pointer(true);
        reg.add_class_to_group(name, "CutElementHandlerFlatTop", tag);
    }
    
    // CutElementHandlerImmersed
    {
        typedef MovingInterfaceBase::CutElementHandlerImmersed<dim> T;
        string name = string("CutElementHandlerImmersed").append(suffix);
        reg.add_class_<T>(name, grp)
        .template add_constructor<void (*)(SmartPtr<MultiGrid> mg, const char*, SmartPtr<MovingInterfaceBase::DiffusionInterfaceProvider<dim> >)>("multigrid, fct names")
        .set_construct_as_smart_pointer(true);
        reg.add_class_to_group(name, "CutElementHandlerImmersed", tag);
    }
    
    //	ParticleProvider
    {
        typedef MovingParticle::ParticleProvider<dim> T;
        string name = string("ParticleProvider").append(suffix);
        reg.add_class_<T>(name, grp)
        .template add_constructor<void (*)( )>("")
        .add_method("print", &T::print)
        .set_construct_as_smart_pointer(true);
        reg.add_class_to_group(name, "ParticleProvider", tag);
    }
    
    
    //	ParticleProviderSphere
    {
        typedef MovingParticle::ParticleProviderSphere<dim> T;
        string name = string("ParticleProviderSphere").append(suffix);
        reg.add_class_<T>(name, grp)
        .template add_constructor<void (*)( )>("")
        .add_method("print", &T::print)
        .add_method("add", &T::add)
        .add_method("add_moving", &T::add_moving)
//      .add_method("get_collision_time", &T::get_collision_time)
        .set_construct_as_smart_pointer(true);
        reg.add_class_to_group(name, "ParticleProviderSphere", tag);
    }
    
    //	ParticleProviderEllipse
    {
        typedef MovingParticle::ParticleProviderEllipse<dim> T;
        string name = string("ParticleProviderEllipse").append(suffix);
        reg.add_class_<T>(name, grp)
        .template add_constructor<void (*)( )>("")
        .add_method("print", &T::print)
        .add_method("add", &T::add)
        .add_method("add_moving", &T::add_moving)
        .set_construct_as_smart_pointer(true);
        reg.add_class_to_group(name, "ParticleProviderEllipse", tag);
    }
    
    // DiffusionInterfaceProvider
    {
        typedef MovingInterfaceBase::DiffusionInterfaceProvider<dim> T;
        string name = string("DiffusionInterfaceProvider").append(suffix);
        reg.add_class_<T>(name, grp)
        .template add_constructor<void (*)( )>("")
        .add_method("print", &T::print)
        .add_method("add", &T::add)
        .set_construct_as_smart_pointer(true);
        reg.add_class_to_group(name, "DiffusionInterfaceProvider", tag);
    }
    
    
}


static void Common(Registry& reg, string grp)
{
	// write numbers into file
	/* TODO: obsolete / replace by lua code */
	{
		//reg.add_function("writeNumbers", static_cast<void (*)(std::string filename,const size_t step,const number t,const number data)>(&writeNumbers), grp);
	}
	// write numbers into file
	/* TODO: obsolete / replace by lua code */
	{
		//reg.add_function("clearFile", static_cast<void (*)(std::string filename)>(&clearFile), grp);
	}
}

}; // end Functionality
} // end namespace ParticleLadenFlow


/**
 * This function is called when the plugin is loaded.
 */
void Init___ParticleLadenFlow(Registry* reg, string grp)
{
	grp.append("SpatialDisc/NavierStokes/");
	typedef ParticleLadenFlow::FunctionalityIncomp Functionality;

	try{
        RegisterDimension2d3dDependent<Functionality>(*reg,grp);
        RegisterDomain2d3dDependent<Functionality>(*reg,grp);
        RegisterDomain2d3dAlgebraDependent<Functionality>(*reg,grp);

//      RegisterAlgebraDependent<Functionality>(*reg,grp);

	}
	UG_REGISTRY_CATCH_THROW(grp);
}



}// namespace ug


