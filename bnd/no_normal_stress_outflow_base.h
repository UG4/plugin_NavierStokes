/*
 * no_normal_stress_outflow_base.h
 *
 *  Created on: 27.03.2012
 *  D. Logashenko, A. Vogel
 */

#ifndef __H__UG__PLUGINS__NAVIER_STOKES__BND__NO_NORMAL_STRESS_OUTFLOW_BASE__
#define __H__UG__PLUGINS__NAVIER_STOKES__BND__NO_NORMAL_STRESS_OUTFLOW_BASE__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"

#include "../navier_stokes_base.h"

namespace ug{
namespace NavierStokes{


/// \ingroup lib_disc_elem_disc
/// @{

/// The zero-stress (neutral) outflow boundary condition for the incompressible NS equation
/**
 * This class implements the so-called neutral boundary condition for outflow
 * boundaries. Note that this class can be used only with the stabilized
 * vertex-centered discretization of the Navier-Stokes equations.
 *
 * This boundary condition imposes two equations on the unknown functions
 * at the outflow boundary \f$ F \f$:
 * <ul>
 * <li> \f$ \int_F p ds = 0 \f$, and
 * <li> \f$ \sigma \vec{n} \cdot \vec{n} = 0 \f$ on \f$ F \f$.
 * </ul>
 * where \f$ \vec{n} \f$ is the outer normal at \f$ F \f$ and
 * \f$ \sigma = \mu (\nabla \vec{u} + (\nabla \vec{u})^T) \f$ the stress tensor.
 */
template<	typename TDomain>
class NavierStokesNoNormalStressOutflowBase
	: public IElemDisc<TDomain>
{
	protected:
	///	Base class type
		typedef IElemDisc<TDomain> base_type;

	///	own type
		typedef NavierStokesNoNormalStressOutflowBase<TDomain> this_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	public:
	///	Constructor (setting default values)
		NavierStokesNoNormalStressOutflowBase(SmartPtr< NavierStokesBase<TDomain> > spMaster);
	
	///	adds a boundary segment
		void add(const char* subsets);
	
	protected:
	///	sets the kinematic viscosity
		virtual void set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > data) = 0;

	///	sets the density
		virtual void set_density(SmartPtr<CplUserData<number, dim> > data) = 0;

	public:
	///	returns if local time series is needed
		virtual bool requests_local_time_series() {return true;}

	protected:
	/// The master discretization:
		SmartPtr< NavierStokesBase<TDomain> > m_spMaster;
	
	/// The boundary subsets:
		std::vector<std::string> m_vScheduledBndSubSets; // names
		std::vector<int> m_vBndSubSetIndex; // indices
	
		void extract_scheduled_data(); // convert m_vScheduledBndSubSets -> m_vBndSubSetIndex
};

/// @}

} // namespace NavierStokes
} // end namespace ug

#endif /*__H__UG__PLUGINS__NAVIER_STOKES__BND__NO_NORMAL_STRESS_OUTFLOW_BASE__*/
