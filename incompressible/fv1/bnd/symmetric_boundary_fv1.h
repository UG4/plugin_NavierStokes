/*
 * symmetric_boundary.h
 *
 *  Created on: 18.03.2013
 *      Author: Christian Wehner
 */

#ifndef __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FV1__BND__SYMMETRIC___
#define __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FV1__BND__SYMMETRIC___

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"

#include "../../incompressible_navier_stokes_base.h"
#include "../stabilization.h"


namespace ug{
namespace NavierStokes{
	
/// symmetric boundary condition
/// description, Naegele p. 55
template<	typename TDomain>
class NavierStokesSymBCFV1
	: public IElemDisc<TDomain>
{
	private:
	///	Base class type
		typedef IElemDisc<TDomain> base_type;

	///	own type
		typedef NavierStokesSymBCFV1<TDomain> this_type;

	public:
	///	Domain type
		typedef typename base_type::domain_type domain_type;

	///	World dimension
		static const int dim = base_type::dim;

	///	Position type
		typedef typename base_type::position_type position_type;

	public:
	///	Constructor (setting default values)
		NavierStokesSymBCFV1(SmartPtr< IncompressibleNavierStokesBase<TDomain> > spMaster);

	///	adds a boundary segment
		void add(const char* subsets);

	///	sets the kinematic viscosity
	/**
	 * This method sets the kinematic viscosity parameter.
	 *
	 * \param[in]	data		kinematic Viscosity
	 */
		void set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > data)
			{m_imKinViscosity.set_data(data);}

	///	sets the density
	/**
	 * This method sets the density parameter.
	 *
	 * \param[in]	data		Density
	 */
		void set_density(SmartPtr<CplUserData<number, dim> > data)
			{m_imDensity.set_data(data);}

	public:
	///	type of trial space for each function used
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

	public:
	///	returns if local time series is needed
		virtual bool requests_local_time_series() {return true;}

	///	prepares the element loop
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);

	///	prepares the element for evaluation
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		void prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[]);

	///	finishes the element loop
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		void fsh_elem_loop();

	///	adds the stiffness part to the local jacobian
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	adds the stiffness part to the local defect
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		void add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]) {}
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		void add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]) {}
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		void add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[]) {}

	private:
	/// The master discretization:
		SmartPtr< IncompressibleNavierStokesBase<TDomain> > m_spMaster;

	/// The boundary subsets:
		std::vector<std::string> m_vScheduledBndSubSets; // names
		std::vector<int> m_vBndSubSetIndex; // indices

		void extract_scheduled_data(); // convert m_vScheduledBndSubSets -> m_vBndSubSetIndex

	///	Data import for kinematic viscosity
		DataImport<number, dim> m_imKinViscosity;
	/// Data import for density
		DataImport<number, dim> m_imDensity;
	/// Boundary integration points of the viscosity and the density
		std::vector<MathVector<dim> > m_vLocIP;
		std::vector<MathVector<dim> > m_vGloIP;
		
	/// Data import for source
		DataImport<MathVector<dim>, dim> m_imSourceSCVF;

	/// abbreviation for pressure
		static const size_t _P_ = dim;

	private:
		void register_all_funcs(bool bHang);

		template <template <class Elem, int WorldDim> class TFVGeom>
		struct RegisterFV1 {
				RegisterFV1(this_type* pThis) : m_pThis(pThis){}
				this_type* m_pThis;
				template< typename TElem > void operator()(TElem&)
				{m_pThis->register_func<TElem, TFVGeom>();}
		};

		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		void register_func();
};

} // namespace NavierStokes
} // end namespace ug

#endif /* __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FV1__BND__SYMMETRIC___ */
