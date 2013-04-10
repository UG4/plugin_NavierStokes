/*
 * symmetric_boundary_fv1.h
 *
 *  Created on: 18.03.2013
 *      Author: Christian Wehner
 */

#ifndef __H__UG__PLUGINS__NAVIER_STOKES__SYMMETRIC__FV1__
#define __H__UG__PLUGINS__NAVIER_STOKES__SYMMETRIC__FV1__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"

#include "../../navier_stokes_base.h"
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
		NavierStokesSymBCFV1(SmartPtr< NavierStokesBase<TDomain> > spMaster);

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
		virtual bool request_finite_element_id(const std::vector<LFEID>& vLfeID);

	///	switches between non-regular and regular grids
		virtual bool request_non_regular_grid(bool bNonRegular);

	public:
	///	returns if local time series is needed
		virtual bool requests_local_time_series() {return true;}

	///	prepares the element loop
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		void prep_elem_loop_fv1(const ReferenceObjectID roid, const int si);

	///	prepares the element for evaluation
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		void prep_elem_fv1(TElem* elem, const LocalVector& u);

	///	finishes the element loop
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		void fsh_elem_loop_fv1();

	///	adds the stiffness part to the local jacobian
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		void add_JA_elem_fv1(LocalMatrix& J, const LocalVector& u);

	///	adds the stiffness part to the local defect
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		void add_dA_elem_fv1(LocalVector& d, const LocalVector& u);

		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		void add_JM_elem(LocalMatrix& J, const LocalVector& u) {}
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		void add_dM_elem(LocalVector& d, const LocalVector& u) {}
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		void add_rhs_elem(LocalVector& d) {}

	private:
	/// The master discretization:
		SmartPtr< NavierStokesBase<TDomain> > m_spMaster;

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

	/// position access
		const position_type* m_vCornerCoords;

	/// abbreviation for pressure
		static const size_t _P_ = dim;

	private:
		void register_all_fv1_funcs(bool bHang);

		template <template <class Elem, int WorldDim> class TFVGeom>
		struct RegisterFV1 {
				RegisterFV1(this_type* pThis) : m_pThis(pThis){}
				this_type* m_pThis;
				template< typename TElem > void operator()(TElem&)
				{m_pThis->register_fv1_func<TElem, TFVGeom>();}
		};

		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		void register_fv1_func();
};

} // namespace NavierStokes
} // end namespace ug

#endif /* __H__UG__PLUGINS__NAVIER_STOKES__SYMMETRIC__FV1__ */
