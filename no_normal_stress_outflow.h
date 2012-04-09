/*
 * no_normal_stress_outflow.h
 *
 *  Created on: 27.03.2011
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__NO_NORMAL_STRESS_OUTFLOW__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__NO_NORMAL_STRESS_OUTFLOW__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/ip_data/data_import_export.h"

#include "upwind.h"
#include "stabilization.h"

#include "navier_stokes.h"

namespace ug{


/// \ingroup lib_disc_elem_disc
/// @{

template<	typename TDomain>
class FVNavierStokesNoNormalStressOutflow
	: public IDomainElemDisc<TDomain>
{
	private:
	///	Base class type
		typedef IDomainElemDisc<TDomain> base_type;

	///	own type
		typedef FVNavierStokesNoNormalStressOutflow<TDomain> this_type;

	public:
	///	Domain type
		typedef typename base_type::domain_type domain_type;

	///	World dimension
		static const int dim = base_type::dim;

	///	Position type
		typedef typename base_type::position_type position_type;

	public:
	///	Constructor (setting default values)
		FVNavierStokesNoNormalStressOutflow(SmartPtr< FV1NavierStokes<TDomain> > spMaster);
	
	///	adds a boundary segment
		void add(const char* subsets);
	
	///	sets the kinematic viscosity
	/**
	 * This method sets the kinematic viscosity value.
	 *
	 * \param[in]	data		kinematic Viscosity
	 */
		void set_kinematic_viscosity(SmartPtr<IPData<number, dim> > data)
			{m_imKinViscosity.set_data(data);}

	public:
	///	type of trial space for each function used
		virtual bool request_finite_element_id(const std::vector<LFEID>& vLfeID)
		{
		//	check number
			if(vLfeID.size() != dim+1) return false;

		//	check that Lagrange 1st order
			for(size_t i = 0; i < vLfeID.size(); ++i)
				if(vLfeID[i] != LFEID(LFEID::LAGRANGE, 1)) return false;
			return true;
		}

	///	switches between non-regular and regular grids
	/**
	 * \param[in]	bNonRegular		flag if non-regular grid needed.
	 */
		virtual bool request_non_regular_grid(bool bNonRegular)
		{
		//	switch, which assemble functions to use.
			if(bNonRegular)
			{
				UG_LOG("ERROR in 'FV1NavierStokes::request_non_regular_grid':"
						" Non-regular grid not implemented.\n");
				return false;
			}

		//	this disc supports regular grids
			return true;
		}

	public:
	///	returns if local time series is needed
		virtual bool requests_local_time_series() {return true;}

	///	prepares the element loop
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		void prepare_element_loop();

	///	prepares the element for evaluation
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		void prepare_element(TElem* elem, const LocalVector& u);

	///	finishes the element loop
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		void finish_element_loop();

	///	adds the stiffness part to the local jacobian
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		void ass_JA_elem(LocalMatrix& J, const LocalVector& u);

	///	adds the stiffness part to the local defect
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		void ass_dA_elem(LocalVector& d, const LocalVector& u);

	///	adds the mass part to the local jacobian
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		void ass_JM_elem(LocalMatrix& J, const LocalVector& u);

	///	adds the mass part to the local defect
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		void ass_dM_elem(LocalVector& d, const LocalVector& u);

	///	adds the source part to the local defect
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		void ass_rhs_elem(LocalVector& d);

	private:
	/// The master discretization:
		SmartPtr< FV1NavierStokes<TDomain> > m_spMaster;
	
	/// The boundary subsets:
		std::vector<std::string> m_vScheduledBndSubSets; // names
		std::vector<int> m_vBndSubSetIndex; // indices
	
		void extract_scheduled_data(); // converst m_vScheduledBndSubSets -> m_vBndSubSetIndex
		virtual void approximation_space_changed() {extract_scheduled_data();}
		
	///	Data import for kinematic viscosity
		DataImport<number, dim> m_imKinViscosity;

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

/// @}

} // end namespace ug

#endif /*__H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__NO_NORMAL_STRESS_OUTFLOW__*/
