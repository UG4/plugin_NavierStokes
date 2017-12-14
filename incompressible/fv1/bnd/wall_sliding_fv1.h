#ifndef __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FV1__BND__WALL_SLIDING___
#define __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FV1__BND__WALL_SLIDING___

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"

#include "../../incompressible_navier_stokes_base.h"
//#include "../stabilization.h"

namespace ug{
namespace NavierStokes{

template<	typename TDomain>
class NavierStokesWSBCFV1
	: public IElemDisc<TDomain>
{
	private:
	///	Base class type
		typedef IElemDisc<TDomain> base_type;

	///	own type
		typedef NavierStokesWSBCFV1<TDomain> this_type;

	public:
	///	Domain type
		typedef typename base_type::domain_type domain_type;

	///	World dimension
		static const int dim = base_type::dim;

	///	Position type
		typedef typename base_type::position_type position_type;

	public:
	///	Constructor (setting default values)
		NavierStokesWSBCFV1(SmartPtr< IncompressibleNavierStokesBase<TDomain> > spMaster);

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

	///	sets the wall sliding factor
	/**
	 * This method sets the wall sliding factor parameter.
	 *
	 * \param[in]	data		wall sliding factor
	 */
		void set_sliding_factor(number val);
#ifdef UG_FOR_LUA
		void set_sliding_factor(const char* fctName);
#endif

	///	sets the sliding limit
	/**
	 * This method sets the sliding limit parameter.
	 *
	 * \param[in]	data		sliding limit
	 */
		void set_sliding_limit(number val);
#ifdef UG_FOR_LUA
		void set_sliding_limit(const char* fctName);
#endif

	public:
	///	type of trial space for each function used
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

	public:
	///	returns if local time series is needed
		virtual bool requests_local_time_series() {return true;}

	///	prepares the element loop
		template <typename TElem, typename TFVGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);

	///	prepares the element for evaluation
		template <typename TElem, typename TFVGeom>
		void prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[]);

	///	finishes the element loop
		template <typename TElem, typename TFVGeom>
		void fsh_elem_loop();

	///	adds the stiffness part to the local jacobian
		template <typename TElem, typename TFVGeom>
		void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	adds the stiffness part to the local defect
		template <typename TElem, typename TFVGeom>
		void add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		template <typename TElem, typename TFVGeom>
		void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]) {}
		template <typename TElem, typename TFVGeom>
		void add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]) {}
		template <typename TElem, typename TFVGeom>
		void add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[]) {}

	private:
	/// The master discretization:
		SmartPtr< IncompressibleNavierStokesBase<TDomain> > m_spMaster;

	/// The boundary subsets:
		std::vector<std::string> m_vScheduledBndSubSets; // names
		std::vector<int> m_vBndSubSetIndex; // indices

		void extract_scheduled_data(); // convert m_vScheduledBndSubSets -> m_vBndSubSetIndex

	/// Data import for kinematic viscosity
		DataImport<number, dim> m_imKinViscosity;
	/// Data import for density
		DataImport<number, dim> m_imDensity;
	/// Data import for sliding factor
		DataImport<number, dim> m_imSlidingFactor;
	/// Data import for sliding limit
		DataImport<number, dim> m_imSlidingLimit;
	/// Boundary integration points of the viscosity and the density
		std::vector<MathVector<dim> > m_vLocIP;
		std::vector<MathVector<dim> > m_vGloIP;

	/// Data import for source
		DataImport<MathVector<dim>, dim> m_imSourceSCVF;

	/// abbreviation for pressure
		static const size_t _P_ = dim;

	protected:
		void register_all_funcs(bool bHang);
		template<typename TElem, typename TFVGeom>
		void register_func();
};

} // namespace NavierStokes
} // end namespace ug

#endif /* __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FV1__BND__WALL_SLIDING___ */
