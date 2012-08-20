/*
 * navier_stokes_cr_bnd.h
 *
 *  Created on: 23.07.2012
 *      Author: Christian Wehner
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__NAVIER_STOKES_CR_BND__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__NAVIER_STOKES_CR_BND__

#include "lib_disc/spatial_disc/disc_item.h"

#include "lib_disc/spatial_disc/elem_disc/neumann_boundary/neumann_boundary.h"
#include "lib_disc/spatial_disc/constraints/dirichlet_boundary/lagrange_dirichlet_boundary.h"

namespace ug{
namespace NavierStokes{

template <	typename TDomain,
			typename TAlgebra>
class CRNavierStokesInflow
	: public IDiscretizationItem<TDomain, TAlgebra>
{
	private:
		const static int dim = TDomain::dim;

	public:
	///	returns the number of element discs
		virtual size_t num_elem_disc() const {return 0;}

	///	returns the element disc
		virtual SmartPtr<IDomainElemDisc<TDomain> > elem_disc(size_t i) {return NULL;}

	///	returns the number of constraints
		virtual size_t num_constraint() const {return 1;}

	///	returns an element disc
		virtual SmartPtr<IDomainConstraint<TDomain, TAlgebra> > constraint(size_t i) {return m_spDirichletConstraint;}

	public:
	///	Constructor
		CRNavierStokesInflow(const char* functions, const char* subsets);

	///	sets the velocity to a given value
	///	\{
		void add(SmartPtr<UserData<MathVector<dim>, dim> > user, const char* subsetsBND);
		void add(number vel_x, const char* subsetsBND);
		void add(number vel_x, number vel_y, const char* subsetsBND);
		void add(number vel_x, number vel_y, number vel_z, const char* subsetsBND);
#ifdef UG_FOR_LUA
		void add(const char* fctName, const char* subsetsBND);
#endif
	///	\}

	protected:
	///	sets the symbolic names for velocity
		void set_functions(const char* functions);

	///	dirichlet disc for velocity components
		SmartPtr<DirichletBoundary<TDomain,TAlgebra> > m_spDirichletConstraint;

	///	name of velocity components
		std::string m_velNames;

};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <	typename TDomain,
			typename TAlgebra>
class CRNavierStokesWall
	: public IDiscretizationItem<TDomain,TAlgebra>
{
	public:
	///	returns the number of element discs
		virtual size_t num_elem_disc() const {return 0;}

	///	returns the element disc
		virtual SmartPtr<IDomainElemDisc<TDomain> > elem_disc(size_t i) {return NULL;}

	///	returns the number of constraints
		virtual size_t num_constraint() const {return 1;}

	///	returns an element disc
		virtual SmartPtr<IDomainConstraint<TDomain, TAlgebra> > constraint(size_t i) {return m_spDirichletConstraint;}

	///	virtual destructor
		~CRNavierStokesWall() {}

	public:
	///	Constructor
		CRNavierStokesWall(const char* functions);

	///	sets the velocity to a given value
		void add(const char* subsetsBND);

	protected:
	///	sets the symbolic names for velocity
		void set_functions(const char* functions);

	///	dirichlet disc for velocity components
		SmartPtr<DirichletBoundary<TDomain,TAlgebra> > m_spDirichletConstraint;

	///	name of velocity components
		std::vector<std::string> m_velNames;
};



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
	class CRNavierStokesNoNormalStressOutflow
		: public IDomainElemDisc<TDomain>
	{
		private:
		///	Base class type
			typedef IDomainElemDisc<TDomain> base_type;

		///	own type
			typedef CRNavierStokesNoNormalStressOutflow<TDomain> this_type;

		public:
		///	Domain type
			typedef typename base_type::domain_type domain_type;

		///	World dimension
			static const int dim = base_type::dim;

		///	Position type
			typedef typename base_type::position_type position_type;

		public:
		///	Constructor (setting default values)
			CRNavierStokesNoNormalStressOutflow(SmartPtr< NavierStokes<TDomain> > spMaster);

		///	adds a boundary segment
			void add(const char* subsets);

		///	sets the kinematic viscosity
		/**
		 * This method sets the kinematic viscosity parameter.
		 *
		 * \param[in]	data		kinematic Viscosity
		 */
			void set_kinematic_viscosity(SmartPtr<UserData<number, dim> > data)
				{m_imKinViscosity.set_data(data);}

		///	sets the density
		/**
		 * This method sets the density parameter.
		 *
		 * \param[in]	data		Density
		 */
			void set_density(SmartPtr<UserData<number, dim> > data)
				{m_imDensity.set_data(data);}

		public:
		///	type of trial space for each function used
			virtual bool request_finite_element_id(const std::vector<LFEID>& vLfeID)
			{
			//	check number
				if(vLfeID.size() != dim+1) return false;

			//	check that Lagrange 1st order
			//	for(size_t i = 0; i < vLfeID.size(); ++i)
			//		if(vLfeID[i] != LFEID(LFEID::LAGRANGE, 1)) return false;
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
					UG_LOG("ERROR in 'NavierStokes::request_non_regular_grid':"
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
		/// adds the diffusive part of the local Jacobian of the momentum equation
			template <typename BF>
			inline void ass_diffusive_flux_Jac
			(
				const size_t ip,
				const BF& bf,
				LocalMatrix& J,
				const LocalVector& u
			);
		/// adds the diffusive part of the local defect of the momentum equation
			template <typename BF>
			inline void ass_diffusive_flux_defect
			(
				const size_t ip,
				const BF& bf,
				LocalVector& d,
				const LocalVector& u
			);
		/// adds the convective part of the local Jacobian of the momentum equation
			template <typename BF>
			inline void ass_convective_flux_Jac
			(
				const size_t ip,
				const BF& bf,
				LocalMatrix& J,
				const LocalVector& u
			);
		/// adds the convective part of the local defect of the momentum equation
			template <typename BF>
			inline void ass_convective_flux_defect
			(
				const size_t ip,
				const BF& bf,
				LocalVector& d,
				const LocalVector& u
			);

		/// The master discretization:
			SmartPtr< NavierStokes<TDomain> > m_spMaster;

		/// The boundary subsets:
			std::vector<std::string> m_vScheduledBndSubSets; // names
			std::vector<int> m_vBndSubSetIndex; // indices

			void extract_scheduled_data(); // converst m_vScheduledBndSubSets -> m_vBndSubSetIndex
			virtual void approximation_space_changed() {extract_scheduled_data();}

		///	Data import for kinematic viscosity
			DataImport<number, dim> m_imKinViscosity;
		/// Data import for density
			DataImport<number, dim> m_imDensity;
		/// Boundary integration points of the viscosity and the density
			std::vector<MathVector<dim> > m_vLocIP;
			std::vector<MathVector<dim> > m_vGloIP;

		/// position access
			const position_type* m_vCornerCoords;

		/// abbreviation for pressure
			static const size_t _P_ = dim;

		private:
			void register_all_cr_funcs(bool bHang);

			template <template <class Elem, int WorldDim> class TFVGeom>
			struct RegisterCR {
					RegisterCR(this_type* pThis) : m_pThis(pThis){}
					this_type* m_pThis;
					template< typename TElem > void operator()(TElem&)
					{m_pThis->register_cr_func<TElem, TFVGeom>();}
			};

			template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
			void register_cr_func();
	};

} // namespace NavierStokes
} // end namespace ug

#include "navier_stokes_cr_bnd_impl.h"

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__NAVIER_STOKES_CR_BND__ */
