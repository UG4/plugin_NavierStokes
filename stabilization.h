/*
 * stabilization.h
 *
 *  Created on: 10.03.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__FV__STABILIZATION__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__FV__STABILIZATION__

//#define UG_NSSTAB_ASSERT(cond, exp)
// include define below to assert arrays used in stabilization
#define UG_NSSTAB_ASSERT(cond, exp) UG_ASSERT((cond), (exp))

#include "upwind.h"

namespace ug{

template <int dim>
class INavierStokesStabilization
{
	protected:
	///	used traits
		typedef fv1_dim_traits<dim, dim> traits;

	public:
	///	number of SubControlVolumes
		static const size_t maxNumSCV = traits::maxNumSCV;

	///	max number of SubControlVolumeFaces
		static const size_t maxNumSCVF = traits::maxNumSCVF;

	/// max number of shape functions
		static const size_t maxNumSH = traits::maxNSH;

	public:
	/// Abbreviation for own type
		typedef INavierStokesStabilization<dim> this_type;

		enum DiffusionLength
		{
		    RAW = 0,
		    FIVEPOINT,
		    COR
		};

	public:
	///	constructor
		INavierStokesStabilization()
			:  m_spUpwind(NULL), m_spConstUpwind(NULL),
			   m_numScvf(0), m_numSh(0)
		{
			m_vUpdateFunc.clear();

		//	default setup
			set_diffusion_length("RAW");
		}

	///	sets the type of diff length used for evaluation
		void set_diffusion_length(std::string diffLength);

	///	sets the upwind method
		void set_upwind(SmartPtr<INavierStokesUpwind<dim> > spUpwind)
		{
			m_spUpwind = spUpwind;
			m_spConstUpwind = spUpwind;
		}

	///	returns the upwind
		ConstSmartPtr<INavierStokesUpwind<dim> > upwind() const {return m_spConstUpwind;}

	///	diff length
		number diff_length_sq_inv(size_t scvf) const
		{
			UG_NSSTAB_ASSERT(scvf < m_numScvf, "Invalid index.");
			return m_vDiffLengthSqInv[scvf];
		}

	/// stabilized velocity
		const MathVector<dim>& stab_vel(size_t scvf) const
		{
			UG_NSSTAB_ASSERT(scvf < m_numScvf, "Invalid index.");
			return m_vStabVel[scvf];
		}

	///	returns if stab velocity comp depends on other vel components
		bool vel_comp_connected() const {return m_bVelCompConnected;}

	/// computed stab shape for velocity. This is: The stab_vel derivative
	/// w.r.t velocity unknowns in the corner for each component
		number stab_shape_vel(size_t scvf, size_t compOut, size_t compIn, size_t sh) const
		{
			UG_NSSTAB_ASSERT(scvf < m_numScvf, "Invalid index.");
			UG_NSSTAB_ASSERT(compOut < dim, "Invalid index.");
			UG_NSSTAB_ASSERT(compIn < dim, "Invalid index.");
			UG_NSSTAB_ASSERT(sh < m_numSh, "Invalid index.");
			return m_vvvvStabShapeVel[scvf][compOut][compIn][sh];
		}

	///	computed stab shape for pressure.
		number stab_shape_p(size_t scvf, size_t compOut, size_t sh) const
		{
			UG_NSSTAB_ASSERT(scvf < m_numScvf, "Invalid index.");
			UG_NSSTAB_ASSERT(compOut < dim, "Invalid index.");
			UG_NSSTAB_ASSERT(sh < m_numSh, "Invalid index.");
			return m_vvvvStabShapePressure[scvf][compOut][sh];
		}


	///	compute values for new geometry and corner velocities
		void update(const FVGeometryBase* geo,
		            const LocalVector& vCornerValue,
		            const bool bStokes,
		            const DataImport<number, dim>& kinVisco,
		            const DataImport<MathVector<dim>, dim>* pSource,
		            const LocalVector* pvCornerValueOldTime, number dt)
			{(this->*(m_vUpdateFunc[m_id]))(geo, vCornerValue,
											bStokes, kinVisco, pSource,
											pvCornerValueOldTime, dt);}

	/////////////////////////////////////////
	// forward methods of Upwind Velocity
	/////////////////////////////////////////

	///	Convection Length
		number upwind_conv_length(size_t scvf) const
		{
			UG_NSSTAB_ASSERT(m_spConstUpwind.valid(), "No upwind object");
			return m_spConstUpwind->upwind_conv_length(scvf);
		}

	///	Convection Length
		number downwind_conv_length(size_t scvf) const
		{
			UG_NSSTAB_ASSERT(m_spConstUpwind.valid(), "No upwind object");
			return m_spConstUpwind->downwind_conv_length(scvf);
		}

	///	upwind shape for corner vel
		number upwind_shape_sh(size_t scvf, size_t sh) const
		{
			UG_NSSTAB_ASSERT(m_spConstUpwind.valid(), "No upwind object");
			return m_spConstUpwind->upwind_shape_sh(scvf, sh);
		}

	///	upwind shape for corner vel
		number downwind_shape_sh(size_t scvf, size_t sh) const
		{
			UG_NSSTAB_ASSERT(m_spConstUpwind.valid(), "No upwind object");
			return m_spConstUpwind->downwind_shape_sh(scvf, sh);
		}

	///	returns if upwind shape w.r.t. ip vel is non-zero
		bool non_zero_shape_ip() const
		{
			UG_NSSTAB_ASSERT(m_spConstUpwind.valid(), "No upwind object");
			return m_spConstUpwind->non_zero_shape_ip();
		}

	///	upwind shapes for ip vel
		number upwind_shape_ip(size_t scvf, size_t scvf2) const
		{
			UG_NSSTAB_ASSERT(m_spConstUpwind.valid(), "No upwind object");
			return m_spConstUpwind->upwind_shape_ip(scvf, scvf2);
		}

	///	upwind shapes for ip vel
		number downwind_shape_ip(size_t scvf, size_t scvf2) const
		{
			UG_NSSTAB_ASSERT(m_spConstUpwind.valid(), "No upwind object");
			return m_spConstUpwind->downwind_shape_ip(scvf, scvf2);
		}

	//////////////////////////
	// internal handling
	//////////////////////////

	protected:
	///	computes the diffusion length
		template <typename TFVGeom>
		void compute_upwind(const TFVGeom& geo, const LocalVector& vCornerValue);

	///	computes the diffusion length
		template <typename TFVGeom>
		void compute_downwind(const TFVGeom& geo);

	///	computes the diffusion length
		template <typename TFVGeom>
		void compute_diff_length(const TFVGeom& geo);

	///	sets the vel comp connected flag
		void set_vel_comp_connected(bool bVelCompConnected) {m_bVelCompConnected = bVelCompConnected;}

	protected:
	/// stabilized velocity
		MathVector<dim>& stab_vel(size_t scvf)
		{
			UG_NSSTAB_ASSERT(scvf < m_numScvf, "Invalid index.");
			return m_vStabVel[scvf];
		}

	/// computed stab shape for velocity. This is: The stab_vel derivative
	/// w.r.t velocity unknowns in the corner for each component
		number& stab_shape_vel(size_t scvf, size_t compOut, size_t compIn, size_t sh)
		{
			UG_NSSTAB_ASSERT(scvf < m_numScvf, "Invalid index.");
			UG_NSSTAB_ASSERT(compOut < dim, "Invalid index.");
			UG_NSSTAB_ASSERT(compIn < dim, "Invalid index.");
			UG_NSSTAB_ASSERT(sh < m_numSh, "Invalid index.");
			return m_vvvvStabShapeVel[scvf][compOut][compIn][sh];
		}

	///	computed stab shape for pressure.
		number& stab_shape_p(size_t scvf, size_t compOut, size_t sh)
		{
			UG_NSSTAB_ASSERT(scvf < m_numScvf, "Invalid index.");
			UG_NSSTAB_ASSERT(compOut < dim, "Invalid index.");
			UG_NSSTAB_ASSERT(sh < m_numSh, "Invalid index.");
			return m_vvvvStabShapePressure[scvf][compOut][sh];
		}

	protected:
	///	Upwind values
		SmartPtr<INavierStokesUpwind<dim> > m_spUpwind;
		ConstSmartPtr<INavierStokesUpwind<dim> > m_spConstUpwind;

	///	number of current scvf
		size_t m_numScvf;

	///	number of current shape functions (usually in corners)
		size_t m_numSh;

	///	type of diffusion length computation
		DiffusionLength m_diffLengthType;

	///	vector holding diffusion Length squared and inverted
		number m_vDiffLengthSqInv[maxNumSCVF];

	///	values of stabilized velocity at ip
		MathVector<dim> m_vStabVel[maxNumSCVF];

	///	flag if velocity components are interconnected
		bool m_bVelCompConnected;

	///	stab shapes w.r.t vel
		number m_vvvvStabShapeVel[maxNumSCVF][dim][dim][maxNumSH];

	///	stab shapes w.r.t pressure
		number m_vvvvStabShapePressure[maxNumSCVF][dim][maxNumSH];

	//////////////////////////
	// registering process
	//////////////////////////

	protected:
	///	type of update function
		typedef void (this_type::*UpdateFunc)(	const FVGeometryBase* geo,
												const LocalVector& vCornerValue,
												const bool bStokes,
												const DataImport<number, dim>& kinVisco,
												const DataImport<MathVector<dim>, dim>* pSource,
												const LocalVector* pvCornerValueOldTime, number dt);

	public:
	///	register a update function for a Geometry
		template <typename TFVGeom, typename TAssFunc>
		void register_update_func(TAssFunc func);

	///	set the Geometry type to use for next updates
		template <typename TFVGeom>
		void set_geometry_type();

	protected:
	///	Vector holding all update functions
		std::vector<UpdateFunc> m_vUpdateFunc;

	///	id of current geometry type
		int m_id;
};



/////////////////////////////////////////////////////////////////////////////
// FIELDS
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
class NavierStokesFIELDSStabilization
	: public INavierStokesStabilization<TDim>
{
	public:
	///	Base class
		typedef INavierStokesStabilization<TDim> base_type;

	///	This class
		typedef NavierStokesFIELDSStabilization<TDim> this_type;

	///	Dimension
		static const int dim = TDim;

	protected:
	//	explicitly forward some function
		using base_type::register_update_func;
		using base_type::diff_length_sq_inv;
		using base_type::stab_shape_vel;
		using base_type::stab_shape_p;
		using base_type::stab_vel;
		using base_type::set_vel_comp_connected;

	//	functions from upwind
		using base_type::upwind_conv_length;
		using base_type::downwind_conv_length;
		using base_type::upwind_shape_sh;
		using base_type::downwind_shape_sh;
		using base_type::non_zero_shape_ip;
		using base_type::upwind_shape_ip;
		using base_type::downwind_shape_ip;

	public:
	///	constructor
		NavierStokesFIELDSStabilization()
		{
		//	vel comp not interconnected
			set_vel_comp_connected(false);

		//	register evaluation function
			register_func(Int2Type<dim>());
		}

	///	update of values for FV1Geometry
		template <typename TElem>
		void update(const FV1Geometry<TElem, dim>* geo,
		            const LocalVector& vCornerValue,
		            const bool bStokes,
		            const DataImport<number, dim>& kinVisco,
		            const DataImport<MathVector<dim>, dim>* pSource,
		            const LocalVector* pvCornerValueOldTime, number dt);

	private:
		void register_func(Int2Type<1>)
		{register_func<Edge>();}

		void register_func(Int2Type<2>)
		{	register_func(Int2Type<1>());
			register_func<Triangle>();
			register_func<Quadrilateral>();}

		void register_func(Int2Type<3>)
		{	register_func(Int2Type<2>());
			register_func<Tetrahedron>();
			register_func<Pyramid>();
			register_func<Prism>();
			register_func<Hexahedron>();}

		template <typename TElem>
		void register_func()
		{
			typedef FV1Geometry<TElem, dim> TGeom;
			typedef void (this_type::*TFunc)(const TGeom* geo,
											 const LocalVector& vCornerValue,
											 const bool bStokes,
											 const DataImport<number, dim>& kinVisco,
											 const DataImport<MathVector<dim>, dim>* pSource,
											 const LocalVector* pvCornerValueOldTime, number dt);

			this->template register_update_func<TGeom, TFunc>(&this_type::template update<TElem>);
		}
};


/////////////////////////////////////////////////////////////////////////////
// FLOW
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
class NavierStokesFLOWStabilization
	: public INavierStokesStabilization<TDim>
{
	public:
	///	Base class
		typedef INavierStokesStabilization<TDim> base_type;

	///	This class
		typedef NavierStokesFLOWStabilization<TDim> this_type;

	///	Dimension
		static const int dim = TDim;

	protected:
	//	explicitly forward some function
		using base_type::register_update_func;
		using base_type::diff_length_sq_inv;
		using base_type::stab_shape_vel;
		using base_type::stab_shape_p;
		using base_type::stab_vel;
		using base_type::set_vel_comp_connected;

	//	functions from upwind
		using base_type::upwind_conv_length;
		using base_type::downwind_conv_length;
		using base_type::upwind_shape_sh;
		using base_type::downwind_shape_sh;
		using base_type::non_zero_shape_ip;
		using base_type::upwind_shape_ip;
		using base_type::downwind_shape_ip;

	public:
	///	constructor
		NavierStokesFLOWStabilization()
		{
		//	vel comp not interconnected
			set_vel_comp_connected(true);

		//	register evaluation function
			register_func(Int2Type<dim>());
		}

	///	update of values for FV1Geometry
		template <typename TElem>
		void update(const FV1Geometry<TElem, dim>* geo,
					const LocalVector& vCornerValue,
					const bool bStokes,
					const DataImport<number, dim>& kinVisco,
					const DataImport<MathVector<dim>, dim>* pSource,
					const LocalVector* pvCornerValueOldTime, number dt);

	private:
		void register_func(Int2Type<1>)
		{register_func<Edge>();}

		void register_func(Int2Type<2>)
		{	register_func(Int2Type<1>());
			register_func<Triangle>();
			register_func<Quadrilateral>();}

		void register_func(Int2Type<3>)
		{	register_func(Int2Type<2>());
			register_func<Tetrahedron>();
			register_func<Pyramid>();
			register_func<Prism>();
			register_func<Hexahedron>();}

		template <typename TElem>
		void register_func()
		{
			typedef FV1Geometry<TElem, dim> TGeom;
			typedef void (this_type::*TFunc)(const TGeom* geo,
											 const LocalVector& vCornerValue,
											 const bool bStokes,
											 const DataImport<number, dim>& kinVisco,
											 const DataImport<MathVector<dim>, dim>* pSource,
											 const LocalVector* pvCornerValueOldTime, number dt);

			this->template register_update_func<TGeom, TFunc>(&this_type::template update<TElem>);
		}
};


} // end namespace ug

// include implementation
#include "stabilization_impl.h"

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__FV__STABILIZATION__ */
