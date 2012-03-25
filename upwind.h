/*
 * upwind.h
 *
 *  Created on: 10.03.2011
 *      Author: andreasvogel
 */

#ifndef NEW_STABILIZATION_IMPL_H___H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__FV__UPWIND__
#define NEW_STABILIZATION_IMPL_H___H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__FV__UPWIND__

//#define UG_NSUPWIND_ASSERT(cond, exp)
// include define below to assert arrays used in stabilization
#define UG_NSUPWIND_ASSERT(cond, exp) UG_ASSERT((cond), (exp))

#include <vector>

#include "lib_disc/common/local_algebra.h"
#include "lib_disc/spatial_disc/disc_util/finite_volume_util.h"
#include "lib_disc/spatial_disc/disc_util/finite_volume_geometry.h"

namespace ug{

/////////////////////////////////////////////////////////////////////////////
// Interface for Upwinds
/////////////////////////////////////////////////////////////////////////////


template <int dim>
class INavierStokesUpwind
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
		typedef INavierStokesUpwind<dim> this_type;

	public:
	///	constructor
		INavierStokesUpwind()
			: m_pCornerValue(NULL),
			  m_numScvf(0), m_numSh(0),
			  m_bNonZeroShapeIp(true)
		{
			m_vComputeFunc.clear();
			m_vUpdateIPVelFunc.clear();

		//	register evaluation function
			register_func(Int2Type<dim>());
		}

	///	returns number of shapes
		size_t num_sh() const {return m_numSh;}

	///	returns number of sub control volume faces
		size_t num_scvf() const {return m_numScvf;}

	///	Convection Length
		number upwind_conv_length(size_t scvf) const
		{
			UG_NSUPWIND_ASSERT(scvf < m_numScvf, "Invalid index");
			return m_vUpConvLength[scvf];
		}

	///	Convection Length
		number downwind_conv_length(size_t scvf) const
		{
			UG_NSUPWIND_ASSERT(scvf < m_numScvf, "Invalid index");
			return m_vDownConvLength[scvf];
		}

	///	ip velocity (i.e. interpolated velocity at ip)
		const MathVector<dim>& ip_vel(size_t scvf) const
		{
			UG_NSUPWIND_ASSERT(scvf < m_numScvf, "Invalid index");
			return m_vIPVel[scvf];
		}

	/// returns the upwind velocity
		MathVector<dim> upwind_vel(size_t scvf) const;

	///	upwind shape for corner vel
		number upwind_shape_sh(size_t scvf, size_t sh) const
		{
			UG_NSUPWIND_ASSERT(scvf < m_numScvf, "Invalid index");
			UG_NSUPWIND_ASSERT(sh < m_numSh, "Invalid index");
			return m_vvUpShapeSh[scvf][sh];
		}

	///	upwind shape for corner vel
		number downwind_shape_sh(size_t scvf, size_t sh) const
		{
			UG_NSUPWIND_ASSERT(scvf < m_numScvf, "Invalid index");
			UG_NSUPWIND_ASSERT(sh < m_numSh, "Invalid index");
			return m_vvDownShapeSh[scvf][sh];
		}

	///	returns if upwind shape w.r.t. ip vel is non-zero
		bool non_zero_shape_ip() const {return m_bNonZeroShapeIp;}

	///	upwind shapes for ip vel
		number upwind_shape_ip(size_t scvf, size_t scvf2) const
		{
			UG_NSUPWIND_ASSERT(scvf < m_numScvf, "Invalid index");
			UG_NSUPWIND_ASSERT(scvf2 < m_numScvf, "Invalid index");
			return m_vvUpShapeIp[scvf][scvf2];
		}

	///	upwind shapes for ip vel
		number downwind_shape_ip(size_t scvf, size_t scvf2) const
		{
			UG_NSUPWIND_ASSERT(scvf < m_numScvf, "Invalid index");
			UG_NSUPWIND_ASSERT(scvf2 < m_numScvf, "Invalid index");
			return m_vvDownShapeIp[scvf][scvf2];
		}

	///	compute values for new geometry and corner velocities
		bool update(const FVGeometryBase* geo, const LocalVector& vCornerValue)
		{
			m_pCornerValue = &vCornerValue;
			bool bRet = true;
			bRet &= update_ip_vel(geo, vCornerValue);
			bRet &= update_upwind(geo);
			return bRet;
		}

	///	compute values for new geometry and corner velocities
		bool update_downwind(const FVGeometryBase* geo)
		{
			MathVector<dim> vDownIPVel[maxNumSCVF];
			for(size_t ip = 0; ip < m_numScvf; ++ip)
				VecScale(vDownIPVel[ip], m_vIPVel[ip], -1.0);

			return compute(geo, vDownIPVel, m_vvDownShapeSh, m_vvDownShapeIp, m_vDownConvLength);
		}

	///	compute values for new geometry and corner velocities
		bool update_upwind(const FVGeometryBase* geo)
		{
			return compute(geo, m_vIPVel, m_vvUpShapeSh, m_vvUpShapeIp, m_vUpConvLength);
		}

	///	compute values for new geometry and corner velocities
		bool update_ip_vel(const FVGeometryBase* geo, const LocalVector& vCornerValue)
		{
			return (this->*(m_vUpdateIPVelFunc[m_id]))(geo, vCornerValue);
		}

	//////////////////////////
	// internal handling
	//////////////////////////

	protected:
	///	sets the shape ip flag
		void set_shape_ip_flag(bool flag) {m_bNonZeroShapeIp = flag;}

	/// non-const access to ip velocity (i.e. interpolated velocity at ip)
		MathVector<dim>& ip_vel(size_t scvf)
		{
			UG_NSUPWIND_ASSERT(scvf < m_numScvf, "Invalid index");
			return m_vIPVel[scvf];
		}

	///	non-const access to upwind shapes for corner vel
		number& upwind_shape_sh(size_t scvf, size_t sh)
		{
			UG_NSUPWIND_ASSERT(scvf < m_numScvf, "Invalid index");
			UG_NSUPWIND_ASSERT(sh < m_numSh, "Invalid index");
			return m_vvUpShapeSh[scvf][sh];
		}

	///	non-const access to upwind shapes for corner vel
		number& downwind_shape_sh(size_t scvf, size_t sh)
		{
			UG_NSUPWIND_ASSERT(scvf < m_numScvf, "Invalid index");
			UG_NSUPWIND_ASSERT(sh < m_numSh, "Invalid index");
			return m_vvDownShapeSh[scvf][sh];
		}

	///	non-const access to upwind shapes for ip vel
		number& upwind_shape_ip(size_t scvf, size_t scvf2)
		{
			UG_NSUPWIND_ASSERT(scvf < m_numScvf, "Invalid index");
			UG_NSUPWIND_ASSERT(scvf2 < m_numScvf, "Invalid index");
			return m_vvUpShapeIp[scvf][scvf2];
		}

	///	non-const access to upwind shapes for ip vel
		number& downwind_shape_ip(size_t scvf, size_t scvf2)
		{
			UG_NSUPWIND_ASSERT(scvf < m_numScvf, "Invalid index");
			UG_NSUPWIND_ASSERT(scvf2 < m_numScvf, "Invalid index");
			return m_vvDownShapeIp[scvf][scvf2];
		}

	///	non-const access to Convection Length
		number& upwind_conv_length(size_t scvf)
		{
			UG_NSUPWIND_ASSERT(scvf < m_numScvf, "Invalid index");
			return m_vUpConvLength[scvf];
		}

	///	non-const access to Convection Length
		number& down_upwind_conv_length(size_t scvf)
		{
			UG_NSUPWIND_ASSERT(scvf < m_numScvf, "Invalid index");
			return m_vDownConvLength[scvf];
		}

	///	pointer to currently used values
		const LocalVector* m_pCornerValue;

	///	number of current scvf
		size_t m_numScvf;

	///	number of current shape functions (usually in corners)
		size_t m_numSh;

	///	interpolated value at ip
		MathVector<dim> m_vIPVel[maxNumSCVF];

	///	convection length
		number m_vUpConvLength[maxNumSCVF];
		number m_vDownConvLength[maxNumSCVF];

	///	upwind shapes for corners shape functions
		number m_vvUpShapeSh[maxNumSCVF][maxNumSH];
		number m_vvDownShapeSh[maxNumSCVF][maxNumSH];

	///	flag if ip shapes are non-zero
		bool m_bNonZeroShapeIp;

	///	upwind shapes for ip vels
		number m_vvUpShapeIp[maxNumSCVF][maxNumSCVF];
		number m_vvDownShapeIp[maxNumSCVF][maxNumSCVF];

	///	compute values for new geometry and corner velocities
		bool compute(const FVGeometryBase* geo,
		             const MathVector<dim> vIPVel[maxNumSCVF],
		             number vUpShapeSh[maxNumSCVF][maxNumSH],
		             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
		             number vConvLength[maxNumSCVF])
		{
			return (this->*(m_vComputeFunc[m_id]))(geo, vIPVel, vUpShapeSh, vUpShapeIp, vConvLength);
		}

	protected:
	///	update of values for FV1Geometry
		template <typename TElem>
		bool update_ip_vel(const FV1Geometry<TElem, dim>* geo, const LocalVector& vCornerValue);

	//////////////////////////
	// registering process
	//////////////////////////
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
			typedef bool (this_type::*TFunc)(const TGeom* geo, const LocalVector& vCornerValue);

			this->template register_update_ip_vel_func<TGeom, TFunc>(&this_type::template update_ip_vel<TElem>);
		}

	public:
	///	register a update function for a Geometry
		template <typename TFVGeom, typename TAssFunc>
		void register_update_func(TAssFunc func);

	///	register a update function for a Geometry
		template <typename TFVGeom, typename TAssFunc>
		void register_update_ip_vel_func(TAssFunc func);

	///	set the Geometry type to use for next updates
		template <typename TFVGeom>
		void set_geometry_type();

	protected:

	///	type of update function
		typedef bool (this_type::*ComputeFunc)(
								const FVGeometryBase* obj,
								const MathVector<dim> vIPVel[maxNumSCVF],
								number vUpShapeSh[maxNumSCVF][maxNumSH],
								number vUpShapeIp[maxNumSCVF][maxNumSCVF],
								number vConvLength[maxNumSCVF]);

	///	type of update function
		typedef bool (this_type::*UpdateIPVelFunc)(const FVGeometryBase* obj,
													const LocalVector& vCornerVels);

	///	Vector holding all update functions
		std::vector<ComputeFunc> m_vComputeFunc;

	///	Vector holding all update functions
		std::vector<UpdateIPVelFunc> m_vUpdateIPVelFunc;

	///	id of current geometry type
		int m_id;
};


/////////////////////////////////////////////////////////////////////////////
// No Upwind
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
class NavierStokesNoUpwind
	: public INavierStokesUpwind<TDim>
{
	public:
	///	Base class
		typedef INavierStokesUpwind<TDim> base_type;

	///	This class
		typedef NavierStokesNoUpwind<TDim> this_type;

	///	Dimension
		static const int dim = TDim;

	protected:
	//	explicitly forward some function
		using base_type::set_shape_ip_flag;
		using base_type::register_update_func;

		static const size_t maxNumSCV = base_type::maxNumSCV;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesNoUpwind()
		{
		//	shapes for ip vels are zero (no dependency between ip shapes, only to corners)
		//	Note: during resize, values are initialized to zero, thus values
		//		  for shapes depending on ip values are always correct (i.e. zero)
			set_shape_ip_flag(false);

		//	register evaluation function
			register_func(Int2Type<dim>());
		}

	///	update of values for FV1Geometry
		template <typename TElem>
		bool compute(const FV1Geometry<TElem, dim>* geo,
		             const MathVector<dim> vIPVel[maxNumSCVF],
		             number vUpShapeSh[maxNumSCVF][maxNumSH],
		             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
		             number vConvLength[maxNumSCVF]);

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
			typedef bool (this_type::*TFunc)(
									const TGeom* obj,
						             const MathVector<dim> vIPVel[maxNumSCVF],
						             number vUpShapeSh[maxNumSCVF][maxNumSH],
						             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
						             number vConvLength[maxNumSCVF]);

			this->template register_update_func<TGeom, TFunc>(&this_type::template compute<TElem>);
		}
};

/////////////////////////////////////////////////////////////////////////////
// Full Upwind
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
class NavierStokesFullUpwind
	: public INavierStokesUpwind<TDim>
{
	public:
	///	Base class
		typedef INavierStokesUpwind<TDim> base_type;

	///	This class
		typedef NavierStokesFullUpwind<TDim> this_type;

	///	Dimension
		static const int dim = TDim;

	protected:
	//	explicitly forward some function
		using base_type::set_shape_ip_flag;
		using base_type::register_update_func;

		static const size_t maxNumSCV = base_type::maxNumSCV;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesFullUpwind()
		{
		//	shapes for ip vels are zero (no dependency between ip shapes, only to corners)
		//	Note: during resize, values are initialized to zero, thus values
		//		  for shapes depending on ip values are always correct (i.e. zero)
			set_shape_ip_flag(false);

		//	register evaluation function
			register_func(Int2Type<dim>());
		}

	///	update of values for FV1Geometry
		template <typename TElem>
		bool compute(const FV1Geometry<TElem, dim>* geo,
		             const MathVector<dim> vIPVel[maxNumSCVF],
		             number vUpShapeSh[maxNumSCVF][maxNumSH],
		             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
		             number vConvLength[maxNumSCVF]);

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
			typedef bool (this_type::*TFunc)(
									const TGeom* obj,
						             const MathVector<dim> vIPVel[maxNumSCVF],
						             number vUpShapeSh[maxNumSCVF][maxNumSH],
						             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
						             number vConvLength[maxNumSCVF]);

			this->template register_update_func<TGeom, TFunc>(&this_type::template compute<TElem>);
		}
};

/////////////////////////////////////////////////////////////////////////////
// Skewed Upwind
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
class NavierStokesSkewedUpwind
	: public INavierStokesUpwind<TDim>
{
	public:
	///	Base class
		typedef INavierStokesUpwind<TDim> base_type;

	///	This class
		typedef NavierStokesSkewedUpwind<TDim> this_type;

	///	Dimension
		static const int dim = TDim;

	protected:
	//	explicitly forward some function
		using base_type::set_shape_ip_flag;
		using base_type::register_update_func;

		static const size_t maxNumSCV = base_type::maxNumSCV;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesSkewedUpwind()
		{
		//	shapes for ip vels are zero (no dependency between ip shapes, only to corners)
		//	Note: during resize, values are initialized to zero, thus values
		//		  for shapes depending on ip values are always correct (i.e. zero)
			set_shape_ip_flag(false);

		//	register evaluation function
			register_func(Int2Type<dim>());
		}

	///	update of values for FV1Geometry
		template <typename TElem>
		bool compute(const FV1Geometry<TElem, dim>* geo,
		             const MathVector<dim> vIPVel[maxNumSCVF],
		             number vUpShapeSh[maxNumSCVF][maxNumSH],
		             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
		             number vConvLength[maxNumSCVF]);

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
			typedef bool (this_type::*TFunc)(
									const TGeom* obj,
						             const MathVector<dim> vIPVel[maxNumSCVF],
						             number vUpShapeSh[maxNumSCVF][maxNumSH],
						             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
						             number vConvLength[maxNumSCVF]);

			this->template register_update_func<TGeom, TFunc>(&this_type::template compute<TElem>);
		}
};

/////////////////////////////////////////////////////////////////////////////
// Linear Profile Skewed Upwind
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
class NavierStokesLinearProfileSkewedUpwind
	: public INavierStokesUpwind<TDim>
{
	public:
	///	Base class
		typedef INavierStokesUpwind<TDim> base_type;

	///	This class
		typedef NavierStokesLinearProfileSkewedUpwind<TDim> this_type;

	///	Dimension
		static const int dim = TDim;

	protected:
	//	explicitly forward some function
		using base_type::set_shape_ip_flag;
		using base_type::register_update_func;

		static const size_t maxNumSCV = base_type::maxNumSCV;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesLinearProfileSkewedUpwind()
		{
		//	shapes for ip vels are zero (no dependency between ip shapes, only to corners)
		//	Note: during resize, values are initialized to zero, thus values
		//		  for shapes depending on ip values are always correct (i.e. zero)
			set_shape_ip_flag(false);

		//	register evaluation function
			register_func(Int2Type<dim>());
		}

	///	update of values for FV1Geometry
		template <typename TElem>
		bool compute(const FV1Geometry<TElem, dim>* geo,
		             const MathVector<dim> vIPVel[maxNumSCVF],
		             number vUpShapeSh[maxNumSCVF][maxNumSH],
		             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
		             number vConvLength[maxNumSCVF]);

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
			typedef bool (this_type::*TFunc)(
									const TGeom* obj,
						             const MathVector<dim> vIPVel[maxNumSCVF],
						             number vUpShapeSh[maxNumSCVF][maxNumSH],
						             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
						             number vConvLength[maxNumSCVF]);

			this->template register_update_func<TGeom, TFunc>(&this_type::template compute<TElem>);
		}
};


/////////////////////////////////////////////////////////////////////////////
// Positive Upwind
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
class NavierStokesPositiveUpwind
	: public INavierStokesUpwind<TDim>
{
	public:
	///	Base class
		typedef INavierStokesUpwind<TDim> base_type;

	///	This class
		typedef NavierStokesPositiveUpwind<TDim> this_type;

	///	Dimension
		static const int dim = TDim;

	protected:
	//	explicitly forward some function
		using base_type::set_shape_ip_flag;
		using base_type::register_update_func;

		static const size_t maxNumSCV = base_type::maxNumSCV;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesPositiveUpwind()
		{
		//	shapes for ip vels are non-zero (dependency between ip shapes)
			set_shape_ip_flag(true);

		//	register evaluation function
			register_func(Int2Type<dim>());
		}

	///	update of values for FV1Geometry
		template <typename TElem>
		bool compute(const FV1Geometry<TElem, dim>* geo,
		             const MathVector<dim> vIPVel[maxNumSCVF],
		             number vUpShapeSh[maxNumSCVF][maxNumSH],
		             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
		             number vConvLength[maxNumSCVF]);

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
			typedef bool (this_type::*TFunc)(
									const TGeom* obj,
						             const MathVector<dim> vIPVel[maxNumSCVF],
						             number vUpShapeSh[maxNumSCVF][maxNumSH],
						             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
						             number vConvLength[maxNumSCVF]);

			this->template register_update_func<TGeom, TFunc>(&this_type::template compute<TElem>);
		}
};

} // end namespace ug

// include implementation
#include "upwind_impl.h"

#endif /* NEW_STABILIZATION_IMPL_H___H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__FV__UPWIND__ */
