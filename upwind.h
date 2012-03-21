/*
 * upwind.h
 *
 *  Created on: 10.03.2011
 *      Author: andreasvogel
 */

#ifndef NEW_STABILIZATION_IMPL_H___H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__FV__UPWIND__
#define NEW_STABILIZATION_IMPL_H___H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__FV__UPWIND__

#include <vector>

#include "lib_disc/common/local_algebra.h"
#include "lib_disc/spatial_disc/disc_util/finite_volume_geometry.h"

namespace ug{

/////////////////////////////////////////////////////////////////////////////
// Interface for Upwinds
/////////////////////////////////////////////////////////////////////////////


template <int dim>
class INavierStokesUpwind
{
	public:
	/// Abbreviation for own type
		typedef INavierStokesUpwind<dim> this_type;

	public:
	///	constructor
		INavierStokesUpwind()
			: m_pCornerValue(NULL), m_numScvf(0), m_numSh(0), m_bNonZeroShapeIp(true)
		{
			m_vComputeFunc.clear();
			m_vUpdateIPVelFunc.clear();

			m_vIPVel.clear();

			m_vUpConvLength.clear();
			m_vDownConvLength.clear();
			m_vUpShapeSh.clear();
			m_vDownShapeSh.clear();
			m_vUpShapeIp.clear();
			m_vDownShapeIp.clear();

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
			UG_ASSERT(scvf < m_vUpConvLength.size(), "Invalid index");
			return m_vUpConvLength[scvf];
		}

	///	Convection Length
		number downwind_conv_length(size_t scvf) const
		{
			UG_ASSERT(scvf < m_vDownConvLength.size(), "Invalid index");
			return m_vDownConvLength[scvf];
		}

	///	ip velocity (i.e. interpolated velocity at ip)
		const MathVector<dim>& ip_vel(size_t scvf) const
		{
			UG_ASSERT(scvf < m_vIPVel.size(), "Invalid index");
			return m_vIPVel[scvf];
		}

	/// returns the upwind velocity
		MathVector<dim> upwind_vel(size_t scvf) const;

	///	upwind shape for corner vel
		number upwind_shape_sh(size_t scvf, size_t sh) const
		{
			UG_ASSERT(scvf < m_vUpShapeSh.size(), "Invalid index");
			UG_ASSERT(sh < m_vUpShapeSh[scvf].size(), "Invalid index");
			return m_vUpShapeSh[scvf][sh];
		}

	///	upwind shape for corner vel
		number downwind_shape_sh(size_t scvf, size_t sh) const
		{
			UG_ASSERT(scvf < m_vDownShapeSh.size(), "Invalid index");
			UG_ASSERT(sh < m_vDownShapeSh[scvf].size(), "Invalid index");
			return m_vDownShapeSh[scvf][sh];
		}

	///	returns if upwind shape w.r.t. ip vel is non-zero
		bool non_zero_shape_ip() const {return m_bNonZeroShapeIp;}

	///	upwind shapes for ip vel
		number upwind_shape_ip(size_t scvf, size_t scvf2) const
		{
			UG_ASSERT(scvf < m_vUpShapeIp.size(), "Invalid index");
			UG_ASSERT(scvf2 < m_vUpShapeIp[scvf].size(), "Invalid index");
			return m_vUpShapeIp[scvf][scvf2];
		}

	///	upwind shapes for ip vel
		number downwind_shape_ip(size_t scvf, size_t scvf2) const
		{
			UG_ASSERT(scvf < m_vDownShapeIp.size(), "Invalid index");
			UG_ASSERT(scvf2 < m_vDownShapeIp[scvf].size(), "Invalid index");
			return m_vDownShapeIp[scvf][scvf2];
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
			std::vector<MathVector<dim> > vDownIPVel(m_vIPVel.size());
			for(size_t ip = 0; ip < vDownIPVel.size(); ++ip)
				VecScale(vDownIPVel[ip], m_vIPVel[ip], -1.0);

			return compute(geo, vDownIPVel, m_vDownShapeSh, m_vDownShapeIp, m_vDownConvLength);
		}

	///	compute values for new geometry and corner velocities
		bool update_upwind(const FVGeometryBase* geo)
		{
			return compute(geo, m_vIPVel, m_vUpShapeSh, m_vUpShapeIp, m_vUpConvLength);
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
	///	resize the data arrays
		void set_sizes(size_t numScvf, size_t numSh);

	///	sets the shape ip flag
		void set_shape_ip_flag(bool flag) {m_bNonZeroShapeIp = flag;}

	/// non-const access to ip velocity (i.e. interpolated velocity at ip)
		MathVector<dim>& ip_vel(size_t scvf)
		{
			UG_ASSERT(scvf < m_vIPVel.size(), "Invalid index");
			return m_vIPVel[scvf];
		}

	///	non-const access to upwind shapes for corner vel
		number& upwind_shape_sh(size_t scvf, size_t sh)
		{
			UG_ASSERT(scvf < m_vUpShapeSh.size(), "Invalid index");
			UG_ASSERT(sh < m_vUpShapeSh[scvf].size(), "Invalid index");
			return m_vUpShapeSh[scvf][sh];
		}

	///	non-const access to upwind shapes for corner vel
		number& downwind_shape_sh(size_t scvf, size_t sh)
		{
			UG_ASSERT(scvf < m_vDownShapeSh.size(), "Invalid index");
			UG_ASSERT(sh < m_vDownShapeSh[scvf].size(), "Invalid index");
			return m_vDownShapeSh[scvf][sh];
		}

	///	non-const access to upwind shapes for ip vel
		number& upwind_shape_ip(size_t scvf, size_t scvf2)
		{
			UG_ASSERT(scvf < m_vUpShapeIp.size(), "Invalid index");
			UG_ASSERT(scvf2 < m_vUpShapeIp[scvf].size(), "Invalid index");
			return m_vUpShapeIp[scvf][scvf2];
		}

	///	non-const access to upwind shapes for ip vel
		number& downwind_shape_ip(size_t scvf, size_t scvf2)
		{
			UG_ASSERT(scvf < m_vDownShapeIp.size(), "Invalid index");
			UG_ASSERT(scvf2 < m_vDownShapeIp[scvf].size(), "Invalid index");
			return m_vDownShapeIp[scvf][scvf2];
		}

	///	non-const access to Convection Length
		number& upwind_conv_length(size_t scvf)
		{
			UG_ASSERT(scvf < m_vUpConvLength.size(), "Invalid index");
			return m_vUpConvLength[scvf];
		}

	///	non-const access to Convection Length
		number& down_upwind_conv_length(size_t scvf)
		{
			UG_ASSERT(scvf < m_vDownConvLength.size(), "Invalid index");
			return m_vDownConvLength[scvf];
		}

	///	pointer to currently used values
		const LocalVector* m_pCornerValue;

	///	interpolated value at ip
		std::vector<MathVector<dim> > m_vIPVel;

	///	number of current scvf
		size_t m_numScvf;

	///	number of current shape functions (usually in corners)
		size_t m_numSh;

	///	convection length
		std::vector<number> m_vUpConvLength;
		std::vector<number> m_vDownConvLength;

	///	upwind shapes for corners shape functions
		std::vector<std::vector<number> > m_vUpShapeSh;
		std::vector<std::vector<number> > m_vDownShapeSh;

	///	flag if ip shapes are non-zero
		bool m_bNonZeroShapeIp;

	///	upwind shapes for ip vels
		std::vector<std::vector<number> > m_vUpShapeIp;
		std::vector<std::vector<number> > m_vDownShapeIp;

	///	compute values for new geometry and corner velocities
		bool compute(const FVGeometryBase* geo,
		             const std::vector<MathVector<dim> >& vIPVel,
		             std::vector<std::vector<number> >& vUpShapeSh,
		             std::vector<std::vector<number> >& vUpShapeIp,
		             std::vector<number>& vConvLength)
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
		bool set_geometry_type();

	protected:

	///	type of update function
		typedef bool (this_type::*ComputeFunc)(
								const FVGeometryBase* obj,
						        const std::vector<MathVector<dim> >& vIPVel,
								std::vector<std::vector<number> >& vUpShapeSh,
								std::vector<std::vector<number> >& vUpShapeIp,
								std::vector<number>& vConvLength);

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
		             const std::vector<MathVector<dim> >& vIPVel,
		             std::vector<std::vector<number> >& vUpShapeSh,
		             std::vector<std::vector<number> >& vUpShapeIp,
		             std::vector<number>& vConvLength);

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
									const std::vector<MathVector<dim> >& vIPVel,
									std::vector<std::vector<number> >& vUpShapeSh,
									std::vector<std::vector<number> >& vUpShapeIp,
									std::vector<number>& vConvLength);

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
		             const std::vector<MathVector<dim> >& vIPVel,
		             std::vector<std::vector<number> >& vUpShapeSh,
		             std::vector<std::vector<number> >& vUpShapeIp,
		             std::vector<number>& vConvLength);

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
									const std::vector<MathVector<dim> >& vIPVel,
									std::vector<std::vector<number> >& vUpShapeSh,
									std::vector<std::vector<number> >& vUpShapeIp,
									std::vector<number>& vConvLength);

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
		             const std::vector<MathVector<dim> >& vIPVel,
		             std::vector<std::vector<number> >& vUpShapeSh,
		             std::vector<std::vector<number> >& vUpShapeIp,
		             std::vector<number>& vConvLength);

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
									const std::vector<MathVector<dim> >& vIPVel,
									std::vector<std::vector<number> >& vUpShapeSh,
									std::vector<std::vector<number> >& vUpShapeIp,
									std::vector<number>& vConvLength);

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
		             const std::vector<MathVector<dim> >& vIPVel,
		             std::vector<std::vector<number> >& vUpShapeSh,
		             std::vector<std::vector<number> >& vUpShapeIp,
		             std::vector<number>& vConvLength);

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
									const std::vector<MathVector<dim> >& vIPVel,
									std::vector<std::vector<number> >& vUpShapeSh,
									std::vector<std::vector<number> >& vUpShapeIp,
									std::vector<number>& vConvLength);

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
		             const std::vector<MathVector<dim> >& vIPVel,
		             std::vector<std::vector<number> >& vUpShapeSh,
		             std::vector<std::vector<number> >& vUpShapeIp,
		             std::vector<number>& vConvLength);

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
									const std::vector<MathVector<dim> >& vIPVel,
									std::vector<std::vector<number> >& vUpShapeSh,
									std::vector<std::vector<number> >& vUpShapeIp,
									std::vector<number>& vConvLength);

			this->template register_update_func<TGeom, TFunc>(&this_type::template compute<TElem>);
		}
};

} // end namespace ug

// include implementation
#include "upwind_impl.h"

#endif /* NEW_STABILIZATION_IMPL_H___H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__FV__UPWIND__ */
