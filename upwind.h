/*
 * upwind.h
 *
 *  Created on: 10.03.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__NAVIER_STOKES__UPWIND__
#define __H__UG__NAVIER_STOKES__UPWIND__

//#define UG_NSUPWIND_ASSERT(cond, exp)
// include define below to assert arrays used in stabilization
#define UG_NSUPWIND_ASSERT(cond, exp) UG_ASSERT((cond), exp)

#include <vector>

#include "lib_disc/common/local_algebra.h"
#include "lib_disc/spatial_disc/disc_util/finite_volume_util.h"
#include "lib_disc/spatial_disc/disc_util/finite_volume_geometry.h"
#include "lib_disc/spatial_disc/disc_util/cr_finite_volume_geometry.h"

namespace ug{
namespace NavierStokes{

/////////////////////////////////////////////////////////////////////////////
// Interface for Upwinds on collocated grid (stabilization method)
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
			: m_numScvf(0), m_numSh(0),
			  m_bNonZeroShapeIp(true)
		{
			m_vComputeFunc.clear();
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

	/// returns the upwind velocity
		MathVector<dim> upwind_vel(const size_t scvf,
		                           const LocalVector& CornerVel,
		                           const MathVector<dim> vStdVel[]) const;

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
			UG_NSUPWIND_ASSERT(scvf < m_numScvf, "Invalid index: " << scvf);
			UG_NSUPWIND_ASSERT(scvf2 < m_numScvf, "Invalid index2: " << scvf2);
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
		void update(const FVGeometryBase* geo,
		            const MathVector<dim> vStdVel[])
		{
			update_upwind(geo, vStdVel);
		}

	///	compute values for new geometry and corner velocities
		void update_upwind(const FVGeometryBase* geo,
						   const MathVector<dim> vStdVel[])
		{
			compute(geo, vStdVel, m_vvUpShapeSh, m_vvUpShapeIp, m_vUpConvLength);
		}

	///	compute values for new geometry and corner velocities
		void update_downwind(const FVGeometryBase* geo,
		                     const MathVector<dim> vStdVel[])
		{
			MathVector<dim> vDownIPVel[maxNumSCVF];
			for(size_t ip = 0; ip < m_numScvf; ++ip)
				VecScale(vDownIPVel[ip], vStdVel[ip], -1.0);

			compute(geo, vDownIPVel, m_vvDownShapeSh, m_vvDownShapeIp, m_vDownConvLength);
		}

	//////////////////////////
	// internal handling
	//////////////////////////

	protected:
	///	sets the shape ip flag
		void set_shape_ip_flag(bool flag) {m_bNonZeroShapeIp = flag;}

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

	///	number of current scvf
		size_t m_numScvf;

	///	number of current shape functions (usually in corners)
		size_t m_numSh;

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
		void compute(const FVGeometryBase* geo,
		             const MathVector<dim> vIPVel[maxNumSCVF],
		             number vUpShapeSh[maxNumSCVF][maxNumSH],
		             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
		             number vConvLength[maxNumSCVF])
		{
			(this->*(m_vComputeFunc[m_id]))(geo, vIPVel, vUpShapeSh, vUpShapeIp, vConvLength);
		}

	//////////////////////////
	// registering process
	//////////////////////////
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
		typedef void (this_type::*ComputeFunc)(
								const FVGeometryBase* obj,
								const MathVector<dim> vIPVel[maxNumSCVF],
								number vUpShapeSh[maxNumSCVF][maxNumSH],
								number vUpShapeIp[maxNumSCVF][maxNumSCVF],
								number vConvLength[maxNumSCVF]);

	///	Vector holding all update functions
		std::vector<ComputeFunc> m_vComputeFunc;

	///	id of current geometry type
		int m_id;
};

//	register a update function for a Geometry
template <int dim>
template <typename TFVGeom, typename TAssFunc>
void
INavierStokesUpwind<dim>::
register_update_func(TAssFunc func)
{
//	get unique geometry id
	size_t id = GetUniqueFVGeomID<TFVGeom>();

//	make sure that there is enough space
	if((size_t)id >= m_vComputeFunc.size())
		m_vComputeFunc.resize(id+1, NULL);

//	set pointer
	m_vComputeFunc[id] = (ComputeFunc)func;
}


//	set the Geometry type to use for next updates
template <int dim>
template <typename TFVGeom>
void
INavierStokesUpwind<dim>::
set_geometry_type()
{
//	get unique geometry id
	size_t id = GetUniqueFVGeomID<TFVGeom>();

//	check that function exists
	if(id >= m_vComputeFunc.size() || m_vComputeFunc[id] == NULL)
		UG_THROW("No update function registered for Geometry "<<id);

//	set current geometry
	m_id = id;

//	set sizes
	TFVGeom& geo = Provider<TFVGeom>::get();
	m_numScvf = geo.num_scvf();
	m_numSh = geo.num_sh();
	UG_NSUPWIND_ASSERT(m_numScvf <= maxNumSCVF, "Invalid index");
	UG_NSUPWIND_ASSERT(m_numSh <= maxNumSH, "Invalid index");
}


///	upwind velocity
template <int dim>
MathVector<dim>
INavierStokesUpwind<dim>::
upwind_vel(const size_t scvf,
           const LocalVector& CornerVel,
           const MathVector<dim> vStdVel[]) const
{
//	reset result
	MathVector<dim> vel; VecSet(vel, 0.0);

//	add corner shapes
	for(size_t sh = 0; sh < num_sh(); ++sh)
		for(int d = 0; d < dim; ++d)
			vel[d] += upwind_shape_sh(scvf, sh) * CornerVel(d, sh);

//	done if only depending on shapes
	if(!non_zero_shape_ip()) return vel;

//	compute ip vel
	for(size_t scvf2 = 0; scvf2 < num_scvf(); ++scvf2)
		VecScaleAppend(vel, upwind_shape_ip(scvf, scvf2), vStdVel[scvf2]);

//	return value
	return vel;
}

template <int dim, typename TImpl>
class NavierStokesUpwindBase : public INavierStokesUpwind<dim>
{
	public:
	///	Base class
		typedef INavierStokesUpwind<dim> base_type;

	///	This class
		typedef NavierStokesUpwindBase<dim, TImpl> this_type;

	protected:
		static const size_t maxNumSCV = base_type::maxNumSCV;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesUpwindBase() {register_func(Int2Type<dim>());}

	///	update of values for FV1Geometry
		template <typename TElem>
		void compute(const FV1Geometry<TElem, dim>* geo,
		             const MathVector<dim> vIPVel[maxNumSCVF],
		             number vUpShapeSh[maxNumSCVF][maxNumSH],
		             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
		             number vConvLength[maxNumSCVF])
		{
			getImpl().template compute<TElem>(geo, vIPVel, vUpShapeSh, vUpShapeIp, vConvLength);
		}

	private:
		void register_func(Int2Type<1>)
		{
			register_func<Edge>();
		}

		void register_func(Int2Type<2>)
		{
			register_func<Triangle>();
			register_func<Quadrilateral>();
		}

		void register_func(Int2Type<3>)
		{
			register_func<Tetrahedron>();
			register_func<Pyramid>();
			register_func<Prism>();
			register_func<Hexahedron>();
		}

		template <typename TElem>
		void register_func()
		{
			typedef FV1Geometry<TElem, dim> TGeom;
			typedef void (this_type::*TFunc)(
									const TGeom* obj,
						             const MathVector<dim> vIPVel[maxNumSCVF],
						             number vUpShapeSh[maxNumSCVF][maxNumSH],
						             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
						             number vConvLength[maxNumSCVF]);

			this->template register_update_func<TGeom, TFunc>(&this_type::template compute<TElem>);
		}

	protected:
	///	access to implementation
		TImpl& getImpl() {return static_cast<TImpl&>(*this);}

	///	const access to implementation
		const TImpl& getImpl() const {return static_cast<const TImpl&>(*this);}
};

/////////////////////////////////////////////////////////////////////////////
// No Upwind
/////////////////////////////////////////////////////////////////////////////

template <int dim>
class NavierStokesNoUpwind
	: public NavierStokesUpwindBase<dim, NavierStokesNoUpwind<dim> >
{
	public:
		typedef INavierStokesUpwind<dim> base_type;
		static const size_t maxNumSCV = base_type::maxNumSCV;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesNoUpwind(){this->set_shape_ip_flag(false);}

	///	update of values for FV1Geometry
		template <typename TElem>
		void compute(const FV1Geometry<TElem, dim>* geo,
		             const MathVector<dim> vIPVel[maxNumSCVF],
		             number vUpShapeSh[maxNumSCVF][maxNumSH],
		             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
		             number vConvLength[maxNumSCVF]);
};

/////////////////////////////////////////////////////////////////////////////
// Full Upwind
/////////////////////////////////////////////////////////////////////////////

template <int dim>
class NavierStokesFullUpwind
: public NavierStokesUpwindBase<dim, NavierStokesFullUpwind<dim> >
{
	public:
		typedef INavierStokesUpwind<dim> base_type;
		static const size_t maxNumSCV = base_type::maxNumSCV;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesFullUpwind(){this->set_shape_ip_flag(false);}

	///	update of values for FV1Geometry
		template <typename TElem>
		void compute(const FV1Geometry<TElem, dim>* geo,
		             const MathVector<dim> vIPVel[maxNumSCVF],
		             number vUpShapeSh[maxNumSCVF][maxNumSH],
		             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
		             number vConvLength[maxNumSCVF]);
};

/////////////////////////////////////////////////////////////////////////////
// Skewed Upwind
/////////////////////////////////////////////////////////////////////////////

template <int dim>
class NavierStokesSkewedUpwind
: public NavierStokesUpwindBase<dim, NavierStokesSkewedUpwind<dim> >
{
	public:
		typedef INavierStokesUpwind<dim> base_type;
		static const size_t maxNumSCV = base_type::maxNumSCV;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesSkewedUpwind(){this->set_shape_ip_flag(false);}

	///	update of values for FV1Geometry
		template <typename TElem>
		void compute(const FV1Geometry<TElem, dim>* geo,
		             const MathVector<dim> vIPVel[maxNumSCVF],
		             number vUpShapeSh[maxNumSCVF][maxNumSH],
		             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
		             number vConvLength[maxNumSCVF]);
};

/////////////////////////////////////////////////////////////////////////////
// Linear Profile Skewed Upwind
/////////////////////////////////////////////////////////////////////////////

template <int dim>
class NavierStokesLinearProfileSkewedUpwind
: public NavierStokesUpwindBase<dim, NavierStokesLinearProfileSkewedUpwind<dim> >
{
	public:
		typedef INavierStokesUpwind<dim> base_type;
		static const size_t maxNumSCV = base_type::maxNumSCV;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesLinearProfileSkewedUpwind(){this->set_shape_ip_flag(false);}

	///	update of values for FV1Geometry
		template <typename TElem>
		void compute(const FV1Geometry<TElem, dim>* geo,
		             const MathVector<dim> vIPVel[maxNumSCVF],
		             number vUpShapeSh[maxNumSCVF][maxNumSH],
		             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
		             number vConvLength[maxNumSCVF]);
};


/////////////////////////////////////////////////////////////////////////////
// Positive Upwind
/////////////////////////////////////////////////////////////////////////////

template <int dim>
class NavierStokesPositiveUpwind
: public NavierStokesUpwindBase<dim, NavierStokesPositiveUpwind<dim> >
{
	public:
		typedef INavierStokesUpwind<dim> base_type;
		static const size_t maxNumSCV = base_type::maxNumSCV;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesPositiveUpwind(){this->set_shape_ip_flag(true);}

	///	update of values for FV1Geometry
		template <typename TElem>
		void compute(const FV1Geometry<TElem, dim>* geo,
		             const MathVector<dim> vIPVel[maxNumSCVF],
		             number vUpShapeSh[maxNumSCVF][maxNumSH],
		             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
		             number vConvLength[maxNumSCVF]);
};

/////////////////////////////////////////////////////////////////////////////
// Regular Upwind
/////////////////////////////////////////////////////////////////////////////

template <int dim>
class NavierStokesRegularUpwind
: public NavierStokesUpwindBase<dim, NavierStokesRegularUpwind<dim> >
{
	public:
		typedef INavierStokesUpwind<dim> base_type;
		static const size_t maxNumSCV = base_type::maxNumSCV;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesRegularUpwind(){this->set_shape_ip_flag(true);}

	///	update of values for FV1Geometry
		template <typename TElem>
		void compute(const FV1Geometry<TElem, dim>* geo,
					 const MathVector<dim> vIPVel[maxNumSCVF],
					 number vUpShapeSh[maxNumSCVF][maxNumSH],
					 number vUpShapeIp[maxNumSCVF][maxNumSCVF],
					 number vConvLength[maxNumSCVF]);
};


template <int dim, typename TImpl>
class NavierStokesCRUpwindBase : public INavierStokesUpwind<dim>
{
	public:
	///	Base class
		typedef INavierStokesUpwind<dim> base_type;

	///	This class
		typedef NavierStokesCRUpwindBase<dim, TImpl> this_type;

	protected:
		static const size_t maxNumSCV = base_type::maxNumSCV;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesCRUpwindBase() {register_func(Int2Type<dim>());}

	///	update of values for FV1Geometry
		template <typename TElem>
		void compute(const CRFVGeometry<TElem, dim>* geo,
		             const MathVector<dim> vIPVel[maxNumSCVF],
		             number vUpShapeSh[maxNumSCVF][maxNumSH],
		             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
		             number vConvLength[maxNumSCVF])
		{
			getImpl().template compute<TElem>(geo, vIPVel, vUpShapeSh, vUpShapeIp, vConvLength);
		}

	private:
		void register_func(Int2Type<1>)
		{
			register_func<Edge>();
		}

		void register_func(Int2Type<2>)
		{
			register_func<Triangle>();
			register_func<Quadrilateral>();
		}

		void register_func(Int2Type<3>)
		{
			register_func<Tetrahedron>();
			register_func<Pyramid>();
			register_func<Prism>();
			register_func<Hexahedron>();
		}

		template <typename TElem>
		void register_func()
		{
			typedef CRFVGeometry<TElem, dim> TGeom;
			typedef void (this_type::*TFunc)(
									const TGeom* obj,
						             const MathVector<dim> vIPVel[maxNumSCVF],
						             number vUpShapeSh[maxNumSCVF][maxNumSH],
						             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
						             number vConvLength[maxNumSCVF]);

			this->template register_update_func<TGeom, TFunc>(&this_type::template compute<TElem>);
		}

	protected:
	///	access to implementation
		TImpl& getImpl() {return static_cast<TImpl&>(*this);}

	///	const access to implementation
		const TImpl& getImpl() const {return static_cast<const TImpl&>(*this);}
};

/////////////////////////////////////////////////////////////////////////////
// No Upwind
/////////////////////////////////////////////////////////////////////////////

template <int dim>
class NavierStokesCRNoUpwind
: public NavierStokesCRUpwindBase<dim, NavierStokesCRNoUpwind<dim> >
{
	public:
		typedef INavierStokesUpwind<dim> base_type;
		static const size_t maxNumSCV = base_type::maxNumSCV;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesCRNoUpwind(){this->set_shape_ip_flag(false);}

	///	update of values for CRFVGeometry
		template <typename TElem>
		void compute(const CRFVGeometry<TElem, dim>* geo,
					 const MathVector<dim> vIPVel[maxNumSCVF],
					 number vUpShapeSh[maxNumSCVF][maxNumSH],
					 number vUpShapeIp[maxNumSCVF][maxNumSCVF],
					 number vConvLength[maxNumSCVF]);
};

/////////////////////////////////////////////////////////////////////////////
// Full Upwind
/////////////////////////////////////////////////////////////////////////////

template <int dim>
class NavierStokesCRFullUpwind
: public NavierStokesCRUpwindBase<dim, NavierStokesCRFullUpwind<dim> >
{
	public:
		typedef INavierStokesUpwind<dim> base_type;
		static const size_t maxNumSCV = base_type::maxNumSCV;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesCRFullUpwind(){this->set_shape_ip_flag(false);}

	///	update of values for CRFVGeometry
		template <typename TElem>
		void compute(const CRFVGeometry<TElem, dim>* geo,
					 const MathVector<dim> vIPVel[maxNumSCVF],
					 number vUpShapeSh[maxNumSCVF][maxNumSH],
					 number vUpShapeIp[maxNumSCVF][maxNumSCVF],
					 number vConvLength[maxNumSCVF]);
};

////////////////////////////////////////////////////////////////////////////////////
// Weighted Upwind
// upwinding between full and no upwind
// shapes computed as m_weight*no_upwind_shape + (1-m_weight)*full_upwind_shape
////////////////////////////////////////////////////////////////////////////////////

template <int dim>
class NavierStokesCRWeightedUpwind
: public NavierStokesCRUpwindBase<dim, NavierStokesCRWeightedUpwind<dim> >
{
	public:
		typedef INavierStokesUpwind<dim> base_type;
		static const size_t maxNumSCV = base_type::maxNumSCV;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	protected:
		number m_weight;

	public:
	///	constructor
		NavierStokesCRWeightedUpwind(number weight)
		{
			this->set_shape_ip_flag(false);
			m_weight = weight;
		}

	///	update of values for CRFVGeometry
		template <typename TElem>
		void compute(const CRFVGeometry<TElem, dim>* geo,
					 const MathVector<dim> vIPVel[maxNumSCVF],
					 number vUpShapeSh[maxNumSCVF][maxNumSH],
					 number vUpShapeIp[maxNumSCVF][maxNumSCVF],
					 number vConvLength[maxNumSCVF]);
};

template <int dim>
class NavierStokesCRLinearProfileSkewedUpwind
: public NavierStokesCRUpwindBase<dim, NavierStokesCRLinearProfileSkewedUpwind<dim> >
{
	public:
		typedef INavierStokesUpwind<dim> base_type;
		static const size_t maxNumSCV = base_type::maxNumSCV;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesCRLinearProfileSkewedUpwind(){this->set_shape_ip_flag(false);}

	///	update of values for FV1Geometry
		template <typename TElem>
		void compute(const CRFVGeometry<TElem, dim>* geo,
					 const MathVector<dim> vIPVel[maxNumSCVF],
					 number vUpShapeSh[maxNumSCVF][maxNumSH],
					 number vUpShapeIp[maxNumSCVF][maxNumSCVF],
					 number vConvLength[maxNumSCVF]);
};

template <int dim>
class NavierStokesCRSkewedUpwind
: public NavierStokesCRUpwindBase<dim, NavierStokesCRSkewedUpwind<dim> >
{
	public:
		typedef INavierStokesUpwind<dim> base_type;
		static const size_t maxNumSCV = base_type::maxNumSCV;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesCRSkewedUpwind(){this->set_shape_ip_flag(false);}

	///	update of values for FV1Geometry
		template <typename TElem>
		void compute(const CRFVGeometry<TElem, dim>* geo,
					 const MathVector<dim> vIPVel[maxNumSCVF],
					 number vUpShapeSh[maxNumSCVF][maxNumSH],
					 number vUpShapeIp[maxNumSCVF][maxNumSCVF],
					 number vConvLength[maxNumSCVF]);
};


} // namespace NavierStokes
} // end namespace ug

#endif /* __H__UG__NAVIER_STOKES__UPWIND__ */
