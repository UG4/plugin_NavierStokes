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
#include <string>

#include "upwind_interface.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/fvcr_geom.h"
#include "lib_disc/spatial_disc/disc_util/hfvcr_geom.h"

namespace ug{
namespace NavierStokes{

/////////////////////////////////////////////////////////////////////////////
// No Upwind
/////////////////////////////////////////////////////////////////////////////

template <int dim>
class NavierStokesNoUpwind
: public INavierStokesUpwind<dim>,
  public NavierStokesUpwindRegister<FV1Geometry, dim, NavierStokesNoUpwind<dim> >,
  public NavierStokesUpwindRegister<CRFVGeometry, dim, NavierStokesNoUpwind<dim> >,
  public NavierStokesUpwindRegister<HCRFVGeometry, dim, NavierStokesNoUpwind<dim> >
{
	public:
		typedef INavierStokesUpwind<dim> base_type;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesNoUpwind(){this->set_shape_ip_flag(false);}

	///	update of values for FV1Geometry
		template <typename TElem>
		static void compute(const FV1Geometry<TElem, dim>* geo,
		             const MathVector<dim> vIPVel[maxNumSCVF],
		             number vUpShapeSh[maxNumSCVF][maxNumSH],
		             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
		             number vConvLength[maxNumSCVF]);

	///	update of values for CRFVGeometry
		template <typename TElem>
		static void compute(const CRFVGeometry<TElem, dim>* geo,
					 const MathVector<dim> vIPVel[maxNumSCVF],
					 number vUpShapeSh[maxNumSCVF][maxNumSH],
					 number vUpShapeIp[maxNumSCVF][maxNumSCVF],
					 number vConvLength[maxNumSCVF]);

		///	update of values for CRFVGeometry
		template <typename TElem>
		static void compute(const HCRFVGeometry<TElem, dim>* geo,
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
: public INavierStokesUpwind<dim>,
  public NavierStokesUpwindRegister<FV1Geometry, dim, NavierStokesFullUpwind<dim> >,
  public NavierStokesUpwindRegister<CRFVGeometry, dim, NavierStokesFullUpwind<dim> >,
  public NavierStokesUpwindRegister<HCRFVGeometry, dim, NavierStokesFullUpwind<dim> >
{
	public:
		typedef INavierStokesUpwind<dim> base_type;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesFullUpwind(){this->set_shape_ip_flag(false);}

	///	update of values for FV1Geometry
		template <typename TElem>
		static void compute(const FV1Geometry<TElem, dim>* geo,
		             const MathVector<dim> vIPVel[maxNumSCVF],
		             number vUpShapeSh[maxNumSCVF][maxNumSH],
		             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
		             number vConvLength[maxNumSCVF]);

	///	update of values for CRFVGeometry
		template <typename TElem>
		static void compute(const CRFVGeometry<TElem, dim>* geo,
					 const MathVector<dim> vIPVel[maxNumSCVF],
					 number vUpShapeSh[maxNumSCVF][maxNumSH],
					 number vUpShapeIp[maxNumSCVF][maxNumSCVF],
					 number vConvLength[maxNumSCVF]);

	///	update of values for CRFVGeometry
		template <typename TElem>
		static void compute(const HCRFVGeometry<TElem, dim>* geo,
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
/*
template <int dim>
class NavierStokesWeightedUpwind
: public INavierStokesUpwind<dim>,
  public NavierStokesUpwindRegister<CRFVGeometry, dim, NavierStokesWeightedUpwind<dim> >
{
	public:
		typedef INavierStokesUpwind<dim> base_type;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesWeightedUpwind(number weight = 0.5) : m_weight(weight)
		{
			this->set_shape_ip_flag(false);
		}

		void set_weight(number weight) {m_weight = weight;}

	///	update of values for CRFVGeometry
		template <typename TElem>
		static void compute(const CRFVGeometry<TElem, dim>* geo,
					 const MathVector<dim> vIPVel[maxNumSCVF],
					 number vUpShapeSh[maxNumSCVF][maxNumSH],
					 number vUpShapeIp[maxNumSCVF][maxNumSCVF],
					 number vConvLength[maxNumSCVF]);

	protected:
		number m_weight;
};
*/
/////////////////////////////////////////////////////////////////////////////
// Skewed Upwind
/////////////////////////////////////////////////////////////////////////////

template <int dim>
class NavierStokesSkewedUpwind
: public INavierStokesUpwind<dim>,
  public NavierStokesUpwindRegister<FV1Geometry, dim, NavierStokesSkewedUpwind<dim> >,
  public NavierStokesUpwindRegister<CRFVGeometry, dim, NavierStokesSkewedUpwind<dim> >
{
	public:
		typedef INavierStokesUpwind<dim> base_type;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesSkewedUpwind(){this->set_shape_ip_flag(false);}

	///	update of values for FV1Geometry
		template <typename TElem>
		static void compute(const FV1Geometry<TElem, dim>* geo,
		             const MathVector<dim> vIPVel[maxNumSCVF],
		             number vUpShapeSh[maxNumSCVF][maxNumSH],
		             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
		             number vConvLength[maxNumSCVF]);

	///	update of values for CRFVGeometry
		template <typename TElem>
		static void compute(const CRFVGeometry<TElem, dim>* geo,
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
: public INavierStokesUpwind<dim>,
  public NavierStokesUpwindRegister<FV1Geometry, dim, NavierStokesLinearProfileSkewedUpwind<dim> >,
  public NavierStokesUpwindRegister<CRFVGeometry, dim, NavierStokesLinearProfileSkewedUpwind<dim> >
{
	public:
		typedef INavierStokesUpwind<dim> base_type;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesLinearProfileSkewedUpwind(){this->set_shape_ip_flag(false);}

	///	update of values for FV1Geometry
		template <typename TElem>
		static void compute(const FV1Geometry<TElem, dim>* geo,
		             const MathVector<dim> vIPVel[maxNumSCVF],
		             number vUpShapeSh[maxNumSCVF][maxNumSH],
		             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
		             number vConvLength[maxNumSCVF]);

	///	update of values for CRFVGeometry
		template <typename TElem>
		static void compute(const CRFVGeometry<TElem, dim>* geo,
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
: public INavierStokesUpwind<dim>,
  public NavierStokesUpwindRegister<FV1Geometry,dim, NavierStokesPositiveUpwind<dim> >
{
	public:
		typedef INavierStokesUpwind<dim> base_type;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesPositiveUpwind(){this->set_shape_ip_flag(true);}

	///	update of values for FV1Geometry
		template <typename TElem>
		static void compute(const FV1Geometry<TElem, dim>* geo,
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
: public INavierStokesUpwind<dim>,
  public NavierStokesUpwindRegister<FV1Geometry, dim, NavierStokesRegularUpwind<dim> >
{
	public:
		typedef INavierStokesUpwind<dim> base_type;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesRegularUpwind(){this->set_shape_ip_flag(true);}

	///	update of values for FV1Geometry
		template <typename TElem>
		static void compute(const FV1Geometry<TElem, dim>* geo,
					 const MathVector<dim> vIPVel[maxNumSCVF],
					 number vUpShapeSh[maxNumSCVF][maxNumSH],
					 number vUpShapeIp[maxNumSCVF][maxNumSCVF],
					 number vConvLength[maxNumSCVF]);
};


} // namespace NavierStokes
} // end namespace ug

#endif /* __H__UG__NAVIER_STOKES__UPWIND__ */
