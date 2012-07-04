/*
 * upwind_impl.cpp
 *
 *  Created on: 10.03.2011
 *      Author: andreasvogel
 */

#ifndef NEW_STABILIZATION_IMPL_H___H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__FV__UPWIND_IMPL__
#define NEW_STABILIZATION_IMPL_H___H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__FV__UPWIND_IMPL__

// for minimum
#include <limits>
#include <algorithm>

// function space, reference element
#include "lib_disc/common/geometry_util.h"
#include "lib_disc/local_finite_element/local_shape_function_set.h"
#include "common/util/provider.h"

#include "upwind.h"
#include "lib_algebra/algebra_type.h"

namespace ug {
namespace NavierStokes{

/////////////////////////////////////////////////////////////////////////////
// Interface for FV1 collocated grid upwinds
/////////////////////////////////////////////////////////////////////////////

//	register a update function for a Geometry
template <int dim>
template <typename TFVGeom, typename TAssFunc>
void
INavierStokesFV1Upwind<dim>::
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
INavierStokesFV1Upwind<dim>::
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
	m_numSh = geo.num_scv();
	UG_NSUPWIND_ASSERT(m_numScvf <= maxNumSCVF, "Invalid index");
	UG_NSUPWIND_ASSERT(m_numSh <= maxNumSH, "Invalid index");
}

///	upwind velocity
template <int dim>
MathVector<dim>
INavierStokesFV1Upwind<dim>::
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

/////////////////////////////////////////////////////////////////////////////
// No Upwind
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
template <typename TElem>
void
NavierStokesNoUpwind<TDim>::
compute(const FV1Geometry<TElem, dim>* geo,
        const MathVector<dim> vIPVel[maxNumSCVF],
        number vUpShapeSh[maxNumSCVF][maxNumSH],
        number vUpShapeIp[maxNumSCVF][maxNumSCVF],
        number vConvLength[maxNumSCVF])
{
//	set shapes
	for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
	//	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
		//	set upwind shape
			vUpShapeSh[ip][sh] = scvf.shape(sh);
		}

	//	compute convection length
	//  \todo: (optional) A convection length is not really defined for no upwind.
	//	       but in the computation of a stabilization the term cancels, so
	//   	   we only have to ensure that the conv_lengh is non-zero
        vConvLength[ip] = 1.0;
	}
}

/////////////////////////////////////////////////////////////////////////////
// Full Upwind
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
template <typename TElem>
void
NavierStokesFullUpwind<TDim>::
compute(const FV1Geometry<TElem, dim>* geo,
        const MathVector<dim> vIPVel[maxNumSCVF],
        number vUpShapeSh[maxNumSCVF][maxNumSH],
        number vUpShapeIp[maxNumSCVF][maxNumSCVF],
        number vConvLength[maxNumSCVF])
{
//	two help vectors
	MathVector<dim> dist;

// 	get corners of elem
    const MathVector<dim>* corners = geo->corners();

// 	set shapes
	for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
    //	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

    //  reset shapes to zero for all IPs
        for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
        	vUpShapeSh[ip][sh]=0.0;

    // 	switch upwind
        const number flux = VecDot(scvf.normal(), vIPVel[ip]);
        if(flux > 0.0)
        {
        	vUpShapeSh[ip][scvf.from()] = 1.0;
        	vConvLength[ip] = VecDistance(scvf.global_ip(), corners[scvf.from()]);
        }
        else
        {
        	vUpShapeSh[ip][scvf.to()] = 1.0;
        	vConvLength[ip] = VecDistance(scvf.global_ip(), corners[scvf.to()]);
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
// Skewed Upwind and Linear Profile Skewed Upwind
/////////////////////////////////////////////////////////////////////////////

/// computes the closest node to a elem side ray intersection
template <typename TRefElem, int TWorldDim>
void GetNodeNextToCut(size_t& coOut,
                      const MathVector<TWorldDim>& IP,
                      const MathVector<TWorldDim>& IPVel,
                      const MathVector<TWorldDim>* vCornerCoords)
{
//	help variables
	size_t side = 0;
	MathVector<TWorldDim> globalIntersection;
	MathVector<TRefElem::dim> localIntersection;

//	compute intersection of ray in direction of ip velocity with elem side
//	we search the ray only in upwind direction
	if(!ElementSideRayIntersection<TRefElem, TWorldDim>
		(	side, globalIntersection, localIntersection,
			IP, IPVel, false /* i.e. search upwind */, vCornerCoords))
		UG_THROW("GetNodeNextToCut: Cannot find cut side.");

//	get reference element
	static const TRefElem& rRefElem = Provider<TRefElem>::get();
	const int dim = TRefElem::dim;

// 	reset minimum
	number min = std::numeric_limits<number>::max();

// 	loop corners of side
	for(size_t i = 0; i < rRefElem.num(dim-1, side, 0); ++i)
	{
	// 	get corner
		const size_t co = rRefElem.id(dim-1, side, 0, i);

	// 	Compute Distance to intersection
		number dist = VecDistanceSq(globalIntersection, vCornerCoords[co]);

	// 	if distance is smaller, choose this node
		if(dist < min)
		{
			min = dist;
			coOut = co;
		}
	}
}

template <int TDim>
template <typename TElem>
void
NavierStokesSkewedUpwind<TDim>::
compute(const FV1Geometry<TElem, dim>* geo,
        const MathVector<dim> vIPVel[maxNumSCVF],
        number vUpShapeSh[maxNumSCVF][maxNumSH],
        number vUpShapeIp[maxNumSCVF][maxNumSCVF],
        number vConvLength[maxNumSCVF])
{
// 	corners of geometry
	const MathVector<dim>* vCornerCoords = geo->corners();

//	loop all scvf
	for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
    //	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

	// 	reset shapes to zero
 		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
 			vUpShapeSh[ip][sh] = 0.0;

	// 	upwind corner
		size_t sh = 0;

	// 	find upwind node
		try{
			GetNodeNextToCut<typename FV1Geometry<TElem, dim>::ref_elem_type, dim>
					(sh, scvf.global_ip(), vIPVel[ip], vCornerCoords);
		}UG_CATCH_THROW("GetSkewedUpwindShapes: Cannot find upwind node.");

	// 	set upwind corner
		vUpShapeSh[ip][sh] = 1.0;

	//	compute convection length
	    vConvLength[ip] = VecDistance(scvf.global_ip(), vCornerCoords[sh]);
	}
}

template <int TDim>
template <typename TElem>
void
NavierStokesLinearProfileSkewedUpwind<TDim>::
compute(const FV1Geometry<TElem, dim>* geo,
        const MathVector<dim> vIPVel[maxNumSCVF],
        number vUpShapeSh[maxNumSCVF][maxNumSH],
        number vUpShapeIp[maxNumSCVF][maxNumSCVF],
        number vConvLength[maxNumSCVF])
{
// 	corners of geometry
	const MathVector<dim>* vCornerCoords = geo->corners();

//	loop all scvf
	for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
    //	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

	// 	reset shapes to zero
 		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
 			vUpShapeSh[ip][sh] = 0.0;

 	//	if the velocity is zero, there will be no possibility to find the
 	//	cutted side. In this case we have no velocity and therefore there is
 	//	no convection. We set all upwind shapes to zero.
 		if(VecTwoNorm(vIPVel[ip]) == 0.0) continue;

 	// 	side and intersection vectors
 		static const int refDim = FV1Geometry<TElem, dim>::dim;
 		size_t side = 0;
 		MathVector<dim> globalIntersection;
 		MathVector<refDim> localIntersection;

 	// 	find local intersection and side
 		try{
 			ElementSideRayIntersection<typename FV1Geometry<TElem, dim>::ref_elem_type, dim>
 			(	side, globalIntersection, localIntersection,
 				scvf.global_ip(), vIPVel[ip], false /* search upwind */, vCornerCoords);
 		}UG_CATCH_THROW("GetLinearProfileSkewedUpwindShapes: Cannot find cut side.");

 	// 	get linear trial space
 		static const ReferenceObjectID roid = reference_element_traits<TElem>::reference_element_type::REFERENCE_OBJECT_ID;
 		const LocalShapeFunctionSet<refDim>& TrialSpace =
 				LocalShapeFunctionSetProvider::get<refDim>(roid, LFEID(LFEID::LAGRANGE, 1));

 	// 	get Reference Element
 		typedef typename FV1Geometry<TElem, dim>::ref_elem_type ref_elem_type;
 		static const ref_elem_type& rRefElem
 			= Provider<ref_elem_type>::get();

 	// 	loop corners of side
 		for(size_t j = 0; j < rRefElem.num(dim-1, side, 0); ++j)
 		{
 		// 	get corner
 			const size_t co = rRefElem.id(dim-1, side, 0, j);

 		//	evaluate trial space
 			vUpShapeSh[ip][co] = TrialSpace.shape(co, localIntersection);
 		}

 	//	compute conv length
		vConvLength[ip] = VecDistance(scvf.global_ip(), globalIntersection);
	}
}



template <int TDim>
template <typename TElem>
void
NavierStokesPositiveUpwind<TDim>::
compute(const FV1Geometry<TElem, dim>* geo,
        const MathVector<dim> vIPVel[maxNumSCVF],
        number vUpShapeSh[maxNumSCVF][maxNumSH],
        number vUpShapeIp[maxNumSCVF][maxNumSCVF],
        number vConvLength[maxNumSCVF])
{

//	1. Reset values and compute ip velocities and Compute mass fluxes at ip's

//	vector for flux values
	std::vector<number> vMassFlux(geo->num_scvf(), 0.0);
	std::vector<bool> vHasFlux(geo->num_scvf(), true);

//	loop all scvf
	const number eps = std::numeric_limits<number>::epsilon() * 10;
	size_t bNumNoFlux = 0;
	for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
	//	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

	// 	reset shapes w.r.t corner value to zero
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			vUpShapeSh[ip][sh] = 0.0;

	// 	reset shapes w.r.t. ip value to zero and extract ip vel
		for(size_t j = 0; j < geo->num_scvf(); ++j)
			vUpShapeIp[ip][j] = 0.0;

		const number normSq = VecTwoNormSq(vIPVel[ip]);
		if(fabs(normSq) <= eps)
		{
			vUpShapeSh[ip][scvf.from()] = 0.5;
			vUpShapeSh[ip][scvf.to()] = 0.5;
			vHasFlux[ip] = false;
			bNumNoFlux++;
			continue;
		}

		vMassFlux[ip] = VecProd(vIPVel[ip], scvf.normal());
		const number vel = std::sqrt(normSq);
		const number len = VecTwoNorm(scvf.normal());
		if(fabs(vMassFlux[ip] / std::sqrt(vel*len))  <= eps)
		{
			vUpShapeSh[ip][scvf.from()] = 0.5;
			vUpShapeSh[ip][scvf.to()] = 0.5;
			vHasFlux[ip] = false;
			bNumNoFlux++;
			continue;
		}
	}


//	2. Handle each SCV separately
	if(bNumNoFlux != geo->num_scvf())
	{
		for(size_t sh = 0; sh < this->num_sh(); ++sh)
		{
		//	reset inflow, outflow
			number m_in = 0, m_out = 0;
			std::vector<size_t> vSCVIP;
			std::vector<number> vFlux;

		//	loop subcontrol volume faces
			for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
			{
			//	if no flux skip
				if(!vHasFlux[ip]) continue;

			//	get SubControlVolumeFace
				const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

			//	if scvf is part of the scv, add fluxes
				if(scvf.from() == sh)
				{
				//	normal directed outwards
					vSCVIP.push_back(ip);
					vFlux.push_back( vMassFlux[ip] );
					m_in += -1.0 * std::min(vMassFlux[ip], 0.0);
					m_out += std::max(vMassFlux[ip], 0.0);
				}
				else if (scvf.to() == sh)
				{
				//	normal directed inwards
					vSCVIP.push_back(ip);
					vFlux.push_back( -1.0 * vMassFlux[ip] );
					m_in += -1.0 * std::min(-1.0 *  vMassFlux[ip], 0.0);
					m_out += std::max(-1.0 * vMassFlux[ip], 0.0);
				}
			}

		//	compute F
			number F = std::max(m_in, m_out);

		//	set shapes
			for(size_t i = 0; i < vSCVIP.size(); ++i)
			{
				if(vFlux[i] > 0)
				{
					number sum = 0.0;
					for(size_t j = 0; j < vSCVIP.size(); ++j)
					{
						if(vFlux[j] < 0)
						{
						//	set ip shapes
							sum += vUpShapeIp[vSCVIP[i]][vSCVIP[j]] = -1.0 * vFlux[j] / F;
						}
					}
				//	set nodal shapes
					vUpShapeSh[vSCVIP[i]][sh] = 1.0 - sum;
				}
			}
		}
	}

//	3. compute convection length

// 	corners of geometry
	const MathVector<dim>* vCornerCoords = geo->corners();

//	compute upwind point
	MathVector<dim> upPos;
	for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
	//	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

	//	reset upwind point
		VecSet(upPos, 0.0);

	//	sum up contributions
        for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
        	VecScaleAppend(upPos, vUpShapeSh[ip][sh], vCornerCoords[sh]);
        for (size_t j = 0; j < geo->num_scvf(); ++j)
        	VecScaleAppend(upPos, vUpShapeIp[ip][j], geo->scvf(j).global_ip());

    //	save convection length
		vConvLength[ip] = VecDistance(scvf.global_ip(), upPos);
	}
}

template <int TDim>
template <typename TElem>
void
NavierStokesRegularUpwind<TDim>::
compute(const FV1Geometry<TElem, dim>* geo,
        const MathVector<dim> vIPVel[maxNumSCVF],
        number vUpShapeSh[maxNumSCVF][maxNumSH],
        number vUpShapeIp[maxNumSCVF][maxNumSCVF],
        number vConvLength[maxNumSCVF])
{
	const number eps = std::numeric_limits<number>::epsilon() * 10;
	const int nc = geo->num_scv();

	// ONLY 2D IMPLEMENTED
	if(dim != 2) UG_THROW("RegularUpwind only implemented for 2d.");

//	loop all scvf
	for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
	//	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

	// 	reset shapes to zero
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			vUpShapeSh[ip][sh] = 0.0;
		for(size_t ip2 = 0; ip2 < scvf.num_sh(); ++ip2)
			vUpShapeIp[ip][ip2] = 0.0;

	//	if the velocity is zero, there will be no possibility to find the
	//	cutted side. In this case we have no velocity and therefore there is
	//	no convection. We set all upwind shapes to zero.
		const number normSq = VecTwoNormSq(vIPVel[ip]);
		if(fabs(normSq) <= eps) continue;

	//	flux over scvf
		number flux = VecProd(vIPVel[ip], scvf.normal());

	//	check if flux is very small (e.g. due to orthogonality of normal and
	//	flux direction)
	//	todo: Think about the small constants.
		if(fabs(flux) <= 100*eps)
		{
			// TODO: THIS IS 2D ONLY !!!!
			// the convection is parallel to the subcontrol volume surface
			flux = vIPVel[ip][0]*scvf.normal()[1] - vIPVel[ip][1]*scvf.normal()[0];
			if (flux>0)
			{
				// the velocity is pointing to the element boundary
				// take lin comb of pred and succ ip
				vUpShapeIp[ip][(ip+nc-1)%nc] = vUpShapeIp[ip][(ip+1)%nc] = 0.5;
			}
			else
			{
				// the velocity is pointing to the element midpoint
				// take lin comb of pred and succ node
				vUpShapeSh[ip][ip] = vUpShapeSh[ip][(ip+1)%nc] = 0.5;
			}
			continue;
	 	}

	//	get upwind scv
		int upwindSCV = -1;
		if(flux > 0.0) upwindSCV = scvf.from();
		else upwindSCV = scvf.to();

	//	get corresponding scv
		const typename FV1Geometry<TElem, dim>::SCV& scv = geo->scv(upwindSCV);

	// 	side and intersection vectors
		static const int refDim = FV1Geometry<TElem, dim>::dim;
		size_t side = 0;
		number lambda = 0.0;
		MathVector<dim> globalIntersection;
		MathVector<refDim> localIntersection;

		// TODO: THIS IS 2D ONLY !!!!
		//a) check for scvf intersection
		if(SCVFofSCVRayIntersection<refDim, dim>(side, lambda, globalIntersection, localIntersection,
		                            scvf.global_ip(), vIPVel[ip], false, scv.global_corners()))
		{
			int ipIntersect = (ip + 1) % nc;
			int ipIntersectOppose = (ip + nc - 1) % nc;
			if(flux > 0.0)
			{
				ipIntersect = (ip + nc - 1) % nc;
				ipIntersectOppose = (ip + 1) % nc;
			}

		// 	a.1) between two ip points
			if(lambda <= 0.5)
			{
				vUpShapeIp[ip][ipIntersectOppose] = lambda - 0.5;
				vUpShapeIp[ip][ipIntersect] = 1.0-(lambda-0.5);
			}

		// 	a.2) between ip point and elem side
			else
			{
			//	interpolation between corners and ip
				vUpShapeSh[ip][scvf.from()] = 0.5*2*(lambda-0.5);
				vUpShapeSh[ip][scvf.to()] = 0.5*2*(lambda-0.5);
				vUpShapeIp[ip][ipIntersect] = 1.0-2*(lambda-0.5);
			}

			continue;
		}
		//b) if not, on elem side intersection
		else
		{
			if(flux > 0.0)
			{
				switch (side)
				{
					case 0:
						// take linear profile between nodes on element side ip+1
						vUpShapeSh[ip][ip]		  = 1.0-0.5*lambda;
						vUpShapeSh[ip][(ip+1)%nc] = 0.5*lambda;
						break;
					case 3:
						// take linear profile between nodes on element side ip
						vUpShapeSh[ip][(ip+nc-1)%nc] = 1.0-0.5*(lambda+1.0);
						vUpShapeSh[ip][ip]			 = 0.5*(lambda+1.0);
						break;
					default:
						UG_THROW("This should not happen.");
				}
			}
			else
			{
				switch (side)
				{
					case 0:
						// take linear profile between nodes on element side ip+1
						vUpShapeSh[ip][(ip+1)%nc] = 1.0-0.5*lambda;
						vUpShapeSh[ip][(ip+2)%nc] = 0.5*lambda;
						break;
					case 3:
						// take linear profile between nodes on element side ip
						vUpShapeSh[ip][ip]		  = 1.0-0.5*(lambda+1.0);
						vUpShapeSh[ip][(ip+1)%nc] = 0.5*(lambda+1.0);
						break;
					default:
						UG_THROW("This should not happen.");
				}
			}
		}
	}

//	compute convection length

// 	corners of geometry
	const MathVector<dim>* vCornerCoords = geo->corners();

//	compute upwind point
	MathVector<dim> upPos;
	for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
	//	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

	//	reset upwind point
		VecSet(upPos, 0.0);

	//	sum up contributions
		for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
			VecScaleAppend(upPos, vUpShapeSh[ip][sh], vCornerCoords[sh]);
		for (size_t j = 0; j < geo->num_scvf(); ++j)
			VecScaleAppend(upPos, vUpShapeIp[ip][j], geo->scvf(j).global_ip());

	//	save convection length
		vConvLength[ip] = VecDistance(scvf.global_ip(), upPos);
	}
}

/////////////////////////////////////////////////////////////////////////////
// Interface for staggered CR type grid upwinds
/////////////////////////////////////////////////////////////////////////////

//	register a update function for a Geometry
template <int dim>
template <typename TFVGeom, typename TAssFunc>
void
INavierStokesCRUpwind<dim>::
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
INavierStokesCRUpwind<dim>::
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
	m_numSh = geo.num_scv();
	UG_NSUPWIND_ASSERT(m_numScvf <= maxNumSCVF, "Invalid index");
	UG_NSUPWIND_ASSERT(m_numSh <= maxNumSH, "Invalid index");
}

///	upwind velocity
template <int dim>
MathVector<dim>
INavierStokesCRUpwind<dim>::
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

/////////////////////////////////////////////////////////////////////////////
// No Upwind
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
template <typename TElem>
void
NavierStokesCRNoUpwind<TDim>::
compute(const CRFVGeometry<TElem, dim>* geo,
        const MathVector<dim> vIPVel[maxNumSCVF],
        number vUpShapeSh[maxNumSCVF][maxNumSH],
        number vUpShapeIp[maxNumSCVF][maxNumSCVF],
        number vConvLength[maxNumSCVF])
{
//	set shapes
	for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
	//	get SubControlVolumeFace
		const typename CRFVGeometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
		//	set upwind shape
			vUpShapeSh[ip][sh] = scvf.shape(sh);
		}

	//	compute convection length
	//  \todo: (optional) A convection length is not really defined for no upwind.
	//	       but in the computation of a stabilization the term cancels, so
	//   	   we only have to ensure that the conv_lengh is non-zero
        vConvLength[ip] = 1.0;
	}
}

/////////////////////////////////////////////////////////////////////////////
// Full Upwind
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
template <typename TElem>
void
NavierStokesCRFullUpwind<TDim>::
compute(const CRFVGeometry<TElem, dim>* geo,
        const MathVector<dim> vIPVel[maxNumSCVF],
        number vUpShapeSh[maxNumSCVF][maxNumSH],
        number vUpShapeIp[maxNumSCVF][maxNumSCVF],
        number vConvLength[maxNumSCVF])
{
//	two help vectors
	MathVector<dim> dist;

// 	get corners of elem
    const MathVector<dim>* corners = geo->corners();

// 	set shapes
	for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
    //	get SubControlVolumeFace
		const typename CRFVGeometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

    //  reset shapes to zero for all IPs
        for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
        	vUpShapeSh[ip][sh]=0.0;

    // 	switch upwind
        const number flux = VecDot(scvf.normal(), vIPVel[ip]);
        if(flux > 0.0)
        {
        	vUpShapeSh[ip][scvf.from()] = 1.0;
        	vConvLength[ip] = VecDistance(scvf.global_ip(), corners[scvf.from()]);
        }
        else
        {
        	vUpShapeSh[ip][scvf.to()] = 1.0;
        	vConvLength[ip] = VecDistance(scvf.global_ip(), corners[scvf.to()]);
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
// Skewed Upwind and Linear Profile Skewed Upwind
/////////////////////////////////////////////////////////////////////////////

/// computes the closest node to a elem side ray intersection
template <typename TRefElem, int TWorldDim>
void GetNodeNextToCutCR(size_t& coOut,
                      const MathVector<TWorldDim>& IP,
                      const MathVector<TWorldDim>& IPVel,
                      const MathVector<TWorldDim>* vCornerCoords)
{
//	help variables
	size_t side = 0;
	MathVector<TWorldDim> globalIntersection;
	MathVector<TRefElem::dim> localIntersection;

//	compute intersection of ray in direction of ip velocity with elem side
//	we search the ray only in upwind direction
	if(!ElementSideRayIntersection<TRefElem, TWorldDim>
		(	side, globalIntersection, localIntersection,
			IP, IPVel, false /* i.e. search upwind */, vCornerCoords))
		UG_THROW("GetNodeNextToCut: Cannot find cut side.");

//	get reference element
	static const TRefElem& rRefElem = Provider<TRefElem>::get();
	const int dim = TRefElem::dim;

// 	reset minimum
	number min = std::numeric_limits<number>::max();

// 	loop corners of side
	for(size_t i = 0; i < rRefElem.num(dim-1, side, 0); ++i)
	{
	// 	get corner
		const size_t co = rRefElem.id(dim-1, side, 0, i);

	// 	Compute Distance to intersection
		number dist = VecDistanceSq(globalIntersection, vCornerCoords[co]);

	// 	if distance is smaller, choose this node
		if(dist < min)
		{
			min = dist;
			coOut = co;
		}
	}
}

} // namespace NavierStokes
} // end namespace ug

#endif /* NEW_STABILIZATION_IMPL_H___H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__FV__UPWIND_IMPL__ */
