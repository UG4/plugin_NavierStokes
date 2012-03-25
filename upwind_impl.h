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
#include "lib_disc/local_finite_element/local_shape_function_set.h"
#include "common/util/provider.h"

#include "upwind.h"
#include "lib_algebra/algebra_type.h"

namespace ug {

/////////////////////////////////////////////////////////////////////////////
// Interface for Upwinds
/////////////////////////////////////////////////////////////////////////////

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

template <int dim>
template <typename TFVGeom, typename TAssFunc>
void
INavierStokesUpwind<dim>::
register_update_ip_vel_func(TAssFunc func)
{
//	get unique geometry id
	size_t id = GetUniqueFVGeomID<TFVGeom>();

//	make sure that there is enough space
	if((size_t)id >= m_vUpdateIPVelFunc.size())
		m_vUpdateIPVelFunc.resize(id+1, NULL);

//	set pointer
	m_vUpdateIPVelFunc[id] = (UpdateIPVelFunc)func;
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
		UG_THROW_FATAL("No update function registered for Geometry "<<id);

	if(id >= m_vUpdateIPVelFunc.size() || m_vUpdateIPVelFunc[id] == NULL)
		UG_THROW_FATAL("No update ip vel function registered for Geometry "<<id);

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
INavierStokesUpwind<dim>::
upwind_vel(size_t scvf) const
{
	UG_NSUPWIND_ASSERT(m_pCornerValue != NULL, "corner vals not set.");

//	reset result
	MathVector<dim> vel; VecSet(vel, 0.0);

//	add corner shapes
	for(size_t sh = 0; sh < num_sh(); ++sh)
		for(size_t d = 0; d < (size_t)dim; ++d)
			vel[d] += upwind_shape_sh(scvf, sh) * (*m_pCornerValue)(d, sh);

//	done if only depending on shapes
	if(!non_zero_shape_ip()) return vel;

//	compute ip vel
	for(size_t scvf2 = 0; scvf2 < num_scvf(); ++scvf)
		VecScaleAppend(vel, upwind_shape_ip(scvf, scvf2), ip_vel(scvf2));

//	return value
	return vel;
}

template <int dim>
template <typename TElem>
void
INavierStokesUpwind<dim>::
update_ip_vel(const FV1Geometry<TElem, dim>* geo, const LocalVector& vCornerValue)
{
//	set shapes
	for(size_t i = 0; i < geo->num_scvf(); ++i)
	{
	//	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(i);

	//  reset IP velocity values to zero
	    VecSet(ip_vel(i), 0.0);

	// 	Compute the Velocity in IPs
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			for(size_t d = 0; d < (size_t)dim; ++d)
				ip_vel(i)[d] += scvf.shape(sh) * vCornerValue(d, sh);
	}
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
	for(size_t i = 0; i < geo->num_scvf(); ++i)
	{
	//	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(i);

		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
		//	set upwind shape
			vUpShapeSh[i][sh] = scvf.shape(sh);
		}

	//	compute convection length
	//  \todo: (optional) A convection length is not really defined for no upwind.
	//	       but in the computation of a stabilization the term cancels, so
	//   	   we only have to ensure that the conv_lengh is non-zero
        vConvLength[i] = 1.0;
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
	for(size_t i = 0; i < geo->num_scvf(); ++i)
	{
    //	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(i);

    //  reset shapes to zero for all IPs
        for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
        	vUpShapeSh[i][sh]=0.0;

    // 	switch upwind
        const number flux = VecDot(scvf.normal(), vIPVel[i]);
        if(flux > 0.0)
        {
        	vUpShapeSh[i][scvf.from()] = 1.0;
            VecSubtract(dist, scvf.global_ip(), corners[scvf.from()]);
        }
        else
        {
        	vUpShapeSh[i][scvf.to()] = 1.0;
            VecSubtract(dist, scvf.global_ip(), corners[scvf.to()]);
        }

     // compute convection length as distance between upwind point and
     // integration point
        vConvLength[i] = VecTwoNorm(dist);
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
		UG_THROW_FATAL("GetNodeNextToCut: Cannot find cut side.");

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
	for(size_t i = 0; i < geo->num_scvf(); ++i)
	{
    //	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(i);

	// 	reset shapes to zero
 		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
 			vUpShapeSh[i][sh] = 0.0;

	// 	upwind corner
		size_t sh = 0;

	// 	find upwind node
		try{
			GetNodeNextToCut<typename FV1Geometry<TElem, dim>::ref_elem_type, dim>
					(sh, scvf.global_ip(), vIPVel[i], vCornerCoords);
		}UG_CATCH_THROW("GetSkewedUpwindShapes: Cannot find upwind node.");

	// 	set upwind corner
		vUpShapeSh[i][sh] = 1.0;

	//	compute convection length
		MathVector<dim> dist;
	    VecSubtract(dist, scvf.global_ip(), vCornerCoords[sh]);
	    vConvLength[i] = VecTwoNorm(dist);
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
	for(size_t i = 0; i < geo->num_scvf(); ++i)
	{
    //	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(i);

	// 	reset shapes to zero
 		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
 			vUpShapeSh[i][sh] = 0.0;

 		if(VecTwoNorm(vIPVel[i]) == 0.0)
 		{
 		//	no upwind -> central differences
 			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
 				vUpShapeSh[i][sh] = scvf.shape(sh);

 		//	compute convection length
 		//  \todo: (optional) A convection length is not really defined for no upwind.
 		//	       but in the computation of a stabilization the term cancels, so
 		//   	   we only have to ensure that the conv_lengh is non-zero
 			vConvLength[i] = 1.0;

 	    //	next ip
 			continue;
 		}

 	// 	side and intersection vectors
 		static const int refDim = FV1Geometry<TElem, dim>::dim;
 		size_t side = 0;
 		MathVector<dim> globalIntersection;
 		MathVector<refDim> localIntersection;

 	// 	find local intersection and side
 		try{
 			ElementSideRayIntersection<typename FV1Geometry<TElem, dim>::ref_elem_type, dim>
 			(	side, globalIntersection, localIntersection,
 				scvf.global_ip(), vIPVel[i], false /* search upwind */, vCornerCoords);
 		}UG_CATCH_THROW("GetLinearProfileSkewedUpwindShapes: Cannot find cut side.");

 	// 	get linear trial space
 		const LocalShapeFunctionSet<typename FV1Geometry<TElem, dim>::ref_elem_type>& TrialSpace =
 				LocalShapeFunctionSetProvider::
 					get<typename FV1Geometry<TElem, dim>::ref_elem_type>
 					(LFEID(LFEID::LAGRANGE, 1));

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
 			vUpShapeSh[j][co] = TrialSpace.shape(co, localIntersection);
 		}

 	//	compute conv length
		MathVector<dim> dist;
	    VecSubtract(dist, scvf.global_ip(), globalIntersection);
        vConvLength[i] = VecTwoNorm(dist);
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

//	loop all scvf
	for(size_t i = 0; i < geo->num_scvf(); ++i)
	{
	//	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(i);

	// 	reset shapes w.r.t corner value to zero
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			vUpShapeSh[i][sh] = 0.0;

	// 	reset shapes w.r.t. ip value to zero and extract ip vel
		for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
			vUpShapeIp[i][ip] = 0.0;

	//	compute flux
		vMassFlux[i] = VecProd(vIPVel[i], scvf.normal());
	}


//	2. Handle each SCV separately

	for(size_t sh = 0; sh < this->num_sh(); ++sh)
	{
	//	reset inflow, outflow
		number m_in = 0, m_out = 0;
		std::vector<size_t> vSCVIP;
		std::vector<number> vFlux;

	//	loop subcontrol volume faces
		for(size_t i = 0; i < geo->num_scvf(); ++i)
		{
		//	get SubControlVolumeFace
			const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(i);

		//	if scvf is part of the scv, add fluxes
			if(scvf.from() == sh)
			{
			//	normal directed outwards
				vSCVIP.push_back(i);
				vFlux.push_back( vMassFlux[i] );
				m_in += -1.0 * std::min(vMassFlux[i], 0.0);
				m_out += std::max(vMassFlux[i], 0.0);
			}
			else if (scvf.to() == sh)
			{
			//	normal directed inwards
				vSCVIP.push_back(i);
				vFlux.push_back( -1.0 * vMassFlux[i] );
				m_in += -1.0 * std::min(-1.0 *  vMassFlux[i], 0.0);
				m_out += std::max(-1.0 * vMassFlux[i], 0.0);
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
						vUpShapeIp[i][j] = -1.0 * vFlux[j] / F;
						sum += vUpShapeIp[i][j];
					}
				}
			//	set nodal shapes
				vUpShapeIp[i][sh] = 1.0 - sum;
			}
		}
	}

//	3. compute convection length

// 	corners of geometry
	const MathVector<dim>* vCornerCoords = geo->corners();

//	compute upwind point
	MathVector<dim> upPos; VecSet(upPos, 0.0);
	for(size_t i = 0; i < geo->num_scvf(); ++i)
	{
	//	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(i);

	//	sum up contributions
        for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
        	VecScaleAppend(upPos, vUpShapeSh[i][sh], vCornerCoords[sh]);
        for (size_t j = 0; j < geo->num_scvf(); ++j)
        	VecScaleAppend(upPos, vUpShapeIp[i][j], geo->scvf(j).global_ip());

    //	save convection length
		MathVector<dim> dist;
	    VecSubtract(dist, scvf.global_ip(), upPos);
        vConvLength[i] = VecTwoNorm(dist);
	}
}

} // end namespace ug

#endif /* NEW_STABILIZATION_IMPL_H___H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NAVIER_STOKES__FV__UPWIND_IMPL__ */
