/*
 * moving_particle_impl.h
 *
 *  Created on: 20.01.2015
 *      Author: suze
 */

#ifndef PARTICLE_MAPPER_IMPL_H_
#define PARTICLE_MAPPER_IMPL_H_

namespace ug {
namespace NavierStokes {

///////////////////////////////////////////////////////////
// Implementation of the methods class
// 	 		'ParticleHandlerGlobal'
///////////////////////////////////////////////////////////

// see particle_handler_global_impl.h/_tools.h


///////////////////////////////////////////////////////////
// Implementation of the methods class
// 	 		'ParticleMapper'
///////////////////////////////////////////////////////////


template<typename TDomain, typename TAlgebra>
number ParticleMapper<TDomain, TAlgebra>::
compute_volume(int levIndex, size_t prtIndex)
{
	bool output = false;

	if (dim == 2 && m_spInterfaceHandlerLocal->m_vBF.size() != 2
			&& m_spInterfaceHandlerLocal->m_vBF.size() != 0) {
		UG_LOG("m_spInterfaceHandlerLocal->m_vBF.size() = " << m_spInterfaceHandlerLocal->m_vBF.size() << "\n");
		UG_THROW("oha, this->m_vBF.size() != 2 && != 0:\n");
	}

	const bool useResized = true;
	number dist = 0.0;
	number baseArea = 0.0;
	number height = 0.0;
	number volume = 0.0;
	MathVector<dim> midPoint(0.0);
	std::vector<MathVector<dim> > vCorners;

	const MathVector<dim> center = m_spCutElementHandler->get_center(
			prtIndex);
 	std::vector < interfaceBF > &vBF =
			m_spInterfaceHandlerLocal->get_boundary_faces();

// loop bf
	for (size_t i = 0; i < vBF.size(); ++i) {
		interfaceBF bf = vBF[i];

		size_t nodeID = bf.node_id();
		if (!useResized)
			nodeID = m_spInterfaceHandlerLocal->corner_orig(nodeID);

		vCorners.push_back(m_spInterfaceHandlerLocal->corner(nodeID));
		midPoint += m_spInterfaceHandlerLocal->corner(nodeID);

		baseArea += VecLength(bf.normal());

		if (output) {
			UG_LOG("i = " << i << "\n");
			UG_LOG("nodeID = " << nodeID << "\n");
			UG_LOG(
					"this->corner(" << nodeID << ") = " << m_spInterfaceHandlerLocal->corner(nodeID) << "\n");
			UG_LOG("vCorners[" << i << "] = " << vCorners[i] << "\n");
			UG_LOG("midPoint = " << midPoint << "\n");
			UG_LOG(
					"VecLength(bf.normal()) = " << VecLength(bf.normal()) << "\n");
			UG_LOG("baseArea = " << baseArea << "\n");
		}
	}

	if (dim == 2 && vCorners.size() != 2 && vCorners.size() != 0)
		UG_THROW("vCorners.size() != 2:" << vCorners.size() << "\n");

	midPoint *= 0.5;
	height = VecDistance(midPoint, center);
	dist = VecDistance(vCorners[0], vCorners[1]);

	if (fabs(baseArea - dist) > 1e-9)
		UG_THROW(
				"baseArea != dist(corners): " << baseArea << " != " << dist << "\n");

	volume = 0.5 * baseArea * height;

	if (output) {
		UG_LOG("midPoint = " << midPoint << "\n");
		UG_LOG("height = " << height << "\n");
		UG_LOG("dist = " << dist << "\n");
		UG_LOG("return:   volume = " << volume << "\n");
	}

// update volume:
	m_volume += volume;

	if (output) {
		UG_LOG("nachher (levIndex = " << levIndex << "):  m_volume = " << m_volume << "\n");
		UG_LOG("ref Volume = " << Volume(levIndex,prtIndex) << "\n");
	}
    
	return volume;

}

template<typename TDomain, typename TAlgebra>
int ParticleMapper<TDomain, TAlgebra>::
getMappingModus(size_t fct1, size_t dof1, size_t fct2, size_t dof2)
{

// if 'dof1' lies on interface, it is a particle from-DoF, i.e. momEq for particle
// if 'dof2' lies on interface, it is a particle to-DoF
	bool isPrtNode = m_spInterfaceHandlerLocal->remapped_fromInterface(dof1);
	bool isPrtConn = m_spInterfaceHandlerLocal->remapped_fromInterface(dof2);

// fluidDoF -> fluidDoF:
	if (!isPrtNode && !isPrtConn)
		return 0;
// fluidDoF -> prtDoF:
	if (!isPrtNode && isPrtConn) {
    // fluidDoF -> interfacePressure:
		if (fct2 == dim)
			return 0;  // = connection to pressure on interface!
		else
			return 1;
	}
// prtDoF -> fluidDoF:
	if (isPrtNode && !isPrtConn) {
    // prtPressure -> fluidDoF
		if (fct1 == dim)
			return 0; // = ContEq on interface!
		else
			return 2;
	}
// prtDoF -> prtDoF:
	if (isPrtNode && isPrtConn) {
      //  interfacePressure -> interfacePressure
		if (fct1 == dim && fct2 == dim)
			return 0;
		else if (fct1 == dim)
			return 1; // prtDoF in ContEq	   => = fluidDoF -> prtDoF
		else if (fct2 == dim) {
			if (m_spInterfaceHandlerLocal->elementModus() != OUTSIDE_DOM)
				return 2; // pressure in prt-ImpEq => prtDoF -> fluidDoF
			else
				return 4;
		} else
			return 3; // trans-rot/rot-trans
	}

	return -1;

}

template<typename TDomain, typename TAlgebra>
int ParticleMapper<TDomain, TAlgebra>::
getMappingModus_for2(size_t fct1, size_t dof1, size_t fct2, size_t dof2)
{
// fluid->...
	if (fct1 == _P_ || dof1 == 4) {
		// fluid->fluid:
		if (fct2 == _P_ || dof2 == 4)
			return 0;
		// fluid->particle:
		else
			return 1;
	}
// particle->...
	else {
		// particle->fluid:
		if (fct2 == _P_ || dof2 == 4)
			return 2;
		// particle->particle:
		else
			return 3;
	}

	return -1;

}

template<typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
set_extraSol(vector_type& vec, std::vector<DoFIndex> transInd, std::vector<DoFIndex> rotInd,
		const int timeIndex, const int prtIndex)
{
	for (int d = 0; d < dim; ++d)
    {
    // A. write linear velocity to data storage in 'ParticleProvider:m_vvLinearVelocity'
    //   --> if the linear velocity is given by the user, it does NOT need to be set
		if (!m_spCutElementHandler->get_DoF_modus_linear(prtIndex))
			m_spCutElementHandler->set_extraSolTrans(
					DoFRef(vec, transInd[d]), prtIndex, timeIndex, d);
        
    // B. write angular velocity to data storage in 'ParticleProvider:m_vvAngularVelocity'
    //   --> if the angular velocity is given by the user, it does NOT need to be set
		if (!m_spCutElementHandler->get_DoF_modus_anguar(prtIndex))
			m_spCutElementHandler->set_extraSolRot(
                    DoFRef(vec, rotInd[d]), prtIndex, timeIndex, d);
	}
}

template<typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
modify_GlobalSol(
		SmartPtr<VectorTimeSeries<vector_type> > vSolMod,
		ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		ConstSmartPtr<DoFDistribution> dd)
{
    bool output = false;
    
//	some checks
	if (vSol->size() != 2)
		UG_THROW("ParticleMapper::modify_GlobalSol: method needs exactly two time points.");
	if (!is_time_dependent())
		set_time_dependent(true);

	const int levIndex = get_Index(dd->grid_level(), dd);

	for (size_t p = 0; p < m_spCutElementHandler->num_particles(); ++p) {
		m_bFlagInertial[p] = true;

		if (gravitation_force()) {
			m_bFlagGravity[p] = true;
		}

    // iFF linear AND angular velocity are given by the user: no storing necessary
		if (m_spCutElementHandler->get_DoF_modus_linear(p)
				&& m_spCutElementHandler->get_DoF_modus_angular(p))
			continue;

#ifdef UG_PARALLEL
    // use size of member 'CutElementHandler_FlatTop::m_vvvElemListCut' in order to
    // indicate, whether a particle lies on a processor or not
        std::vector<grid_base_object*> ElemList =
        m_spCutElementHandler->m_vvvElemListCut[levIndex][p];
		UG_LOG("1 NavierStokes::modify_GlobalSol: ElemList.size(): " << ElemList.size() << "\n");
		if (ElemList.size() == 0) {
			UG_LOG("2 NavierStokes::modify_GlobalSol: ElemList.size(): "
                   << ElemList.size() << " => skip assembling! \n");
			continue;
		}
#endif

    // REMARK: 'get_transInd()' returns NOT YET updated indices (updated during 'movingParticle:update()' vie lua)
		std::vector < DoFIndex > transInd = m_spCutElementHandler->get_transInd(levIndex, p);
		std::vector < DoFIndex > rotInd   = m_spCutElementHandler->get_rotInd(levIndex, p);
        
        if ( output ){
		UG_LOG("in modify_GlobalSol: transInd: " << transInd[0] << "\n");
		UG_LOG("in modify_GlobalSol: rotInd: " << rotInd[0] << "\n");

		UG_LOG("VORHER: modify_GlobalSol transSol(0) = " << m_spCutElementHandler->get_transSol(p,0) << "\n");
		UG_LOG("VORHER: modify_GlobalSol transSol(1) = " << m_spCutElementHandler->get_transSol(p, 1) << "\n");
		UG_LOG("VORHER: modify_GlobalSol rotSol(0) = " << m_spCutElementHandler->get_rotSol(p,0) << "\n");
		UG_LOG("VORHER: modify_GlobalSol rot1Sol(1) = " << m_spCutElementHandler->get_rotSol(p, 1) << "\n");
        }
        
    //	loop all time points and assemble them
		for (int i = vSol->size() - 1; i >= 0; --i) {
            number solution;
    
            for (int d = 0; d < dim; ++d)
            {
            // A. if the linear velocity is given by the user, it needs NOT to be set
                if (!m_spCutElementHandler->get_DoF_modus_linear(p))
                {
                    solution = DoFRef(*vSol->solution(i), transInd[d]);
                    m_spCutElementHandler->set_extraSolTrans(solution, p, i, d);
                }
                
            // B. if the angular velocity is given by the user, is needs NOT to be set
                if (!m_spCutElementHandler->get_DoF_modus_angular(p))
                {
                    solution = DoFRef(*vSol->solution(i), rotInd[d]);
                    m_spCutElementHandler->set_extraSolRot(solution, p, i, d);
                }
            }
            
		} // end time series loop
        
        if ( output ){
		UG_LOG("NACHHER: modify_GlobalSol transSol(0) = " << m_spCutElementHandler->get_transSol(p,0) << "\n");
		UG_LOG("NACHHER: modify_GlobalSol transSol(1) = " << m_spCutElementHandler->get_transSol(p, 1) << "\n");
		UG_LOG("NACHHER: modify_GlobalSol rotSol(0) = " << m_spCutElementHandler->get_rotSol(p,0) << "\n");
		UG_LOG("NACHHER: modify_GlobalSol rot1Sol(1) = " << m_spCutElementHandler->get_rotSol(p, 1) << "\n");
        }
    } // end particle loop

}

template<typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
modify_GlobalSol(vector_type& uMod,
		const vector_type& u, ConstSmartPtr<DoFDistribution> dd)
{
    bool output = false;
    
	if (is_time_dependent())
		set_time_dependent(false);

	const int levIndex = get_Index(dd->grid_level(), dd);

	for (size_t p = 0; p < m_spCutElementHandler->num_particles(); ++p)
    {
		if (gravitation_force()) {
			m_bFlagGravity[p] = true;
		}

    // IFF linear AND angular velocity are given by the user: no need for storing
		if (m_spCutElementHandler->get_DoF_modus_linear(p)
         && m_spCutElementHandler->get_DoF_modus_angular(p)) {
			continue;
		}


#ifdef UG_PARALLEL
    // use size of member 'CutElementHandler_FlatTop::m_vvvElemListCut' in order to
    // indicate, whether a particle lies on a processor or not
        std::vector<grid_base_object*> ElemList =
        m_spCutElementHandler->m_vvvElemListCut[levIndex][p];
		UG_LOG("1 NavierStokes::modify_GlobalSol: ElemList.size(): " << ElemList.size() << "\n");
		if (ElemList.size() == 0) {
			UG_LOG("2 NavierStokes::modify_GlobalSol: ElemList.size(): "
                   << ElemList.size() << " => skip assembling! \n");
			continue;
		}
#endif

    // REMARK: 'get_transInd()' returns NOT YET updated indices (updated during 'movingParticle:update()' vie lua)
		std::vector < DoFIndex > transInd = m_spCutElementHandler->get_transInd(levIndex, p);
		std::vector < DoFIndex > rotInd   = m_spCutElementHandler->get_rotInd(levIndex, p);
        
        if (output ){
		UG_LOG("VORHER: modify_GlobalSol transSol(0) = " << m_spCutElementHandler->get_transSol(p,0) << "\n");
		UG_LOG("VORHER: modify_GlobalSol transSol(1) = " << m_spCutElementHandler->get_transSol(p, 1) << "\n");
		UG_LOG("VORHER: modify_GlobalSol rotSol(0) = " << m_spCutElementHandler->get_rotSol(p,0) << "\n");
		UG_LOG("VORHER: modify_GlobalSol rot1Sol(1) = " << m_spCutElementHandler->get_rotSol(p, 1) << "\n");
        }
    //	set_extraSol(u, transInd, rotInd, 0, p);
		number solution;
		for (int d = 0; d < dim; ++d)
        {
        // A. if the linear velocity is given by the user, is needs NOT to be set
			if (!m_spCutElementHandler->get_DoF_modus_linear(p)) {
				solution = DoFRef(u, transInd[d]);
				m_spCutElementHandler->set_extraSolTrans(solution, p, 0, d);
 			}

        // B. if the angular velocity is given by the user, is needs NOT to be set
			if (!m_spCutElementHandler->get_DoF_modus_angular(p)) {
				solution = DoFRef(u, rotInd[d]);
				m_spCutElementHandler->set_extraSolRot(solution, p, 0, d);
 			}

		}
        if (output ){
		UG_LOG("NACHHER: modify_GlobalSol transSol(0) = " << m_spCutElementHandler->get_transSol(p,0) << "\n");
		UG_LOG("NACHHER: modify_GlobalSol transSol(1) = " << m_spCutElementHandler->get_transSol(p, 1) << "\n");
		UG_LOG("NACHHER: modify_GlobalSol rotSol(0) = " << m_spCutElementHandler->get_rotSol(p,0) << "\n");
		UG_LOG("NACHHER: modify_GlobalSol rot1Sol(1) = " << m_spCutElementHandler->get_rotSol(p, 1) << "\n");
        }
	} // end particle loop

}

template<typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
map_local_data(LocalVector& d)
{
	LocalIndices ind = d.get_indices();

// initialize tmp data
	LocalVector tmpLocD;
	tmpLocD.resize(d.get_indices());

// operator= exists only for single value! => set to zero and add
	tmpLocD = 0.0;
	tmpLocD += d;

// copy old local vector to new local algebra
	for (size_t fct = 0; fct < ind.num_fct(); ++fct)
		for (size_t dof = 0; dof < ind.num_dof(fct); ++dof) {
			size_t _dof = m_spInterfaceHandlerLocal->corner_orig(dof);
			d(fct, dof) = tmpLocD(fct, _dof);
		}

}

template<typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
modify_LocalData(LocalMatrix& locJ, LocalVector& locU, ConstSmartPtr<DoFDistribution> dd)
{
// 	INSIDE_DOM: do nothing
	if (m_spInterfaceHandlerLocal->elementModus() == INSIDE_DOM)
		return;

// resize local indices:
	if (m_spInterfaceHandlerLocal->elementModus() == CUT_BY_2_INTERFACE
			&& m_spInterfaceHandlerLocal->interface_id_all().size() == 2)
	{
		if (!m_spInterfaceHandlerLocal->StdFV_assembling())
			resize_local_indices(locU, 5);
		else
			resize_local_indices(locU);
	}
	else
		resize_local_indices(locU);

// resize local data:
	const LocalIndices& ind = m_spInterfaceHandlerLocal->get_local_indices();
	locU.resize(ind);
	locJ.resize(ind);

	return;
}

template<typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
modify_LocalData(LocalVectorTimeSeries& uT, LocalMatrix& locJ, LocalVector& locU,
		ConstSmartPtr<DoFDistribution> dd)
{
// 	INSIDE_DOM: do nothing
	if (m_spInterfaceHandlerLocal->elementModus() == INSIDE_DOM)
		return;

	for (int i = uT.size() - 1; i >= 0; --i)
    {
    // resize local indices:
		if (m_spInterfaceHandlerLocal->elementModus() == CUT_BY_2_INTERFACE
				&& m_spInterfaceHandlerLocal->interface_id_all().size() == 2) {
			if (!m_spInterfaceHandlerLocal->StdFV_assembling())
				resize_local_indices(locU, 5);
			else
				resize_local_indices(locU);
		} else
			resize_local_indices(locU);

    // resize local data:
		const LocalIndices& ind =
				m_spInterfaceHandlerLocal->get_local_indices();
		locU.resize(ind);
		locJ.resize(ind);
		uT.solution(i).resize(ind);

	}

	return;
}

template<typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
set_QuadriSol(LocalVector& locU, LocalVector& locD)
{

	if (locU.num_all_dof(0) != 5)
		UG_THROW(
				"ParticleMapper:set_QuadriSol(): number of DoFs should be 5, but is " << locU.num_all_dof(0) << "\n");

// 3. get velocities of the 2 particles:
	MathVector<dim> transSol1 = m_spCutElementHandler->get_transSol(0, 0);
	MathVector<dim> rotSol1 = m_spCutElementHandler->get_rotSol(0, 0);
	MathVector<dim> transSol2 = m_spCutElementHandler->get_transSol(1, 0);
	MathVector<dim> rotSol2 = m_spCutElementHandler->get_rotSol(1, 0);

// resize local data as done during 'modify_LocalData'
//	 --> but there with num_co = 3 for the Triangle with inside node!
	for (size_t fct = 0; fct < locU.num_fct() - 1; ++fct)
		for (size_t dof = 0; dof < locU.num_all_dof(fct); ++dof) {
			MathMatrix<dim, dim> rotationMatCo =
					m_spCutElementHandler->get_rotationMat(
							m_spInterfaceHandlerLocal->radial_at_co(dof));

			UG_LOG("radial(" << dof << ": " << m_spInterfaceHandlerLocal->radial_at_co(dof) << "\n");

			// write solution of particle with prtIndex = 0:
			if (dof < 2) {
				locU(fct, dof) = transSol1[fct];
				for (int d = 0; d < dim; ++d)
					locU(fct, dof) += rotationMatCo[fct][d] * rotSol1[d];
			}
			// write solution of particle with prtIndex = 1:
			else {
				locU(fct, dof) = transSol2[fct];
				for (int d = 0; d < dim; ++d)
					locU(fct, dof) += rotationMatCo[fct][d] * rotSol2[d];
			}
		}

	m_spInterfaceHandlerLocal->set_QuadriSol(locD, locU);

	UG_LOG("after set_QuadriSol: locU = \n" << locU << "\n");
	UG_LOG("after set_QuadriSol: locD = \n" << locD << "\n");

	UG_LOG("after set_QuadriSol: m_quadriLocU = \n" << m_spInterfaceHandlerLocal->m_quadriLocU << "\n");
	UG_LOG("after set_QuadriSol: m_quadriLocD = \n" << m_spInterfaceHandlerLocal->m_quadriLocD << "\n");

}

template<typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
modify_LocalData(LocalVector& locD,
		LocalVector& tmpLocD, LocalVector& locU,
		ConstSmartPtr<DoFDistribution> dd)
{
///////////////////////////////////////////////////////////////////////////////
// 	A. INSIDE_DOM: do nothing
	if (m_spInterfaceHandlerLocal->elementModus() == INSIDE_DOM)
		return;

///////////////////////////////////////////////////////////////////////////////
// B. Resize local data: for QUADRILATERALS (2d) or HEXA, PRISM, PYRAMID (3d)
    
// resize local indices:
	if (m_spInterfaceHandlerLocal->elementModus() == CUT_BY_2_INTERFACE
			&& m_spInterfaceHandlerLocal->interface_id_all().size() == 2) {
		if (!m_spInterfaceHandlerLocal->StdFV_assembling())
			resize_local_indices(locU, 5);
		else
			resize_local_indices(locU);
	} else
		resize_local_indices(locU);

// resize local data:
	const LocalIndices& ind = m_spInterfaceHandlerLocal->get_local_indices();
	locU.resize(ind);
	locD.resize(ind);
	tmpLocD.resize(ind);

///////////////////////////////////////////////////////////////////////////////
// C. Write the particle rigid body velocity into resized solution 'locU'
    
	int prtIndex = get_prtIndex();
	const MathVector<dim> center = m_spCutElementHandler->get_center(prtIndex);

	MathVector<dim> transSol = m_spCutElementHandler->get_transSol(prtIndex,0);
	MathVector<dim> rotSol   = m_spCutElementHandler->get_rotSol(prtIndex, 0);

//  OUTSIDE_DOM: write solution for ALL corners
	if (m_spInterfaceHandlerLocal->elementModus() == OUTSIDE_DOM) {
		for (size_t fct = 0; fct < locU.num_fct() - 1; ++fct)
			for (size_t dof = 0; dof < locU.num_all_dof(fct); ++dof) {

				MathVector<dim> radialCo;
				VecSubtract(radialCo, m_spInterfaceHandlerLocal->corner(dof),
						center);
				MathMatrix<dim, dim> rotationMatCo =
						m_spCutElementHandler->get_rotationMat(radialCo);

				//	set solution
				locU(fct, dof) = transSol[fct];
				for (int d = 0; d < dim; ++d)
					locU(fct, dof) += rotationMatCo[fct][d] * rotSol[d];

			}
	}

// 	CUT_BY_INTERFACE: write solution for interface corners
	if (m_spInterfaceHandlerLocal->elementModus() == CUT_BY_INTERFACE
			|| (m_spInterfaceHandlerLocal->elementModus() == CUT_BY_2_INTERFACE
					&& m_spInterfaceHandlerLocal->StdFV_assembling())) {
		map_local_data(locU);

		for (size_t fct = 0; fct < locU.num_fct() - 1; ++fct) {
			for (size_t i = 0;
					i < m_spInterfaceHandlerLocal->interface_id_all().size();
					++i) {
				size_t interfaceID = m_spInterfaceHandlerLocal->interface_id(i);

				MathMatrix<dim, dim> rotationMatCo =
						m_spCutElementHandler->get_rotationMat(
								m_spInterfaceHandlerLocal->radial_at_co(
										interfaceID));

            //	set solution
				locU(fct, interfaceID) = transSol[fct];

				for (int d = 0; d < dim; ++d)
					locU(fct, interfaceID) += rotationMatCo[fct][d] * rotSol[d];

			}
		}
	}

// CUT_BY_2_INTERFACE:
	if (m_spInterfaceHandlerLocal->elementModus() == CUT_BY_2_INTERFACE) {
		if (!m_spInterfaceHandlerLocal->StdFV_assembling())
			set_QuadriSol(locU, locD);
	}

}

template<typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
modify_LocalData(
		LocalVectorTimeSeries& uT, LocalVector& locD, LocalVector& tmpLocD,
		LocalVector& locU, ConstSmartPtr<DoFDistribution> dd, size_t t)
{

// 	'modify_LocalData()' is called twice (for t=0, t=1) in elem_disc_assemble_util;
//   BUT: 'map_local_data()' may be called only once!
	if (t == 1)
		return;

///////////////////////////////////////////////////////////////////////////////
// 	A. INSIDE_DOM: do nothing
    
	if (m_spInterfaceHandlerLocal->elementModus() == INSIDE_DOM)
		return;
    
///////////////////////////////////////////////////////////////////////////////
// B. Resize local data: for QUADRILATERALS (2d) or HEXA, PRISM, PYRAMID (3d)
    
// resize local indices:
	if (m_spInterfaceHandlerLocal->elementModus() == CUT_BY_2_INTERFACE
			&& m_spInterfaceHandlerLocal->interface_id_all().size() == 2) {
		if (!m_spInterfaceHandlerLocal->StdFV_assembling())
			resize_local_indices(locU, 5);
		else
			resize_local_indices(locU);
	} else
		resize_local_indices(locU);

// resize local data:
	const LocalIndices& ind = m_spInterfaceHandlerLocal->get_local_indices();
	locU.resize(ind);
	locD.resize(ind);
	tmpLocD.resize(ind);

///////////////////////////////////////////////////////////////////////////////
// C. Write the particle rigid body velocity into resized solution 'locU'
    
    int prtIndex = get_prtIndex();
    const MathVector<dim> center = m_spCutElementHandler->get_center(prtIndex);

//	loop all time points and assemble them
	for (int i = uT.size() - 1; i >= 0; --i) {

		uT.solution(i).resize(ind);

		LocalVector* vecMod = &uT.solution(i);

		MathVector<dim> transSol = m_spCutElementHandler->get_transSol(prtIndex, i);
		MathVector<dim> rotSol   = m_spCutElementHandler->get_rotSol(prtIndex,i);

    // OUTSIDE_DOM: write solution for ALL corners
		if (m_spInterfaceHandlerLocal->elementModus() == OUTSIDE_DOM) {
			for (size_t fct = 0; fct < locU.num_fct() - 1; ++fct)
				for (size_t dof = 0; dof < locU.num_all_dof(fct); ++dof) {
					MathVector<dim> radialCo;
					VecSubtract(radialCo,
							m_spInterfaceHandlerLocal->corner(dof), center);
					MathMatrix<dim, dim> rotationMatCo =
							m_spCutElementHandler->get_rotationMat(
									radialCo);

                //	set solution
					(*vecMod)(fct, dof) = transSol[fct];
					for (int d = 0; d < dim; ++d)
						(*vecMod)(fct, dof) += rotationMatCo[fct][d]
								* rotSol[d];

				}
		}

    // CUT_BY_INTERFACE: write solution for interface corners
		if (m_spInterfaceHandlerLocal->elementModus() == CUT_BY_INTERFACE
				|| (m_spInterfaceHandlerLocal->elementModus()
						== CUT_BY_2_INTERFACE
						&& m_spInterfaceHandlerLocal->StdFV_assembling())) {
			map_local_data(*vecMod);

			for (size_t fct = 0; fct < locU.num_fct() - 1; ++fct) {
				for (size_t i = 0;
						i < m_spInterfaceHandlerLocal->interface_id_all().size();
						++i) {
					size_t interfaceID =
							m_spInterfaceHandlerLocal->interface_id(i);
					MathMatrix<dim, dim> rotationMatCo =
							m_spCutElementHandler->get_rotationMat(
									m_spInterfaceHandlerLocal->radial_at_co(
											interfaceID));

                //	set solution
					(*vecMod)(fct, interfaceID) = transSol[fct];
					for (int d = 0; d < dim; ++d)
						(*vecMod)(fct, interfaceID) += rotationMatCo[fct][d]
								* rotSol[d];
				}
			}

		}

    // CUT_BY_2_INTERFACE:
		if (m_spInterfaceHandlerLocal->elementModus() == CUT_BY_2_INTERFACE) {
			if (!m_spInterfaceHandlerLocal->StdFV_assembling())
				UG_THROW("modify_LocalData(): CUT_BY_2_INTERFACE\n"); //set_QuadriSol(locU, locD);
		}

	} // end time loop

}


template<typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
set_identity_mat_constraint(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd)
{
// get levIndex
	const int levIndex = get_Index(dd->grid_level(), dd);

// get local indices
	const LocalIndices& rowInd = lmat.get_row_indices();

	for (size_t fct1 = 0; fct1 < lmat.num_all_row_fct(); ++fct1)
		for (size_t dof1 = 0; dof1 < lmat.num_all_row_dof(fct1); ++dof1) {
			const size_t _dof1 = m_spInterfaceHandlerLocal->corner_orig(dof1);
			const size_t rowIndex = rowInd.index(fct1, _dof1);
			const size_t rowComp = rowInd.comp(fct1, _dof1);

        // 	m_spCutElementHandler->is_extraDoF returns false, if particle velocity is given by the model
			bool is_extraDoF = m_spCutElementHandler->is_extraDoF(
					DoFIndex(rowIndex, rowComp), levIndex);

        // IFF node is no transInd/rotInd-node and no interface node!
			bool is_interfaceDoF = true;
        // on interface, only pressure is real DoF:
			if (fct1 != dim)
				is_interfaceDoF = false;
			else if (m_spInterfaceHandlerLocal->StdFV_assembling()) {

				is_interfaceDoF =
						m_spInterfaceHandlerLocal->remapped_fromInterface(
								_dof1);
				//is_interfaceDoF = true; //m_spInterfaceHandlerLocal->lies_onInterface(_dof1);

			} else
				is_interfaceDoF =
						m_spInterfaceHandlerLocal->remapped_fromInterface(
								_dof1);

        // DoFs on interface OR for translational/angular velocity need to remain free
			if (is_interfaceDoF || is_extraDoF)
				continue;
			else
				BlockRef(mat(rowIndex, rowIndex), rowComp, rowComp) = 1.0;
		}

}

template<typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
add_local_mat_to_global(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd)
{
// reset print modus for cut element data to 'false'
//  => only printed during assemling of defect, i.e. ONE loop over all element
    m_spInterfaceHandlerLocal->set_print_cutElemData(false);
    
	ElementModus modus = m_spInterfaceHandlerLocal->elementModus();
	LocalMatrix locJ = lmat;

	switch (modus) {
	case OUTSIDE_DOM: // no DoFs: set identity-rows
 		set_identity_mat_constraint(mat, lmat, dd); //IFF element does NOT include transInd/rotInd-node!
		break;
	case INSIDE_DOM: // call usual local-to-global-mapping
		AddLocalMatrixToGlobal(mat, lmat);
		break;
	case CUT_BY_INTERFACE: // call adapted local-to-global-mapping
		add_local_mat_to_global_interface(mat, lmat, dd);
		break;
	case CUT_BY_2_INTERFACE: // call adapted local-to-global-mapping
		add_local_mat_to_global_interface_for2(mat, lmat, dd);
		break;
	default:
		throw(UGError("Error in IInterfaceMapper::add_local_mat_to_global()!"));

	}

}

template<typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
add_local_vec_to_global(vector_type& vec, const LocalVector& lvec,
		ConstSmartPtr<DoFDistribution> dd)
{
 
	ElementModus modus = m_spInterfaceHandlerLocal->elementModus();

	switch (modus) {
	case OUTSIDE_DOM: // no DoFs: no contribution to defect
		break;
	case INSIDE_DOM: // call usual local-to-global-mapping
		AddLocalVector(vec, lvec);
		break;
	case CUT_BY_INTERFACE: // call adapted local-to-global-mapping
		add_local_vec_to_global_interface(vec, lvec, dd);
		break;
	case CUT_BY_2_INTERFACE: // call adapted local-to-global-mapping
		add_local_vec_to_global_interface_for2(vec, lvec, dd);
		break;
	default:
		throw(UGError("Error in IInterfaceMapper::add_local_vec_to_global()!"));

	}

}


template<typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
add_mass_part_jac(matrix_type& mat,	std::vector<DoFIndex> transInd, std::vector<DoFIndex> rotInd,
		const int levIndex, const int prtIndex)
{
	number volume;

	if (m_bVolumeCompExact)
		volume = Volume(levIndex, prtIndex);
	else
		volume = compute_volume(levIndex, prtIndex);

// add mass term to global matrix (instead of using IConstraint::adjust_jacobian() )
	for (int d = 0; d < dim; ++d) {
		DoFRef(mat, transInd[d], transInd[d]) += m_spInterfaceHandlerLocal->get_density(prtIndex)*volume;
    // for dim = 2: the second DoF of rotInd is unused!
		if (dim == 3 || d == 0)
			DoFRef(mat, rotInd[d], rotInd[d]) += MomOfInertia(levIndex, prtIndex, volume);
	}

// adding mass/inertia terms only ONCE!!
    if (m_bVolumeCompExact)
        m_bFlagInertial[prtIndex] = false;
    
}



// call for 'CUT_BY_2_INTERFACE'
template<typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
add_local_mat_to_global_interface_for2(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd)
{
// REMARK: not finally tested!
    UG_THROW("ParticleMapper::add_local_mat_to_global_interface_for2(): not finally tested!!!!\n");

// get data
	const int levIndex = get_Index(dd->grid_level(), dd);
	std::vector < DoFIndex > transInd1, rotInd1, transInd2, rotInd2;
	m_spCutElementHandler->get_global_indices(transInd1, rotInd1, levIndex,0);
	m_spCutElementHandler->get_global_indices(transInd2, rotInd2, levIndex,1);


	if (!UsualAss())
    {
		if (m_spInterfaceHandlerLocal->StdFV_assembling())
			add_local_mat_to_global_FT_for2_StdFV(mat, lmat, transInd1, rotInd1,
					transInd2, rotInd2);
		else
			add_local_mat_to_global_FT_for2(mat, lmat, transInd1, rotInd1,
					transInd2, rotInd2);
    }
    else {
		UG_THROW("add_local_mat_to_global_interface_for2(): -> not implemented!\n");
    }
}

// call for 'CUT_BY_INTERFACE'
template<typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
add_local_mat_to_global_interface(matrix_type& mat, const LocalMatrix& lmat,
		ConstSmartPtr<DoFDistribution> dd)
{
// get data
	const int levIndex = get_Index(dd->grid_level(), dd);
    int prtIndex = get_prtIndex();

  	std::vector < DoFIndex > transInd, rotInd;
	m_spCutElementHandler->get_global_indices(transInd, rotInd, levIndex, prtIndex);

// A. add mass part to jacobian
	if (m_bFlagInertial[prtIndex])
		add_mass_part_jac(mat, transInd, rotInd, levIndex, prtIndex);

// B. add local to global
    if (!UsualAss()){
		add_local_mat_to_global_FT(mat, lmat, transInd, rotInd);
    }
    else {
        UG_THROW("ParticleMapper::add_local_mat_to_global_usual(): -> not implemented!\n");
    }
    
}

    
    
template<typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
add_glowinski_repulsive_force_rhs(vector_type& vec, std::vector<DoFIndex> transInd,
                std::vector<DoFIndex> rotInd,const int levIndex, const int prtIndex)
{
////////////////////////////////////////////////////////////////
// compute repulsive force:
////////////////////////////////////////////////////////////////
        
    for (size_t i = 0; i < num_particles(); ++i) {
        if ((int)i == prtIndex) {continue;}
    // Get particle values
        MathVector<dim> center_i = m_spCutElementHandler->get_center(i);
        MathVector<dim> center_prtIndex = m_spCutElementHandler->get_center(prtIndex);
        double Mass_i = Mass(levIndex,i);
        double Mass_prtIndex = Mass(levIndex,prtIndex);
        
    // Calculate distance vector between center points
        double dist = 0;
        for (size_t d = 0; d < dim; ++d){
            dist += (center_i[d]-center_prtIndex[d])*(center_i[d]-center_prtIndex[d]);
        }
        dist = sqrt(dist);
        MathVector<dim> distanceVec = center_prtIndex;
        distanceVec -=  center_i;
        distanceVec /= dist;
            
    // Calculate force
        double force = (Mass_i+Mass_prtIndex)/2*981/(0.5*m_epsilon);
        double max1 = std::max(0.0,-(dist-m_spCutElementHandler->get_radius(i)-m_spCutElementHandler->get_radius(prtIndex)-m_rho)/(m_rho));
        force *= pow(max1,2);
            
    // Scale with time step
        number timeScale = 1.0;
        if (is_time_dependent()) {
            timeScale = m_dt;
            force *= timeScale;
        }
            
        //			UG_LOG("Max/original force is " << force << ".\n");
            
        UG_LOG("Unscaled Max/original repulsive force is "<<force<<".\n");
            
            
        DoFRef(vec, transInd[0]) -= force*distanceVec[0];
        DoFRef(vec, transInd[1]) -= force*distanceVec[1];
        if (dim == 3)
            DoFRef(vec, transInd[2]) -= force*distanceVec[2];
        //		UG_LOG("Max/original force = (" << Mass_prtIndex*force*distanceVec[0] << ", " << Mass_prtIndex*force*distanceVec[1] << "), " << Mass_prtIndex*force << "\n");
            
    }
}
    
template<typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
add_minimum_correction_force_rhs(vector_type& vec,std::vector<DoFIndex> transInd,
            std::vector<DoFIndex> rotInd, const int levIndex, const int prtIndex)
{
////////////////////////////////////////////////////////////////
// compute repulsive force:
////////////////////////////////////////////////////////////////
    for (size_t i = 0; i < num_particles(); ++i) {
        if (!m_bGlowRepulsiveForce && (i != prtIndex)) {

            if (i == prtIndex) {continue;}
            // Get particle values
            MathVector<dim> center_i = m_spCutElementHandler->get_center(i);
            MathVector<dim> center_prtIndex = m_spCutElementHandler->get_center(prtIndex);
            double Mass_i = Mass(levIndex,i);
            double Mass_prtIndex = Mass(levIndex,prtIndex);

            // Calculate distance vector between center points
            double dist = 0;
            for (size_t d = 0; d < dim; ++d){
                dist += (center_i[d]-center_prtIndex[d])*(center_i[d]-center_prtIndex[d]);
            }
            dist = sqrt(dist);
            MathVector<dim> distanceVec = center_prtIndex;
            distanceVec -=  center_i;
            distanceVec /= dist;
            if (m_bForceLog) {
                if (dim == 2) {
                    UG_LOG("distanceVec is ("<< distanceVec[0] <<","<< distanceVec[1] <<").\n");
                } else {
                    UG_LOG("distanceVec is ("<< distanceVec[0] <<","<< distanceVec[1] <<","<< distanceVec[2] <<").\n");
                }
            }

            // Calculate force
            double force = (Mass_i+Mass_prtIndex)/2*m_gravityConst/(m_epsilon);
            double max1 = std::max(0.0,-(dist-m_spCutElementHandler->get_radius(i)-m_spCutElementHandler->get_radius(prtIndex)-m_rho)/(m_rho));
            force *= pow(max1,2);

            // Scale with time step
            number timeScale = 1.0;
            if (is_time_dependent()) {
                timeScale = m_dt;
                force *= timeScale;
            }

            if (m_bForceLog) {
                UG_LOG("Unscaled Max/original repulsive force is "<<force<<".\n");
            }
        }
    }
    if (m_bMinimumCorrectionForce){
        double Mass_prtIndex = Mass(levIndex,prtIndex);
        MathVector<dim> Delta_r = m_repulsiveDistance[prtIndex];
        MathVector<dim> force;
        VecScale(force, Delta_r, 1/(m_dt*m_dt));

        DoFRef(vec, transInd[0]) -= Mass_prtIndex*force[0];
        DoFRef(vec, transInd[1]) -= Mass_prtIndex*force[1];
        if (dim == 3)
            DoFRef(vec, transInd[2]) -= Mass_prtIndex*force[2];
        if (m_bForceLog) {
            UG_LOG("Unscaled equivalent repulsive force is "<<Mass_prtIndex << "*" << force <<".\n");
        }
    }
    m_bForceLog = false;
}

    
template<typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
add_repulsive_force_rhs(vector_type& vec, std::vector<DoFIndex> transInd,
    std::vector<DoFIndex> rotInd, const int levIndex, const int prtIndex)
{
    ////////////////////////////////////////////////////////////////
    // compute repulsive force:
    ////////////////////////////////////////////////////////////////
        
    double force = m_repulsiveForce;
        
    for (size_t i = 0; i < num_particles(); ++i) {
        if ((int)i == prtIndex) {continue;}
    // Get particle values
        MathVector<dim> center_i = m_spCutElementHandler->get_center(i);
        MathVector<dim> center_prtIndex = m_spCutElementHandler->get_center(prtIndex);
        double Mass_prtIndex = Mass(levIndex,prtIndex);
        
    // Calculate distance vector between center points
        double dist = 0;
        for (size_t d = 0; d < dim; ++d){
            dist += (center_i[d]-center_prtIndex[d])*(center_i[d]-center_prtIndex[d]);
        }
        dist = sqrt(dist);
        MathVector<dim> distanceVec = center_prtIndex;
        distanceVec -=  center_i;
        distanceVec /= dist;
            
        UG_LOG("Unscaled equivalent repulsive force is "<<force<<".\n");
    // Scale with time step
        number timeScale = 1.0;
        if (is_time_dependent()) {
            timeScale = m_dt;
            force *= timeScale;
        }
        
        DoFRef(vec, transInd[0]) -= Mass_prtIndex*force*distanceVec[0];
        DoFRef(vec, transInd[1]) -= Mass_prtIndex*force*distanceVec[1];
        if (dim == 3)
            DoFRef(vec, transInd[2]) -= Mass_prtIndex*force*distanceVec[2];
            
        if (m_bForceLog) {
            UG_LOG("Unscaled force is set to "<< force << ".\n");
        }
    }
        
}

    
template<typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
set_gravitational_rhs(vector_type& vec, std::vector<DoFIndex> transInd, std::vector<DoFIndex> rotInd,
		const int levIndex, const int prtIndex)
{
    bool logGravity = false;

// some initial checks: Parameters need to be set via lua-script!
    if (m_gravityConst == 0.0)
        UG_THROW("m_gravityConst not defined!\n");
    if (is_time_dependent() && m_dt == 0.0)
        UG_THROW("time step 'm_dt' not defined! Needs to be set via lua-script: 'set_time_step(dt)'...\n");
 
	if (logGravity)
    {
		UG_LOG("//////////////////////////////// - log_Gravity (" << prtIndex << ") - ///////////////////////////////\n");
    
 		UG_LOG("VORHER: defect(trans,0): " << DoFRef(vec, transInd[0]) << "\n");
 		UG_LOG("VORHER: defect(trans,1): " << DoFRef(vec, transInd[1]) << "\n");
        if (dim == 3)
            UG_LOG("VORHER: defect(trans,2): " << DoFRef(vec, transInd[2]) << "\n");

 		UG_LOG("VORHER: defect(rot,0): " << DoFRef(vec, rotInd[0]) << "\n");
 		UG_LOG("VORHER: defect(rot,1): " << DoFRef(vec, rotInd[1]) << "\n");
        if (dim == 3)
            UG_LOG("VORHER: defect(rot,2): " << DoFRef(vec, rotInd[2]) << "\n");
        
        UG_LOG("m_gravityConst = " << m_gravityConst << "\n");

    }
    
    ////////////////////////////////////////////////////////////////
    // compute gravitational force to rhs:
    ////////////////////////////////////////////////////////////////
    
	number volume;
	if (m_bVolumeCompExact)
		volume = Volume(levIndex, prtIndex);
	else
		volume = compute_volume(levIndex, prtIndex);

	const number gravitationalMass = Mass(levIndex, prtIndex, volume);

	number gravityForce = m_gravityConst * gravitationalMass;
	if (dim == 3)
		gravityForce = m_gravityConst * gravitationalMass;
   
    number timeScale = 1.0;
    if (is_time_dependent())
        timeScale = m_dt;
    
   ////////////////////////////////////////////////////////////////
   // add gravitational force to rhs:
   ////////////////////////////////////////////////////////////////
	if (dim == 2)
		DoFRef(vec, transInd[0]) -= timeScale * gravityForce;
	if (dim == 3)
		DoFRef(vec, transInd[2]) -= timeScale * gravityForce;

    
	if (logGravity)
    {
 		UG_LOG("gravitationalMass: " << gravitationalMass << "\n");
 		UG_LOG("gravForce added: " << gravityForce << "dt: " << m_dt <<"\n");
 		UG_LOG("//////////////////////////////// - log_Gravity - ///////////////////////////////\n");
 		UG_LOG("NACHHER: defect(trans,0): " << DoFRef(vec, transInd[0]) << "\n");
 		UG_LOG("NACHHER: defect(trans,1): " << DoFRef(vec, transInd[1]) << "\n");
        if (dim == 3)
            UG_LOG("NACCHHER: defect(trans,2): " << DoFRef(vec, transInd[2]) << "\n");
 		UG_LOG("NACHHER: defect(rot,0): " << DoFRef(vec, rotInd[0]) << "\n");
 		UG_LOG("NACHHER: defect(rot,1): " << DoFRef(vec, rotInd[1]) << "\n");
        if (dim == 3)
            UG_LOG("NACHHER: defect(rot,2): " << DoFRef(vec, rotInd[2]) << "\n");
    }
    
}


template<typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
add_mass_part_def(vector_type& vec,std::vector<DoFIndex> transInd, std::vector<DoFIndex> rotInd,
		const int levIndex, const int prtIndex)
{
	bool logGravity = false;
	if (logGravity)
    {
		UG_LOG("//////////////////////////////// - log momentum - ///////////////////////////////\n");
 		UG_LOG("VORHER: defect(trans,0): " << DoFRef(vec, transInd[0]) << "\n");
 		UG_LOG("VORHER: defect(trans,1): " << DoFRef(vec, transInd[1]) << "\n");
 		UG_LOG("VORHER: defect(rot,0): " << DoFRef(vec, rotInd[0]) << "\n");
 		UG_LOG("VORHER: defect(rot,1): " << DoFRef(vec, rotInd[1]) << "\n");
    }
    
	number volume;
	if (m_bVolumeCompExact) // compute volume analytically
		volume = Volume(levIndex, prtIndex);
	else // compute volume by adding volumes of elements, covered by the prtIndex-th particle
		volume = compute_volume(levIndex, prtIndex);

	const number mass = m_spInterfaceHandlerLocal->get_density(prtIndex)*volume;
	const number momOfInertia = MomOfInertia(levIndex, prtIndex, volume);
	MathVector < 2 > vScaleMass;
	vScaleMass[0] = 1.0;
	vScaleMass[1] = -1.0;

//	loop all time points and add mass parts to global defect
	for (int i = 1; i >= 0; --i)
    {
		if (logGravity)
			UG_LOG("i = " << i << ":\n");
        
		for (int d = 0; d < dim; ++d) {
			double transSol = m_spCutElementHandler->get_transSol(prtIndex,i)[d];
			double rotSol   = m_spCutElementHandler->get_rotSol(prtIndex, i)[d];
            
			DoFRef(vec, transInd[d]) += vScaleMass[i] * mass * transSol;
			if (dim == 3 || d == 0)
				DoFRef(vec, rotInd[d]) += vScaleMass[i] * momOfInertia * rotSol;
            
            if (logGravity)
            {
                UG_LOG("transSol(" << i << ") = " << transSol << "\n");
                UG_LOG("rotSol(" << i << ") = " << rotSol << "\n");
                UG_LOG("+= vScaleMass[" << i << "]*mass*transSol = " << vScaleMass[i]*mass*transSol << "\n");
                UG_LOG("+= vScaleMass[" << i << "]*momOfInertia*rotSol = " << vScaleMass[i]*momOfInertia*rotSol << "\n");
            }
		}
	}

	if (logGravity)
    {
		UG_LOG("NACHHER: defect(trans,0): " << DoFRef(vec, transInd[0]) << "\n");
 		UG_LOG("NACCHHER: defect(trans,1): " << DoFRef(vec, transInd[1]) << "\n");
 		UG_LOG("NACHHER: defect(rot,0): " << DoFRef(vec, rotInd[0]) << "\n");
 		UG_LOG("NACHHER: defect(rot,1): " << DoFRef(vec, rotInd[1]) << "\n");
    }
    
// call only ONCE!! (during 'modify_GlobalSol(): m_bFlagInertial := true)
    if (m_bVolumeCompExact)
        m_bFlagInertial[prtIndex] = false;

    
    return;

}

template<typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
add_rhs(vector_type& vec, std::vector<DoFIndex> transInd, std::vector<DoFIndex> rotInd,
		const int levIndex, const int prtIndex)
{
    
// call assembling of gravitational force acting on prtIndex-th particle
    if (m_bFlagGravity[prtIndex]){
        set_gravitational_rhs(vec, transInd, rotInd, levIndex, prtIndex);
    }
    
// call assembling of repuslive force acting on prtIndex-th particle
    if (m_bRepulsiveForce && m_bFlagGravity[prtIndex]){
        add_repulsive_force_rhs(vec, transInd, rotInd, levIndex, prtIndex);
    }
    
// call assembling of glowinski repuslive force acting on prtIndex-th particle
    if (m_bGlowRepulsiveForce && m_bFlagGravity[prtIndex]){
        add_glowinski_repulsive_force_rhs(vec, transInd, rotInd, levIndex, prtIndex);
    }
    
// no additional if here because of calculations and a UG_LOG call inside the function that should be printed everytime.
// if clause is inside the function and force is added only if m_bMaxRepulsiveForce
    if (m_bFlagGravity[prtIndex]) {
        add_minimum_correction_force_rhs(vec, transInd, rotInd, levIndex, prtIndex);
    }


// call this method only ONCE!! (during 'modify_GlobalSol(): m_bFlagGravity := true)
    if (m_bVolumeCompExact)
        m_bFlagGravity[prtIndex] = false;

}

// call for 'CUT_BY_2_INTERFACE'
template<typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
add_local_vec_to_global_interface_for2(vector_type& vec, const LocalVector& lvec,
		ConstSmartPtr<DoFDistribution> dd)
{
    
// REMARK: not finally tested!
    UG_THROW("ParticleMapper::add_local_vec_to_global_interface_for2(): not finally tested!!!!\n");

// get data
	const int levIndex = get_Index(dd->grid_level(), dd);
	std::vector < DoFIndex > transInd1, rotInd1, transInd2, rotInd2;
	m_spCutElementHandler->get_global_indices(transInd1, rotInd1, levIndex,0);
	m_spCutElementHandler->get_global_indices(transInd2, rotInd2, levIndex,1);

// A. add local to global
	if (!UsualAss())
		if (m_spInterfaceHandlerLocal->StdFV_assembling())
			add_local_vec_to_global_FT_for2_StdFV(vec, lvec, transInd1, rotInd1,
					transInd2, rotInd2);
		else
			add_local_vec_to_global_FT_for2(vec, lvec, transInd1, rotInd1,
					transInd2, rotInd2);
	else
		UG_THROW("add_local_vec_to_global_interface_for2(): -> not implemented!\n");

}

// call for 'CUT_BY_INTERFACE'
template<typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
add_local_vec_to_global_interface(vector_type& vec, const LocalVector& lvec,
		ConstSmartPtr<DoFDistribution> dd)
{
// get data
	const int levIndex = get_Index(dd->grid_level(), dd);
	const int prtIndex = get_prtIndex();
	std::vector < DoFIndex > transInd, rotInd;
	m_spCutElementHandler->get_global_indices(transInd, rotInd, levIndex,
			prtIndex);

// A. add local to global
    if (!UsualAss()){
		add_local_vec_to_global_FT(vec, lvec, transInd, rotInd);
    }
    else {
        UG_THROW("in ParticleMapper::add_local_vec_to_global_interface(): !UsualAss() not implemented!\n");
    }
    
// B. add mass part to defect
	if (m_bFlagInertial[prtIndex])
		add_mass_part_def(vec, transInd, rotInd, levIndex, prtIndex);

// C. add rhs to defect
	if (m_bFlagGravity[prtIndex])
		add_rhs(vec, transInd, rotInd, levIndex, prtIndex);

}


template<typename TDomain, typename TAlgebra>
size_t ParticleMapper<TDomain, TAlgebra>::
map_for2(size_t dof)
{
    UG_THROW("ParticleMapper::map_for2(): not finally tested!\n");
    
	if (dof == 4)
		return m_spInterfaceHandlerLocal->m_vNOInterfaceID[0];
	else if (dof < 2)
		return m_spInterfaceHandlerLocal->m_vQuadriOrigID[dof];
	else if (dof >= 2)
		return m_spInterfaceHandlerLocal->m_vQuadriOrigID[dof];
	else
		UG_THROW("ParticleMapper::map_for2(): error: dof = "
                 << dof << " should lie between 0 and 4!\n");

}

template<typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
add_local_mat_to_global_FT_for2_StdFV(
		matrix_type& mat, const LocalMatrix& lmat,
		std::vector<DoFIndex> transInd1, std::vector<DoFIndex> rotInd1,
		std::vector<DoFIndex> transInd2, std::vector<DoFIndex> rotInd2)
{
// REMARK: not finally tested!
    UG_THROW("ParticleMapper::add_local_mat_to_global_FT_for2_StdFV(): not finally tested!!!!\n");

	MathMatrix<dim, dim> rotationMatCo, rotationMatIP, rotationMatIP_transposed;
	MathVector<dim> radialVector;

	const LocalIndices& rowInd = lmat.get_row_indices();
	const LocalIndices& colInd = lmat.get_col_indices();

	bool remap = true;

	for (size_t fct1 = 0; fct1 < lmat.num_all_row_fct(); ++fct1)
		for (size_t dof1 = 0; dof1 < lmat.num_all_row_dof(fct1); ++dof1) {
			size_t _dof1 = dof1;
			if (remap)
				_dof1 = m_spInterfaceHandlerLocal->corner_orig(dof1);
			const DoFIndex indexRow = DoFIndex(rowInd.index(fct1, _dof1),
					rowInd.comp(fct1, _dof1));

			for (size_t fct2 = 0; fct2 < lmat.num_all_col_fct(); ++fct2)
				for (size_t dof2 = 0; dof2 < lmat.num_all_col_dof(fct2);
						++dof2) {
					size_t _dof2 = dof2;
					if (remap)
						_dof2 = m_spInterfaceHandlerLocal->corner_orig(dof2);
					const DoFIndex indexCol = DoFIndex(
							colInd.index(fct2, _dof2),
							colInd.comp(fct2, _dof2));


					size_t modus = getMappingModus(fct1, _dof1, fct2, _dof2);

					rotationMatCo = m_spCutElementHandler->get_rotationMat(
							m_spInterfaceHandlerLocal->radial_at_co(dof2));
					rotationMatIP = m_spCutElementHandler->get_rotationMat(
							m_spInterfaceHandlerLocal->radial_at_ip(dof1));

					const int prtIndex1 = getPrtIndex(_dof1);
					const int prtIndex2 = getPrtIndex(_dof2);

					switch (modus) {
					case 0: {
						DoFRef(mat, indexRow, indexCol) += lmat.value(fct1,
								dof1, fct2, dof2);
					}
						break;
					case 1: {
						std::vector < DoFIndex > transInd_to, rotInd_to;
                    // coupling TO particle1
						if (prtIndex2 == 0) {
							transInd_to = transInd1;
							rotInd_to = rotInd1;
                        }
                    // coupling TO particle2
						else {
							transInd_to = transInd2;
							rotInd_to = rotInd2;
						}

						DoFRef(mat, indexRow, transInd_to[fct2]) += lmat.value(
								fct1, dof1, fct2, dof2);

						for (size_t cmp = 0; cmp < dim; ++cmp)
							DoFRef(mat, indexRow, rotInd_to[cmp]) += lmat.value(
									fct1, dof1, fct2, dof2)
									* rotationMatCo[fct2][cmp];

					}
						break;
					case 2: {
						std::vector < DoFIndex > transInd_from, rotInd_from;
                    // coupling FROM particle1
						if (prtIndex1 == 0) {
							transInd_from = transInd1;
							rotInd_from = rotInd1;
						}
                    // coupling FROM particle2
						else {
							transInd_from = transInd2;
							rotInd_from = rotInd2;
						}

						DoFRef(mat, transInd_from[fct1], indexCol) +=
								lmat.value(fct1, dof1, fct2, dof2);

						for (size_t cmp = 0; cmp < dim; ++cmp)
							DoFRef(mat, rotInd_from[cmp], indexCol) +=
									lmat.value(fct1, dof1, fct2, dof2)
											* rotationMatIP[fct1][cmp];

					}
						break;
					case 3: {
						std::vector < DoFIndex > transInd_from, rotInd_from;
						std::vector < DoFIndex > transInd_to, rotInd_to;
                    // coupling TO particle1
						if (prtIndex2 == 0) {
							transInd_to = transInd1;
							rotInd_to = rotInd1;
						}
                    // coupling TO particle2
						else {
							transInd_to = transInd2;
							rotInd_to = rotInd2;
						}

                    // coupling FROM particle1
						if (prtIndex1 == 0) {
							transInd_from = transInd1;
							rotInd_from = rotInd1;
						}
                    // coupling FROM particle2
						else {
							transInd_from = transInd2;
							rotInd_from = rotInd2;
						}

                    // trans-trans
						DoFRef(mat, transInd_from[fct1], transInd_to[fct2]) +=
								lmat.value(fct1, dof1, fct2, dof2);

                    // trans-rot ( w x r = w x r_co = R_co^T*w); (R^T*... <=> [fct1][cmp])
						for (size_t cmp = 0; cmp < dim; ++cmp)
							DoFRef(mat, transInd_from[fct1], rotInd_to[cmp]) +=
									lmat.value(fct1, dof1, fct2, dof2)
											* rotationMatCo[fct2][cmp];

						for (size_t cmp = 0; cmp < dim; ++cmp)
							DoFRef(mat, rotInd_from[cmp], transInd_to[fct2]) +=
									lmat.value(fct1, dof1, fct2, dof2)
											* rotationMatIP[fct1][cmp];

						for (size_t cmp1 = 0; cmp1 < dim; ++cmp1)
							for (size_t cmp2 = 0; cmp2 < dim; ++cmp2)
								DoFRef(mat, rotInd_from[cmp1], rotInd_to[cmp2]) +=
										lmat.value(fct1, dof1, fct2, dof2)
												* rotationMatCo[fct2][cmp2]
												* rotationMatIP[fct1][cmp1];

					}
						break;
					case 4: { 		// to pressure outside: sum = e-20 => skip!
					}
						break;

					default:
						throw(UGError("Error in LocToGlob!"));

					} // end switch-cases
				} // end dof2-loop

		}

}

template<typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
add_local_mat_to_global_FT_for2(
		matrix_type& mat, const LocalMatrix& lmat,
		std::vector<DoFIndex> transInd1, std::vector<DoFIndex> rotInd1,
		std::vector<DoFIndex> transInd2, std::vector<DoFIndex> rotInd2)
{
// REMARK: not finally tested!
    UG_THROW("ParticleMapper::add_local_mat_to_global_FT_for2(): not finally tested!!!!\n");

	MathMatrix<dim, dim> rotationMatCo, rotationMatIP, rotationMatIP_transposed;
	MathVector<dim> radialVector;

	const LocalIndices& rowInd = lmat.get_row_indices();
	const LocalIndices& colInd = lmat.get_col_indices();

	for (size_t fct1 = 0; fct1 < lmat.num_all_row_fct(); ++fct1)
		for (size_t dof1 = 0; dof1 < lmat.num_all_row_dof(fct1); ++dof1) {
			const size_t dofMap1 = map_for2(dof1);
			const DoFIndex indexRow = DoFIndex(rowInd.index(fct1, dofMap1),
					rowInd.comp(fct1, dofMap1));

			for (size_t fct2 = 0; fct2 < lmat.num_all_col_fct(); ++fct2)
				for (size_t dof2 = 0; dof2 < lmat.num_all_col_dof(fct2);
						++dof2) {
					const size_t dofMap2 = map_for2(dof2);
					const DoFIndex indexCol = DoFIndex(
							colInd.index(fct2, dofMap2),
							colInd.comp(fct2, dofMap2));

 					rotationMatCo = m_spCutElementHandler->get_rotationMat(
							m_spInterfaceHandlerLocal->radial_at_co(dof2));
					rotationMatIP = m_spCutElementHandler->get_rotationMat(
							m_spInterfaceHandlerLocal->radial_at_ip(dof1));
 
					size_t modus = getMappingModus_for2(fct1, dof1, fct2, dof2);

					switch (modus) {
					case 0: {
						DoFRef(mat, indexRow, indexCol) += lmat.value(fct1,
								dof1, fct2, dof2);
					}
						break;
					case 1: {
						std::vector < DoFIndex > transInd_to, rotInd_to;
                    // coupling TO particle1
						if (dof2 < 2) {
							transInd_to = transInd1;
							rotInd_to = rotInd1;
						}
                    // coupling TO particle2
						else {
							transInd_to = transInd2;
							rotInd_to = rotInd2;
						}
 
						DoFRef(mat, indexRow, transInd_to[fct2]) += lmat.value(
								fct1, dof1, fct2, dof2);

						for (size_t cmp = 0; cmp < dim; ++cmp)
							DoFRef(mat, indexRow, rotInd_to[cmp]) += lmat.value(
									fct1, dof1, fct2, dof2)
									* rotationMatCo[fct2][cmp];

					}
						break;
					case 2: {
						std::vector < DoFIndex > transInd_from, rotInd_from;
                    // coupling FROM particle1
						if (dof1 < 2) {
							transInd_from = transInd1;
							rotInd_from = rotInd1;
						}
                    // coupling FROM particle2
						else {
							transInd_from = transInd2;
							rotInd_from = rotInd2;
						}
						UG_LOG("transInd_from: " << transInd_from[0] << "\n");
						UG_LOG("rotInd_from: " << rotInd_from[0] << "\n");

						DoFRef(mat, transInd_from[fct1], indexCol) +=
								lmat.value(fct1, dof1, fct2, dof2);

						for (size_t cmp = 0; cmp < dim; ++cmp)
							DoFRef(mat, rotInd_from[cmp], indexCol) +=
									lmat.value(fct1, dof1, fct2, dof2)
											* rotationMatIP[fct1][cmp];

					}
						break;
					case 3: {
						std::vector < DoFIndex > transInd_from, rotInd_from;
						std::vector < DoFIndex > transInd_to, rotInd_to;
                    // coupling TO particle1
						if (dof2 < 2) {
							transInd_to = transInd1;
							rotInd_to = rotInd1;
						}
                    // coupling TO particle2
						else {
							transInd_to = transInd2;
							rotInd_to = rotInd2;
						}
						UG_LOG("transInd_to: " << transInd_to[0] << "\n");
						UG_LOG("rotInd_to: " << rotInd_to[0] << "\n");

                    // coupling FROM particle1
						if (dof1 < 2) {
							transInd_from = transInd1;
							rotInd_from = rotInd1;
						}
                    // coupling FROM particle2
						else {
							transInd_from = transInd2;
							rotInd_from = rotInd2;
						}
						UG_LOG("transInd_from: " << transInd_from[0] << "\n");
						UG_LOG("rotInd_from: " << rotInd_from[0] << "\n");

                    // trans-trans
						DoFRef(mat, transInd_from[fct1], transInd_to[fct2]) +=
								lmat.value(fct1, dof1, fct2, dof2);

                    // trans-rot ( w x r = w x r_co = R_co^T*w); (R^T*... <=> [fct1][cmp])
						for (size_t cmp = 0; cmp < dim; ++cmp)
							DoFRef(mat, transInd_from[fct1], rotInd_to[cmp]) +=
									lmat.value(fct1, dof1, fct2, dof2)
											* rotationMatCo[fct2][cmp];

						for (size_t cmp = 0; cmp < dim; ++cmp)
							DoFRef(mat, rotInd_from[cmp], transInd_to[fct2]) +=
									lmat.value(fct1, dof1, fct2, dof2)
											* rotationMatIP[fct1][cmp];

						for (size_t cmp1 = 0; cmp1 < dim; ++cmp1)
							for (size_t cmp2 = 0; cmp2 < dim; ++cmp2)
								DoFRef(mat, rotInd_from[cmp1], rotInd_to[cmp2]) +=
										lmat.value(fct1, dof1, fct2, dof2)
												* rotationMatCo[fct2][cmp2]
												* rotationMatIP[fct1][cmp1];

					}
						break;
					case 4: { // to pressure outside: sum = e-20 => skip!
					}
						break;

					default:
						throw(UGError("Error in LocToGlob!"));

					} // end switch-cases

				} // end dof2-loop

				////////////////////////////////////////////////////////////////////////////////////////////////
				// finally set dirichlet row for velocity DoFs on interface:
				// 	---> also done during 'set_identity_mat_constraint()', BUT only for OUTSIDE_DOM-elements!!
				////////////////////////////////////////////////////////////////////////////////////////////////

				// dof1 != 4 => DoF lies on an interface of one of the two particles
			if (dof1 != 4) {
				const size_t indexDiagInd = rowInd.index(0, dofMap1);
				if (indexDiagInd != transInd1[0][0]
						&& indexDiagInd != rotInd1[0][0]
						&& indexDiagInd != transInd2[0][0]
						&& indexDiagInd != rotInd2[0][0]) {
					for (size_t d = 0; d < dim; ++d) {
						const DoFIndex indexDiag = DoFIndex(
								rowInd.index(d, dofMap1),
								rowInd.comp(d, dofMap1));
						DoFRef(mat, indexDiag, indexDiag) = 1.0;
					}
				}
			}
		} // end dof1-loop


}

template<typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
add_local_mat_to_global_FT(matrix_type& mat, const LocalMatrix& lmat,
		std::vector<DoFIndex> transInd, std::vector<DoFIndex> rotInd)
{
	bool bOld = true;
	bool remap = true;
	bool transposed = false;
    
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ATTENTION: lmat has the algebra-structure and size of m_locJ (i.e. bigger than standard, in case of QUADIR)
// BUT: the indices-infos are still the original one for the associated global indices:
//          --> lmat.m_pRowIndex, lmat.m_pColIndex
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	LocalMatrix locMatOrig = lmat;
	LocalMatrix locMatTrans = lmat;
	LocalMatrix locMatRot = lmat;
	locMatOrig = 0.0;
	locMatTrans = 0.0;
	locMatRot = 0.0;

// cutElem = true: => re-associate to trans/rot-DoFs
	////////////////////////////////////////////////////////////////////////////////////////////////////
	// for LocToGLob-mapping only the following geo-data are needed:
	// 		geo.vOriginalCornerID, geo.vInterfaceID, geo.m_bUseOriginalCornerRadial,
	// 		geo.m_vRadialAtCo, geo.m_vRadialAtIP
	////////////////////////////////////////////////////////////////////////////////////////////////////

	MathMatrix<dim, dim> rotationMatCo, rotationMatIP;
	MathVector<dim> radialVector;

	const LocalIndices& rowInd = lmat.get_row_indices();
	const LocalIndices& colInd = lmat.get_col_indices();


	const int prtIndex = get_prtIndex();

	if (prtIndex == -1)
		UG_THROW("ParticleMapper::add_local_mat_to_global_FT(): prtIndex invalid!: prtIndex = " << prtIndex << "\n");

	for (size_t fct1 = 0; fct1 < lmat.num_all_row_fct(); ++fct1)
		for (size_t dof1 = 0; dof1 < lmat.num_all_row_dof(fct1); ++dof1) {
			size_t _dof1 = dof1;
			if (remap)
				_dof1 = m_spInterfaceHandlerLocal->corner_orig(dof1);

			const DoFIndex indexRow = DoFIndex(rowInd.index(fct1, _dof1),
					rowInd.comp(fct1, _dof1));

			for (size_t fct2 = 0; fct2 < lmat.num_all_col_fct(); ++fct2)
				for (size_t dof2 = 0; dof2 < lmat.num_all_col_dof(fct2);
						++dof2)
                {

					size_t _dof2 = dof2;
					if (remap)
						_dof2 = m_spInterfaceHandlerLocal->corner_orig(dof2);
					const DoFIndex indexCol = DoFIndex(
							colInd.index(fct2, _dof2),
							colInd.comp(fct2, _dof2));

					size_t modus = getMappingModus(fct1, _dof1, fct2, _dof2);
                   
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                // compute according rotational matrices:
                //  --> the couplings to/from the angular velocity need to be multiplied by
                //      the matrix resutling from the cross product
					if (1) //or: if ( m_spInterfaceHandlerLocal->m_bUseOriginalCornerRadial )
					{
						rotationMatCo =
								m_spCutElementHandler->get_rotationMat(
										m_spInterfaceHandlerLocal->radial_at_co(
												dof2));
						rotationMatIP =
								m_spCutElementHandler->get_rotationMat(
										m_spInterfaceHandlerLocal->radial_at_ip(
												dof1));
					} else {
						rotationMatCo =
								m_spCutElementHandler->get_rotationMat(
										m_spInterfaceHandlerLocal->radial_at_co(
												dof2));
						rotationMatIP =
								m_spCutElementHandler->get_rotationMat(
										m_spInterfaceHandlerLocal->radial_at_ip(
												dof1));
					}

					MathMatrix<dim, dim> rotationMatIP_transposed;
					Transpose(rotationMatIP_transposed, rotationMatIP);

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // for explanations on the mapping modus, see 'loc_to_glob_mapper_particle.h'
					switch (modus) {
					case 0: {
						DoFRef(mat, indexRow, indexCol) += lmat.value(fct1, dof1, fct2, dof2);
					}
						break;
					case 1: {
						DoFRef(mat, indexRow, transInd[fct2]) += lmat.value(fct1, dof1, fct2, dof2);

						for (size_t cmp = 0; cmp < dim; ++cmp) {
							DoFRef(mat, indexRow, rotInd[cmp]) += lmat.value(fct1, dof1, fct2, dof2)
									* rotationMatCo[fct2][cmp];
						}

					}
						break;
					case 2: {
						if (!m_spCutElementHandler->get_DoF_modus_linear(prtIndex))
							DoFRef(mat, transInd[fct1], indexCol) += lmat.value(
                                                        fct1, dof1, fct2, dof2);

                    // assemble entries, EVEN THOUGH it is: \int (r x p*n) = 0
                    // --> fct2 == dim NECESSARY, since in case of irregular meshes, the entries can get inconsistent!
                    // --> this effect was tested and observed for some numerical test simulations
						for (size_t cmp = 0; cmp < dim; ++cmp) {
							if (transposed) {
								if (!m_spCutElementHandler->get_DoF_modus_angular(prtIndex))
									DoFRef(mat, rotInd[cmp], indexCol) +=
											lmat.value(fct1, dof1, fct2, dof2)
													* rotationMatIP[cmp][fct1];
							} else {
								if (1) //or: if ( fct2 != dim )
								{
									if (!m_spCutElementHandler->get_DoF_modus_angular(prtIndex))
										DoFRef(mat, rotInd[cmp], indexCol) +=
												lmat.value(fct1, dof1, fct2,dof2)
														* rotationMatIP[fct1][cmp];
								}
							}

						}

					}
						break;
					case 3: {

						if (!m_spCutElementHandler->get_DoF_modus_linear(prtIndex))
                        {
                        // trans -> trans
							DoFRef(mat, transInd[fct1], transInd[fct2]) +=
									lmat.value(fct1, dof1, fct2, dof2);

                        // trans -> rot: ( w x r = w x r_co = R_co^T*w); (R^T*... <=> [fct1][cmp])
							for (size_t cmp = 0; cmp < dim; ++cmp)
								DoFRef(mat, transInd[fct1], rotInd[cmp]) +=
										lmat.value(fct1, dof1, fct2, dof2)
												* rotationMatCo[fct2][cmp];
						}

						if (!m_spCutElementHandler->get_DoF_modus_angular(prtIndex))
                        {
							if (bOld) {
								if (transposed) {
									for (size_t cmp = 0; cmp < dim; ++cmp)
										DoFRef(mat, rotInd[cmp], transInd[fct2]) +=
												lmat.value(fct1, dof1, fct2,
														dof2)
														* rotationMatIP[cmp][fct1];
								} else {
									for (size_t cmp = 0; cmp < dim; ++cmp)
										DoFRef(mat, rotInd[cmp], transInd[fct2]) +=
												lmat.value(fct1, dof1, fct2,
														dof2)
														* rotationMatIP[fct1][cmp];
								}

                            // rot -> rot: ( r x (w x r) = r_ip x (w x r_co) = R_ip*R_co^T*w = (R_co*R_ip^T)^T*w) ; ) (...^T*w <=> [fct1][cmp])
								if (transposed) {
									for (size_t cmp1 = 0; cmp1 < dim; ++cmp1)
										for (size_t cmp2 = 0; cmp2 < dim;
												++cmp2)
											DoFRef(mat, rotInd[cmp1],
													rotInd[cmp2]) += lmat.value(
													fct1, dof1, fct2, dof2)
													* rotationMatCo[fct2][cmp2]
													* rotationMatIP[cmp1][fct1];
								} else {
									for (size_t cmp1 = 0; cmp1 < dim; ++cmp1)
										for (size_t cmp2 = 0; cmp2 < dim;
												++cmp2)
											DoFRef(mat, rotInd[cmp1],
													rotInd[cmp2]) += lmat.value(
													fct1, dof1, fct2, dof2)
													* rotationMatCo[fct2][cmp2]
													* rotationMatIP[fct1][cmp1];
								}
							} else {
                            // use local couplings computed by class 'ParticleBndCond'
                            // and add those to global matrix
                                DoFRef(mat, rotInd[fct1], transInd[fct2]) +=
										get_rotJ_ind(fct1, dof1, fct2, dof2);
								DoFRef(mat, rotInd[fct1], rotInd[fct2]) +=
										get_rotJ_rot(fct1, dof1, fct2, dof2);
							}
						} // end if (...DoF_modus_angular())

					}
						break;
					case 4: { // to pressure outside: sum = e-20 => skip!
					}
						break;

					default:
						throw(UGError("Error in LocToGlob!"));

					} // end switch-cases

				} // end dof2-loop

		////////////////////////////////////////////////////////////////////////////////////////////////
		// finally set dirichlet row for velocity DoFs on interface:
		// 	---> also done during 'set_identity_mat_constraint()', BUT only for OUTSIDE_DOM-elements!!
		////////////////////////////////////////////////////////////////////////////////////////////////

			bool bSetDirichletRowsForVelocity =
					m_spInterfaceHandlerLocal->remapped_fromInterface(_dof1);
			if (bSetDirichletRowsForVelocity) {
				const size_t indexDiagInd = rowInd.index(0, _dof1);
				if (indexDiagInd != transInd[0][0]
						&& indexDiagInd != rotInd[0][0]) {
					for (size_t d = 0; d < dim; ++d) {
						const DoFIndex indexDiag = DoFIndex(
								rowInd.index(d, _dof1), rowInd.comp(d, _dof1));
						DoFRef(mat, indexDiag, indexDiag) = 1.0;
					}
				} else //if transInd or rotInd lie on cut element and are given by the user, i.e. NO DoFs!
				{
					// set dirichlet row for translation vel
					if (m_spCutElementHandler->get_DoF_modus_linear(
							prtIndex)) {
						for (size_t d = 0; d < dim; ++d) {
							const DoFIndex indexDiag = DoFIndex(
									rowInd.index(d, _dof1),
									rowInd.comp(d, _dof1));
							DoFRef(mat, indexDiag, indexDiag) = 1.0;
						}
					}
					// set dirichlet row for rotation vel
					if (m_spCutElementHandler->get_DoF_modus_angular(
							prtIndex)) {
						for (size_t d = 0; d < dim; ++d) {
							const DoFIndex indexDiag = DoFIndex(
									rowInd.index(d, _dof1),
									rowInd.comp(d, _dof1));
							DoFRef(mat, indexDiag, indexDiag) = 1.0;
						}
					}
				}
			}
		} // end dof1-loop

}

template<typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
assemble_QuadriCorners(vector_type& vec,
		const LocalVector& lvec, std::vector<DoFIndex> transInd1,
		std::vector<DoFIndex> rotInd1, std::vector<DoFIndex> transInd2,
		std::vector<DoFIndex> rotInd2)
{
// REMARK: not finally tested!
    return;
    
    
	if (dim == 3)
		UG_THROW("in 'assemble_QuadriCorners()': not implemented!\n");

	UG_LOG("transInd1: " << transInd1[0] << "transInd2: " << transInd2[0] << "\n");
	UG_LOG("rotInd1: " << rotInd1[0] << "rotInd2: " << rotInd2[0] << "\n");

	bool remap = true;

	for (size_t fct = 0; fct < lvec.num_all_fct(); ++fct) {
		for (size_t dof = 0; dof < lvec.num_all_dof(fct); ++dof) {
			if (lvec.value(fct, dof) != lvec.value(fct, dof))
				UG_THROW("NAN in 'add_local_vec_to_global_FT()'!...\n");

			size_t _dof = dof;
			if (remap)
				_dof = m_spInterfaceHandlerLocal->corner_orig(dof);


			DoFRef(vec, transInd1[fct]) += lvec.value(fct, dof);

			MathVector<dim> RotVec;
			RotVec[0] = -m_spInterfaceHandlerLocal->radial_at_ip(dof)[1];
			RotVec[1] = m_spInterfaceHandlerLocal->radial_at_ip(dof)[0];

			DoFRef(vec, rotInd1[0]) += lvec.value(fct, dof) * RotVec[fct];

		} // end dof-loop
	}

	UG_THROW(" end assemble_QuadriCorners()\n\n");

}

template<typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
add_local_vec_to_global_FT_for2_StdFV(
		vector_type& vec, const LocalVector& lvec,
		std::vector<DoFIndex> transInd1, std::vector<DoFIndex> rotInd1,
		std::vector<DoFIndex> transInd2, std::vector<DoFIndex> rotInd2)
{
// REMARK: not finally tested!
    UG_THROW("ParticleMapper::add_local_vec_to_global_FT_for2_StdFV(): not finally tested!!!!\n");
    
	if (dim == 3)
		UG_THROW("add_local_vec_to_global_FT_for2_StdFV() not implemented for 3d!!!\n");

	bool remap = true;

	const LocalIndices& ind = lvec.get_indices();

	for (size_t fct = 0; fct < lvec.num_all_fct(); ++fct)
		for (size_t dof = 0; dof < lvec.num_all_dof(fct); ++dof) {
			if (lvec.value(fct, dof) != lvec.value(fct, dof))
				UG_THROW("NAN in 'add_local_vec_to_global_FT()'!...\n");
			if (lvec.value(fct, dof) > 1e+10)
				UG_THROW("-----------------------> 1e+10: " << lvec.value(fct,dof) << "\n");

			size_t _dof = dof;
			if (remap)
				_dof = m_spInterfaceHandlerLocal->corner_orig(dof);

			const DoFIndex multi_index = DoFIndex(ind.index(fct, _dof),ind.comp(fct, _dof));
			const int prtIndex = getPrtIndex(_dof);

			bool isPrtNode = m_spInterfaceHandlerLocal->remapped_fromInterface(_dof);

        // usual assembling for fluid-dof and pressure-fct
			if (!isPrtNode || fct == dim) {
				DoFRef(vec, multi_index) += lvec.value(fct, dof);
				UG_LOG("added  for fct = " << fct << ", dof = " << dof << ": " << lvec.value(fct,dof) << "\n");
			} else {
				MathVector<dim> RotVec;
				RotVec[0] = -m_spInterfaceHandlerLocal->radial_at_ip(dof)[1];
				RotVec[1] = m_spInterfaceHandlerLocal->radial_at_ip(dof)[0];

				DoFIndex transInd, rotInd;

				if (prtIndex == 0) {
					transInd = transInd1[fct];
					rotInd = rotInd1[0];
  				} else {
					if (prtIndex != 1)
						UG_THROW("prtIndex = " << prtIndex << "\n");
                    
					transInd = transInd2[fct];
					rotInd = rotInd2[0];
				}
				// finally add connection between the 2 particles:
				DoFRef(vec, transInd) += lvec.value(fct, dof);
				DoFRef(vec, rotInd) += lvec.value(fct, dof) * RotVec[fct];
			}

		}
}

template<typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
add_local_vec_to_global_FT_for2(
		vector_type& vec, const LocalVector& lvec,
		std::vector<DoFIndex> transInd1, std::vector<DoFIndex> rotInd1,
		std::vector<DoFIndex> transInd2, std::vector<DoFIndex> rotInd2)
{
// REMARK: not finally tested!
    UG_THROW("ParticleMapper::add_local_vec_to_global_FT_for2(): not finally tested!!!!\n");

	const LocalIndices& ind = lvec.get_indices();

	for (size_t fct = 0; fct < lvec.num_all_fct(); ++fct)
		for (size_t dof = 0; dof < lvec.num_all_dof(fct); ++dof) {
			if (lvec.value(fct, dof) != lvec.value(fct, dof))
				UG_THROW("NAN in 'add_local_vec_to_global_FT()'!...\n");

			const size_t dofMap = map_for2(dof);
			const DoFIndex multi_index = DoFIndex(ind.index(fct, dofMap),
					ind.comp(fct, dofMap));

        // usual assembling for fluid-dof and pressure-fct
			if (dof == 4 || fct == dim) {
				DoFRef(vec, multi_index) += lvec.value(fct, dof);
			} else {

				MathVector<dim> RotVec;
				RotVec[0] = -m_spInterfaceHandlerLocal->radial_at_ip(dof)[1];
				RotVec[1] = m_spInterfaceHandlerLocal->radial_at_ip(dof)[0];

				DoFIndex transInd, rotInd;

				if (dof < 2) {
					transInd = transInd1[fct];
					rotInd = rotInd1[0];
 				} else {
					transInd = transInd2[fct];
					rotInd = rotInd2[0];
 				}
				// finally add connection between the 2 particles:
				DoFRef(vec, transInd) += lvec.value(fct, dof);
				DoFRef(vec, rotInd) += lvec.value(fct, dof) * RotVec[fct];

   			}

		}
}

template<typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
add_local_vec_to_global_FT(vector_type& vec, const LocalVector& lvec,
		std::vector<DoFIndex> transInd, std::vector<DoFIndex> rotInd)
{
	const LocalIndices& ind = lvec.get_indices();
	bool remap = true;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ATTENTION: lmat has the algebra-structure and size of m_locJ (i.e. bigger than standard, in case of QUADIR)
// BUT: the indices-infos are still the original one for the associated global indices:
//          --> lmat.m_pRowIndex, lmat.m_pColIndex
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	LocalVector locMatOrig = lvec;
	LocalVector locMatTrans = lvec;
	LocalVector locMatRot = lvec;
	locMatOrig = 0.0;
	locMatTrans = 0.0;
	locMatRot = 0.0;

	const int prtIndex = get_prtIndex();
	if (prtIndex == -1)
		UG_THROW("ParticleMapper::add_local_mat_to_global_FT(): prtIndex invalid!: prtIndex = "
                 << prtIndex << "\n");

	for (size_t fct = 0; fct < lvec.num_all_fct(); ++fct) {
		for (size_t dof = 0; dof < lvec.num_all_dof(fct); ++dof) {
            
        // some security checks:
			if (lvec.value(fct, dof) != lvec.value(fct, dof))
				UG_THROW("NAN in 'add_local_vec_to_global_FT()'!...\n");

			if (lvec.value(fct, dof) > 1e+10)
				UG_THROW("-----------------------> 1e+10: " << lvec.value(fct,dof) << "\n");

        // if remapping is done, the solution space turns into the 'flat top' shape, since ALL
        // DoFs on the interface get associated (mapped) onto the original vertex of the mesh
			size_t _dof = dof;
			if (remap)
				_dof = m_spInterfaceHandlerLocal->corner_orig(dof);

			const DoFIndex multi_index = DoFIndex(ind.index(fct, _dof),
					ind.comp(fct, _dof));

			bool isPrtNode = m_spInterfaceHandlerLocal->remapped_fromInterface(_dof);

		// usual assembling for fluid-dof OR pressure-fct (continuity equation)
			if (!isPrtNode || fct == dim) {
				DoFRef(vec, multi_index) += lvec.value(fct, dof);
			}

			if (!isPrtNode || fct == dim) {
				// do nothing...
			} else {

            // check, if the linear velocity is given by the user: => NO need to be assembled
				if (!m_spCutElementHandler->get_DoF_modus_linear(prtIndex))
					DoFRef(vec, transInd[fct]) += lvec.value(fct, dof);
                
            // check, if the angular velocity is given by the user: => NO need to be assembled
				if (!m_spCutElementHandler->get_DoF_modus_angular(prtIndex))
                {
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                // for the entries of the angular velocity, they need to be multiplied by
                // the matrix resutling from the cross product
					MathVector<dim> RotVec;

					if (dim == 2) {
						RotVec[0] = -m_spInterfaceHandlerLocal->radial_at_ip(dof)[1];
						RotVec[1] =  m_spInterfaceHandlerLocal->radial_at_ip(dof)[0];
					}
					MathMatrix<dim, dim> rotationMatCo, rotationMatIP;
					MathMatrix<dim, dim> rotationMatIP_transposed;
					if (dim == 3) {
						rotationMatIP = m_spCutElementHandler->get_rotationMat(
										m_spInterfaceHandlerLocal->radial_at_ip(dof));
						Transpose(rotationMatIP_transposed, rotationMatIP);
					}
					MathVector<dim> Pos;
					Pos[0] = m_spInterfaceHandlerLocal->radial_at_ip(dof)[0];
					Pos[1] = m_spInterfaceHandlerLocal->radial_at_ip(dof)[1];

					if (dim == 2)
						DoFRef(vec, rotInd[0]) += lvec.value(fct, dof) * RotVec[fct];
					if (dim == 3)
						for (size_t cmp = 0; cmp < dim; ++cmp)
							DoFRef(vec, rotInd[cmp]) += lvec.value(fct, dof)
									* rotationMatIP_transposed[cmp][fct];
                    
				} // end if (...DoF_modus_angular() )
			}

		} // end dof-loop
	}

}



} // end namespace NavierStokes
} // end namespace ug

#endif /* PARTICLE_MAPPER_IMPL_H_ */
