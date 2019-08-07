/*
 * moving_particle_tools.h
 *
 *  Created on: 20.01.2015
 *      Author: suze
 */

#ifndef PARTICLE_MAPPER_TOOLS_H_
#define PARTICLE_MAPPER_TOOLS_H_

namespace ug{
namespace NavierStokes{

/*
template <typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
add_local_mat_to_global_usual(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd)
{
// get level index for 'isFluidDoF()' - level dependent!
	const int Index = m_myParticle->get_Index(dd->grid_level());

	PROFILE_FUNC();


	const LocalIndices& rowInd = lmat.get_row_indices();
	const LocalIndices& colInd = lmat.get_col_indices();


	// index of particle, the element is associated to
	int prtIndex, prtIndex_ref = -1;
	size_t indexList1;    // indexList1 = entry in 'm_vvvIndexRadialPair_List', see 'is_marked()'
							 	 	  // indexList2 = entry in 'm_vvvIndexPosPair_List2', see 'is_borderFluid()'

	bool cutElem = false;
	cutElem = elemIsCut(rowInd, lmat.num_all_row_dof(0), Index, &prtIndex, dd);

	for(size_t fct1=0; fct1 < lmat.num_all_row_fct(); ++fct1)
	{
   		// scv-loop:
		for(size_t dof1=0; dof1 < lmat.num_all_row_dof(fct1); ++dof1)
		{

			const size_t rowIndex = rowInd.index(fct1,dof1);
			const size_t rowComp = rowInd.comp(fct1,dof1);

			DoFIndex multi_index_row = DoFIndex(rowIndex, rowComp);

			// if the corresponding DoF is 'inside the fluid' OR 'pressure-DoF':
			///////////////////////////////////////////////////////////
			// 				REMARK:
			// if m_bExtrapolatePrs = true:
			// is_marked() = true also for pressure-DoFs inside the particle!!
			// 	-> see: particle_tools.h: if ( m_dExtrapolatePrs ) ...
			///////////////////////////////////////////////////////////
			if ( !is_marked(multi_index_row, Index, &prtIndex, &indexList1) )
			{

				for(size_t fct2=0; fct2 < lmat.num_all_col_fct(); ++fct2)
					for(size_t dof2=0; dof2 < lmat.num_all_col_dof(fct2); ++dof2)
					{
						const size_t colIndex = colInd.index(fct2,dof2);
						const size_t colComp = colInd.comp(fct2,dof2);

						DoFIndex multi_index_col = DoFIndex(colIndex, colComp);

						// if the target DoF is 'inside the fluid' OR 'pressure-DoF': do usual assembling...
						///////////////////////////////////////////////////////////
						// 				REMARK:
						// if m_bExtrapolatePrs = true:
						// is_marked() = true also for pressure-DoFs inside the particle!!
						// 	-> see: particle_tools.h: if ( m_dExtrapolatePrs ) ...
						///////////////////////////////////////////////////////////
						if ( !is_marked(multi_index_col, Index, &prtIndex, &indexList1) )
						{
  							DoFRef(mat, multi_index_row, multi_index_col)
										+= lmat.value(fct1,dof1,fct2,dof2);
						}
						// if the target DoF is a velocity DoF in the particle: RECONNECT to Trans/Rot
						else
						{
 						// CHECK for elements, cut by two different particles:
							if ( prtIndex_ref == -1 )
								prtIndex_ref = prtIndex;

 						// get indices for re-Association:
							DoFIndex transInd = m_myParticle->get_transInd_Comp(Index, prtIndex, fct2);
 							if ( !m_myParticle->m_bExtrapolate && (int)fct2 != dim )
							{
 							// CONNECT to Trans:
								DoFRef(mat, multi_index_row, transInd)
											+= lmat.value(fct1,dof1,fct2,dof2);

								//UG_THROW("lmat.value(fct1,dof1,fct2,dof2) =" << lmat.value(fct1,dof1,fct2,dof2) << "\n");

							// CONNECT to Rot:
								if ( 1 ) {//fct1 != dim ){
								MathMatrix<dim,dim> rotationMat;
								get_rotationMat(multi_index_col, Index, prtIndex, rotationMat);
								std::vector<DoFIndex> rotInd_   = m_myParticle->get_rotInd(Index, prtIndex);
								for ( int cmp = 0; cmp < dim; ++cmp )
   									DoFRef(mat, multi_index_row, rotInd_[cmp])
											+= lmat.value(fct1,dof1,fct2,dof2)*rotationMat[fct2][cmp];
								}
							}
							else if (  (!m_myParticle->m_bExtrapolatePrs && (int)fct2 == dim) ) //( m_myParticle->m_ExtrapolationPrsMethod == 6 && (int)fct2 == dim) ||
							{
 								// do usual assembling for pressure DoF inside particle:
	 							DoFRef(mat, multi_index_row, multi_index_col)
											+= lmat.value(fct1,dof1,fct2,dof2);
							}
							else if ( m_myParticle->m_bExtrapolatePrs && (int)fct1 == dim && (int)fct2 != dim )
							{
								UG_THROW("Haehhh...dieser fall sollte schon abgedeckt sein?...:(...\n");
								if ( (int)fct2 != dim )
								{
									// fct2 != dim => die standart coeffs fuer u\cdot n uebernehmen:
								// CONNECT to Trans:
									DoFRef(mat, multi_index_row, transInd)
												+= lmat.value(fct1,dof1,fct2,dof2);

								// CONNECT to Rot:
									MathMatrix<dim,dim> rotationMat;
									get_rotationMat(multi_index_col, Index, prtIndex, rotationMat);
									std::vector<DoFIndex> rotInd_   = m_myParticle->get_rotInd(Index, prtIndex);
									for ( int cmp = 0; cmp < dim; ++cmp )
										DoFRef(mat, multi_index_row, rotInd_[cmp])
												+= lmat.value(fct1,dof1,fct2,dof2)*rotationMat[fct2][cmp];
								}

							}
							else //if ( (int) fct1 != dim && (int) fct2 == dim ) // druck beitraege aus ContEq (also fuer stab) werden einfach weggelassen!
							{

								//////////////////////////////////////////////////////////////////////////
								//				 EXTRAPOLATION-AKTION
								//
								//	'isBorder' and 'coefficients' are computed depending on the
								// chosen 'm_ExtrapolationMethod'
								//////////////////////////////////////////////////////////////////////////
					//			for ( size_t f = 0; f < lmat.num_all_row_dof(0); ++f )
					//				UG_LOG("ind(" << f << ") = " << rowInd.index(0,f) << "\n");

								std::vector<size_t> ind;
								for ( size_t f = 0; f < lmat.num_all_row_dof(0); ++f )
									ind.push_back(rowInd.index(0,f));

								bool isBorder = false; // -> only substitute velocity-DoFs
								vVertexExtrapolData coefficients;
								size_t coeffComp;

								//////////////////////////////////////////////////////////////////////////
								//	'isBorder' and 'coefficients'
								//////////////////////////////////////////////////////////////////////////
								if ( (int)fct2 != dim )
								{
									if ( !m_myParticle->m_bExtrapolate )
									{UG_THROW("Mince!....case wrong handeled!....:  IFF extrapolation of VELOCITY is NOT set, this case may not occur!\n");}

									coeffComp = colComp;
									isBorder = m_myParticle->m_myExtrapolation->getCoefficients_byDoFIndex(&coefficients, multi_index_col, ind, Index, prtIndex);
								}
								else
								{
									if ( !m_myParticle->m_bExtrapolatePrs )
									{UG_THROW("Mince!....case wrong handeled!....:  IFF extrapolation of PRESSURE is NOT set, this case may not occur!\n");}

									// re-set 'colComp', since in coefficients[] for pressure-extrapolation
									// only 1 connection is stored!
									coeffComp = 0;
									isBorder = m_myParticle->m_myExtrapolationPrs->getCoefficients_byDoFIndex(&coefficients, multi_index_col, ind, Index, prtIndex);

								}

								//////////////////////////////////////////////////////////////////////////
								// final assembling with coefficient-weighting
								// identically for BOTH extrapolationMethods!!
								//////////////////////////////////////////////////////////////////////////
								if ( isBorder )
								{

								/// wird fuer Aufruf von 'projectDefect()' in 'particle_constraint_impl.h':
	//								std::pair<DoFIndex, DoFIndex> indexPair = std::make_pair(multi_index_row, multi_index_col);

	//								if ( m_myParticle->m_bExtrapolatePrs && (int)fct2 == dim )
	//									m_myParticle->m_myExtrapolationPrs->set_locJData(lmat.value(fct1,dof1,fct2,dof2), indexPair, coefficients, Index, prtIndex);

									if ( multi_index_row[0] == 1783 ){
									UG_LOG("_____________ multi_index_row: " << multi_index_row << "\n");

									UG_LOG("_____________ coeff.size(): " << coefficients[0].size() << "\n");
									UG_LOG("_____________ lmat.value(): " << lmat.value(fct1,dof1,fct2,dof2) << "\n");
									}
									// write new entries coming from extrapolating equation for outVrt:
									for ( size_t i = 0; i < coefficients[0].size(); ++i )
									{

										const size_t writeIndex = coefficients[coeffComp][i].second[0];
										const size_t writeComp  = coefficients[coeffComp][i].second[1];
										const number lambda = coefficients[coeffComp][i].first;
										DoFIndex writeInd = DoFIndex(writeIndex, writeComp);


										//UG_LOG("ind_row: " << multi_index_row << "Ind_col: " << writeInd << "\n");
										//UG_LOG("lambda: " << lambda << "\n");

										DoFRef(mat, multi_index_row, writeInd)
											+= lmat.value(fct1,dof1,fct2,dof2)*lambda;
									}
									//UG_LOG("_____\n");

								} // end 'm_myParticle->m_bExtrapolate '


							} // else ( m_bExtrapolate = true )
						}
					} // end dof2_loop


			} // end !is_marked()


		} // end dof1-loop
	} // end fct1-loop

}

template <typename TDomain, typename TAlgebra>
void ParticleMapper<TDomain, TAlgebra>::
add_local_vec_to_global_usual(vector_type& vec, const LocalVector& lvec, ConstSmartPtr<DoFDistribution> dd)
{

// get level index for 'isFluidDoF()' - level dependent!
	const int Index = m_myParticle->get_Index(dd->grid_level());

	const LocalIndices& ind = lvec.get_indices();

	// index of particle, the element is associated to
	int prtIndex, prtIndex_ref = -1;
	size_t indexList1 = 1013;

	bool cutElem = false;
	if (1 ) // m_myParticle->m_bExtrapolatePrs )
	{
//		bool skipAssemblingForOutsideElem;
//		cutElem = elemIsCut(ind, lvec.num_all_dof(0), skipAssemblingForOutsideElem, Index, &prtIndex);
		cutElem = elemIsCut(ind, lvec.num_all_dof(0), Index, &prtIndex, dd);

		if ( 0 ) //!cutElem )
			if ( skipAssemblingForOutsideElem )
				return;
	}



 	for(size_t fct=0; fct < lvec.num_all_fct(); ++fct)
	{
		// skip assembling of ContEy on cutElem in order to discretize the pure
		// mass flux usind u \dot n via 'add_PureMassFlux' in 'particle_constraint_imol.h'
		if ( m_myParticle->m_bExtrapolatePrs )
			if ( m_myParticle->m_myExtrapolationPrs->computePureMassFlux() && cutElem && (int)fct == dim )
				continue;

		for(size_t dof=0; dof < lvec.num_all_dof(fct); ++dof)
		{
			const size_t VecIndex = ind.index(fct,dof);
			const size_t VecComp = ind.comp(fct,dof);

			DoFIndex multi_index = DoFIndex(VecIndex, VecComp);

			// if the corresponding DoF is 'inside the fluid' OR 'pressure-DoF':
			if ( !is_marked(multi_index, Index, &prtIndex, &indexList1) )
			{

			// throw error, if defect has value 'nan':
				if ( lvec.value(fct,dof) != lvec.value(fct,dof))
				{
					UG_LOG("multiindex: " << multi_index << "\t" << indexList1 << "\n");
					UG_THROW("ParticleMap:add_local_vec_to_global: Matrix entry with error: lvec.value(fct,dof) = " << lvec.value(fct,dof) << "\n");
				}

				DoFRef(vec, multi_index) += lvec.value(fct,dof);

			}
			else if ( !m_myParticle->m_bDefectByHand )
			{
				if ( dim == 3 )
					UG_THROW("in 'add_local_vec_to_global': NOT implemented for 3d!...EXIT!..\n");

			// CHECK for elements, cut by two different particles:
				if ( prtIndex_ref == -1 )
					prtIndex_ref = prtIndex;

			// get ExtraIndices - for WRITING DoFRef(vec,ind):
				DoFIndex transInd = m_myParticle->get_transInd_Comp(Index, prtIndex, fct);
				DoFIndex rotInd   = m_myParticle->get_rotInd_Comp(Index, prtIndex, 0);

				MathVector<dim>  RotVec;
				size_t cmp;


			// FIRST: get radialVector before adding to global vector
				DoFIndex Ind = m_myParticle->m_vvvIndexRadialPair_List[Index][prtIndex][indexList1].first;
				MathVector<dim>  radialVector = m_myParticle->m_vvvIndexRadialPair_List[Index][prtIndex][indexList1].second;
				RotVec[0] = -radialVector[1];
				RotVec[1] = radialVector[0];
				cmp = Ind[1];


			// SECOND: adding to global vector
			// ADD to Trans:
				DoFRef(vec, transInd) += lvec.value(fct,dof);
			// ADD to Rot:
				DoFRef(vec, rotInd) += lvec.value(fct,dof)*RotVec[cmp];

			} // end else


		} // end dof-loop
	} // end fct-loop


}
*/

} // end namespace NavierStokes
} // end namespace ug



#endif /* PARTICLE_MAPPER_TOOLS_H_ */
