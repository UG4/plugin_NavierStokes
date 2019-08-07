/*
 * diffusion_interface_handler_local_impl.h
 *
 *  Created on: 15.01.2015
 *      Author: suze
 */

#ifndef MAPPER_2PF_IMPL_
#define MAPPER_2PF_IMPL_


namespace ug{
namespace NavierStokes{


template<typename TDomain, typename TAlgebra>
void InterfaceMapper2PF<TDomain, TAlgebra>::
set_identity_mat_constraint(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd)
{
 	DoFIndex index;

	for (size_t i = 0; i < 18; ++i)
	{
		index = DoFIndex(85 + i,0);

		BlockRef(mat(index, index), 0, 0) = 1.0;
	}

}



template<typename TDomain, typename TAlgebra>
void InterfaceMapper2PF<TDomain, TAlgebra>::
modify_GlobalSol(vector_type& vecMod, const vector_type& vec, ConstSmartPtr<DoFDistribution> dd)
{

	size_t numDoFs = vecMod.size();
	const size_t numNewDoFs = numDoFs - m_numDoFs;

	UG_LOG("---------------modify_GlobalSol--------------\n");
	UG_LOG(" vecMod.size(): " << numDoFs << "\n");
    UG_LOG(" vec.size(): " << vec.size() << "\n");
    UG_LOG("m_numDoFs: " << m_numDoFs << "\n");
	UG_LOG(" computed numNewDoFs: " << numNewDoFs << "\n");
	UG_LOG(" m_numNewDoFs: " << m_numNewDoFs << "\n");


	DoFIndex index;
    std::vector<double > verticesValue;
    verticesValue.clear();

    // numNewDoFs enthaelt schon fact *fct! --> siehe MovingInterface2PF.init()
    for (size_t i = 0; i < numNewDoFs; ++i)
    {
        index = DoFIndex(m_numDoFs + i, 0);
        double value = DoFRef(vec, index);
        
        verticesValue.push_back(value);
    }
	

	UG_LOG(" computed numNewDoFs: " << numNewDoFs << "\n");
	UG_LOG("length of verticesValue should be " << numNewDoFs << " = " << numNewDoFs << "\n");

	if ( numNewDoFs != verticesValue.size() )
	{
		UG_LOG(" computed numNewDoFs: " << numNewDoFs << "\n");
		UG_THROW("length of verticesValue should be " << numNewDoFs << " = " << numNewDoFs << "\n");
	}

	// call InterfaceHandlerLocal-method:
	// no not doing it! Done locally in add_def: locU_tri and locU_quad:
	write_solution(verticesValue);

}

template<typename TDomain, typename TAlgebra>
void InterfaceMapper2PF<TDomain, TAlgebra>::
add_local_mat_to_global(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd)
{
	bool print = false;

// resize global matrix dynamically:
	const size_t numDoFs = mat.num_rows();

	if ( print ) UG_LOG("*** vorher: vec.size(): " << numDoFs << "\n");

	size_t numAllDoFs = m_numDoFs + m_numNewDoFs;
	if ( m_scaleDoFs )
		numAllDoFs = m_numDoFs + 2*m_numNewDoFs;

	const int diffDoFs = numAllDoFs - numDoFs;

// resize global defect ONCE:
	if ( diffDoFs > 0 )
	{
		mat.resize_and_keep_values(numAllDoFs, numAllDoFs);
		if ( print )
		{
			UG_LOG("*** m_numDoFs: " << m_numDoFs << "\n");
			UG_LOG("*** m_numNewDoFs: " << m_numNewDoFs << "\n");
			UG_LOG("*** numAllDoFs: " << numAllDoFs << "\n");
			UG_LOG("*** nachher: mat.num_rows(): " << mat.num_rows() << "\n");
		}
	}
	else if (  diffDoFs == 0  )
	{   if ( print ) UG_LOG("no resizing!\n");}
	else if ( diffDoFs < 0 )
	{
		if ( print ) UG_LOG("diffDoFs = " << diffDoFs << "\n");
		UG_THROW("error in add_local_mat_to_global: diffDofs < 0\n");
	}


	ElementModus modus = m_spInterfaceHandlerLocal->elementModus();

	switch(modus)
	{
		case OUTSIDE_DOM:
			AddLocalMatrixToGlobal(mat, lmat);
			break;
		case INSIDE_DOM:
			AddLocalMatrixToGlobal(mat, lmat);
			break;
		case CUT_BY_INTERFACE:
			add_local_mat_to_global_interface(mat, lmat, dd);
			break;
 		default:
			throw(UGError("Error in IInterfaceMapper::add_local_mat_to_global()!"));

	}

}


template<typename TDomain, typename TAlgebra>
void InterfaceMapper2PF<TDomain, TAlgebra>::
add_local_vec_to_global(vector_type& vec, const LocalVector& lvec, ConstSmartPtr<DoFDistribution> dd)
{
	bool print = false;

	const size_t numDoFs = vec.size();

	if ( print ) UG_LOG("*** vorher: vec.size(): " << numDoFs << "\n");

	size_t numAllDoFs = m_numDoFs + m_numNewDoFs;
	if ( m_scaleDoFs )
		numAllDoFs = m_numDoFs + 2*m_numNewDoFs;

	const int diffDoFs = numAllDoFs - numDoFs;

// resize global defect ONCE:
	if ( diffDoFs > 0 )
	{
		vec.resize(numAllDoFs, true);
        vec.set(0.0);
/*
        bool bJac = m_spInterfaceHandlerLocal->get_jac_bool();
        if ( !bJac )
        {
           // vec.set(0.0);
            m_spInterfaceHandlerLocal->set_jac_bool(true);
        }
*/
		if ( print )
		{
			UG_LOG("*** m_numDoFs: " << m_numDoFs << "\n");
			UG_LOG("*** m_numNewDoFs: " << m_numNewDoFs << "\n");
			UG_LOG("*** numAllDoFs: " << numAllDoFs << "\n");
			UG_LOG("*** nachher: vec.size(): " << vec.size() << "\n");
		}
	}
	else if (  diffDoFs == 0  )
	{   if ( print ) UG_LOG("no resizing!\n");}
	else if ( diffDoFs < 0 )
	{
		if ( print ) UG_LOG("diffDoFs = " << diffDoFs << "\n");
		UG_THROW("error in add_local_vec_to_global: diffDofs < 0\n");
	}


	ElementModus modus = m_spInterfaceHandlerLocal->elementModus();

	switch(modus)
	{
		case OUTSIDE_DOM:
			AddLocalVector(vec, lvec);
			break;
		case INSIDE_DOM:
			AddLocalVector(vec, lvec);
			break;
		case CUT_BY_INTERFACE:
			add_local_vec_to_global_interface(vec, lvec, dd);
			break;
 		default:
			throw(UGError("Error in IInterfaceMapper::add_local_vec_to_global()!"));

	}
}

//get_real_index(dof):
// if dof NOT on interface: returns orig corner index
// if dof ON interface: returns index of m_vertex-array

template<typename TDomain, typename TAlgebra>
void InterfaceMapper2PF<TDomain, TAlgebra>::
add_local_mat_to_global_interface(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd)
{
	UG_LOG("------------------------ START InterfaceMapper2PF ------------------------\n");

	const LocalIndices& rowInd = lmat.get_row_indices();
	const LocalIndices& colInd = lmat.get_col_indices();

	DoFIndex indexRow, indexCol;

	///////////////////////////////////////////////////////////////
	/// FIRST: add loc to glob for locJ_tri:
	///////////////////////////////////////////////////////////////
	const LocalMatrix& locJ_tri = get_local_jacobian_tri();
 	const bool shift_global_index_tri = m_spInterfaceHandlerLocal->get_index_shift_tri();

 	size_t numAllDoFs = m_numDoFs;

 	if ( shift_global_index_tri )
 		numAllDoFs = m_numDoFs + m_numNewDoFs;

//	UG_LOG("in InterfaceMapper2PF::add_loc_mat(): locJ_tri = " << locJ_tri << "\n");

	for(size_t fct1=0; fct1 < locJ_tri.num_all_row_fct(); ++fct1)
		for (size_t dof1 = 0; dof1 < locJ_tri.num_all_row_dof(fct1); ++dof1)
		{
			const size_t dof1_real = m_spInterfaceHandlerLocal->real_index_tri(dof1);

		// 	if dof1_real is index of m_vertex: compute global index:
			if ( m_spInterfaceHandlerLocal->lies_onInterface_tri(dof1) )
				indexRow = DoFIndex(numAllDoFs + dof1_real, fct1);
			else
				indexRow = DoFIndex(rowInd.index(fct1, dof1_real), rowInd.comp(fct1, dof1_real));

			for(size_t fct2=0; fct2 < locJ_tri.num_all_col_fct(); ++fct2)
				for (size_t dof2 = 0; dof2 < locJ_tri.num_all_col_dof(fct2); ++dof2)
				{
					const size_t dof2_real = m_spInterfaceHandlerLocal->real_index_tri(dof2);

				// if dof1_real is index of m_vertex: compute global index:
					if ( m_spInterfaceHandlerLocal->lies_onInterface_tri(dof2) )
						indexCol = DoFIndex(numAllDoFs + dof2_real, fct2);
					else
						indexCol = DoFIndex(colInd.index(fct2, dof2_real), colInd.comp(fct2, dof2_real));

				// finally add loc to glob:
					DoFRef(mat, indexRow, indexCol) += locJ_tri.value(fct1, dof1, fct2, dof2);

					if ( indexRow[0] == 0 || indexRow[0] == 30 ){
						UG_LOG("------------------------ tri: += " << locJ_tri.value(0, dof1, 0, dof2) << "\n");
						UG_LOG("(indexRow, indexCol) = " << indexRow[0] << "," << indexCol[0] << "\n");
					}
				}
		}

	///////////////////////////////////////////////////////////////
	/// SECOND: add loc to glob for locJ_tri:
	///////////////////////////////////////////////////////////////
	const LocalMatrix& locJ_quad = get_local_jacobian_quad();
  	const bool shift_global_index_quad = m_spInterfaceHandlerLocal->get_index_shift_quad();

  	// reset numAllDoFs!
  	 numAllDoFs = m_numDoFs;

 	if ( shift_global_index_quad )
 		numAllDoFs = m_numDoFs + m_numNewDoFs;

//	UG_LOG("in InterfaceMapper2PF::add_loc_mat(): locJ_quad = " << locJ_quad << "\n");

	for(size_t fct1=0; fct1 < locJ_quad.num_all_row_fct(); ++fct1)
		for (size_t dof1 = 0; dof1 < locJ_quad.num_all_row_dof(fct1); ++dof1)
		{
			size_t dof1_real = m_spInterfaceHandlerLocal->real_index_quad(dof1);

		// if dof1_real is index of m_vertex: compute global index:
			if ( m_spInterfaceHandlerLocal->lies_onInterface_quad(dof1) )
				indexRow = DoFIndex(numAllDoFs + dof1_real, fct1);
			else
				indexRow = DoFIndex(rowInd.index(fct1, dof1_real), rowInd.comp(fct1, dof1_real));

			for(size_t fct2=0; fct2 < locJ_quad.num_all_col_fct(); ++fct2)
				for (size_t dof2 = 0; dof2 < locJ_quad.num_all_col_dof(fct2); ++dof2)
				{
					size_t dof2_real = m_spInterfaceHandlerLocal->real_index_quad(dof2);

				// if dof1_real is index of m_vertex: compute global index:
					if ( m_spInterfaceHandlerLocal->lies_onInterface_quad(dof2) )
						indexCol = DoFIndex(numAllDoFs + dof2_real, fct2);
					else
						indexCol = DoFIndex(colInd.index(fct2, dof2_real), colInd.comp(fct2, dof2_real));

				// finally add loc to glob:
					DoFRef(mat, indexRow, indexCol) += locJ_quad.value(fct1, dof1, fct2, dof2);

					if ( indexRow[0] == 0 || indexRow[0] == 30){
						UG_LOG("------------------------ quad: += " << locJ_quad.value(0, dof1, 0, dof2) << "\n");
						UG_LOG("(indexRow, indexCol) = " << indexRow[0] << "," << indexCol[0] << "\n");
					}
				}
		}

	UG_LOG("------------------------ END InterfaceMapper2PF ------------------------\n");

}

template<typename TDomain, typename TAlgebra>
void InterfaceMapper2PF<TDomain, TAlgebra>::
add_local_vec_to_global_interface(vector_type& vec, const LocalVector& lvec, ConstSmartPtr<DoFDistribution> dd)
{
	DoFIndex index_print;

	const LocalIndices& ind = lvec.get_indices();
	DoFIndex index;

	///////////////////////////////////////////////////////////////
	/// FIRST: add loc to glob for locD_tri:
	///////////////////////////////////////////////////////////////
	const LocalVector& locD_tri = get_local_defect_tri();
	const bool shift_global_index_tri = m_spInterfaceHandlerLocal->get_index_shift_tri();

 	size_t numAllDoFs = m_numDoFs;

 	if ( shift_global_index_tri )
 	{		numAllDoFs = m_numDoFs + m_numNewDoFs;
 			UG_LOG("it is shifted!\n");
 	}
//	UG_LOG("in InterfaceMapper2PF::add_loc_vec(): locD_tri = " << locD_tri << "\n");

 	for(size_t fct=0; fct < locD_tri.num_all_fct(); ++fct)
 		for (size_t dof = 0; dof < locD_tri.num_all_dof(fct); ++dof)
 		{
 			const size_t dof_real = m_spInterfaceHandlerLocal->real_index_tri(dof);

 		// if dof_real is index of m_vertex: compute global index:
 			if ( m_spInterfaceHandlerLocal->lies_onInterface_tri(dof) )
 				index = DoFIndex(numAllDoFs + dof_real,fct);
 			else
 				index = DoFIndex(ind.index(fct, dof_real), ind.comp(fct, dof_real));

 		// finally add loc to glob:
 			DoFRef(vec, index) += locD_tri.value(fct, dof);

 		}


	///////////////////////////////////////////////////////////////
	/// SECOND: add loc to glob for locU_quad:
	///////////////////////////////////////////////////////////////
	const LocalVector& locD_quad = get_local_defect_quad();
	const bool shift_global_index_quad = m_spInterfaceHandlerLocal->get_index_shift_quad();

// reset numAllDoFs!
    numAllDoFs = m_numDoFs;
    
 	if ( shift_global_index_quad )
 		numAllDoFs = m_numDoFs + m_numNewDoFs;

//	UG_LOG("in InterfaceMapper2PF::add_loc_vec(): locD_quad = " << locD_quad << "\n");

 	for(size_t fct=0; fct < locD_quad.num_all_fct(); ++fct)
 		for (size_t dof = 0; dof < locD_quad.num_all_dof(0); ++dof)
 		{
 			size_t dof_real = m_spInterfaceHandlerLocal->real_index_quad(dof);

 		// if dof_real is index of m_vertex: compute global index:
 			if ( m_spInterfaceHandlerLocal->lies_onInterface_quad(dof) )
 				index = DoFIndex(numAllDoFs + dof_real,fct);
 			else
 				index = DoFIndex(ind.index(fct, dof_real), ind.comp(fct, dof_real));

 		// finally add loc to glob:
 			DoFRef(vec, index) += locD_quad.value(fct, dof);
 		}

}

} // namespace NavierStokes
} // end ug namespace

#endif /* MAPPER_2PF_IMPL_ */
