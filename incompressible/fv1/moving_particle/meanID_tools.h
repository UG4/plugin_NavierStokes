/*
 * moving_particle_tools.h
 *
 *  Created on: 20.01.2015
 *      Author: suze
 */

#ifndef MEAN_ID_H_
#define MEAN_ID_H_

namespace ug{
namespace NavierStokes{

/*
     
    template <typename TDomain, typename TAlgebra>
    number MovingParticle<TDomain, TAlgebra>::
    compute_functional_fixed(const size_t n, vector_type& u, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace,
                       SmartPtr<NavierStokesFV1<TDomain> > spMaster, const int topLevel)
    {
        ConstSmartPtr<DoFDistribution> dd = spApproxSpace->dof_distribution(GridLevel(topLevel, GridLevel::LEVEL));
        
        //	get data
        typename TDomain::position_accessor_type aaPos = this->m_spApproxSpace->domain()->position_accessor();
        
        //	create Function Group
        FunctionGroup fctGrp(dd->function_pattern(), TokenizeString(m_spParticleHandlerGlobal->m_fctNames));
        
        // get subset handler
        ConstSmartPtr<ISubsetHandler> rSH = dd->subset_handler();
        
        //	get iterators for all elems on subset
        typedef typename DoFDistribution::dim_traits<dim>::grid_base_object grid_base_object;
        
        typename DoFDistribution::traits<grid_base_object>::const_iterator iter, iterEnd;
        iter = dd->begin<grid_base_object>();
        iterEnd = dd->end<grid_base_object>();
        
        size_t numDoFs = 0.5 * (n+1) * (n+2); // = 45; 36 :)
        
        
        number D_kVol = 0.0;
        number D_kSurf = 0.0;
        std::vector<number> gradD(2*numDoFs, 0.0);
        
        size_t counterFixed = 0;
        //	loop elements in order to compute 'U_global' and 'omega_global':
        for( ; iter != iterEnd; ++iter)
        {
            //	get element
            grid_base_object* elem = *iter;
            
            std::vector<Vertex*> vVertex;
            CollectVertices(vVertex, *m_spParticleHandlerGlobal->m_spMG, elem);
            
            //	loop vertices
            for (size_t v = 0; v < vVertex.size(); ++v)
            {
                UG_LOG("vrt_" << v << ": " << aaPos[vVertex[v]][0] << "\t" <<  aaPos[vVertex[v]][1] << "\n");
            }
            
            //////////////////////////////////////////////////
            // compute data:
            std::vector<MathVector<dim> > vCornerCoords;
            CollectCornerCoordinates(vCornerCoords, *elem, aaPos, true);
            
            number areaTria = ElementSize<dim>(ROID_TRIANGLE, &vCornerCoords[0]);
            
            //////////////////////////////////////////////////
            // loop edges and:
            //  -> compute 'baseLineSquared', 'baseLine', 'dist'
            //  -> collect associated vertex 'vrtOut', 'vIndex_Out'
            
            std::vector<Edge*> vEdges;
            CollectEdgesSorted(vEdges, *m_spParticleHandlerGlobal->m_spMG, elem);
            
            // values to be computed:
            number baseLineSquared, baseLine;
            MathVector<dim> dist;
            size_t vIndex_Out;
            Vertex* vrtOut;
            
            for(size_t e = 0; e < vEdges.size(); ++e)
            {
                Edge* edge = vEdges[e];
                
                if ( rSH->get_subset_index(edge) != 3 )
                {
                    UG_LOG("...continue...subset = " << rSH->get_subset_index(edge) << "\n");
                    continue;
                }
                
                UG_LOG("subset = " << rSH->get_subset_index(edge) << "\n");
                
                std::vector<Vertex*> vVertexEdge;
                CollectVertices(vVertexEdge, *m_spParticleHandlerGlobal->m_spMG, edge);
                if ( vVertexEdge.size() != 2 )
                    UG_THROW("error in collecting vertices associated to an edge!....EXIT!...\n");
                
                Vertex* vrt1 = vVertexEdge[0];
                Vertex* vrt2 = vVertexEdge[1];
                
                baseLineSquared = VecDistanceSq(aaPos[vrt1], aaPos[vrt2]);
                baseLine = sqrt(baseLineSquared);
                dist[0] = fabs(aaPos[vrt1][0] - aaPos[vrt2][0]);
                dist[1] = fabs(aaPos[vrt1][1] - aaPos[vrt2][1]);
                
                //	loop vertices and get vrtOut:
                for (size_t v = 0; v < vVertex.size(); ++v)
                {
                    if ( (vVertex[v] != vrt1) && (vVertex[v] != vrt2) )
                    {
                        vrtOut = vVertex[v];
                        vIndex_Out = v;
                    }
                }
                
                UG_LOG("vrt1: " << aaPos[vrt1][0] << "\t" <<  aaPos[vrt1][1] << "\n");
                UG_LOG("vrt2: " << aaPos[vrt2][0] << "\t" <<  aaPos[vrt2][1] << "\n");
                UG_LOG("vrtOut: " << aaPos[vrtOut][0] << "\t" <<  aaPos[vrtOut][1] << "\n");
                
            } // end edge-loop
            UG_LOG("new edge\n");
            
            ////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////
            // (1) functional computations:
            ////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////
            
            //////////////////////////////////////////////////
            // volume computations:
            if ( rSH->get_subset_index(elem) == 0 )
                
            {
                D_kVol += 0.5 * baseLineSquared/areaTria;
            }
            //////////////////////////////////////////////////
            // surface computations:
            else if ( rSH->get_subset_index(elem) == 4 )
            {
                D_kSurf += 0.5 * baseLineSquared/areaTria;
            }
            else {
                UG_LOG("elem subset index = " << rSH->get_subset_index(elem) << "\n");
                
                UG_THROW("why more than 2 subsets possible?\n");
            }
            
            
            ////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////
            // (2) gradient computations:
            ////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////
            
            
            //////////////////////////////////////////////////
            // (1) get global indices of vrt1, vrt2, vrtOut for final adding up to DoF-vector in (4):
            
            std::vector<DoFIndex> Ind;
            dd->dof_indices(elem, 1, Ind);  // #instead-ToDo: dd->dof_indices(elem, cmp, Ind);
            
            for (size_t v = 0; v < vVertex.size(); ++v)
                UG_LOG("2 Ind = " << Ind[v] << "\n");
            
            
            
            ////////////////////////////////////////////////////////////////////////////////
            // #instead-ToDo: loop components 'cmp' = x- and y-direction:
            for (size_t cmp = 0; cmp < dim; ++cmp)
            {
                for (size_t v = 0; v < vVertex.size(); ++v)
                {

                    std::vector<DoFIndex>  vInd;
                    if(dd->inner_dof_indices(vVertex[v], cmp, vInd) != 1)
                        UG_THROW("1: error in inner_dof_indices operatio!\n");
                    
                    size_t Index = Ind[v][0];
                    
                    
                    UG_LOG("Index = " << Index << "\n");
                    UG_LOG("vInd_0[0] = " << vInd[0] << "\n");
                    
                    
        // boundary nodes are fixed and therefore no DoFs:
        if ( rSH->get_subset_index(vVertex[v]) == 1 || rSH->get_subset_index(vVertex[v]) == 2 || rSH->get_subset_index(vVertex[v]) == 3 || rSH->get_subset_index(vVertex[v]) == 5 )
        {
            DoFRef(u, vInd[0]) = 0.0;
            
            ++counterFixed;
            continue;
        }
                    
                    //////////////////////////////////////////////////
                    // (2) data to be computed for the gradient computations
                    number gradD_elem;
                    number gradA_elem;
                    
                    // get indices modulo 2
                    size_t index_1 = (v+1)%3;
                    size_t index_2 = (v+2)%3;
                    
                    gradA_elem = fabs(aaPos[vVertex[index_1]][(cmp+1)%2] - aaPos[vVertex[index_2]][(cmp+1)%2]);
                    
                    
                    // ToDo: Vorzeichen für Ableitung von Area A(x1, x2, x3) checken: wirkliche fabs()?
                    
                    if ( (vVertex[v] != vrtOut) )
                    {
                        if ( v == vIndex_Out ) UG_THROW("1: error in vIndex_Out computation!\n");
                        
                        gradD_elem = 2*dist[cmp];
                    }
                    else
                    {
                        if ( v != vIndex_Out ) UG_THROW("2: error in vIndex_Out computation!\n");
                        
                        gradD_elem = 0.0;
                    }
                    
                    
                    
                    //////////////////////////////////////////////////
                    // (3) gradient computations:
                    
                    
                    //////////////////////////////////////////////////
                    // (3.1) volume computations:
                    if ( rSH->get_subset_index(elem) == 0 )
                        
                    {
                        UG_LOG("elem subset index = 0:" << rSH->get_subset_index(elem) << "\n");
                        
                        number gradD_kVol = (areaTria*gradD_elem - baseLineSquared*gradA_elem)/(2*areaTria*areaTria);
 
                        //////////////////////////////////////////////////
                        // (4) add up all computations to entry of DoF-vector:
                        //      gradD[Index]           +=  gradD_kVol_x/(n*n); // #instead-ToDo:: DoFRef(u, Ind[v]) += ...
                        //     gradD[numDoFs + Index] +=  gradD_kVol_y/(n*n); // ToDo!
                        DoFRef(u, vInd[0]) += (-1.0) * gradD_kVol/(n*n);
                    }
                    
                    //////////////////////////////////////////////////
                    // (3.2) surface computations:
                    else if ( rSH->get_subset_index(elem) == 4 )
                    {
                        UG_LOG("elem subset index = 4:" << rSH->get_subset_index(elem) << "\n");
                        
                        number gradD_kSurf = (areaTria*gradD_elem - baseLineSquared*gradA_elem)/(2*areaTria*areaTria);
                        
                         //////////////////////////////////////////////////////////
                        // (4) add up all computations to entry of DoF-vector:
                        //     gradD[Index]           -= gradD_kSurf_x*(n-1)/(n*n); // #instead-ToDo:: DoFRef(u, Ind[v]) += ...
                        //    gradD[numDoFs + Index] -= gradD_kSurf_y*(n-1)/(n*n); // ToDo!
                        DoFRef(u, vInd[0]) -= (-1.0) * gradD_kSurf*(n-1)/(n*n);
                    }
                    else
                    {
                        UG_LOG("elem subset index = " << rSH->get_subset_index(elem) << "\n");
                        UG_THROW("--> why more than 2 subsets possible?\n");
                    }
                    
                    
                    
                    UG_LOG("\n ---> v = " << v << " and Ind = " << Ind[v] << "\n");
                    
                }// end vrt-loop
                
            }// end cmp-loop
            
        } // end elem-loop
        
        number functional = D_kVol/(n*n) - D_kSurf*(n-1)/(n*n);
        
        UG_LOG("counterFixed = " << counterFixed << "\n");
        
        return functional;
        
    }
*/
/*
    template <typename TDomain, typename TAlgebra>
    number MovingParticle<TDomain, TAlgebra>::
    compute_functional_combined(const size_t n, vector_type& u, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace,
                       SmartPtr<NavierStokesFV1<TDomain> > spMaster, const int topLevel)
    {
        ConstSmartPtr<DoFDistribution> dd = spApproxSpace->dof_distribution(GridLevel(topLevel, GridLevel::LEVEL));
        
        //	get data
        typename TDomain::position_accessor_type aaPos = this->m_spApproxSpace->domain()->position_accessor();
        
        //	create Function Group
        FunctionGroup fctGrp(dd->function_pattern(), TokenizeString(m_spParticleHandlerGlobal->m_fctNames));
        
        // get subset handler
        ConstSmartPtr<ISubsetHandler> rSH = dd->subset_handler();
        
        //	get iterators for all elems on subset
        typedef typename DoFDistribution::dim_traits<dim>::grid_base_object grid_base_object;
        
        typename DoFDistribution::traits<grid_base_object>::const_iterator iter, iterEnd;
        iter = dd->begin<grid_base_object>();
        iterEnd = dd->end<grid_base_object>();
        
        size_t numDoFs = 0.5 * (n+1) * (n+2); // = 45; 36 :)
        
        
        number D_kVol = 0.0;
        number D_kSurf = 0.0;
        std::vector<number> gradD(2*numDoFs, 0.0);
        
        //	loop elements in order to compute 'U_global' and 'omega_global':
        for( ; iter != iterEnd; ++iter)
        {
            //	get element
            grid_base_object* elem = *iter;
            
            std::vector<Vertex*> vVertex;
            CollectVertices(vVertex, *m_spParticleHandlerGlobal->m_spMG, elem);
            
            //	loop vertices
            for (size_t v = 0; v < vVertex.size(); ++v)
            {
                UG_LOG("vrt_" << v << ": " << aaPos[vVertex[v]][0] << "\t" <<  aaPos[vVertex[v]][1] << "\n");
            }
            
            //////////////////////////////////////////////////
            // compute data:
            std::vector<MathVector<dim> > vCornerCoords;
            CollectCornerCoordinates(vCornerCoords, *elem, aaPos, true);
            
            number areaTria = ElementSize<dim>(ROID_TRIANGLE, &vCornerCoords[0]);
            
            //////////////////////////////////////////////////
            // loop edges and:
            //  -> compute 'baseLineSquared', 'baseLine', 'dist'
            //  -> collect associated vertex 'vrtOut', 'vIndex_Out'
            
            std::vector<Edge*> vEdges;
            CollectEdgesSorted(vEdges, *m_spParticleHandlerGlobal->m_spMG, elem);
            
            // values to be computed:
            number baseLineSquared, baseLine;
            MathVector<dim> dist;
            size_t vIndex_Out;
            Vertex* vrtOut;
            
            for(size_t e = 0; e < vEdges.size(); ++e)
            {
                Edge* edge = vEdges[e];
                
                if ( rSH->get_subset_index(edge) != 3 )
                {
                    UG_LOG("...continue...subset = " << rSH->get_subset_index(edge) << "\n");
                    continue;
                }
                
                UG_LOG("subset = " << rSH->get_subset_index(edge) << "\n");
                
                std::vector<Vertex*> vVertexEdge;
                CollectVertices(vVertexEdge, *m_spParticleHandlerGlobal->m_spMG, edge);
                if ( vVertexEdge.size() != 2 )
                    UG_THROW("error in collecting vertices associated to an edge!....EXIT!...\n");
                
                Vertex* vrt1 = vVertexEdge[0];
                Vertex* vrt2 = vVertexEdge[1];
                
                baseLineSquared = VecDistanceSq(aaPos[vrt1], aaPos[vrt2]);
                baseLine = sqrt(baseLineSquared);
                dist[0] = fabs(aaPos[vrt1][0] - aaPos[vrt2][0]);
                dist[1] = fabs(aaPos[vrt1][1] - aaPos[vrt2][1]);
                
                //	loop vertices and get vrtOut:
                for (size_t v = 0; v < vVertex.size(); ++v)
                {
                    if ( (vVertex[v] != vrt1) && (vVertex[v] != vrt2) )
                    {
                        vrtOut = vVertex[v];
                        vIndex_Out = v;
                    }
                }
                
                UG_LOG("vrt1: " << aaPos[vrt1][0] << "\t" <<  aaPos[vrt1][1] << "\n");
                UG_LOG("vrt2: " << aaPos[vrt2][0] << "\t" <<  aaPos[vrt2][1] << "\n");
                UG_LOG("vrtOut: " << aaPos[vrtOut][0] << "\t" <<  aaPos[vrtOut][1] << "\n");
                
            } // end edge-loop
            UG_LOG("new edge\n");
            
            ////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////
            // (1) functional computations:
            ////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////
            
            //////////////////////////////////////////////////
            // volume computations:
            if ( rSH->get_subset_index(elem) == 0 )
                
            {
                D_kVol += 0.5 * baseLineSquared/areaTria;
            }
            //////////////////////////////////////////////////
            // surface computations:
            else if ( rSH->get_subset_index(elem) == 4 )
            {
                D_kSurf += 0.5 * baseLineSquared/areaTria;
            }
            else {
                UG_LOG("elem subset index = " << rSH->get_subset_index(elem) << "\n");
                
                UG_THROW("why more than 2 subsets possible?\n");
            }
            
            
            ////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////
            // (2) gradient computations:
            ////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////
            
            
            //////////////////////////////////////////////////
            // (1) get global indices of vrt1, vrt2, vrtOut for final adding up to DoF-vector in (4):
            
            std::vector<DoFIndex> Ind;
            dd->dof_indices(elem, 1, Ind);  // #instead-ToDo: dd->dof_indices(elem, cmp, Ind);
            
            for (size_t v = 0; v < vVertex.size(); ++v)
                UG_LOG("2 Ind = " << Ind[v] << "\n");
            
            
            number _A_ = (aaPos[vVertex[1]][1]-aaPos[vVertex[0]][1])*(aaPos[vVertex[2]][0]-aaPos[vVertex[0]][0])
            -(aaPos[vVertex[1]][0]-aaPos[vVertex[0]][0])*(aaPos[vVertex[2]][1]-aaPos[vVertex[0]][1]);
            
            number gradA_scale = 0.5 * _A_ / fabs(_A_);
            
            //        if ( gradA_scale > 0 )
            //           UG_LOG("gradA_scale = " << gradA_scale << "\n");
            
            ////////////////////////////////////////////////////////////////////////////////
            // #instead-ToDo: loop components 'cmp' = x- and y-direction:
            for (size_t cmp = 0; cmp < dim; ++cmp)
            {
                for (size_t v = 0; v < vVertex.size(); ++v)
                {
                    std::vector<DoFIndex>  vInd;
                    if(dd->inner_dof_indices(vVertex[v], cmp, vInd) != 1)
                        UG_THROW("1: error in inner_dof_indices operatio!\n");
                    
                    size_t Index = Ind[v][0];
                    
                    
                    UG_LOG("Index = " << Index << "\n");
                    UG_LOG("vInd_0[0] = " << vInd[0] << "\n");
                    
                    //////////////////////////////////////////////////
                    // (2) data to be computed for the gradient computations
                    number gradD_elem;
                    number gradA_elem;
                    
                    // get indices modulo 2
                    size_t ind_1 = (v+1)%3;
                    size_t ind_2 = (v+2)%3;
                    size_t cmp_shift = (cmp+1)%2;
                    
                    UG_LOG("-------> cmp_shift = " << cmp_shift << "\n");
                    
                    
                    if ( cmp_shift == 0 ) // ind_2 - ind_1
                        gradA_elem = 0.5 * aaPos[vVertex[ind_2]][cmp_shift] - aaPos[vVertex[ind_1]][cmp_shift];
                    else if ( cmp_shift == 1 ) // ind_1 - ind_2
                        gradA_elem = 0.5 * aaPos[vVertex[ind_1]][cmp_shift] - aaPos[vVertex[ind_2]][cmp_shift];
                    else UG_THROW("cmp_shift not valid: " << cmp_shift << "\n");
                    
                    //gradA_elem = 0.5 * aaPos[vVertex[ind_1]][(cmp+1)%2] - aaPos[vVertex[ind_2]][(cmp+1)%2]);
                    UG_LOG("-------> gradA_elem = " << gradA_elem << "\n");
                    
                    //gradA_elem *= gradA_scale;
                    
                    
                    // ToDo: Vorzeichen für Ableitung von Area A(x1, x2, x3) checken: wirkliche fabs()?
                    
                    if ( (vVertex[v] != vrtOut) )
                    {
                        if ( v == vIndex_Out ) UG_THROW("1: error in vIndex_Out computation!\n");
                        
                        gradD_elem = 2*dist[cmp];
                    }
                    else
                    {
                        if ( v != vIndex_Out ) UG_THROW("2: error in vIndex_Out computation!\n");
                        
                        gradD_elem = 0.0;
                    }
                    
                    
                    
                    //////////////////////////////////////////////////
                    // (3) gradient computations:
                    
                    
                    //////////////////////////////////////////////////
                    // (3.1) volume computations:
                    if ( rSH->get_subset_index(elem) == 0 )
                        
                    {
                        UG_LOG("elem subset index = 0:" << rSH->get_subset_index(elem) << "\n");
                        
                        number gradD_kVol = (areaTria*gradD_elem - baseLineSquared*gradA_elem)/(2*areaTria*areaTria);
                        
                        //////////////////////////////////////////////////
                        // (4) add up all computations to entry of DoF-vector:
                        //      gradD[Index]           +=  gradD_kVol_x/(n*n); // #instead-ToDo:: DoFRef(u, Ind[v]) += ...
                        //     gradD[numDoFs + Index] +=  gradD_kVol_y/(n*n); // ToDo!
                        DoFRef(u, vInd[0]) += (-1.0) * gradD_kVol/(n*n);
                        if ( vInd[0][0] == 18 && cmp == 0 )
                        {
                            UG_LOG("gradD_elem = " << gradD_elem << "\n");
                            UG_LOG("gradA_elem = " << gradA_elem << "\n");
                            UG_LOG("gradA_scale = " << gradA_scale << "\n");
                            UG_LOG("areaTria = " << areaTria << "\n");
                            UG_LOG("baseLineSquared = " << baseLineSquared << "\n");
                            
                            UG_LOG("gradD_kVol = " << gradD_kVol << "\n");
                            UG_LOG("DoFRef(u, vInd[0][" << cmp << "]) = " << DoFRef(u, vInd[0]) << "\n");
                        }
                    }
                    
                    //////////////////////////////////////////////////
                    // (3.2) surface computations:
                    else if ( rSH->get_subset_index(elem) == 4 )
                    {
                        UG_LOG("elem subset index = 4:" << rSH->get_subset_index(elem) << "\n");
                        
                        number gradD_kSurf = (areaTria*gradD_elem - baseLineSquared*gradA_elem)/(2*areaTria*areaTria);
                        
                        //////////////////////////////////////////////////////////
                        // (4) add up all computations to entry of DoF-vector:
                        //     gradD[Index]           -= gradD_kSurf_x*(n-1)/(n*n); // #instead-ToDo:: DoFRef(u, Ind[v]) += ...
                        //    gradD[numDoFs + Index] -= gradD_kSurf_y*(n-1)/(n*n); // ToDo!
                        DoFRef(u, vInd[0]) -= (-1.0) * gradD_kSurf*(n-1)/(n*n);
                        if ( vInd[0][0] == 18 && cmp == 0 )
                        {
                            UG_LOG("gradD_elem = " << gradD_elem << "\n");
                            UG_LOG("gradA_elem = " << gradA_elem << "\n");
                            UG_LOG("gradA_scale = " << gradA_scale << "\n");
                            UG_LOG("areaTria = " << areaTria << "\n");
                            UG_LOG("baseLineSquared = " << baseLineSquared << "\n");
                            
                            UG_LOG("gradD_kSurf = " << gradD_kSurf << "\n");
                            UG_LOG("DoFRef(u, vInd[0][" << cmp << "]) = " << DoFRef(u, vInd[0]) << "\n");
                        }
                    }
                    else
                    {
                        UG_LOG("elem subset index = " << rSH->get_subset_index(elem) << "\n");
                        UG_THROW("--> why more than 2 subsets possible?\n");
                    }
                    
                    
                    
                    UG_LOG("\n ---> v = " << v << " and Ind = " << Ind[v] << "\n");
                    
                }// end vrt-loop
                
            }// end cmp-loop
            
        } // end elem-loop
        
        number functional = D_kVol/(n*n) - D_kSurf*(n-1)/(n*n);
        functional *= functional;
        
        return functional;
        
    }
 */

template <typename TDomain, typename TAlgebra>
number MovingParticle<TDomain, TAlgebra>::
compute_functional_all(const size_t n, vector_type& u, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace,
                       SmartPtr<NavierStokesFV1<TDomain> > spMaster, const int topLevel)
{
    ConstSmartPtr<DoFDistribution> dd = spApproxSpace->dof_distribution(GridLevel(topLevel, GridLevel::LEVEL));
        
    //	get data
    typename TDomain::position_accessor_type aaPos = this->m_spApproxSpace->domain()->position_accessor();
        
    //	create Function Group
    FunctionGroup fctGrp(dd->function_pattern(), TokenizeString(m_spParticleHandlerGlobal->m_fctNames));
    
    // get subset handler
    ConstSmartPtr<ISubsetHandler> rSH = dd->subset_handler();
    
    //	get iterators for all elems on subset
    typedef typename DoFDistribution::dim_traits<dim>::grid_base_object grid_base_object;
    
    typename DoFDistribution::traits<grid_base_object>::const_iterator iter, iterEnd;
    iter = dd->begin<grid_base_object>();
    iterEnd = dd->end<grid_base_object>();
    
    size_t numDoFs = 0.5 * (n+1) * (n+2); // = 45; 36 :)
    
    number D_kVol = 0.0;
    number D_kSurf = 0.0;
    std::vector<number> gradD(2*numDoFs, 0.0);
    
    int counter1, counter2, counter3, counter0, counter5;
    counter1 = counter2 = counter3 = counter5 = counter0 = 0;
    
    
    //	loop elements in order to compute 'U_global' and 'omega_global':
    for( ; iter != iterEnd; ++iter)
    {
        counter0 += 1;
        //	get element
        grid_base_object* elem = *iter;
        
        std::vector<Vertex*> vVertex;
        CollectVertices(vVertex, *m_spParticleHandlerGlobal->m_spMG, elem);
        for (size_t v = 0; v < vVertex.size(); ++v)
            UG_LOG("vrt_" << v << ": " << aaPos[vVertex[v]][0] << "\t" <<  aaPos[vVertex[v]][1] << "\n");
        
            
        //////////////////////////////////////////////////
        // compute data:
        std::vector<MathVector<dim> > vCornerCoords;
        CollectCornerCoordinates(vCornerCoords, *elem, aaPos, true);
        
        number areaTria = ElementSize<dim>(ROID_TRIANGLE, &vCornerCoords[0]);
        
        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        // (1) functional computations:
        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        
        //////////////////////////////////////////////////
        // volume computations ---> ALWAYS!!
        
        std::vector<Edge*> vEdges;
        CollectEdgesSorted(vEdges, *m_spParticleHandlerGlobal->m_spMG, elem);
        
        // loop edges
        for(size_t e = 0; e < vEdges.size(); ++e)
        {
            Edge* edge = vEdges[e];

            std::vector<Vertex*> vVertexEdge;
            CollectVertices(vVertexEdge, *m_spParticleHandlerGlobal->m_spMG, edge);
            if ( vVertexEdge.size() != 2 )
                UG_THROW("error in collecting vertices associated to an edge!....EXIT!...\n");
            
            number baseLineSquared = VecDistanceSq(aaPos[vVertexEdge[0]], aaPos[vVertexEdge[1]]);
 
            D_kVol += baseLineSquared/(2*areaTria);

        }
        
        
        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        // surface computations - again:
        
        ////////////////////////////////////////////////////////////////////////////////
        // FIRST: get subSet of boundary element (exclude inner elements (0) )
        int subSet = rSH->get_subset_index(elem);

        ////////////////////////////////////////////////////////////////////////////////
        // data filled for subSet = 1,2,3 and 5
        // AND needed for later computations!
        size_t indexA, indexB, indexC;
        std::vector<Edge*> vEdgesS(3);
        std::vector<number> vBaseLine(3);
        number sSq = 0.0;
        number sSq_all = 0.0;
        
        std::vector<Edge*> vEdges_buffer;
        CollectEdgesSorted(vEdges_buffer, *m_spParticleHandlerGlobal->m_spMG, elem);
        
        if ( subSet == 5 )
        {
            counter5 += 1;
        //	collect all vertices of the element
        std::vector<Vertex*> vVertexS_buffer;
        std::vector<Vertex*> vVertexS(3);
        CollectVertices(vVertexS_buffer, *m_spParticleHandlerGlobal->m_spMG, elem);
            
        // loop all vertices; treat them as vertex A as in the case below; add up and weight with 1/3:
        for(size_t countA = 0; countA < 3; ++countA)
        {
            // now store the other indices in order:
            indexA = countA;
            indexC = (countA+1)%3;
            indexB = (countA+2)%3;
            
            ////////////////////////////////////////////////////////////////////////////////
            // (1) collect 'vVertexS':
            vVertexS[0] = vVertexS_buffer[indexA];
            vVertexS[1] = vVertexS_buffer[indexC];
            vVertexS[2] = vVertexS_buffer[indexB];
            
            UG_LOG("---------> vVertexS: " << aaPos[vVertexS[0]][0] << "\t" <<  aaPos[vVertexS[0]][1] << "\n");
            UG_LOG("---------> vVertexS: " << aaPos[vVertexS[1]][0] << "\t" <<  aaPos[vVertexS[1]][1] << "\n");
            UG_LOG("---------> vVertexS: " << aaPos[vVertexS[2]][0] << "\t" <<  aaPos[vVertexS[2]][1] << "\n");
            
            ////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////
            // (B) loop edges ---> for boundary elements with subset 'subSet':
            
            // loop edges
            for(size_t e = 0; e < vEdges_buffer.size(); ++e)
            {
                Edge* edge = vEdges_buffer[e];
                
                std::vector<Vertex*> vVertexEdge;
                CollectVertices(vVertexEdge, *m_spParticleHandlerGlobal->m_spMG, edge);
                if ( vVertexEdge.size() != 2 )
                    UG_THROW("error in collecting vertices associated to an edge!....EXIT!...\n");
 
            // choose the edge, which does NOT contain the vertex A:
                if ( (vVertexEdge[0] != vVertexS[0]) && (vVertexEdge[1] != vVertexS[0]) )
                {
                    vEdgesS[0] = edge;
                    indexA = e;
                    break;
                }
            }
            // now store the other indices in order:
            indexB = (indexA+1)%3;
            indexC = (indexA+2)%3;
            
            //ToDo: egal, was indexB und indexC ist, solange NICHT indexA = subSet??
            
            ////////////////////////////////////////////////////////////////////////////////
            // (2) collect 'vEdgesS':
            vEdgesS[1] = vEdges_buffer[indexB];
            vEdgesS[2] = vEdges_buffer[indexC];
            
            ////////////////////////////////////////////////////////////////////////////////
            // (3) collect 'vBaseLine':
            for(size_t i = 0; i < 3; ++i)
            {
                std::vector<Vertex*> vVertexEdge;
                CollectVertices(vVertexEdge, *m_spParticleHandlerGlobal->m_spMG, vEdgesS[i]);
                if ( vVertexEdge.size() != 2 )
                    UG_THROW("error in collecting vertices associated to an edge!....EXIT!...\n");
                
                UG_LOG("---------> vVertexEdge1: " << aaPos[vVertexEdge[0]][0] << "\t" <<  aaPos[vVertexEdge[0]][1] << "\n");
                UG_LOG("---------> vVertexEdge2: " << aaPos[vVertexEdge[1]][0] << "\t" <<  aaPos[vVertexEdge[1]][1] << "\n");
                
                number baseLineSquared = VecDistanceSq(aaPos[vVertexEdge[0]], aaPos[vVertexEdge[1]]);
                vBaseLine[i] = sqrt(baseLineSquared);
            }
            
            ////////////////////////////////////////////////////////////////////////////////
            // THIRD: compute data for surface integrals ---> s, alphaB, alphaC
            
            number a = vBaseLine[0];
            number b = vBaseLine[1];
            number c = vBaseLine[2];
            number aSq = a*a;
            number bSq = b*b;
            number cSq = c*c;
            
            number s = 0.5*sqrt(2 * (bSq+cSq) - aSq);
            sSq = s*s;
            sSq_all += sSq;
            UG_LOG("---------> sSq: " << sSq << "\n");
            UG_LOG("---------> sSq_all: " << sSq_all << "\n");

            number cosAlphaB = (bSq + sSq - 0.25*aSq)/(2*b*s);
            number cosAlphaC = (cSq + sSq - 0.25*aSq)/(2*c*s);
            
            ////////////////////////////////////////////////////////////////////////////////
            // FINALLY: compute fluxes
            
            number fluxB = 0.5*s*b*cosAlphaB/areaTria ;
            number fluxC = 0.5*s*c*cosAlphaC/areaTria ;
            
            D_kSurf += (fluxB + fluxC)/3.0;
            
            number added = (fluxB + fluxC)/3.0;;
            UG_LOG("subSet = " << subSet << "added = " << added << "\n");
            
        } // end indexA-loop
        } // end case subSet == 5
        
        
        ////////////////////////////////////////////////////////////////////////////////
        // SECOND: loop vertices AND edges --> write data according to specified indices 'indexA, indexB, indexC':
        
        if ( subSet == 1 )  counter1 += 1;
        if ( subSet == 2 )  counter2 += 1;
        if ( subSet == 3 )  counter3 += 1;
        
        ////////////////////////////////////////////////////////////////////////////////
   
        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        // loop edges ---> for boundary elements with subset 'subSet':
        
        ////////////////////////////////////////////////////////////////////////////////
        // THIRD: compute data for surface integrals ---> s, alphaB, alphaC
        
        if ( subSet == 1 || subSet == 2 || subSet == 3 )
        {
            
        // loop edges
        for(size_t e = 0; e < vEdges_buffer.size(); ++e)
        {
            if ( rSH->get_subset_index(vEdges_buffer[e]) == subSet )
            {
                vEdgesS[0] = vEdges_buffer[e];
                indexA = e;
                break;
            }
        }
        // now store the other indices in order:
        indexB = (indexA+1)%3;
        indexC = (indexA+2)%3;
        
        //ToDo: egal, was indexB und indexC ist, solange NICHT indexA = subSet??
        
        ////////////////////////////////////////////////////////////////////////////////
        // (2) collect 'vEdgesS':
        vEdgesS[1] = vEdges_buffer[indexB];
        vEdgesS[2] = vEdges_buffer[indexC];
        
        ////////////////////////////////////////////////////////////////////////////////
        // (3) collect 'vBaseLine':
        for(size_t i = 0; i < 3; ++i)
        {
            std::vector<Vertex*> vVertexEdge;
            CollectVertices(vVertexEdge, *m_spParticleHandlerGlobal->m_spMG, vEdgesS[i]);
            if ( vVertexEdge.size() != 2 )
                UG_THROW("error in collecting vertices associated to an edge!....EXIT!...\n");
            
            UG_LOG("---------> vVertexEdge1: " << aaPos[vVertexEdge[0]][0] << "\t" <<  aaPos[vVertexEdge[0]][1] << "\n");
            UG_LOG("---------> vVertexEdge2: " << aaPos[vVertexEdge[1]][0] << "\t" <<  aaPos[vVertexEdge[1]][1] << "\n");

            number BaseLineSquared = VecDistanceSq(aaPos[vVertexEdge[0]], aaPos[vVertexEdge[1]]);
            vBaseLine[i]        = sqrt(BaseLineSquared);

        }
            
        number a = vBaseLine[0];
        number b = vBaseLine[1];
        number c = vBaseLine[2];
        number aSq = a*a;
        number bSq = b*b;
        number cSq = c*c;
        
        number s = 0.5*sqrt(2 * (bSq+cSq) - aSq);
        sSq = s*s;
        number cosAlphaB = (bSq + sSq - 0.25*aSq)/(2*b*s);
        number cosAlphaC = (cSq + sSq - 0.25*aSq)/(2*c*s);
        
        ////////////////////////////////////////////////////////////////////////////////
        // FINALLY: compute fluxes
        
        number fluxB_ = 0.5*s*b*cosAlphaB/areaTria ;
        number fluxC_ = 0.5*s*c*cosAlphaC/areaTria ;
        
        number fluxB = 0.5*0.25*(3*bSq + cSq - aSq)/areaTria ;
        number fluxC = 0.5*0.25*(3*cSq + bSq - aSq)/areaTria ;
        
        if ( fabs(fluxB - fluxB_) > 0.00001 )
            UG_THROW(" fluxB: " << fluxB << " and " << " fluxB_: " << fluxB_ << "\n");
        if ( fabs(fluxC - fluxC_) > 0.00001 )
            UG_THROW(" fluxC: " << fluxC << " and " << " fluxC_: " << fluxC_ << "\n");
  
        D_kSurf += fluxB + fluxC;
            number added = fluxB + fluxC;
            UG_LOG("subSet = " << subSet << "added = " << added << "\n");

        
        } // END: functional computations
        
        
        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        // compute gradients
        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////
        // compute scaling for gradA_elem:
        
        number _A_ =(aaPos[vVertex[1]][0]-aaPos[vVertex[0]][0])*(aaPos[vVertex[2]][1]-aaPos[vVertex[0]][1])
        -(aaPos[vVertex[2]][0]-aaPos[vVertex[0]][0])*(aaPos[vVertex[1]][1]-aaPos[vVertex[0]][1]);
        
        number gradA_scale = _A_ / fabs(_A_);
        
        ////////////////////////////////////////////////////////////////////////////////
        // loop cmp and vertices to compute 'gradS_elem' and
        //    ADD 'gradA_elem'/'gradA_elem' to gradient:
        
        number gradA_elem;
        number gradD_elem_all;
        number gradS_elem;
        std::vector<Vertex*> vVertexEdge;

        for (size_t cmp = 0; cmp < dim; ++cmp)
        {
            
            for (size_t v = 0; v < vVertex.size(); ++v)
            {
                std::vector<DoFIndex>  vInd;
                if(dd->inner_dof_indices(vVertex[v], cmp, vInd) != 1)
                    UG_THROW("1: error in inner_dof_indices operatio!\n");
                
                UG_LOG("vInd[0] = " << vInd[0] << "\n");
                
            ////////////////////////////////////////////////////////////////////////////////
            // initialize gradient for the denominator and add all the contributions from the edges
                gradS_elem = 0.0;
                gradD_elem_all = 0.0;
                
            ////////////////////////////////////////////////////////////////////////////////
            // compute 'gradA_elem':
                
                // get indices modulo 2
                size_t ind_1 = (v+1)%3;
                size_t ind_2 = (v+2)%3;
                size_t cmp_shift = (cmp+1)%2;
                
                // ind_2 - ind_1
                if ( cmp_shift == 0 )
                    gradA_elem = 0.5 * (aaPos[vVertex[ind_2]][cmp_shift] - aaPos[vVertex[ind_1]][cmp_shift]);
                // ind_1 - ind_2
                else if ( cmp_shift == 1 )
                    gradA_elem = 0.5 * (aaPos[vVertex[ind_1]][cmp_shift] - aaPos[vVertex[ind_2]][cmp_shift]);
                else
                    UG_THROW("cmp_shift not valid: " << cmp_shift << "\n");
                
                gradA_elem *= gradA_scale;
                
                
            ////////////////////////////////////////////////////////////////////////////////
            // compute 'gradS_elem' and 'gradD_elem_all':
                    
            //////////////////////////////////////////////////
            // edge contributions:
                
            number gradS_elem_scale = 1.0;
            std::vector<number> vBaseLineSquared(3);

            if ( subSet == 5 )
                { UG_LOG("---> elem subset index = :" << rSH->get_subset_index(elem) << "\n");}
                
            for(size_t e = 0; e < 3; ++e)
            {
                CollectVertices(vVertexEdge, *m_spParticleHandlerGlobal->m_spMG, vEdges_buffer[e]);
                if ( vVertexEdge.size() != 2 )
                    UG_THROW("---> error in collecting vertices associated to an edge!....EXIT!...\n");

                vBaseLineSquared[e] = VecDistanceSq(aaPos[vVertexEdge[0]], aaPos[vVertexEdge[1]]);

                UG_LOG("for 'gradS_elem_scale': vVertexEdge1: " << aaPos[vVertexEdge[0]][0] << "\t" <<  aaPos[vVertexEdge[0]][1] << "\n");
                UG_LOG("for 'gradS_elem_scale': vVertexEdge2: " << aaPos[vVertexEdge[1]][0] << "\t" <<  aaPos[vVertexEdge[1]][1] << "\n");
                

                // adapt scale factor for edge = baseLine:
                if ( rSH->get_subset_index(vEdges_buffer[e]) == 5 )
                    UG_THROW("---> hmmmm: EDGE subset index can NOT be 5! ---> index = :" << rSH->get_subset_index(vEdges_buffer[e]) << "\n");

                if ( rSH->get_subset_index(vEdges_buffer[e]) == subSet )
                    gradS_elem_scale = -0.5;
                if ( subSet == 5 )
                    gradS_elem_scale = 1.5;

                if ( vVertex[v] == vVertexEdge[0] )
                {
                    gradS_elem += gradS_elem_scale * 2.0 * (aaPos[vVertexEdge[0]][cmp] - aaPos[vVertexEdge[1]][cmp]);
                    gradD_elem_all += 2.0 * (aaPos[vVertexEdge[0]][cmp] - aaPos[vVertexEdge[1]][cmp]);
                }
                else if ( vVertex[v] == vVertexEdge[1] )
                {
                    gradS_elem += gradS_elem_scale * 2.0 * (aaPos[vVertexEdge[1]][cmp] - aaPos[vVertexEdge[0]][cmp]);
                    gradD_elem_all += 2.0 * (aaPos[vVertexEdge[0]][cmp] - aaPos[vVertexEdge[1]][cmp]);
                }
                // else: add 0.0 ==> do nothing...
                
            } // end edge-loop
            
            number vBaseLineSquared_all = vBaseLineSquared[0] + vBaseLineSquared[1] + vBaseLineSquared[2];
            
                
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            // ADD up all v-th derivatives to 'cmp' of vector 'u':
             
            //////////////////////////////////////////////////
            // (3.1) volume computations:
                
                number gradD_kVol = (areaTria*gradD_elem_all - vBaseLineSquared_all*gradA_elem)/(2*areaTria*areaTria);
                
                DoFRef(u, vInd[0]) += (-1.0) * gradD_kVol/(n*n);

                if (vInd[0][0] == 18 )
                    UG_LOG("18: added: " << (-1.0) * gradD_kVol/(n*n) << "\t value: " << DoFRef(u, vInd[0]) << "\n");
                
                if (vInd[0][0] == 19 )
                    UG_LOG("19: added: " << (-1.0) * gradD_kVol/(n*n) << "\t value: " << DoFRef(u, vInd[0]) << "\n");
                
                if (vInd[0][0] == 20 )
                    UG_LOG("20: added: " << (-1.0) * gradD_kVol/(n*n) << "\t value: " << DoFRef(u, vInd[0]) << "\n");
                
                
                if (vInd[0][0] == 21 )
                    UG_LOG("21: added: " << (-1.0) * gradD_kVol/(n*n) << "\t value: " << DoFRef(u, vInd[0]) << "\n");
                
                
                if (vInd[0][0] == 24 )
                    UG_LOG("24: added: " << (-1.0) * gradD_kVol/(n*n) << "\t value: " << DoFRef(u, vInd[0]) << "\n");
                
                if (vInd[0][0] == 16 )
                    UG_LOG("16: added: " << (-1.0) * gradD_kVol/(n*n) << "\t value: " << DoFRef(u, vInd[0]) << "\n");
                
                
                if ( vInd[0][0] == 18 && cmp == 0 )
                {
                    UG_LOG("gradD_elem_all = " << gradD_elem_all << "\n");
                    UG_LOG("gradA_elem = " << gradA_elem << "\n");
                    UG_LOG("gradA_scale = " << gradA_scale << "\n");
                    UG_LOG("areaTria = " << areaTria << "\n");
                    UG_LOG("vBaseLineSquared_all = " << vBaseLineSquared_all << "\n");
                    
                    UG_LOG("gradD_kVol = " << gradD_kVol << "\n");
                    UG_LOG("DoFRef(u, vInd[0][" << cmp << "]) = " << DoFRef(u, vInd[0]) << "\n");
                }
            
            //////////////////////////////////////////////////
            // (3.2) surface computations:
                if ( 0 ) //rSH->get_subset_index(elem) != 0 )
                {
                    if ( subSet != 1 && subSet != 2 && subSet != 3 )
                    {
                        if ( subSet == 5 )
                        {
                            sSq = sSq_all/3.0;
                            gradS_elem = gradS_elem/3.0;
                            UG_LOG("elem subset index = :" << rSH->get_subset_index(elem) << "\n");
                        }
                        else
                        { UG_THROW("elem subset index = :" << rSH->get_subset_index(elem) << "\n");}
                    }
                    
                    number gradS_kSurf = 0.5 * (areaTria*gradS_elem - sSq*gradA_elem)/(areaTria*areaTria);
                    
                    DoFRef(u, vInd[0]) -= (-1.0) * gradS_kSurf/n;
                    if ( vInd[0][0] == 18 && cmp == 0 )
                    {
                        UG_LOG("gradS_elem = " << gradS_elem << "\n");
                        UG_LOG("gradA_elem = " << gradA_elem << "\n");
                        UG_LOG("gradA_scale = " << gradA_scale << "\n");
                        UG_LOG("areaTria = " << areaTria << "\n");
                        
                        UG_LOG("gradS_kSurf = " << gradS_kSurf << "\n");
                        UG_LOG("DoFRef(u, vInd[0][" << cmp << "]) = " << DoFRef(u, vInd[0]) << "\n");
                    }
                }
                
                
            }// end vrt-loop
            
        }// end cmp-loop
        
    } // end elem-loop
    
    number functional = 0.5 * D_kVol/(n*n) - D_kSurf/n;
    functional *= functional;
    
    UG_LOG(" counter0: " << counter0 << "\n");
    UG_LOG(" counter1: " << counter1 << "\n");
    UG_LOG(" counter2: " << counter2 << "\n");
    UG_LOG(" counter3: " << counter3 << "\n");
    UG_LOG(" counter5: " << counter5 << "\n");

    
    return functional;
    
}

    
template <typename TDomain, typename TAlgebra>
number MovingParticle<TDomain, TAlgebra>::
compute_functional(const size_t n, vector_type& u, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace,
                      SmartPtr<NavierStokesFV1<TDomain> > spMaster, const int topLevel)
{
        ConstSmartPtr<DoFDistribution> dd = spApproxSpace->dof_distribution(GridLevel(topLevel, GridLevel::LEVEL));

        //	get data
        typename TDomain::position_accessor_type aaPos = this->m_spApproxSpace->domain()->position_accessor();
        
        //	create Function Group
        FunctionGroup fctGrp(dd->function_pattern(), TokenizeString(m_spParticleHandlerGlobal->m_fctNames));
    
        // get subset handler
        ConstSmartPtr<ISubsetHandler> rSH = dd->subset_handler();

        //	get iterators for all elems on subset
        typedef typename DoFDistribution::dim_traits<dim>::grid_base_object grid_base_object;
        
        typename DoFDistribution::traits<grid_base_object>::const_iterator iter, iterEnd;
        iter = dd->begin<grid_base_object>();
        iterEnd = dd->end<grid_base_object>();
        
        size_t numDoFs = 0.5 * (n+1) * (n+2); // = 45; 36 :)
    

        number D_kVol = 0.0;
        number D_kSurf = 0.0;
        std::vector<number> gradD(2*numDoFs, 0.0);

        //	loop elements in order to compute 'U_global' and 'omega_global':
        for( ; iter != iterEnd; ++iter)
        {
            //	get element
            grid_base_object* elem = *iter;
            
            std::vector<Vertex*> vVertex;
            CollectVertices(vVertex, *m_spParticleHandlerGlobal->m_spMG, elem);

            //	loop vertices
            for (size_t v = 0; v < vVertex.size(); ++v)
            {
                UG_LOG("vrt_" << v << ": " << aaPos[vVertex[v]][0] << "\t" <<  aaPos[vVertex[v]][1] << "\n");
            }
            
        //////////////////////////////////////////////////
        // compute data:
            std::vector<MathVector<dim> > vCornerCoords;
            CollectCornerCoordinates(vCornerCoords, *elem, aaPos, true);
            
            number areaTria = ElementSize<dim>(ROID_TRIANGLE, &vCornerCoords[0]);
         
        //////////////////////////////////////////////////
        // loop edges and:
        //  -> compute 'baseLineSquared', 'baseLine', 'dist'
        //  -> collect associated vertex 'vrtOut', 'vIndex_Out'
            
            std::vector<Edge*> vEdges;
            CollectEdgesSorted(vEdges, *m_spParticleHandlerGlobal->m_spMG, elem);
            
            // values to be computed:
            number baseLineSquared, baseLine;
            MathVector<dim> dist;
            size_t vIndex_Out, vIndex_1, vIndex_2;
            Vertex* vrtOut;
            Vertex* vrt1;
            Vertex* vrt2;
            
            for(size_t e = 0; e < vEdges.size(); ++e)
            {
                Edge* edge = vEdges[e];
                
                if ( rSH->get_subset_index(edge) != 3 )
                {
                    UG_LOG("...continue...subset = " << rSH->get_subset_index(edge) << "\n");
                    continue;
                }
                
                UG_LOG("subset = " << rSH->get_subset_index(edge) << "\n");
                
                std::vector<Vertex*> vVertexEdge;
                CollectVertices(vVertexEdge, *m_spParticleHandlerGlobal->m_spMG, edge);
                if ( vVertexEdge.size() != 2 )
                    UG_THROW("error in collecting vertices associated to an edge!....EXIT!...\n");
                
                vrt1 = vVertexEdge[0];
                vrt2 = vVertexEdge[1];
                
                baseLineSquared = VecDistanceSq(aaPos[vrt1], aaPos[vrt2]);
                baseLine = sqrt(baseLineSquared);
                dist[0] = fabs(aaPos[vrt1][0] - aaPos[vrt2][0]);
                dist[1] = fabs(aaPos[vrt1][1] - aaPos[vrt2][1]);

                //	loop vertices and get vrtOut:
                for (size_t v = 0; v < vVertex.size(); ++v)
                {
                    if ( (vVertex[v] != vrt1) && (vVertex[v] != vrt2) )
                    {
                        vrtOut = vVertex[v];
                        vIndex_Out = v;
                    }
                    if ( vVertex[v] == vrt1 )
                        vIndex_1 = v;
                    if ( vVertex[v] == vrt2 )
                        vIndex_2 = v;
                }
                
                UG_LOG("vrt1: " << aaPos[vrt1][0] << "\t" <<  aaPos[vrt1][1] << "\n");
                UG_LOG("vrt2: " << aaPos[vrt2][0] << "\t" <<  aaPos[vrt2][1] << "\n");
                UG_LOG("vrtOut: " << aaPos[vrtOut][0] << "\t" <<  aaPos[vrtOut][1] << "\n");
                
            } // end edge-loop
            UG_LOG("new edge\n");
            
        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        // (1) functional computations:
        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////

        //////////////////////////////////////////////////
        // volume computations ---> ALWAYS!!
            if ( 1 ) //rSH->get_subset_index(elem) == 0 )
                
            {
            // Attention: the factor 1/2 is removed from D_kVol AND D_kSurf, since it arises in both summands!
                D_kVol += baseLineSquared/(2.0*areaTria);
             }
        //////////////////////////////////////////////////
        // surface computations:
            if ( rSH->get_subset_index(elem) == 4 )
            {
            // Attention: the factor 1/2 is removed from D_kVol AND D_kSurf, since it arises in both summands!
                D_kSurf += baseLineSquared/(2.0*areaTria);
            }

            
        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        // (2) gradient computations:
        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////

            
        //////////////////////////////////////////////////
        // (1) get global indices of vrt1, vrt2, vrtOut for final adding up to DoF-vector in (4):
            
            std::vector<DoFIndex> Ind;
            dd->dof_indices(elem, 1, Ind);  // #instead-ToDo: dd->dof_indices(elem, cmp, Ind);
            
            for (size_t v = 0; v < vVertex.size(); ++v)
                UG_LOG("2 Ind = " << Ind[v] << "\n");
    /*
     template <>
     inline number ElementSize<ReferenceTriangle, 2>(const MathVector<2>* vCornerCoords)
     {
        return(0.5*fabs((vCornerCoords[1][1]-vCornerCoords[0][1])*(vCornerCoords[2][0]-vCornerCoords[0][0])
                       -(vCornerCoords[1][0]-vCornerCoords[0][0])*(vCornerCoords[2][1]-vCornerCoords[0][1])));
     }
     */

            number _A_ =(aaPos[vVertex[1]][0]-aaPos[vVertex[0]][0])*(aaPos[vVertex[2]][1]-aaPos[vVertex[0]][1])
                       -(aaPos[vVertex[2]][0]-aaPos[vVertex[0]][0])*(aaPos[vVertex[1]][1]-aaPos[vVertex[0]][1]);
            
            number gradA_scale = _A_ / fabs(_A_);
            
    //        if ( gradA_scale > 0 )
     //           UG_LOG("gradA_scale = " << gradA_scale << "\n");
            
        ////////////////////////////////////////////////////////////////////////////////
        // #instead-ToDo: loop components 'cmp' = x- and y-direction:
        for (size_t cmp = 0; cmp < dim; ++cmp)
        {
            for (size_t v = 0; v < vVertex.size(); ++v)
            {
                std::vector<DoFIndex>  vInd;
                if(dd->inner_dof_indices(vVertex[v], cmp, vInd) != 1)
                    UG_THROW("1: error in inner_dof_indices operatio!\n");
 
                size_t Index = Ind[v][0];

                
                UG_LOG("Index = " << Index << "\n");
                UG_LOG("vInd_0[0] = " << vInd[0] << "\n");
                
                //////////////////////////////////////////////////
                // (2) data to be computed for the gradient computations
                number gradD_elem;
                number gradA_elem;
                
                // get indices modulo 2
                size_t ind_1 = (v+1)%3;
                size_t ind_2 = (v+2)%3;
                size_t cmp_shift = (cmp+1)%2;

                UG_LOG("-------> cmp_shift = " << cmp_shift << "\n");

                
                if ( cmp_shift == 0 ) // ind_2 - ind_1
                    gradA_elem = 0.5 * (aaPos[vVertex[ind_2]][cmp_shift] - aaPos[vVertex[ind_1]][cmp_shift]);
                else if ( cmp_shift == 1 ) // ind_1 - ind_2
                    gradA_elem = 0.5 * (aaPos[vVertex[ind_1]][cmp_shift] - aaPos[vVertex[ind_2]][cmp_shift]);
                else UG_THROW("cmp_shift not valid: " << cmp_shift << "\n");
                
                gradA_elem *= gradA_scale;
                
                
                // ToDo: Vorzeichen für Ableitung von Area A(x1, x2, x3) checken: wirkliche fabs()?
                
                if ( vVertex[v] != vrtOut )
                {
                    if ( v == vIndex_Out ) UG_THROW("1: error in vIndex_Out computation!\n");

                    gradD_elem = 2.0 * (aaPos[vrt1][cmp] - aaPos[vrt2][cmp]);

                    if ( vVertex[v] == vrt2 )
                    {
                        gradD_elem *= -1.0;
                        UG_LOG("ind_1: " << ind_1 << "\n");
                        UG_LOG("vIndex_1: " << vIndex_1 << "\n");
                        UG_LOG("ind_2: " << ind_2 << "\n");
                        UG_LOG("vIndex_2: " << vIndex_2 << "\n");
                    }

                }
                else
                {
                    if ( v != vIndex_Out ) UG_THROW("2: error in vIndex_Out computation!\n");

                    gradD_elem = 0.0;
                }
                

                
            //////////////////////////////////////////////////
            // (3) gradient computations:

                
            //////////////////////////////////////////////////
            // (3.1) volume computations:
                if ( 1 ) //rSH->get_subset_index(elem) == 0 )
                    
                {
                    UG_LOG("elem subset index = 0:" << rSH->get_subset_index(elem) << "\n");
                    
                    number gradD_kVol = (areaTria*gradD_elem - baseLineSquared*gradA_elem)/(2*areaTria*areaTria);
                    
                    //////////////////////////////////////////////////
                    // (4) add up all computations to entry of DoF-vector:
              //      gradD[Index]           +=  gradD_kVol_x/(n*n); // #instead-ToDo:: DoFRef(u, Ind[v]) += ...
               //     gradD[numDoFs + Index] +=  gradD_kVol_y/(n*n); // ToDo!
                    DoFRef(u, vInd[0]) += (-1.0) * gradD_kVol/(n*n);
                    if ( vInd[0][0] == 18 && cmp == 0 )
                    {
                        UG_LOG("gradD_elem = " << gradD_elem << "\n");
                        UG_LOG("gradA_elem = " << gradA_elem << "\n");
                        UG_LOG("gradA_scale = " << gradA_scale << "\n");
                        UG_LOG("areaTria = " << areaTria << "\n");
                        UG_LOG("baseLineSquared = " << baseLineSquared << "\n");

                        UG_LOG("gradD_kVol = " << gradD_kVol << "\n");
                        UG_LOG("DoFRef(u, vInd[0][" << cmp << "]) = " << DoFRef(u, vInd[0]) << "\n");
                    }
                }
                
            //////////////////////////////////////////////////
            // (3.2) surface computations:
                if ( rSH->get_subset_index(elem) == 4 )
                {
                    UG_LOG("elem subset index = 4:" << rSH->get_subset_index(elem) << "\n");

                    number gradD_kSurf = (areaTria*gradD_elem - baseLineSquared*gradA_elem)/(2*areaTria*areaTria);
                    
                    //////////////////////////////////////////////////////////
                    // (4) add up all computations to entry of DoF-vector:
               //     gradD[Index]           -= gradD_kSurf_x*(n-1)/(n*n); // #instead-ToDo:: DoFRef(u, Ind[v]) += ...
                //    gradD[numDoFs + Index] -= gradD_kSurf_y*(n-1)/(n*n); // ToDo!
                    DoFRef(u, vInd[0]) -= (-1.0) * gradD_kSurf/n;
                    if ( vInd[0][0] == 18 && cmp == 0 )
                    {
                        UG_LOG("gradD_elem = " << gradD_elem << "\n");
                        UG_LOG("gradA_elem = " << gradA_elem << "\n");
                        UG_LOG("gradA_scale = " << gradA_scale << "\n");
                        UG_LOG("areaTria = " << areaTria << "\n");
                        UG_LOG("baseLineSquared = " << baseLineSquared << "\n");
                        
                        UG_LOG("gradD_kSurf = " << gradD_kSurf << "\n");
                        UG_LOG("DoFRef(u, vInd[0][" << cmp << "]) = " << DoFRef(u, vInd[0]) << "\n");
                    }
                }
                


                UG_LOG("\n ---> v = " << v << " and Ind = " << Ind[v] << "\n");

            }// end vrt-loop
            
        }// end cmp-loop
            
    } // end elem-loop
    
    // Attention: the factor 1/2 is removed from D_kVol AND D_kSurf, since it arises in both summands!
    number functional = D_kVol/(n*n) - D_kSurf/n;
    functional *= functional;
    
    return functional;
    
}

template <typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::
rescale_directions(const number functional, vector_type& u, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace,
                       SmartPtr<NavierStokesFV1<TDomain> > spMaster, const int topLevel)
{
    /////////////////////////////////////////////////
    // loop vertices and project directions
    //  of 'u' on outer boundary edges:
        
    ConstSmartPtr<DoFDistribution> dd = spApproxSpace->dof_distribution(GridLevel(topLevel, GridLevel::LEVEL));
    DoFDistribution::traits<Vertex>::const_iterator iterBegin 	= dd->begin<Vertex>();
    DoFDistribution::traits<Vertex>::const_iterator iterEnd 	= dd->end<Vertex>();
        
    // get subset handler
    ConstSmartPtr<ISubsetHandler> rSH = dd->subset_handler();
 
    number scaleFactor = sqrt(functional);
    
    //	loop vertices
    for(DoFDistribution::traits<Vertex>::const_iterator iter = iterBegin; iter != iterEnd; iter++)
    {
        //	get vertex
        Vertex* vrt = *iter;
        
        for (size_t cmp = 0; cmp < dim; ++cmp)
        {
            std::vector<DoFIndex>  vInd;
            if(dd->inner_dof_indices(vrt, cmp, vInd) != 1)
                UG_THROW("1: error in inner_dof_indices operatio!\n");
            
            DoFRef(u, vInd[0]) *= 2.0 * scaleFactor;
        }

            
    } // end vertices-loop
        
}
    
template <typename TDomain, typename TAlgebra>
number MovingParticle<TDomain, TAlgebra>::
compute_max_step_size(SmartPtr<ApproximationSpace<TDomain> > spApproxSpace,
                     SmartPtr<NavierStokesFV1<TDomain> > spMaster, const int topLevel)
{
    bool output = false;
    
    ConstSmartPtr<DoFDistribution> dd = spApproxSpace->dof_distribution(GridLevel(topLevel, GridLevel::LEVEL));
    
    //	get data
    typename TDomain::position_accessor_type aaPos = this->m_spApproxSpace->domain()->position_accessor();
    
    //	create Function Group
    FunctionGroup fctGrp(dd->function_pattern(), TokenizeString(m_spParticleHandlerGlobal->m_fctNames));
    
    // get subset handler
    ConstSmartPtr<ISubsetHandler> rSH = dd->subset_handler();
    
    //	get iterators for all elems on subset
    typedef typename DoFDistribution::dim_traits<dim>::grid_base_object grid_base_object;
    
    typename DoFDistribution::traits<grid_base_object>::const_iterator iter, iterEnd;
    iter = dd->begin<grid_base_object>();
    iterEnd = dd->end<grid_base_object>();
    
    //	loop elements in order to compute minimal edge size:
    number minEdgeSize = 1000.0;
    int N = 0;
    for( ; iter != iterEnd; ++iter)
    {
        //	get element
        grid_base_object* elem = *iter;
        
        //////////////////////////////////////////////////
        // loop edges and collect associated vertices 'vrt1' and 'vrt2' and 'vrtOut:
        std::vector<Edge*> vEdges;
        CollectEdgesSorted(vEdges, *m_spParticleHandlerGlobal->m_spMG, elem);
        for(size_t e = 0; e < vEdges.size(); ++e)
        {
            Edge* edge = vEdges[e];
            
            std::vector<Vertex*> vVertexEdge;
            CollectVertices(vVertexEdge, *m_spParticleHandlerGlobal->m_spMG, edge);
            if ( vVertexEdge.size() != 2 )
                UG_THROW("error in collecting vertices associated to an edge!....EXIT!...\n");
            
            number edgeSize = VecDistance(aaPos[vVertexEdge[0]], aaPos[vVertexEdge[1]]);
            
            std::vector<DoFIndex>  vInd1;
            if(dd->inner_dof_indices(vVertexEdge[0], 0, vInd1) != 1)
                UG_THROW("1: error in inner_dof_indices operatio!\n");
            std::vector<DoFIndex>  vInd2;
            if(dd->inner_dof_indices(vVertexEdge[1], 0, vInd2) != 1)
                UG_THROW("1: error in inner_dof_indices operatio!\n");
            
            if ( vInd1[0][0] == 11 || vInd2[0][0] == 11 )
            {
                if ( vInd1[0][0] == 4 || vInd2[0][0] == 4 )
                {
                    //UG_LOG("minEdgeSize = " << minEdgeSize << "\n");
                    //UG_LOG("aaPos[vVertexEdge[0]] = " << aaPos[vVertexEdge[0]][0] << "\t" << aaPos[vVertexEdge[0]][1] << "\n");
                    //UG_LOG("aaPos[vVertexEdge[1]] = " << aaPos[vVertexEdge[1]][0] << "\t" << aaPos[vVertexEdge[1]][1] << "\n");

                }
            }
            
            
            if ( minEdgeSize > edgeSize )
            {
                minEdgeSize = edgeSize;
            
                if ( output ) UG_LOG("minEdgeSize = " << minEdgeSize << "\n");
            }
            else
            {
                if ( output ) UG_LOG("--> edgeSize = " << edgeSize << "\n");

            }
            
        }// end edge-loop
        N = N+1;

        if ( output ) UG_LOG("--> new element: " << N << "\n");

    }// end elem-loop
    
    UG_LOG("1 -----------> minEdgeSize = " << minEdgeSize << "\n");
    
    return minEdgeSize; 
}


template <typename TDomain, typename TAlgebra>
void MovingParticle<TDomain, TAlgebra>::
project_directions(const number functional, vector_type& u, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace,
                     SmartPtr<NavierStokesFV1<TDomain> > spMaster, const int topLevel)
{
    /////////////////////////////////////////////////
    // loop vertices and project directions
    //  of 'u' on outer boundary edges:
    
    ConstSmartPtr<DoFDistribution> dd = spApproxSpace->dof_distribution(GridLevel(topLevel, GridLevel::LEVEL));
    DoFDistribution::traits<Vertex>::const_iterator iterBegin 	= dd->begin<Vertex>();
    DoFDistribution::traits<Vertex>::const_iterator iterEnd 	= dd->end<Vertex>();
    
    // get subset handler
    ConstSmartPtr<ISubsetHandler> rSH = dd->subset_handler();
    
    size_t counter1 = 0;
    size_t counter2 = 0;
    size_t counter3 = 0;
    size_t counter4 = 0;
    size_t counter5 = 0;
    size_t counter6 = 0;

    //	loop vertices
    for(DoFDistribution::traits<Vertex>::const_iterator iter = iterBegin; iter != iterEnd; iter++)
    {
        //	get vertex
        Vertex* vrt = *iter;
         
        /////////////////////////////////////////////////
        // (1) get direction to be projected:
        MathVector<dim> project;

        for (size_t cmp = 0; cmp < dim; ++cmp)
        {
            std::vector<DoFIndex>  vInd;
            if(dd->inner_dof_indices(vrt, cmp, vInd) != 1)
                UG_THROW("1: error in inner_dof_indices operatio!\n");
        
            project[cmp] = DoFRef(u, vInd[0]);
        }
        
        /////////////////////////////////////////////////
        // (2) project according to projection edge:
        MathVector<dim> edge;
        // (2.1) project onto x-coordinate: (1,0)
        if ( rSH->get_subset_index(vrt) == 3 )
        {
            edge[0] = 1.0; edge[1] = 0.0;
            number scalar = VecProd(project, edge);
            
            VecScale(project, edge, scalar);
            ++counter1;
            //project[1] = 0.0;
        }
        // (2.2) project onto left-coordinate: (0.5,1)
        else if ( rSH->get_subset_index(vrt) == 2 )
        {
            edge[0] = 0.5; edge[1] = 1.0;
            number scalar = VecProd(project, edge);

            VecScale(project, edge, scalar);
            ++counter2;
        }
        
        // (2.3) project onto x-coordinate: (-0.5,1)
        else if ( rSH->get_subset_index(vrt) == 1 )
        {
            edge[0] = -0.5; edge[1] = 1.0;
            number scalar = VecProd(project, edge);
 
            VecScale(project, edge, scalar);
            ++counter3;
        }
        // (2.4) corner may not moved at all: (0,0)
        else if ( rSH->get_subset_index(vrt) == 5 )
        {
            project[0] = project[1] = 0.0;
            ++counter4;
        }
        else
        {
            ++counter5;
            continue;
        }
        
        /////////////////////////////////////////////////
        // (3) write projected direction back to data 'u':
        
        // REMARK: this part of code is only reached for vertices on the boundary!!!
        //   --> see counter5: continue!
        for (size_t cmp = 0; cmp < dim; ++cmp)
        {
            std::vector<DoFIndex>  vInd;
            if(dd->inner_dof_indices(vrt, cmp, vInd) != 1)
                UG_THROW("1: error in inner_dof_indices operatio!\n");
            
            DoFRef(u, vInd[0]) = project[cmp];
            ++counter6;


        }
        
    } // end vertices-loop

    UG_LOG("counter1 = " << counter1 << "\n");
    UG_LOG("counter2 = " << counter2 << "\n");
    UG_LOG("counter3 = " << counter3 << "\n");
    UG_LOG("counter4 = " << counter4 << "\n");
    UG_LOG("counter5 = " << counter5 << "\n");
    UG_LOG("counter6 = " << counter6 << "\n");

}
    
template <typename TDomain, typename TAlgebra>
number MovingParticle<TDomain, TAlgebra>::
gradient_descent(const size_t n, const number functional, const number bound, const number scaleAlpha, vector_type& u, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, SmartPtr<NavierStokesFV1<TDomain> > spMaster, const int topLevel)
{
    // set dummy value in order to pass second while-loop at least once:
    number functionalNew = functional;
    
    UG_LOG("-----------> functionalNew = " << functionalNew << "\n");

    if ( !(fabs(functionalNew) > bound) )
        UG_THROW("attention: first time this condition neads to be satisfied!\n");
    
    //while ( fabs(functionalNew) > bound )
    if ( 1 )
    {
        number alpha = scaleAlpha*compute_max_step_size(spApproxSpace, spMaster, topLevel);
        UG_LOG("2 -----------> alpha = " << alpha << "\n");
        
        bool iterate = true;
        // while ( iterate )
        if ( 1 )
        {
            /////////////////////////////////////////////////
            // loop vertices and update coordinates:

            ConstSmartPtr<DoFDistribution> dd = spApproxSpace->dof_distribution(GridLevel(topLevel, GridLevel::LEVEL));
            DoFDistribution::traits<Vertex>::const_iterator iterBegin 	= dd->begin<Vertex>();
            DoFDistribution::traits<Vertex>::const_iterator iterEnd 	= dd->end<Vertex>();
 
            //	create MultiindexProjector
            std::vector<DoFIndex> multInd;
        
            typename TDomain::position_accessor_type aaPos = this->m_spApproxSpace->domain()->position_accessor();

            //	loop vertices
            for(DoFDistribution::traits<Vertex>::const_iterator iter = iterBegin; iter != iterEnd; iter++)
            {
                //	get vertex
                Vertex* vrt = *iter;
                MathVector<dim> position = aaPos[vrt];

                for (size_t cmp = 0; cmp < dim; ++cmp)
                {
                    std::vector<DoFIndex>  vInd;
                    if(dd->inner_dof_indices(vrt, cmp, vInd) != 1)
                        UG_THROW("1: error in inner_dof_indices operatio!\n");
                
                    /////////////////////////////////////////////////
                    // (1) get gradient = update direction:
                    MathVector<dim> update;
                    update[cmp] = DoFRef(u, vInd[0]);
            
                    /////////////////////////////////////////////////
                    // (2) update position coordinates due to
                    //     gradient direction and step size alpha:
                    aaPos[vrt][cmp] = position[cmp] + alpha * update[cmp];
                
                    UG_LOG("aaPos[vrt][" << cmp << "]: " << aaPos[vrt][cmp] << "\n");
                    
                }// end cmp-loop
        
            
            } // end vertices-loop
        
            UG_LOG("vorher -----------> functionalNew = " << functionalNew << "\n");

            functionalNew = compute_functional(n, u, spApproxSpace, spMaster, topLevel);
            alpha = 0.5*alpha;

            if ( functionalNew < 0.0 )  iterate = true;
            else                        iterate = false;
            
            UG_LOG("nachher -----------> functionalNew = " << functionalNew << "\n");

            
        }// end while-loop

        // ToDo: wie Verschiebumg von aaPos wieder Rückgängig machen bzw. altes Gitter zwischenspeichern? => auf Clone updaten?
        
    }// end while-loop
    
}

    
}// end namespace MovingParticle
} // end namespace ug



#endif /* MEAN_ID_H_ */
