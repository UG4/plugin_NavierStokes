/*
 * CGradUpwind_constraint.h
 *
 *  Created on: 21.06.2013
 *      Author: Christian Wehner
 */

#ifndef CGRAD_CONSTRAINT_H_
#define CGRAD_CONSTRAINT_H_

#include "lib_disc/spatial_disc/constraints/constraint_interface.h"
#include "lib_disc/spatial_disc/disc_util/fvcr_geom.h"

namespace ug{

template <int dim> struct face_type_traits
{
    typedef void face_type0;
	typedef void face_type1;
};

template <> struct face_type_traits<1>
{
    typedef ReferenceVertex face_type0;
	typedef ReferenceVertex face_type1;
};

template <> struct face_type_traits<2>
{
    typedef ReferenceEdge face_type0;
	typedef ReferenceEdge face_type1;
};

template <> struct face_type_traits<3>
{
    typedef ReferenceTriangle face_type0;
	typedef ReferenceQuadrilateral face_type1;
};

template <typename TGridFunction>
class DiscConstraintFVCR: public IDomainConstraint<typename TGridFunction::domain_type, typename TGridFunction::algebra_type>
{
	public:
		typedef typename TGridFunction::domain_type TDomain;
		typedef typename TGridFunction::algebra_type TAlgebra;

	///	world Dimension
		static const int dim = TDomain::dim;

	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	///	Type of Domain
		typedef TDomain domain_type;

	/// blockSize of used algebra
		static const int blockSize = algebra_type::blockSize;

	///	grid type
		typedef typename domain_type::grid_type grid_type;

	/// element type
		typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object elem_type;

	/// side type
		typedef typename elem_type::side side_type;

	/// element iterator
		typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// side iterator
		typedef typename TGridFunction::template traits<side_type>::const_iterator SideIterator;

		static const size_t _P_ = dim;

	///	Type of geometric base object
		typedef typename domain_traits<TDomain::dim>::geometric_base_object geometric_base_object;

	/// position accessor
		typedef typename domain_type::position_accessor_type position_accessor_type;
	
		static const size_t maxShapeSize = 2*DimCRFVGeometry<dim>::maxNumSCV-1;

		typedef std::vector<std::pair<MultiIndex<2>, MathVector<dim> > > vIndexPosPair;
		typedef std::vector<std::vector<std::pair<MultiIndex<2>, MathVector<dim> > > > vvIndexPosPair;
		typedef std::pair<MathVector<dim>, MathVector<dim> > MathVector_Pair;
		
		typedef MathMatrix<dim,dim> dimMat;
		typedef Attachment<dimMat> AMathDimMat;
		typedef PeriodicAttachmentAccessor<side_type,AMathDimMat > aSideDimMat;
		typedef PeriodicAttachmentAccessor<side_type,ANumber > aSideNumber;
	
		typedef Attachment<std::vector< MathVector<dim> > > ANumberArray;
		typedef Attachment<std::vector< MultiIndex<2> > > ASizetArray;
		typedef PeriodicAttachmentAccessor<side_type,ANumberArray> aSideNumberArray;
		typedef PeriodicAttachmentAccessor<side_type,ASizetArray> aSideSizetArray;
	
		aSideDimMat acGrad;
		aSideNumber acVol;
		aSideNumberArray acGradSh;
		aSideSizetArray acGradShInd;
		
		AMathDimMat aGrad;
		ANumber aVol;
		ANumberArray aGradSh;
		ASizetArray aGradShInd;
		
		typedef typename face_type_traits<dim>::face_type0 face_type0;
		typedef typename face_type_traits<dim>::face_type1 face_type1;

	private:
		SmartPtr<TGridFunction> m_u;
		grid_type* m_grid;
		bool m_bAdaptive;
		bool m_bLinPressureDefect;
		bool m_bLinPressureJacobian;
		bool m_bLinUpConvDefect;
		bool m_bLinUpConvJacobian;
	public:
		void init(SmartPtr<TGridFunction> u,bool bLinUpConvDefect,bool bLinUpConvJacobian,bool bLinPressureDefect,bool bLinPressureJacobian,bool bAdaptive){
			m_u = u;
			domain_type& domain = *m_u->domain().get();
			grid_type& grid = *domain.grid();
			m_grid = &grid;
			m_bLinPressureDefect = bLinPressureDefect;
			m_bLinPressureJacobian = bLinPressureJacobian;
			m_bLinUpConvDefect = bLinUpConvDefect;
			m_bLinUpConvJacobian = bLinUpConvJacobian;
			m_bAdaptive=bAdaptive;
			if (m_bLinUpConvDefect==true){
				// attach
				grid.template attach_to<side_type>(aGrad);
				grid.template attach_to<side_type>(aVol);
				// access
				acGrad.access(grid,aGrad);
				acVol.access(grid,aVol);
			}
			if (m_bLinUpConvJacobian==true){
				grid.template attach_to<side_type>(aGradSh);
				grid.template attach_to<side_type>(aGradShInd);
				acGradSh.access(grid,aGradSh);
				acGradShInd.access(grid,aGradShInd);
			}
			if (m_bLinUpConvJacobian==true) if (bAdaptive==false) compute_grad_shapes();
		}
	/// constructor
		DiscConstraintFVCR(SmartPtr<TGridFunction> u){
			init(u,true,false,true,false,false);
		};

		DiscConstraintFVCR(SmartPtr<TGridFunction> u,bool bLinUpConvDefect,bool bLinUpConvJacobian,bool bLinPressureDefect,bool bLinPressureJacobian,bool bAdaptive){
			init(u,bLinUpConvDefect,bLinUpConvJacobian,bLinPressureDefect,bLinPressureJacobian,bAdaptive);
		};
		
	///	destructor
		~DiscConstraintFVCR() {};

		///	adapts jacobian to enforce constraints
			/// \{
		void compute_grad_shapes(){
			domain_type& domain = *m_u->domain().get();
			
			//	create Multiindex
			std::vector<MultiIndex<2> > multInd;
			
			position_accessor_type aaPos = m_u->domain()->position_accessor();
			
			typename grid_type::template traits<side_type>::secure_container sides;
			std::vector<MathVector<dim> > vCorner;
			std::vector<MultiIndex<2> > ind;
			
			//	create a FV Geometry for the dimension
			DimCRFVGeometry<dim> geo;
			
			SetAttachmentValues(acVol, m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);
			
			// resize shape indices
			for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
			{
				SideIterator iter = m_u->template begin<side_type>(si);
				SideIterator iterEnd = m_u->template end<side_type>(si);
				for(  ;iter !=iterEnd; ++iter)
				{
					side_type* side = *iter;
					acGradSh[side].clear();
				}
			}
			
			// compute gradient in sides
			for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
			{
				ElemIterator iter = m_u->template begin<elem_type>(si);
				ElemIterator iterEnd = m_u->template end<elem_type>(si);
				
				//	loop elements of dimension
				for(  ;iter !=iterEnd; ++iter)
				{
					//	get Elem
					elem_type* elem = *iter;
					
					//  get sides of element
					m_grid->template associated_elements_sorted(sides, elem );
					
					//	reference object type
					//  ReferenceObjectID roid = elem->reference_object_id();
					
					//	get corners of element
					CollectCornerCoordinates(vCorner, *elem, aaPos);
					
					//	evaluate finite volume geometry
					geo.update(elem, &(vCorner[0]), domain.subset_handler().get());
					
					static const size_t maxNumSCVF = DimCRFVGeometry<dim>::maxNumSCVF;
					
					MathMatrix<maxNumSCVF,dim> uValue;
					
					typename grid_type::template traits<side_type>::secure_container sides;
					
					UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Only elements of type elem_type are currently supported");
					
					domain.grid()->associated_elements_sorted(sides, static_cast<elem_type*>(elem) );
					
					size_t nofsides = geo.num_scv();
					
					m_u->multi_indices(elem,0,ind);
					
					//	loop sides
					for (size_t s=0;s < nofsides;s++)
					{
						//	get scv for sh
						const typename DimCRFVGeometry<dim>::SCV& scv = geo.scv(s);
						
						side_type* localside = sides[s];
						
						size_t currentInd = acGradSh[localside].size();
						if (currentInd==0){ 
							currentInd=1;
							acGradSh[localside].resize(nofsides);
							acGradShInd[localside].resize(nofsides);
							acGradSh[localside][0]=0; 
							acGradShInd[localside][0]=ind[s];
						}
						else {
							acGradSh[localside].resize(currentInd+nofsides-1); 
							acGradShInd[localside].resize(currentInd+nofsides-1);
						}
						
						//	volume of scv
						number vol = scv.volume();
						
						for (size_t sh=0;sh < nofsides;sh++)
						{
							if (sh==s){
								for (int d0=0;d0<dim;d0++){
									acGradSh[localside][0][d0] += scv.global_grad(sh)[d0]*vol;
									//UG_LOG("# " << d0 << " " << scv.global_grad(sh)[d0] << "\n");
								}
							} else {
								for (int d0=0;d0<dim;d0++){
									acGradSh[localside][currentInd][d0] = scv.global_grad(sh)[d0]*vol;
									//UG_LOG("# " << d0 << " " << scv.global_grad(sh)[d0] << "\n");
								}
								acGradShInd[localside][currentInd] = ind[sh];
								currentInd++;
							}
						}
						acVol[sides[s]] += vol;
						//UG_LOG("uGrad = " << acUGrad[sides[s]] << " vGrad = " << acVGrad[sides[s]] << " vol = " << acVol[sides[s]] << " div err = " << std::abs(acUGrad[sides[s]][0]+acVGrad[sides[s]][1]) << "\n");
					}
					
				}
			}
			PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
			// complete computation by averaging
			for(int si = 0; si <  domain.subset_handler()->num_subsets(); ++si)
			{
				SideIterator sideIter = m_u->template begin<side_type>(si);
				SideIterator sideIterEnd = m_u->template end<side_type>(si);
				for(  ;sideIter !=sideIterEnd; sideIter++)
				{
					side_type* side = *sideIter;
					if (pbm && pbm->is_slave(side)){
						continue;
					}
					number vol = acVol[side];
					for (size_t i=0;i<acGradSh[side].size();i++){ 
						acGradSh[side][i]/=(number)vol;
						//UG_LOG(i << " (" << acGradSh[side][i][0] << "," << acGradSh[side][i][1] << ")\n");
					}
				}
			}
		}
	
		virtual void adjust_jacobian(matrix_type& J, const vector_type& u,
				                             ConstSmartPtr<DoFDistribution> dd, number time = 0.0,
				                             ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = NULL,const number s_a0 = 1.0){
			if ((m_bLinUpConvJacobian==false)&&(m_bLinPressureJacobian==false)) return;
			if ((m_bAdaptive==true)&&(m_bLinUpConvJacobian==true)) compute_grad_shapes();
			
			domain_type& domain = *m_u->domain().get();

			//	create Multiindex
			std::vector<MultiIndex<2> > multInd;
			
			position_accessor_type aaPos = m_u->domain()->position_accessor();
			
			typename grid_type::template traits<side_type>::secure_container sides;
			std::vector<MathVector<dim> > vCorner;
			std::vector<MultiIndex<2> > ind;
			
			MultiIndex<2> localInd;
			localInd[1]=0;
			MultiIndex<2> shapeInd;
			shapeInd[1]=0;

			//	create a FV Geometry for the dimension
			DimCRFVGeometry<dim> geo;

			for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
			{
				//UG_LOG("si = " << si << "\n");
				//UG_LOG("______________________________________________________________\n");
			//	get iterators
				ElemIterator iter = dd->template begin<elem_type>(si);
				ElemIterator iterEnd = dd->template end<elem_type>(si);

			//	loop elements of dimension
				for(  ;iter !=iterEnd; ++iter)
				{
					//	get Elem
					elem_type* elem = *iter;
					
					MathVector<dim> vGlobalGrad=0;
					
					//  get sides of element
					m_grid->template associated_elements_sorted(sides, elem );
					
					//	get corners of element
					CollectCornerCoordinates(vCorner, *elem, aaPos);
					
					//	evaluate finite volume geometry
					geo.update(elem, &(vCorner[0]), domain.subset_handler().get());
					
					/// handle convection
					if (m_bLinUpConvJacobian==true){
					
						static const size_t maxNumSCVF = DimCRFVGeometry<dim>::maxNumSCVF;
						MathVector<dim> StdVel[maxNumSCVF];
					
						m_u->multi_indices(elem,0,ind);

						for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
						{
							const typename DimCRFVGeometry<dim>::SCVF& scvf = geo.scvf(ip);
							VecSet(StdVel[ip], 0.0);
							for(int d1 = 0; d1 < dim; ++d1){
								dd->multi_indices(elem,0,ind);
								for(size_t sh = 0; sh < sides.size(); ++sh){
									localInd[0]=ind[sh][0]+d1;
									StdVel[ip][d1] += DoFRef(u,localInd) * scvf.shape(sh);
								}
							}
						}
						for(size_t ip = 0; ip < geo.num_scvf(); ++ip){
							const typename DimCRFVGeometry<dim>::SCVF& scvf = geo.scvf(ip);
							//	compute product of standard velocity and normal
							//  const number prod = VecProd(StdVel[ip], scvf.normal()) * m_imDensitySCVF[ip];
							const number prod = s_a0 * VecProd(StdVel[ip], scvf.normal());
							size_t base;
							if (prod>0)
								base = scvf.from();
							else
								base = scvf.to();
							side_type* baseSide = sides[base];
							MathVector<dim> distVec;
							VecSubtract(distVec, scvf.global_ip(),geo.scv(base).global_ip());
							for	(int d0=0;d0<dim;d0++){
								size_t shapeSize = acGradShInd[baseSide].size();
								for (size_t sh=0;sh<shapeSize;sh++){
									for (int d1=0;d1<dim;d1++){
										number flux = prod * distVec[d1] * acGradSh[baseSide][sh][d1];
										localInd[0]=ind[scvf.from()][0]+d0;
										shapeInd[0]=acGradShInd[baseSide][sh][0]+d0;
										DoFRef(J,localInd,shapeInd)+= flux;
										localInd[0]=ind[scvf.to()][0]+d0;
										DoFRef(J,localInd,shapeInd)-= flux;
									}
								}	
							}
						} 
					}
					
					/// handle pressure
					if (m_bLinPressureJacobian==true){
					
						MathVector<dim> vGlobalGrad=0;
					
						//  get sides of element
						m_grid->template associated_elements_sorted(sides, elem );
					
						//	reference object type
						ReferenceObjectID roid = elem->reference_object_id();
					
						//	compute size (volume) of element
						const number elemSize = ElementSize<dim>(roid, &vCorner[0]);
					
						typename grid_type::template traits<elem_type>::secure_container assoElements;
					
						std::vector<MultiIndex<2> > elemInd(sides.size()+1);
					
						dd->inner_multi_indices(elem,_P_,ind);
						elemInd[0] = ind[0];
					
						//UG_LOG("0 " << elemInd[0] << "\n");

						size_t gradShapesSize = 1;
					
						//UG_LOG(elem << "\n");

						MathMatrix<2*dim+1,dim> gradShapes;
						for (int sh=0;sh<2*dim+1;sh++) for (int d0=0;d0<dim;d0++) gradShapes[sh][d0]=0;

						for (size_t s=0;s<sides.size();s++){
							m_grid->template associated_elements(assoElements,sides[s]);
							// face value is average of associated elements
							size_t numOfAsso = assoElements.size();
							const typename DimCRFVGeometry<dim>::SCV& scv = geo.scv(s);
							if (numOfAsso==1){
								for (int d=0;d<dim;d++) gradShapes[0][d]+=scv.normal()[d];
								continue;
							}
							for (size_t i=0;i<numOfAsso;i++){
								dd->inner_multi_indices(assoElements[i],_P_,ind);
								//UG_LOG(assoElements[i] << "\n");
								if (assoElements[i]!=elem) break;
							}
							elemInd[gradShapesSize] = ind[0];
							//UG_LOG(gradShapesSize << " " << elemInd[gradShapesSize] << "\n");
							for (int d=0;d<dim;d++){
								gradShapes[0][d]+=0.5*scv.normal()[d];
								gradShapes[gradShapesSize][d]+=0.5*scv.normal()[d];
							}
							gradShapesSize++;
						}
						gradShapes/=(number)elemSize;
						//UG_LOG("gradShapesSize=" << gradShapesSize << "\n");
						//UG_LOG(gradShapes << "\n");
						// for debug set grad shapes to trivial
						//for (int d2=0;d2<dim;d2++) gradShapes[0][d2]=1;
						//for (size_t sh=0;sh<gradShapesSize;sh++)for (int d2=0;d2<dim;d2++) gradShapes[sh][d2]=0;
					

						for (int d1=0;d1<dim;d1++){
							dd->multi_indices(elem, d1 , multInd);
							for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
							{
								// 	get current SCVF
								const typename DimCRFVGeometry<dim>::SCVF& scvf = geo.scvf(ip);
								MathVector<dim> distVec;
								VecSubtract(distVec,scvf.global_ip(),geo.global_bary());

								for (size_t sh=0;sh<gradShapesSize;sh++){
									for (int d2=0;d2<dim;d2++){
										number flux = s_a0 * distVec[d2]*gradShapes[sh][d2];
										//UG_LOG("from = " << multInd[scvf.from()] << " p = " << elemInd[sh] << "\n");
										DoFRef(J,multInd[scvf.from()],elemInd[sh])+= flux *  scvf.normal()[d1];
										//UG_LOG("to = " << multInd[scvf.to()] << " p = " << elemInd[sh] << "\n");
										DoFRef(J,multInd[scvf.to()],elemInd[sh])-= flux *  scvf.normal()[d1];
									}
								}
							}
						}
					}
				}
			}
		};
			/// \}
			
		virtual void add_pressure_defect(vector_type& d, const vector_type& u,
				                           ConstSmartPtr<DoFDistribution> dd,const number time = 0.0,const number s_a = 1.0){
			domain_type& domain = *m_u->domain().get();

			//	create Multiindex
			std::vector<MultiIndex<2> > multInd;
			
			position_accessor_type aaPos = m_u->domain()->position_accessor();
			
			typename grid_type::template traits<side_type>::secure_container sides;
			std::vector<MathVector<dim> > vCorner;
			std::vector<MultiIndex<2> > ind;

			//	create a FV Geometry for the dimension
			DimCRFVGeometry<dim> geo;

			for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
			{
				ElemIterator iter = dd->template begin<elem_type>(si);
				ElemIterator iterEnd = dd->template end<elem_type>(si);

			//	loop elements of dimension
				for(  ;iter !=iterEnd; ++iter)
				{
					//	get Elem
					elem_type* elem = *iter;
					
					MathVector<dim> vGlobalGrad=0;
					
					//  get sides of element
					m_grid->template associated_elements_sorted(sides, elem );
					
					//	reference object type
					ReferenceObjectID roid = elem->reference_object_id();
					
					//	get corners of element
					CollectCornerCoordinates(vCorner, *elem, aaPos);
					
					//	evaluate finite volume geometry
					geo.update_geometric_data(elem, &(vCorner[0]), domain.subset_handler().get());
					
					//	compute size (volume) of element
					const number elemSize = ElementSize<dim>(roid, &vCorner[0]);
					
					typename grid_type::template traits<elem_type>::secure_container assoElements;
					
					dd->inner_multi_indices(elem,_P_,ind);
					
					number elemValue = DoFRef(u,ind[0]);
					
					MathVector<dim> grad;
					grad *= 0;

					for (size_t s=0;s<sides.size();s++){
						m_grid->template associated_elements(assoElements,sides[s]);
						// face value is average of associated elements
						size_t numOfAsso = assoElements.size();
						const typename DimCRFVGeometry<dim>::SCV& scv = geo.scv(s);
						if (numOfAsso==1){
							for (int d=0;d<dim;d++) grad[d]+=scv.normal()[d]*elemValue;
							continue;
						}
						for (size_t i=0;i<numOfAsso;i++){
							dd->inner_multi_indices(assoElements[i],_P_,ind);
							//UG_LOG(assoElements[i] << "\n");
							if (assoElements[i]!=elem) break;
						}
						for (int d=0;d<dim;d++) grad[d]+=0.5*(elemValue+DoFRef(u,ind[0]))*scv.normal()[d];
					}
					grad/=(number)elemSize;
					for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
					{
						// 	get current SCVF
						const typename DimCRFVGeometry<dim>::SCVF& scvf = geo.scvf(ip);
						MathVector<dim> distVec;
						VecSubtract(distVec,scvf.global_ip(),geo.global_bary());
						number pressureAddition = 0;
						for (int j=0;j<dim;j++)
							pressureAddition+=grad[j]*distVec[j];
						pressureAddition*=s_a;						
														
						for (int d1=0;d1<dim;d1++){
							dd->multi_indices(elem, d1 , multInd);
							//UG_LOG("from = " << multInd[scvf.from()] << " p = " << elemInd[sh] << "\n");
							DoFRef(d,multInd[scvf.from()]) += pressureAddition *  scvf.normal()[d1];
							//UG_LOG("to = " << multInd[scvf.to()] << " p = " << elemInd[sh] << "\n");
							DoFRef(d,multInd[scvf.to()]) -= pressureAddition *  scvf.normal()[d1];
						}
					}
				}
			}
		};
			
		virtual void add_defect(vector_type& d, const vector_type& u,
				                           ConstSmartPtr<DoFDistribution> dd,const number time = 0.0,const number s_a = 1.0){
				domain_type& domain = *m_u->domain().get();

			//	create Multiindex
			std::vector<MultiIndex<2> > multInd;
			
			position_accessor_type aaPos = m_u->domain()->position_accessor();
			
			typename grid_type::template traits<side_type>::secure_container sides;
			std::vector<MathVector<dim> > vCorner;
			std::vector<MultiIndex<2> > ind;

			//	create a FV Geometry for the dimension
			DimCRFVGeometry<dim> geo;
			
			// initialize attachment value
			SetAttachmentValues(acVol, dd->template begin<side_type>(), dd->template end<side_type>(), 0);
			SetAttachmentValues(acGrad, dd->template begin<side_type>(), dd->template end<side_type>(), 0);

			// compute gradient in sides
			for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
			{
				ElemIterator iter = dd->template begin<elem_type>(si);
				ElemIterator iterEnd = dd->template end<elem_type>(si);

			//	loop elements of dimension
				for(  ;iter !=iterEnd; ++iter)
				{
					//	get Elem
					elem_type* elem = *iter;
					
					MathVector<dim> vGlobalGrad=0;
					
					//  get sides of element
					m_grid->template associated_elements_sorted(sides, elem );
					
					//	reference object type
					//  ReferenceObjectID roid = elem->reference_object_id();
					
					//	get corners of element
					CollectCornerCoordinates(vCorner, *elem, aaPos);
					
					//	evaluate finite volume geometry
					geo.update(elem, &(vCorner[0]), domain.subset_handler().get());
					
					static const size_t maxNumSCVF = DimCRFVGeometry<dim>::maxNumSCVF;

					MathMatrix<maxNumSCVF,dim> uValue;

					typename grid_type::template traits<side_type>::secure_container sides;

					UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Only elements of type elem_type are currently supported");

					domain.grid()->associated_elements_sorted(sides, static_cast<elem_type*>(elem) );

					size_t nofsides = geo.num_scv();

					for (int d0=0;d0<dim;d0++){
						for (size_t s=0;s<nofsides;s++){
							dd->inner_multi_indices(sides[s], d0, multInd);
							uValue[s][d0]=DoFRef(u,multInd[0]);
							//UG_LOG("u(" << s << "," << d0 << ") = " << uValue[s][d0] << "\n");
						}
					}

					//	storage for global gradient
					MathMatrix<dim,dim> globalGrad;
					for (int d=0;d<dim;d++) for (int d1=0;d1<dim;d1++) globalGrad[d][d1]=0;

					//	loop sides
					for (size_t s=0;s < nofsides;s++)
					{
						//	get scv for sh
						const typename DimCRFVGeometry<dim>::SCV& scv = geo.scv(s);

						//	reset global gradient
						for (int d0=0;d0<dim;d0++)
							for (int d1=0;d1<dim;d1++)
								globalGrad[d0][d1] = 0.0;

						for (int d0=0;d0<dim;d0++)
							//	sum up gradients of shape functions in side
							for(size_t sh = 0 ; sh < nofsides; ++sh)
							{
								for (int d1=0;d1<dim;d1++) globalGrad[d0][d1]+=uValue[sh][d0] * scv.global_grad(sh)[d1];
							//	if (d==1) UG_LOG(" " << scv.global_grad(sh) << "\n");
								//if (d==1) UG_LOG(globalGrad[d] << "\n");
								//if (d==1) UG_LOG("uValue(" << sh << "," << d << ") = " << uValue[sh][d] << "\n");
							}

						//	volume of scv
						number vol = scv.volume();

						//	scale gradient by volume
						globalGrad *= vol;

						//	add both values to attachements
						acGrad[sides[s]] += globalGrad;
						acVol[sides[s]] += vol;
						//UG_LOG("uGrad = " << acUGrad[sides[s]] << " vGrad = " << acVGrad[sides[s]] << " vol = " << acVol[sides[s]] << " div err = " << std::abs(acUGrad[sides[s]][0]+acVGrad[sides[s]][1]) << "\n");
					}

				}
			}
			PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
			// complete computation by averaging
			for(int si = 0; si <  domain.subset_handler()->num_subsets(); ++si)
			{
				SideIterator sideIter = dd->template begin<side_type>(si);
				SideIterator sideIterEnd = dd->template end<side_type>(si);
				for(  ;sideIter !=sideIterEnd; sideIter++)
				{
					side_type* side = *sideIter;
					if (pbm && pbm->is_slave(side)){
						continue;
					}
					acGrad[side]/=(number)acVol[side];
					//UG_LOG(acGrad[side] << "\n");
				}
			}
			// compute defect
			for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
			{
				ElemIterator iter = dd->template begin<elem_type>(si);
				ElemIterator iterEnd = dd->template end<elem_type>(si);

			//	loop elements of dimension
				for(  ;iter !=iterEnd; ++iter)
				{
					//	get Elem
					elem_type* elem = *iter;
					
					MathVector<dim> vGlobalGrad=0;
					
					//  get sides of element
					m_grid->template associated_elements_sorted(sides, elem );
					
					//	reference object type
					//  ReferenceObjectID roid = elem->reference_object_id();
					
					//	get corners of element
					CollectCornerCoordinates(vCorner, *elem, aaPos);
					
					//	evaluate finite volume geometry
					geo.update(elem, &(vCorner[0]), domain.subset_handler().get());
					
					// static const size_t MaxNumSidesOfElem = 10;

					typename grid_type::template traits<side_type>::secure_container sides;

					UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Only elements of type elem_type are currently supported");

					domain.grid()->associated_elements_sorted(sides, static_cast<elem_type*>(elem) );

					// size_t nofsides = geo.num_scv();
					
					static const size_t maxNumSCVF = DimCRFVGeometry<dim>::maxNumSCVF;
					MathVector<dim> StdVel[maxNumSCVF];

					for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
					{
						const typename DimCRFVGeometry<dim>::SCVF& scvf = geo.scvf(ip);
						VecSet(StdVel[ip], 0.0);

						for(size_t sh = 0; sh < sides.size(); ++sh){
							for(int d1 = 0; d1 < dim; ++d1){
								dd->inner_multi_indices(sides[sh], d1, multInd);
								StdVel[ip][d1] += DoFRef(u,multInd[0]) * scvf.shape(sh);
							}
						}
					}
					
					number elemPressureValue=0, pressure=0;
					MathVector<dim> pGrad;
					
					//// PRESSURE GRADIENT IN ELEMENT
					if (m_bLinPressureDefect==true){
						//	reference object type
						ReferenceObjectID roid = elem->reference_object_id();
					
						//	compute size (volume) of element
						const number elemSize = ElementSize<dim>(roid, &vCorner[0]);
					
						typename grid_type::template traits<elem_type>::secure_container assoElements;
					
						dd->inner_multi_indices(elem,dim,multInd);
						elemPressureValue = DoFRef(u,multInd[0]);
					
						for (int i=0;i<dim;i++) pGrad[i]=0;
					
						for (size_t s=0;s<sides.size();s++){
							m_grid->template associated_elements(assoElements,sides[s]);
							// face value is average of associated elements
							size_t numOfAsso = assoElements.size();
							const typename DimCRFVGeometry<dim>::SCV& scv = geo.scv(s);
							if (numOfAsso==1){
								for (int d=0;d<dim;d++) pGrad[d]+=scv.normal()[d]*elemPressureValue;
								continue;
							}
							for (size_t i=0;i<numOfAsso;i++){
								dd->inner_multi_indices(assoElements[i],_P_,ind);
								//UG_LOG(assoElements[i] << "\n");
								if (assoElements[i]!=elem) break;
							}
							for (int d=0;d<dim;d++) pGrad[d]+=0.5*(elemPressureValue+DoFRef(u,ind[0]))*scv.normal()[d];
						}
						pGrad/=(number)elemSize;
					}
					
					size_t base;
					MathVector<dim> distVec;
					for(size_t ip = 0; ip < geo.num_scvf(); ++ip){
						const typename DimCRFVGeometry<dim>::SCVF& scvf = geo.scvf(ip);
						const number flux = s_a * VecProd(StdVel[ip], scvf.normal());
						
						if (flux>0)
							base = scvf.from();
						else
							base = scvf.to();
						
						//// PRESSURE FLUX
						if (m_bLinPressureDefect==true){
							MathVector<dim> distVecBary;
							VecSubtract(distVecBary,scvf.global_ip(),geo.global_bary());
							number pressureFlux = 0;// 
							// pressureFlux = elemPressureValue;
							for (int j=0;j<dim;j++)
								pressureFlux+=pGrad[j]*distVecBary[j];
							pressure = s_a * pressureFlux;
						}
						
						for(int d1 = 0; d1 < dim; ++d1)
						{
//							dd->inner_multi_indices(sides[base],d1,multInd);
//							number upwindVel = DoFRef(u,multInd[0]);
							number upwindVel = 0;
							VecSubtract(distVec, scvf.global_ip(),geo.scv(base).global_ip());
					//		UG_LOG("grad = " << acGrad[sides[base]][d1][0] << "," << acGrad[sides[base]][d1][1] << "\n");
					//		UG_LOG("distVec = " << distVec << "\n");
							for (int d2=0;d2<dim;d2++){
								upwindVel+=(acGrad[sides[base]])[d1][d2]*distVec[d2];
							}
						//	UG_LOG("upv(" << ip+1 << "," << d1+1 << ")=" << upwindVel << "\n");
							dd->multi_indices(elem,d1,multInd);
							DoFRef(d,multInd[scvf.from()]) += upwindVel * flux;
							DoFRef(d,multInd[scvf.to()]) -= upwindVel * flux;
							
							//// PRESSURE ADDITION
							if (m_bLinPressureDefect==true){
								DoFRef(d,multInd[scvf.from()]) += pressure * scvf.normal()[d1];
								DoFRef(d,multInd[scvf.to()]) -= pressure * scvf.normal()[d1];
							}
						}
					}
				}
			}

		}

			///	adapts defect to enforce constraints
			/// \{
		virtual void adjust_defect(vector_type& d, const vector_type& u,
				                           ConstSmartPtr<DoFDistribution> dd, number time = 0.0,
				                           ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = NULL,
										   const std::vector<number>* vScaleMass = NULL,
										   const std::vector<number>* vScaleStiff = NULL)
		{
			if (vSol == NULL)
				if (m_bLinUpConvDefect==false)
					add_pressure_defect(d,u,dd);
				else
					add_defect(d,u,dd);
			else {
				//	loop all time points and assemble them
				for(size_t t = 0; t < vScaleStiff->size(); ++t){
					if ((*vScaleStiff)[t]==0) continue;
					if (m_bLinUpConvDefect==false) 
						add_pressure_defect(d,*(vSol->solution(t)),dd,time,(*vScaleStiff)[t]);
					else 
						add_defect(d,*(vSol->solution(t)),dd,time,(*vScaleStiff)[t]);
				}
			}
		};
			/// \}

			///	adapts matrix and rhs (linear case) to enforce constraints
			/// \{
		virtual void adjust_linear(matrix_type& mat, vector_type& rhs,
				                           ConstSmartPtr<DoFDistribution> dd, number time = 0.0){};
			/// \}

			///	adapts a rhs to enforce constraints
			/// \{
		virtual void adjust_rhs(vector_type& rhs, const vector_type& u,
				                        ConstSmartPtr<DoFDistribution> dd, number time = 0.0){};
			/// \}

			///	sets the constraints in a solution vector
			/// \{
		virtual void adjust_solution(vector_type& u, ConstSmartPtr<DoFDistribution> dd,
				                             number time = 0.0){};

	///	returns the type of the constraints 
		virtual int type() const {return CT_CONSTRAINTS;}
};

} // end namespace ug

// #include "CGradUpwind_constraint_impl.h"

#endif /* CGRAD_CONSTRAINT_H_ */
