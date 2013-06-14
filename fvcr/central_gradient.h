/*	
 *  central_gradient.h
 *
 *  Created on: 18.04.2013
 *      Author: Christian Wehner
 */

#ifndef __H__UG__NAVIER_STOKES_CENTRAL_GRADIENT__
#define __H__UG__NAVIER_STOKES_CENTRAL_GRADIENT__

#include "common/common.h"

#include "lib_disc/common/subset_group.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/local_finite_element/local_shape_function_set.h"
#include "lib_disc/spatial_disc/disc_util/fvcr_geom.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "pressure_separation.h"
#include "lib_grid/tools/periodic_boundary_manager.h"
#include "lib_grid/algorithms/attachment_util.h"

#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{
namespace NavierStokes{

template <typename TGridFunction>
class CentralGradient
	: 	public StdUserData<CentralGradient<TGridFunction>, MathMatrix<TGridFunction::dim,TGridFunction::dim>, TGridFunction::dim>,
	virtual public INewtonUpdate
{
		
		static const int dim = TGridFunction::domain_type::dim;

		///	domain type
		typedef typename TGridFunction::domain_type domain_type;

		/// position accessor type
		typedef typename domain_type::position_accessor_type position_accessor_type;

		///	algebra type
		typedef typename TGridFunction::algebra_type algebra_type;

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

		typedef MathVector<dim> vecDim;
		typedef Attachment<vecDim> AMathVectorDim;
		typedef PeriodicAttachmentAccessor<side_type,AMathVectorDim > aSideDimVector;
		typedef PeriodicAttachmentAccessor<side_type,ANumber > aSideNumber;

		aSideDimVector acUGrad;
		aSideDimVector acVGrad;
		aSideDimVector acWGrad;
		aSideDimVector acLaplacian;
		aSideNumber acVol;

		AMathVectorDim aUGrad;
		AMathVectorDim aVGrad;
		AMathVectorDim aWGrad;
		AMathVectorDim aLaplacian;
		ANumber aVol;

		SmartPtr<TGridFunction> m_u;
		grid_type* m_grid;
	public:
		// constructor
		CentralGradient(SmartPtr<TGridFunction> u){
				m_u = u;

				domain_type& domain = *m_u->domain().get();
				grid_type& grid = *domain.grid();

				m_grid = &grid;

				// attach
				grid.template attach_to<side_type>(aUGrad);
				grid.template attach_to<side_type>(aVGrad);
				grid.template attach_to<side_type>(aLaplacian);

				if (dim==3) grid.template attach_to<side_type>(aWGrad);
				grid.template attach_to<side_type>(aVol);
				// access
				acUGrad.access(grid,aUGrad);
				acVGrad.access(grid,aVGrad);
				if (dim==3) acWGrad.access(grid,aWGrad);
				acLaplacian.access(grid,aLaplacian);
				acVol.access(grid,aVol);
		}
	// format of vValue
	// row 0 gradient u, row 1 gradient v, row 2 gradient w
	template <int refDim>
	inline void evaluate(MathMatrix<dim,dim> vValue[],
	                     const MathVector<dim> vGlobIP[],
	                     number time, int si,
	                     GeometricObject* elem,
	                     const MathVector<dim> vCornerCoords[],
	                     const MathVector<refDim> vLocIP[],
	                     const size_t nip,
	                     LocalVector* u,
	                     const MathMatrix<refDim, dim>* vJT = NULL) const
	{
		
		typename grid_type::template traits<side_type>::secure_container sides;
		
		UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Unsupported element type");
		
		m_grid->associated_elements_sorted(sides, static_cast<elem_type*>(elem) );
		
		for (size_t i=0;i<sides.size();i++){
			vValue[i].assign(acUGrad[sides[i]], 0);
			vValue[i].assign(acVGrad[sides[i]], 1);
  			if (dim==3) vValue[i].assign(acWGrad[sides[i]], dim-1);
		}
	}
	
		// compute p gradient in side by going through neighbour elements and interpolate
		inline void compute_p_grad(MathVector<dim>& grad,side_type* side,MathVector<dim> sideCo,size_t order,const position_accessor_type& posAcc){
			size_t depth=order+2;
		       number nodefactor=(number)2.0;
		       MathVector<dim> baryCo;
		       elem_type* element;
		       typename grid_type::template traits<elem_type>::secure_container assoElements;
		       m_grid->associated_elements(assoElements,side);
		       element =  assoElements[0];
		       computeElemBarycenter<elem_type,position_accessor_type,dim>(baryCo,element,posAcc);
		       std::vector<elem_type*> nbrs;
		       std::vector<elem_type*> nbrCandidates;
		       std::vector<size_t> stageStart(depth+2);
		       std::vector<MathVector<dim> > coord;
		       std::vector<number> values;
		       nbrs.push_back(element);
		       stageStart[0]=0;
		       std::vector<MultiIndex<2> > multInd;
		       for (size_t i=1;i<depth+2;i++){
		    	   stageStart[i]=1;
		       }
		       for (size_t stage=0;stage<depth;stage++){
		    	   for (size_t i=stageStart[stage];i<stageStart[stage+1];i++){
		    		   CollectNeighbors(nbrCandidates,nbrs[i], *m_grid,NHT_VERTEX_NEIGHBORS);
		               for (size_t j=0;j<nbrCandidates.size();j++){
		            	   bool newNeighbor=true;
		                   for (size_t k=0;k<nbrs.size();k++){
		                	   if (nbrCandidates[j]==nbrs[k]){
		                		   newNeighbor=false;
		                           break;
		                        }
		                    };
		                    if (newNeighbor==true){
		                    	nbrs.push_back(nbrCandidates[j]);
		                    };
		                };
		            };
		            stageStart[stage+2]=nbrs.size();
		         };
		         size_t nOfPoints = nbrs.size();
		         coord.resize(nOfPoints);
		         std::vector<number> distToBaseP(nOfPoints);
		         std::vector<int> sortedList(nOfPoints);
		         values.resize(nOfPoints);
		         for (size_t i=0;i<nOfPoints;i++){
		        	 computeElemBarycenter<elem_type,position_accessor_type,dim>(coord[i],nbrs[i],posAcc);
		             m_u->inner_multi_indices(nbrs[i], _P_, multInd);
		             values[i]=DoFRef(*m_u,multInd[0]);
		             distToBaseP[i]=VecDistance(coord[i],sideCo);
		         }
		        bubblesort(sortedList,distToBaseP);
		    		size_t ord = order;
		    		std::vector<number> coeffs;
		    		while (ord>=0){
		    			size_t vlength=(size_t)round(0.5*(ord+1)*(ord+2));
						//UG_LOG("vlength=" << vlength << "\n");
						//UG_LOG("ord=" << ord << "\n");
		    			size_t nrOfInterPoints = (size_t)round(vlength*nodefactor);
						if (nrOfInterPoints>nOfPoints){
							//UG_LOG("not enough points for desired order. reduce order\n");
							ord--;
							continue;
						}
						//UG_LOG("nrOfInterPoints=" << nrOfInterPoints << "\n");
		    			std::vector<number> interM(vlength*nrOfInterPoints);
		    			coeffs.resize(vlength);
		    			std::vector<number> interRhs(nrOfInterPoints);
		    			size_t matindex=0;
		    			for (size_t j=0;j<nrOfInterPoints;j++){
		    				size_t i=sortedList[j];
							//UG_LOG(i << "\n");
		    				// interpolation rhs
		    				interRhs[j] = values[i];
							//UG_LOG(interRhs[j] << "\n");
		    				// interpolation matrix
		    			    for (size_t ii=0;ii<=ord;ii++){
		    			       	for (size_t k=0;k<=ii;k++){
		    			       		interM[matindex]=std::pow((number)coord[i][0],(int)(ii-k))*std::pow((number)coord[i][1],(int)k);
		    			       		matindex++;
		    			       	};
		    			    };
		    			};
		    			if (leastSquares(coeffs,interM,interRhs)==true) break;
		    			ord--;
		    			//UG_LOG("Least squares problem had no regular solution. Reduce order to " << ord << ".\n");
		    			if (ord==0){
		    				//UG_LOG("Feasible set of interpolation nodes not found. Set gradient to zero.\n");
		    				break;
		    			}
		    		}
		    		grad=0.0;
					//UG_LOG("*\n");
				 MathVector<dim> locCo;
				 locCo[0]=sideCo[0];
		         locCo[1]=sideCo[1];
		         // avoid 0^0
				 if (locCo[0]==0) locCo[0]=1e-15;
				 if (locCo[1]==0) locCo[1]=1e-15;
	         switch (ord)
		         {
		           	 case 0:
		           		 grad[0]=0;
		           		 grad[1]=0;
		             	 break;
		             case 1:
		              	 grad[0]=coeffs[1];
		                 grad[1]=coeffs[2];
		                 break;
		             case 2:
		              	 grad[0]=coeffs[1]+2*coeffs[3]*locCo[0]+coeffs[4]*locCo[1];
		               	 grad[1]=coeffs[2]+coeffs[4]*locCo[0]+2*coeffs[5]*locCo[1];
		                 break;
		             case 3:
		               	 grad[0]=coeffs[1]+2*coeffs[3]*locCo[0]+coeffs[4]*locCo[1]+3*coeffs[6]*locCo[0]*locCo[0]+2*coeffs[7]*locCo[0]*locCo[1]+coeffs[8]*locCo[1]*locCo[1];
		               	 grad[1]=coeffs[2]+coeffs[4]*locCo[0]+2*coeffs[5]*locCo[1]+coeffs[7]*locCo[0]*locCo[0]+2*coeffs[8]*locCo[0]*locCo[1]+3*coeffs[9]*locCo[1]*locCo[1];
		                 break;
		             default:
		            	 size_t count = 0;
		                 for (size_t i=1;i<=ord;i++)
		                  	 for (size_t k=0;k<=i;k++){
		                   		 grad[0]+=coeffs[count]*(i-k)*std::pow((number)locCo[0],(int)(i-k-1))*std::pow((number)locCo[1],(int)k);
		                   		 grad[1]+=coeffs[count]*k*std::pow((number)locCo[0],(int)(i-k))*std::pow((number)locCo[1],(int)(k-1));
		                   		 count++;
		                     }
		         }
		}

		void update(){
			SmartPtr<TGridFunction> u = m_u;

			domain_type& domain = *u->domain().get();

			//	create Multiindex
			std::vector<MultiIndex<2> > multInd;

			//	create a FV Geometry for the dimension
			DimCRFVGeometry<dim> geo;

			//	get element iterator type
			typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;
			typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object ElemType;

			// initialize attachment value
			SetAttachmentValues(acVol, u->template begin<side_type>(), u->template end<side_type>(), 0);
			SetAttachmentValues(acUGrad, u->template begin<side_type>(), u->template end<side_type>(), 0);
			SetAttachmentValues(acVGrad, u->template begin<side_type>(), u->template end<side_type>(), 0);
			if (dim==3) SetAttachmentValues(acWGrad, u->template begin<side_type>(), u->template end<side_type>(), 0);
			SetAttachmentValues(acLaplacian, u->template begin<side_type>(), u->template end<side_type>(), 0);

			//	coord and vertex array
			MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
			VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

			//	sum up all contributions of the sub control volumes to one vertex in an attachment
			for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
			{
			//	get iterators
				ElemIterator iter = u->template begin<elem_type>(si);
				ElemIterator iterEnd = u->template end<elem_type>(si);

			//	loop elements of dimension
				for(  ;iter !=iterEnd; ++iter)
				{
				//	get Elem
					elem_type* elem = *iter;

				//	get position accessor
					typedef typename domain_type::position_accessor_type position_accessor_type;
					const position_accessor_type& aaPos = domain.position_accessor();

				//	get vertices and extract corner coordinates
					const size_t numVertices = elem->num_vertices();
					for(size_t i = 0; i < numVertices; ++i){
						vVrt[i] = elem->vertex(i);
						coCoord[i] = aaPos[vVrt[i]];
						// UG_LOG("co = " << coCoord[i] << "\n");
					};

				//	evaluate finite volume geometry
					geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

					static const size_t MaxNumSidesOfElem = 10;

					typedef MathVector<dim> MVD;
					std::vector<MVD> uValue(MaxNumSidesOfElem);

					typename grid_type::template traits<side_type>::secure_container sides;

					UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Only elements of type elem_type are currently supported");

					domain.grid()->associated_elements_sorted(sides, static_cast<elem_type*>(elem) );

					size_t nofsides = geo.num_scv();

					for (size_t s=0;s<nofsides;s++){
						for (int d=0;d<dim;d++){
							u->multi_indices(sides[s], d, multInd);
							uValue[s][d]=DoFRef(*u,multInd[0]);
							// UG_LOG("u(" << s << "," << d << ") = " << uValue[s][d] << "\n");
						}
					}

					//	storage for global gradient
					MVD globalGrad[3];

					//	loop sides
					for (size_t s=0;s < nofsides;s++)
					{
						//	get scv for sh
						const typename DimCRFVGeometry<dim>::SCV& scv = geo.scv(s);

						//	reset global gradient
						for (int d=0;d<dim;d++)
							for (int d1=0;d1<dim;d1++)
								globalGrad[d][d1] = 0.0;

						for (int d=0;d<dim;d++)
							//	sum up gradients of shape functions in side
							for(size_t sh = 0 ; sh < nofsides; ++sh)
							{
							//	for (int d1=0;d1<dim;d1++) globalGrad[d][d1]+=uValue[sh][d1] * scv.global_grad(sh)[d1];
								VecScaleAppend(globalGrad[d], uValue[sh][d], scv.global_grad(sh));
							//	if (d==1) UG_LOG(" " << scv.global_grad(sh) << "\n");
								//if (d==1) UG_LOG(globalGrad[d] << "\n");
								//if (d==1) UG_LOG("uValue(" << sh << "," << d << ") = " << uValue[sh][d] << "\n");
							}

						//	volume of scv
						number vol = scv.volume();

						//	scale gradient by volume
						for (int d=0;d<dim;d++)
							globalGrad[d] *= vol;

						//	add both values to attachements
						acUGrad[sides[s]] += globalGrad[0];
						acVGrad[sides[s]] += globalGrad[1];
						if (dim==3) acWGrad[sides[s]] += globalGrad[2];
						acVol[sides[s]] += vol;
						//UG_LOG("uGrad = " << acUGrad[sides[s]] << " vGrad = " << acVGrad[sides[s]] << " vol = " << acVol[sides[s]] << " div err = " << std::abs(acUGrad[sides[s]][0]+acVGrad[sides[s]][1]) << "\n");
					}
					
					
					// assemble laplacian
					for(size_t ip = 0; ip < geo.num_scvf(); ++ip){
						// 	get current SCVF
						const typename DimCRFVGeometry<dim>::SCVF& scvf = geo.scvf(ip);
						for (int d=0;d<dim;d++){
							MathVector<dim> Dgrad_c, grad_c;
							VecSet(grad_c, 0.0);
							for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
								VecScaleAppend(grad_c, uValue[sh][d], scvf.global_grad(sh));
							Dgrad_c = grad_c;
							const number diff_flux = VecDot(Dgrad_c, scvf.normal());
//							UG_LOG("flux = " << diff_flux << "\n");
							acLaplacian[sides[scvf.from()]][d] += diff_flux;
							acLaplacian[sides[scvf.to()]][d] -= diff_flux;
//							UG_LOG(acLaplacian[sides[scvf.from()]] << "\n");
//							UG_LOG(acLaplacian[sides[scvf.to()]] << "\n");
						}
					}
				}
			}
			PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
			// complete computation by averaging
			// for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
			// boundary subsets change later
			for(int si = 1; si <  domain.subset_handler()->num_subsets(); ++si)
			{
				SideIterator sideIter = u->template begin<side_type>(si);
				SideIterator sideIterEnd = u->template end<side_type>(si);
				for(  ;sideIter !=sideIterEnd; sideIter++)
				{
					side_type* side = *sideIter;
					if (pbm && pbm->is_slave(side)){
						continue;
					}
					const number volume = acVol[side];
					acUGrad[side]/=(number)volume;
					acVGrad[side]/=(number)volume;
					if (dim==3) acWGrad[side]/=(number)volume;
				}
			}
			// inner subsets change later
			for(int si = 0; si <  1; ++si)
			{
				SideIterator sideIter = u->template begin<side_type>(si);
				SideIterator sideIterEnd = u->template end<side_type>(si);
				for(  ;sideIter !=sideIterEnd; sideIter++)
				{
					side_type* side = *sideIter;
					if (pbm && pbm->is_slave(side)){
						continue;
					}
					const number volume = acVol[side];
					acUGrad[side]/=(number)volume;
					acVGrad[side]/=(number)volume;
					if (dim==3) acWGrad[side]/=(number)volume;
					continue;
					number unorm=0;
					for (int d=0;d<dim;d++){
						u->multi_indices(side, d, multInd);
						number value = DoFRef(*u,multInd[0]);
						unorm+=value*value;
					}
					if (unorm<1e-10) continue;
				//	UG_LOG("uGrad = " << acUGrad[side] << " vGrad = " << acVGrad[side] << " vol = " << acVol[side] << " div err = " << std::abs(acUGrad[side][0]+acVGrad[side][1]) << "\n");
				//	continue;
					acLaplacian[side]/=(number)volume;
					//UG_LOG("uGrad = " << acUGrad[side] << " vGrad = " << acVGrad[side] << " vol = " << acVol[side] << " div err = " << std::abs(acUGrad[side][0]-acVGrad[side][1]) << "\n");
					typedef typename domain_type::position_accessor_type position_accessor_type;
					const position_accessor_type& aaPos = domain.position_accessor();
					//UG_LOG("laplace = " << acLaplacian[side] << "\n");
					
					MathVector<dim> coord[1];
					MathVector<dim> source[1];
					MathVector<dim> viscosity[1];
					coord[0]*=0;
					for (size_t i=0;i<side->num_vertices();i++){
						coord[0] += aaPos[side->vertex(i)];
					}
					coord[0]/=(number)side->num_vertices();
					// UG_LOG("coord = " << coord[0] << "\n");
					number m_time = 0;
					(*m_imSource)(&source[0], coord, m_time, si, 1);
					(*m_imSource)(&viscosity[0], coord, m_time, si, 1);

					// computation of gradient
					static const size_t N = dim*(dim+1)+1;
					DenseMatrix< FixedArray2<number, N, N> > mat;
					for (size_t i=0;i<N;i++)
						for (size_t j=0;j<N;j++)
							mat(i,j) = 0;

					for (size_t i=0;i<dim*dim;i++){
						mat(i,i)=1;
					}
					for (int i=0;i<dim;i++){
						for (int d=0;d<dim;d++){
							u->multi_indices(side, d, multInd);
							number value = DoFRef(*u,multInd[0]);
							mat(dim*dim+i,2*i+d) = value;
							mat(2*i+d,dim*dim+i) = value;
						}
					}
					for (int i=0;i<dim;i++){
						mat(N-1,i+i*dim) = 1;
						mat(i+i*dim,N-1) = 1;
					}
					DenseVector< FixedArray1<number, N> > rhs,sol;
					for (int i=0;i<dim;i++)	rhs[i] = acUGrad[side][i];
					for (int i=0;i<dim;i++) rhs[dim+i] = acVGrad[side][i];
					if (dim==3) for (int i=0;i<dim;i++) rhs[2*dim+i] = acWGrad[side][i];
					
					MathVector<dim> pgrad;		
					compute_p_grad(pgrad,side,coord[0],2,aaPos);
//UG_LOG(pgrad << "\n");
					for (int d=0;d<dim;d++)
						rhs[dim*dim+d] = source[0][d] + viscosity[0]*acLaplacian[side][d] - pgrad[d];
					rhs[N-1] = 0;
//UG_LOG(mat << "\n");
//UG_LOG(rhs << "\n");
					InverseMatMult(sol,1,mat,rhs);
					number alpha=0;
					for (int i=0;i<dim;i++) acUGrad[side][i] = alpha*acUGrad[side][i] + (1-alpha)*sol[i];
					for (int i=0;i<dim;i++) acVGrad[side][i] = alpha*acVGrad[side][i] + (1-alpha)*sol[dim+i];
					if (dim==3) for (int i=0;i<dim;i++) acWGrad[side][i] = alpha*acWGrad[side][i] + (1-alpha)*sol[2*dim+i];
					//UG_LOG("uGrad = " << acUGrad[side] << " vGrad = " << acVGrad[side] << " pgrad = " << pgrad << " vol = " << acVol[side] << " div err = " << std::abs(acUGrad[side][0]+acVGrad[side][1]) << "\n");
					//acUGrad[side][0]=0;acUGrad[side][1]=2;
					//acVGrad[side][0]=1;acVGrad[side][1]=0;
					//UG_LOG("--------------\n");
				}
			}
			MathVector<dim> q;
			q[4]=9;
		}
		
		public:
	virtual void operator() (MathMatrix<dim,dim>& value,
	                         const MathVector<dim>& globIP,
	                         number time, int si) const
	{
		UG_THROW("LevelSetUserData: Need element.");
	}

	virtual void operator() (MathMatrix<dim,dim> vValue[],
	                         const MathVector<dim> vGlobIP[],
	                         number time, int si, const size_t nip) const
	{
		UG_THROW("LevelSetUserData: Need element.");
	}

	virtual void compute(LocalVector* u, GeometricObject* elem,
	                     const MathVector<dim> vCornerCoords[], bool bDeriv = false)
	{
		const number t = this->time();
		const int si = this->subset();
		for(size_t s = 0; s < this->num_series(); ++s)
			evaluate<dim>(this->values(s), this->ips(s), t, si,
			              elem, NULL, this->template local_ips<dim>(s),
			              this->num_ip(s), u);
	}

	///	returns if provided data is continuous over geometric object boundaries
	virtual bool continuous() const {return false;}

	///	returns if grid function is needed for evaluation
	virtual bool requires_grid_fct() const {return true;}
	
	private:
	///	Data import for source
	SmartPtr<CplUserData<MathVector<dim>,dim> > m_imSource;
	
	///	Data import for kinematic viscosity
	SmartPtr<CplUserData<number,dim> > m_imKinViscosity;
	
public:
	/////////// Source
	
	void set_source(SmartPtr<CplUserData<MathVector<dim>, dim> > data)
	{
		m_imSource = data;
	}
	
	void set_source(number f_x)
	{
		SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
		for (int i=0;i<dim;i++){
			f->set_entry(i, f_x);
		}
		set_source(f);
	}
	
	void set_source(number f_x, number f_y)
	{
		if (dim!=2){
			UG_THROW("NavierStokes: Setting source vector of dimension 2"
					 " to a Discretization for world dim " << dim);
		} else {
			SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
			f->set_entry(0, f_x);
			f->set_entry(1, f_y);
			set_source(f);
		}
	}
	
	void set_source(number f_x, number f_y, number f_z)
	{
		if (dim<3){
			UG_THROW("NavierStokes: Setting source vector of dimension 3"
					 " to a Discretization for world dim " << dim);
		}
		else
		{
			SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
			f->set_entry(0, f_x);
			f->set_entry(1, f_y);
			f->set_entry(2, f_z);
			set_source(f);
		}
	}
	
	#ifdef UG_FOR_LUA
	void set_source(const char* fctName)
	{
		set_source(LuaUserDataFactory<MathVector<dim>, dim>::create(fctName));
	}
	#endif
	
	void set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > user){
		m_imKinViscosity = user;
	}
	void set_kinematic_viscosity(number val){
		set_kinematic_viscosity(CreateSmartPtr(new ConstUserNumber<dim>(val)));
	}
	
	#ifdef UG_FOR_LUA
	void set_kinematic_viscosity(const char* fctName){
		set_kinematic_viscosity(LuaUserDataFactory<number, dim>::create(fctName));
	}
	#endif
};
	
} // namespace NavierStokes
} // end namespace ug


#endif /* __H__UG__NAVIER_STOKES_CENTRAL_GRADIENT__ */
