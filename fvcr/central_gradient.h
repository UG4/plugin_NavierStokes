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

		///	algebra type
		typedef typename TGridFunction::algebra_type algebra_type;

		/// position accessor type
		typedef typename domain_type::position_accessor_type position_accessor_type;

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

		typedef MathVector<dim> vecDim;
		typedef Attachment<vecDim> AMathVectorDim;
		typedef PeriodicAttachmentAccessor<side_type,AMathVectorDim > aSideDimVector;
		typedef PeriodicAttachmentAccessor<side_type,ANumber > aSideNumber;

		aSideDimVector acUGrad;
		aSideDimVector acVGrad;
		aSideDimVector acWGrad;
		aSideNumber acVol;

		AMathVectorDim aUGrad;
		AMathVectorDim aVGrad;
		AMathVectorDim aWGrad;
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

				if (dim==3) grid.template attach_to<side_type>(aWGrad);
				grid.template attach_to<side_type>(aVol);
				// access
				acUGrad.access(grid,aUGrad);
				acVGrad.access(grid,aVGrad);
				if (dim==3) acWGrad.access(grid,aWGrad);
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
							globalGrad[d] *= 0.0;

						for (int d=0;d<dim;d++)
							//	sum up gradients of shape functions in side
							for(size_t sh = 0 ; sh < nofsides; ++sh)
							{
								VecScaleAppend(globalGrad[d], uValue[sh][d], scv.global_grad(sh));
							}

						//	volume of scv
						number vol = scv.volume();

						////UG_LOG("*** global grad " << i << ": " << globalGrad << "\n");

						//	scale gradient by volume
						for (int d=0;d<dim;d++)
							globalGrad[d] *= vol;

						//	add both values to attachements
						acUGrad[sides[s]] += globalGrad[0];
						acVGrad[sides[s]] += globalGrad[1];
						if (dim==3) acWGrad[sides[s]] += globalGrad[2];
						acVol[sides[s]] += vol;
					}
				}
			}
			PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
			// complete computation by averaging
			for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
			{
				SideIterator sideIter = u->template begin<side_type>(si);
				SideIterator sideIterEnd = u->template end<side_type>(si);
				for(  ;sideIter !=sideIterEnd; sideIter++)
				{
					side_type* side = *sideIter;
					//	get position accessor
					typedef typename domain_type::position_accessor_type position_accessor_type;
					// const position_accessor_type& aaPos = domain.position_accessor();
					//side->vertex(0);
					//number z = aaPos[side->vertex(0)][0];
					//		number co0=0.5*(aaPos[side->vertex(0)][0]  + aaPos[side->vertex(1)][0] ) ;
					//		number co1=0.5*(aaPos[side->vertex(0)][1]  + aaPos[side->vertex(1)][1] ) ;
					//UG_LOG("[" << co0 << "," << co1 << "]\n");
					if (pbm && pbm->is_slave(side)){
						//UG_LOG("----\n");
						continue;
					}
					acUGrad[side]/=(number)acVol[side];
					acVGrad[side]/=(number)acVol[side];
					if (dim==3) acWGrad[side]/=(number)acVol[side];
//					UG_LOG("uGrad=" << acUGrad[side] << " vGrad=" << acVGrad[side] << "vol=" << acVol[side] << "\n");
				}
			}
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

};
	
} // namespace NavierStokes
} // end namespace ug


#endif /* __H__UG__NAVIER_STOKES_CENTRAL_GRADIENT__ */
