/*	
 *  pressure_gradient.h
 *
 *  Created on: 12.06.2013
 *      Author: Christian Wehner
 */

#ifndef __H__UG__NAVIER_STOKES_PRESSURE_GRADIENT__
#define __H__UG__NAVIER_STOKES_PRESSURE_GRADIENT__

#include "common/common.h"

#include "lib_disc/common/subset_group.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/local_finite_element/local_shape_function_set.h"
#include "lib_disc/spatial_disc/user_data/user_data.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "lib_disc/operator/non_linear_operator/newton_solver/newton_update_interface.h"
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
class PressureGradient
: 	public StdUserData<PressureGradient<TGridFunction>, MathVector<TGridFunction::dim>, TGridFunction::dim>,
  	virtual public INewtonUpdate
  	{
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;

	///	algebra type
	typedef typename TGridFunction::algebra_type algebra_type;

	/// position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;

	///	world dimension
	static const int dim = domain_type::dim;

	///	grid type
	typedef typename domain_type::grid_type grid_type;

	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object elem_type;

	/// side type
	typedef typename elem_type::side side_type;

	/// MathVector<dim> attachment
	//		typedef MathVector<dim> vecDim;
	//		typedef Attachment<vecDim> AMathVectorDim;

	/// attachment accessor
	typedef PeriodicAttachmentAccessor<VertexBase,ANumber > aVertexNumber;
	typedef Grid::AttachmentAccessor<elem_type,ANumber > aElementNumber;

	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// vertex iterator
	typedef typename TGridFunction::template traits<VertexBase>::const_iterator VertexIterator;

  		private:
	// pressure gradient attachment accessor
	typedef MathVector<dim> vecDim;
	typedef Attachment<vecDim> AMathVectorDim;
	typedef PeriodicAttachmentAccessor<elem_type,AMathVectorDim > aElemDimVector;

	aElemDimVector acPGrad;
	AMathVectorDim aPGrad;
	
	// level set grid function
	SmartPtr<TGridFunction> m_u;

	//	approximation space for level and surface grid
	SmartPtr<ApproximationSpace<domain_type> > m_spApproxSpace;

	//  grid
	grid_type* m_grid;

	//  pressure index
	static const size_t _P_ = dim;

	public:
	/// constructor
	PressureGradient(SmartPtr<ApproximationSpace<domain_type> > approxSpace,SmartPtr<TGridFunction> spGridFct){
		m_u = spGridFct;
		domain_type& domain = *m_u->domain().get();
		grid_type& grid = *domain.grid();
		m_grid = &grid;
		m_spApproxSpace = approxSpace;
		
		grid.template attach_to<elem_type>(aPGrad);
		// access
		acPGrad.access(grid,aPGrad);
//		SetAttachmentValues(aPGrad, m_u->template begin<elem_type>(), m_u->template end<elem_type>(), 0);
		this->update();
	}

	virtual ~PressureGradient(){};

	template <int refDim>
	inline void evaluate(MathVector<dim> vValue[],
	                     const MathVector<dim> vGlobIP[],
	                     number time, int si,
	                     GeometricObject* elem,
	                     const MathVector<dim> vCornerCoords[],
	                     const MathVector<refDim> vLocIP[],
	                     const size_t nip,
	                     LocalVector* u,
	                     const MathMatrix<refDim, dim>* vJT = NULL) const
	{
		vValue[0] = acPGrad[dynamic_cast<elem_type*>(elem)];
	};

	void update(){
		SmartPtr<TGridFunction> u = m_u;
		domain_type& domain = *u->domain().get();
		DimCRFVGeometry<dim> geo;

		position_accessor_type aaPos = m_u->domain()->position_accessor();

		typename grid_type::template traits<side_type>::secure_container sides;
		std::vector<MathVector<dim> > vCorner;
		std::vector<MultiIndex<2> > ind;

		for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
		{
			//	get iterators
			ElemIterator iter = u->template begin<elem_type>(si);
			ElemIterator iterEnd = u->template end<elem_type>(si);
			for(; iter != iterEnd; ++iter)
			{
				//	get the element
				elem_type* elem = *iter;

				MathVector<dim> vGlobalGrad=0;

				//  get sides of element
				m_grid->template associated_elements_sorted(sides, elem );

				//	reference object type
				ReferenceObjectID roid = elem->reference_object_id();

				//	get corners of element
				CollectCornerCoordinates(vCorner, *elem, aaPos);

				//	evaluate finite volume geometry
				geo.update(elem, &(vCorner[0]), domain.subset_handler().get());

				//	compute size (volume) of element
				const number elemSize = ElementSize<dim>(roid, &vCorner[0]);

				typename grid_type::template traits<elem_type>::secure_container assoElements;

				// assemble element-wise finite volume gradient
				for (size_t s=0;s<sides.size();s++){
					m_grid->template associated_elements(assoElements,sides[s]);
					// face value is average of associated elements
					number faceValue = 0;
					size_t numOfAsso = assoElements.size();
					for (size_t i=0;i<numOfAsso;i++){
						m_u->inner_multi_indices(assoElements[i], _P_, ind);
						faceValue+=DoFRef(*m_u,ind[0]);
					}
					faceValue/=(number)numOfAsso;
					const typename DimCRFVGeometry<dim>::SCV& scv = geo.scv(s);
					for (int d=0;d<dim;d++){
						vGlobalGrad[d] += faceValue * scv.normal()[d];
					}
				}
				vGlobalGrad/=(number)elemSize;

				//	write result in array storage
				acPGrad[elem]=vGlobalGrad;
			}

		}
	}

  		public:
	virtual void operator() (MathVector<dim>& value,
	                         const MathVector<dim>& globIP,
	                         number time, int si) const
	{
		UG_THROW("LevelSetUserData: Need element.");
	}

	virtual void operator() (MathVector<dim> vValue[],
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


#endif /* __H__UG__NAVIER_STOKES_PRESSURE_GRADIENT__ */
