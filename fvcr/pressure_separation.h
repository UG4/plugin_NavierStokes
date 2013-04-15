/*	
 *  pressure_separation.h
 *
 *  Created on: 06.02.2013
 *      Author: Christian Wehner
 */

#ifndef __H__UG__NAVIER_STOKES_PRESSURE_SEPARATION__
#define __H__UG__NAVIER_STOKES_PRESSURE_SEPARATION__

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

/**
concept derived from grid_function_user_data.h
 */
template <typename TGridFunction>
class SeparatedPressureSource
: 	public StdUserData<SeparatedPressureSource<TGridFunction>, MathVector<TGridFunction::dim>, TGridFunction::dim>,
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
	// old pressure attachment accessor
	ANumber m_aPOld;
	aElementNumber m_pOld;

	//	pressure attachment accessor (interpolated pressure in vertices)
	ANumber m_aP;
	aVertexNumber m_p;

	//  volume attachment accessor
	ANumber m_aVol;
	aVertexNumber m_vol;

	// level set grid function
	SmartPtr<TGridFunction> m_u;

	//	approximation space for level and surface grid
	SmartPtr<ApproximationSpace<domain_type> > m_spApproxSpace;

	//  grid
	grid_type* m_grid;

	//  pressure index
	static const size_t _P_ = dim;

  		private:

	///	Data import for source
	SmartPtr<CplUserData<MathVector<dim>,dim> > m_imSource;

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

  		public:
	/// constructor
	SeparatedPressureSource(SmartPtr<ApproximationSpace<domain_type> > approxSpace,SmartPtr<TGridFunction> spGridFct){
		m_u = spGridFct;
		domain_type& domain = *m_u->domain().get();
		grid_type& grid = *domain.grid();
		m_grid = &grid;
		m_spApproxSpace = approxSpace;
		set_source(0.0);
		grid.template attach_to<elem_type>(m_aPOld);
		grid.template attach_to<VertexBase>(m_aP);
		grid.template attach_to<VertexBase>(m_aVol);
		m_pOld.access(grid,m_aPOld);
		m_p.access(grid,m_aP);
		m_vol.access(grid,m_aVol);
		// set all values to zero
		SetAttachmentValues(m_vol, m_u->template begin<VertexBase>(), m_u->template end<VertexBase>(), 0);
		SetAttachmentValues(m_p, m_u->template begin<VertexBase>(), m_u->template end<VertexBase>(), 0);
		SetAttachmentValues(m_pOld, m_u->template begin<elem_type>(), m_u->template end<elem_type>(), 0);
		this->update();
	}

	virtual ~SeparatedPressureSource(){};

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
		UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Unsupported element type");
		elem_type* element = static_cast<elem_type*>(elem);

		//	reference object id
		ReferenceObjectID roid = elem->reference_object_id();

		const size_t numVertices = element->num_vertices();
		//    get domain of grid function
		const domain_type& domain = *m_u->domain().get();

		//    get position accessor
		typedef typename domain_type::position_accessor_type position_accessor_type;
		const position_accessor_type& posAcc = domain.position_accessor();

		position_accessor_type aaPos = m_u->domain()->position_accessor();

		// coord and vertex array
		MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
		VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];
		DimCRFVGeometry<dim> crfvgeo;

		MathVector<dim> grad[domain_traits<dim>::MaxNumVerticesOfElem];

		for(size_t i = 0; i < numVertices; ++i){
			vVrt[i] = element->vertex(i);
			coCoord[i] = posAcc[vVrt[i]];
			// UG_LOG("co=[" << coCoord[i][0] << " " << coCoord[i][1] << "]\n");
		};

		for (size_t i=0;i<domain_traits<dim>::MaxNumVerticesOfElem;i++)
			grad[i]*=0.0;

		// evaluate finite volume geometry
		crfvgeo.update(elem, &(coCoord[0]), domain.subset_handler().get());

		// Lagrange 1 trial space
		const LocalShapeFunctionSet<dim>& lagrange1 =
				LocalShapeFunctionSetProvider::get<dim>(roid, LFEID(LFEID::LAGRANGE, 1));

		std::vector<number> shapes;

		for (size_t ip=0;ip<crfvgeo.num_scvf();ip++){
			const typename DimCRFVGeometry<dim>::SCVF& scvf = crfvgeo.scvf(ip);
			number pinter=0;
			lagrange1.shapes(shapes,scvf.local_ip());
			for (size_t sh=0;sh<numVertices;sh++){
				// UG_LOG("corner(" << sh << ")=" << m_p[vVrt[sh]] << "\n");
				pinter += m_p[vVrt[sh]]*shapes[sh];
			}
			// UG_LOG("pinter = " << pinter << "\n");
			for (int d=0;d<dim;d++){
				number flux = pinter*scvf.normal()[d];
				grad[scvf.from()][d]-=flux;
				grad[scvf.to()][d]+=flux;
			}
		}

		for (size_t i=0;i<crfvgeo.num_scv();i++){
			grad[i]/=crfvgeo.scv(i).volume();
		}

		// UG_LOG("gradient values:\n");
		for (size_t i=0;i<nip;i++){
			// UG_LOG(i << " " << grad[i] << "\n");
		}

		// evaluate source data
		(*m_imSource)(vValue,
				vGlobIP,
				time, si,
				elem,
				vCornerCoords,
				vLocIP,
				nip,
				u,
				vJT);
		// UG_LOG("source values:\n");
		for (size_t i=0;i<nip;i++){
			// UG_LOG(i << " " << vValue[i] << "\n");
		}
		for (size_t i=0;i<nip;i++)
			vValue[i]+=grad[i];
	}; // evaluate

	void test(){

	}

	void update(){
		//	get domain
		domain_type& domain = *m_u->domain().get();
		//	create Multiindex
		std::vector<MultiIndex<2> > multInd;
		DimFV1Geometry<dim> geo;

		//	coord and vertex array
		MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
		VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

		//	get position accessor
		typedef typename domain_type::position_accessor_type position_accessor_type;
		const position_accessor_type& posAcc = domain.position_accessor();

		// set volume and p values to zero
		SetAttachmentValues(m_vol, m_u->template begin<VertexBase>(), m_u->template end<VertexBase>(), 0);
		SetAttachmentValues(m_p, m_u->template begin<VertexBase>(), m_u->template end<VertexBase>(), 0);

		// set p^{old} = p^{old} + p^{sep}
		for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
			ElemIterator iter = m_u->template begin<elem_type>(si);
			ElemIterator iterEnd = m_u->template end<elem_type>(si);
			for(  ;iter !=iterEnd; ++iter)
			{
				elem_type* elem = *iter;
				m_u->inner_multi_indices(elem, _P_, multInd);
				// UG_LOG( DoFRef(*m_u,multInd[0]) << " " << m_pOld[elem] << " " << m_pOld[elem]+DoFRef(*m_u,multInd[0]) << "\n");
				m_pOld[elem]+=DoFRef(*m_u,multInd[0]);
			}
		}
		// compute pressure by averaging
		for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
			ElemIterator iter = m_u->template begin<elem_type>(si);
			ElemIterator iterEnd = m_u->template end<elem_type>(si);
			for(  ;iter !=iterEnd; ++iter)
			{
				elem_type* elem = *iter;
				number pValue = m_pOld[elem];
				const size_t numVertices = elem->num_vertices();
				for(size_t i = 0; i < numVertices; ++i){
					vVrt[i] = elem->vertex(i);
					coCoord[i] = posAcc[vVrt[i]];
				};
				geo.update(elem, &(coCoord[0]), domain.subset_handler().get());
				for(size_t i = 0; i < numVertices; ++i){
					number scvVol = geo.scv(i).volume();
					m_vol[vVrt[i]]+=scvVol;
					m_p[vVrt[i]]+=scvVol*pValue;
				}
			}
		}
		PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
		// go over all vertices and average
		for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
			VertexIterator iter = m_u->template begin<VertexBase>(si);
			VertexIterator iterEnd = m_u->template end<VertexBase>(si);
			for(  ;iter !=iterEnd; ++iter)
			{
				VertexBase* vrt = *iter;
				if (pbm && pbm->is_slave(vrt)) continue;
				m_p[vrt]/=m_vol[vrt];
				//UG_LOG(posAcc[vrt] << "\n");
				//UG_LOG(m_p[vrt] <<  " " << m_vol[vrt] << "\n");
			}
		}
	}

  		private:
	static const size_t max_number_of_ips = 20;

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


#endif /* __H__UG__NAVIER_STOKES_PRESSURE_SEPARATION__ */
