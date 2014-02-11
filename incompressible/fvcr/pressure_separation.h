/*	
 *  pressure_separation.h
 *
 *  Created on: 06.02.2013
 *      Author: Christian Wehner
 */

#ifndef __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FVCR__PRESSURE_SEPARATION__
#define __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FVCR__PRESSURE_SEPARATION__

#include "common/common.h"

#include "lib_disc/common/subset_group.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
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

template<typename TElem,typename TPos,int dim>
void computeElemBarycenter(MathVector<dim>& bary,TElem* elem,TPos posA){
	bary = 0;
	size_t noc=elem->num_vertices();
	for (size_t i=0;i<noc;i++){
		bary+=posA[elem->vertex(i)];
	}
	bary/=noc;
}

bool multMatVec(const std::vector<number>& avec, const std::vector<number>& b,std::vector<number>& c,size_t m,size_t n){
   size_t i;
   size_t j;
   size_t count=0;
   for(i = 0; i < m; i = i + 1){
	   c[i]=0.0;
       for(j = 0; j < n; j = j + 1){
            c[i] = c[i] + avec[count] * b[j];
		    count++;
	   }
   }
   return true;
}


bool solveLS(std::vector<number>& x,/* solution */
		const std::vector<number>& matField,/* matrix given as field */
            const std::vector<number>& b /* rhs */){
	size_t n=b.size();
	number P[n][n];
    std::vector<number> Pvec(n*n);
    std::vector<number> b2(n);
    number A[n][n];
	number L[n][n];
	number U[n][n];
	number z[n];
	size_t i,k,j,l;
	bool boolean=false;
	number d,max;
    size_t count=0;
    for (i=0;i<n;i++){
        for (j=0;j<n;j++){
            A[i][j]=matField[count];
            if (i!=j)
               P[i][j]=0;
            else
               P[i][j]=1;
            L[i][j]=0;
            U[i][j]=0;
            count ++;
        }
    }
	for (i=0;i<n;i++){
		j=i;
		// Suche groesstes Element in Spalten
		max=std::abs(A[i][i]);
		for (k=i+1;k<n;k++){
			if (std::abs(A[k][i])>max){
				j=k;
				max=std::abs(A[k][i]);
			};
		};
		//debug
		////UG_LOG(max << "\n");
		if (max<1e-14){
			// A nicht invertierbar
			// printf(".\n");
		    return false;
		};
		if (i!=j){
			// Vertauschen der Zeilen
			if (boolean==false){
				// beim ersten Durchlauf Vertauschen der Spalten der
				// Permutationsmatrix und der Zeilen von A,L.
				boolean=true;
				for (k=0;k<n;k++){
					d=P[k][i];
				    P[k][i]=P[k][j];
				    P[k][j]=d;
					d=A[i][k];
				    A[i][k]=A[j][k];
				    A[j][k]=d;
					d=L[i][k];
				    L[i][k]=L[j][k];
				    L[j][k]=d;
				};
			} else{
				for (k=0;k<n;k++){
				    // Vertauschen der Zeilen von P,A,L
				    d=P[i][k];
				    P[i][k]=P[j][k];
				    P[j][k]=d;
					d=A[i][k];
				    A[i][k]=A[j][k];
				    A[j][k]=d;
					d=L[i][k];
				    L[i][k]=L[j][k];
				    L[j][k]=d;
				};
			};
		};
		// Elimination
		L[i][i]=1;
		for (k=i+1;k<n;k++){
			L[k][i]=A[k][i]/A[i][i];
			for (l=i+1;l<n;l++){
				A[k][l]=A[k][l]-L[k][i]*A[i][l];
			};
		};
	};
	// U ist der obere Dreiecksteil von A
	for (i=0;i<n;i++){
		for (j=i;j<n;j++){
			U[i][j]=A[i][j];
		};
	};
	// A*x = b <-> P*A*x = P*b
    count = 0;
    for (i=0;i<n;i++){
        for (j=0;j<n;j++){
            Pvec[count]=P[i][j];
			count++;
        }
    }
	multMatVec(Pvec,b,b2,n,n);
	// Loese zuerst L*z = P*b
	for (i=0;i<n;i++){
		number s=b2[i];
		for (j=0;j<i;j++){
			s=s-z[j]*L[i][j];
		};
		z[i]=s;
	};
	// Loese U*x = z
	for (i=n-1; ;i--){
		number s=z[i];
		for (j=n-1;j>i;j--){
			s=s-x[j]*U[i][j];
		};
		x[i]=s/U[i][i];
		if (i==0) break;
	};
	return true;
};


bool leastSquares(std::vector<number>& x,const std::vector<number>& mField,const std::vector<number>& b){
	size_t m = b.size();
	size_t n = x.size();
	if (mField.size()!=m*n){
		n = mField.size()/b.size();
	}
	std::vector<number> tmmField(n*n);
	std::vector<number> tmb(n);
	number z;
	// compute A^t * A
	for (size_t i=0;i<n;i++){
		for (size_t j=i;j<n;j++){
			z=0;
			for (size_t k=0;k<m;k++){
				z+=mField[n*k+i]*mField[n*k+j];
			}
			tmmField[j*n+i]=z;
			if (j!=i) tmmField[i*n+j]=z;
		}
		tmb[i]=0;
		for (size_t k=0;k<m;k++) tmb[i]+= mField[n*k+i]*b[k];
	}
	return solveLS(x,tmmField,tmb);
}

int bubblesort(std::vector<int>& list,std::vector<number> array){
	int i,j,mini,k;
	int length=array.size();
	number min;
	for (i=0;i<length;i++)
		list[i]=i;
	for (i=0;i<length-1;i++){
		min = array[i];
		mini= i;
		for (j=i+1;j<length;j++){
			if (array[j]<min){
				mini = j;
				min = array[j];
			};
		};
		array[mini] = array[i];
		k=list[i];
		list[i]     = list[mini];
		list[mini]  = k;
		array[i]    = min;
	};
	return 0;
};

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

//		position_accessor_type aaPos = m_u->domain()->position_accessor();

		// coord and vertex array
		MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
		VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];
		DimCRFVGeometry<dim> crfvgeo;

		MathVector<dim> grad[domain_traits<dim>::MaxNumVerticesOfElem];

		for(size_t i = 0; i < numVertices; ++i){
			vVrt[i] = element->vertex(i);
			coCoord[i] = posAcc[vVrt[i]];
		};

		for (size_t i=0;i<domain_traits<dim>::MaxNumVerticesOfElem;i++)
			grad[i]*=0.0;

		// evaluate finite volume geometry
		crfvgeo.update(elem, &(coCoord[0]), domain.subset_handler().get());

		// Lagrange 1 trial space
		const LocalShapeFunctionSet<dim>& lagrange1 =
				LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::LAGRANGE, dim, 1));

		std::vector<number> shapes;

		for (size_t ip=0;ip<crfvgeo.num_scvf();ip++){
			const typename DimCRFVGeometry<dim>::SCVF& scvf = crfvgeo.scvf(ip);
			number pinter=0;
			lagrange1.shapes(shapes,scvf.local_ip());
			for (size_t sh=0;sh<numVertices;sh++)
				pinter += m_p[vVrt[sh]]*shapes[sh];
			for (int d=0;d<dim;d++){
				number flux = pinter*scvf.normal()[d];
				grad[scvf.from()][d]-=flux;
				grad[scvf.to()][d]+=flux;
			}
		}

		for (size_t i=0;i<crfvgeo.num_scv();i++){
			grad[i]/=crfvgeo.scv(i).volume();
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
		for (size_t i=0;i<nip;i++)
			vValue[i]+=grad[i];
	}; // evaluate

	void update(){
		//	get domain
		domain_type& domain = *m_u->domain().get();
		//	create Multiindex
		std::vector<DoFIndex> multInd;
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
				m_u->inner_dof_indices(elem, _P_, multInd);
				m_pOld[elem]+=DoFRef(*m_u,multInd[0]);
			}
		}
		// compute pressure in vertices by averaging
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

template <typename TGridFunction>
class SeparatedPressureSourceInter
: 	public StdUserData<SeparatedPressureSourceInter<TGridFunction>, MathVector<TGridFunction::dim>, TGridFunction::dim>,
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
	SeparatedPressureSourceInter(SmartPtr<ApproximationSpace<domain_type> > approxSpace,SmartPtr<TGridFunction> spGridFct){
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

	virtual ~SeparatedPressureSourceInter(){};

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
		};

		for (size_t i=0;i<domain_traits<dim>::MaxNumVerticesOfElem;i++)
			grad[i]*=0.0;

		// evaluate finite volume geometry
		crfvgeo.update(elem, &(coCoord[0]), domain.subset_handler().get());

		// Lagrange 1 trial space
		const LocalShapeFunctionSet<dim>& lagrange1 =
				LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::LAGRANGE, dim, 1));

		std::vector<number> shapes;

		for (size_t ip=0;ip<crfvgeo.num_scvf();ip++){
			const typename DimCRFVGeometry<dim>::SCVF& scvf = crfvgeo.scvf(ip);
			number pinter=0;
			lagrange1.shapes(shapes,scvf.local_ip());
			for (size_t sh=0;sh<numVertices;sh++)
				pinter += m_p[vVrt[sh]]*shapes[sh];
			for (int d=0;d<dim;d++){
				number flux = pinter*scvf.normal()[d];
				grad[scvf.from()][d]-=flux;
				grad[scvf.to()][d]+=flux;
			}
		}

		for (size_t i=0;i<crfvgeo.num_scv();i++){
			grad[i]/=crfvgeo.scv(i).volume();
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
		for (size_t i=0;i<nip;i++)
			vValue[i]+=grad[i];
	}; // evaluate

	void update(){
		//	get domain
		domain_type& domain = *m_u->domain().get();
		//	create Multiindex
		std::vector<DoFIndex> multInd;
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
				m_u->inner_dof_indices(elem, _P_, multInd);
				m_pOld[elem]+=DoFRef(*m_u,multInd[0]);
			}
		}
		// compute pressure in vertices by averaging
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


#endif /* __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FVCR__PRESSURE_SEPARATION__ */
