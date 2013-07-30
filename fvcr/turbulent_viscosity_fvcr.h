/*
 * turbulent_viscosity_data.h
 *
 *  Created on: 01.11.2012
 *      Author: Christian Wehner
 */

#ifndef __H__UG__NAVIER_STOKES_TURBULENT_VISCOSITY_DATA__
#define __H__UG__NAVIER_STOKES_TURBULENT_VISCOSITY_DATA__

#include "common/common.h"

#include "lib_disc/common/subset_group.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/spatial_disc/user_data/user_data.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "lib_disc/operator/non_linear_operator/newton_solver/newton_update_interface.h"
#include "lib_disc/spatial_disc/disc_util/fvcr_geom.h"
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
template <typename TData, int dim, typename TImpl,typename TGridFunction>
class StdTurbulentViscosityData
: 	public StdUserData<StdTurbulentViscosityData<TData,dim,TImpl,TGridFunction>, TData,dim>
{
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

	/// attachment accessor types
	typedef MathMatrix<dim,dim> dimMat;
	typedef Attachment<dimMat> ATensor;

	typedef MathVector<dim> vecDim;
	typedef Attachment<vecDim> AMathVectorDim;

	typedef PeriodicAttachmentAccessor<side_type,ANumber > aSideNumber;
	typedef PeriodicAttachmentAccessor<side_type,ATensor > aSideTensor;
	typedef PeriodicAttachmentAccessor<side_type,AMathVectorDim > aSideDimVector;
	public:
	virtual void operator() (TData& value,
	                         const MathVector<dim>& globIP,
	                         number time, int si) const
	{
		UG_THROW("StdTurbulentViscosityData: Need element.");
	}

	virtual void operator() (TData vValue[],
	                         const MathVector<dim> vGlobIP[],
	                         number time, int si, const size_t nip) const
	{
		UG_THROW("StdTurbulentViscosityData: Need element.");
	}

	template <int refDim>
	void evaluate(TData vValue[],
	             const MathVector<dim> vGlobIP[],
	             number time, int si,
	             GeometricObject* elem,
	             const MathVector<dim> vCornerCoords[],
	             const MathVector<refDim> vLocIP[],
	             const size_t nip,
	             LocalVector* u,
	             const MathMatrix<refDim, dim>* vJT = NULL) const
	{
		getImpl().template evaluate<refDim>(vValue,vGlobIP,time,si,elem,
		                                    vCornerCoords,vLocIP,nip,u,vJT);
	}

	virtual void compute(LocalVector* u, GeometricObject* elem,
	                     const MathVector<dim> vCornerCoords[], bool bDeriv = false)
	{
		const number t = this->time();
		const int si = this->subset();
		for(size_t s = 0; s < this->num_series(); ++s)
			getImpl().template evaluate<dim>(this->values(s), this->ips(s), t, si,
			                                 elem, NULL, this->template local_ips<dim>(s),
			                                 this->num_ip(s), u);
	}

	///	returns if provided data is continuous over geometric object boundaries
	virtual bool continuous() const {return false;}

	///	returns if grid function is needed for evaluation
	virtual bool requires_grid_fct() const {return true;}

	void assembleDeformationTensor(aSideTensor& aaDefTensor,aSideNumber& aaVol,SmartPtr<TGridFunction> u);
	void assembleDeformationTensor(aSideTensor& aaDefTensor,aSideNumber& aaVol,aSideDimVector aaU);

	number FNorm(MathMatrix<dim,dim> M);

	void addUiUjTerm(aSideTensor& aaDefTensor,const number factor,aSideDimVector aaU);
	void addUiUjTerm(aSideTensor& aaDefTensor,const number factor,SmartPtr<TGridFunction> u);

	template <typename VType>
	void elementFilter(PeriodicAttachmentAccessor<side_type,Attachment<VType> >& aaUHat,aSideNumber& aaVol,const PeriodicAttachmentAccessor<side_type,Attachment<VType> >& aaU);

	void elementFilter(aSideDimVector& aaUHat,aSideNumber& aaVol,SmartPtr<TGridFunction> u);

	template <typename VType>
	void scvFilter(PeriodicAttachmentAccessor<side_type,Attachment<VType> >& aaUHat,aSideNumber& aaVol,const PeriodicAttachmentAccessor<side_type,Attachment<VType> >& aaU);

	void scvFilter(aSideDimVector& aaUHat,aSideNumber& aaVol,SmartPtr<TGridFunction> u);

	void fillAttachment(aSideDimVector& aaU,SmartPtr<TGridFunction> u);

	void transferToLowerLevels(aSideNumber& aaData,ApproximationSpace<domain_type>& approximationSpace);

	void scaleTensorByNorm(aSideTensor& aaTensor);

	// set non-periodic boundaries so that viscosity can be set to zero there
	void setTurbulenceZeroBoundaries(const char* subsets){
		try{
			m_turbZeroSg = m_uInfo->subset_grp_by_name(subsets);
		}UG_CATCH_THROW("ERROR while parsing Subsets.");
	}

	protected:
	///	access to implementation
	TImpl& getImpl() {return static_cast<TImpl&>(*this);}

	///	const access to implementation
	const TImpl& getImpl() const {return static_cast<const TImpl&>(*this);}

	// grid function
	SmartPtr<TGridFunction> m_uInfo;

	// subset group
	SubsetGroup m_turbZeroSg;
};


template <typename TGridFunction>
class CRSmagorinskyTurbViscData
: public StdTurbulentViscosityData<number, TGridFunction::dim,
  CRSmagorinskyTurbViscData<TGridFunction>,TGridFunction >, virtual public INewtonUpdate
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

	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// side iterator
	typedef typename TGridFunction::template traits<side_type>::const_iterator SideIterator;

	/// attachment accessor types
	typedef MathMatrix<dim,dim> dimMat;
	typedef Attachment<dimMat> ATensor;

	typedef PeriodicAttachmentAccessor<side_type,ANumber > aSideNumber;
	typedef PeriodicAttachmentAccessor<side_type,ATensor > aSideTensor;

	  public:
	/**
	 * This method sets the kinematic viscosity value. Kinematic viscosity is added to turbulent viscosity in evaluation routine.
	 *
	 * \param[in]	data		kinematic Viscosity
	 */
	///	\{
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
	///	\}

	  private:
	///	Data import for kinematic viscosity
	SmartPtr<CplUserData<number,dim> > m_imKinViscosity;

	  private:
	// grid function
	SmartPtr<TGridFunction> m_u;

	grid_type* m_grid;

	//  turbulent viscosity attachment
	aSideNumber m_acTurbulentViscosity;
	ANumber m_aTurbulentViscosity;

	//  volume attachment
	aSideNumber m_acVolume;
	ANumber m_aVolume;

	//	deformation tensor attachment
	aSideTensor m_acDeformation;
	ATensor m_aDeformation;

	//  Smagorinsky model parameter, typical values [0.01 0.1]
	number m_c;

	//	approximation space for level and surface grid
	SmartPtr<ApproximationSpace<domain_type> > m_spApproxSpace;

	//  periodic boundary manager
	PeriodicBoundaryManager* m_pbm;

	  public:
	/// constructor
	CRSmagorinskyTurbViscData(SmartPtr<ApproximationSpace<domain_type> > approxSpace,SmartPtr<TGridFunction> spGridFct,number c = 0.05){
		m_c = c;
		m_u = spGridFct;
		this->m_uInfo = m_u;
		m_spApproxSpace = approxSpace;
		domain_type& domain = *m_u->domain().get();
		grid_type& grid = *domain.grid();
		m_grid = &grid;
		m_pbm = m_grid->periodic_boundary_manager();
		// attachments
		grid.template attach_to<side_type>(m_aTurbulentViscosity);
		grid.template attach_to<side_type>(m_aVolume);
		grid.template attach_to<side_type>(m_aDeformation);
		// accessors
		m_acTurbulentViscosity.access(grid,m_aTurbulentViscosity);
		m_acVolume.access(grid,m_aVolume);
		m_acDeformation.access(grid,m_aDeformation);

	}

	virtual ~CRSmagorinskyTurbViscData() {
		domain_type& domain = *m_u->domain().get();
		grid_type& grid = *domain.grid();
		grid.template detach_from<side_type>(m_aTurbulentViscosity);
		grid.template detach_from<side_type>(m_aVolume);
		grid.template detach_from<side_type>(m_aDeformation);
	};

	void set_model_parameter(number c){
		m_c = c;
	}

	template <int refDim>
	inline void evaluate(number vValue[],
	                     const MathVector<dim> vGlobIP[],
	                     number time, int si,
	                     GeometricObject* elem,
	                     const MathVector<dim> vCornerCoords[],
	                     const MathVector<refDim> vLocIP[],
	                     const size_t nip,
	                     LocalVector* u,
	                     const MathMatrix<refDim, dim>* vJT = NULL) const
	{
		//	reference object id
		ReferenceObjectID roid = elem->reference_object_id();

		typename grid_type::template traits<side_type>::secure_container sides;

		UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Unsupported element type");

		m_grid->associated_elements_sorted(sides, static_cast<elem_type*>(elem) );

		//	get trial space
		try{
			const LocalShapeFunctionSet<refDim>& rTrialSpace =
					LocalFiniteElementProvider::get<refDim>(roid, LFEID(LFEID::CROUZEIX_RAVIART, refDim, 1));

			//	memory for shapes
			std::vector<number> vShape;

			//	loop ips
			for(size_t ip = 0; ip < nip; ++ip)
			{
				//	evaluate at shapes at ip
				rTrialSpace.shapes(vShape, vLocIP[ip]);

				// 	compute solution at integration point
				vValue[ip] = 0.0;
				for(size_t sh = 0; sh < vShape.size(); ++sh)
				{
					const number valSH = m_acTurbulentViscosity[sides[sh]];
					vValue[ip] += valSH * vShape[sh];
				}
				//	UG_LOG(ip << " " << vValue[ip] << "\n");
			}

		}
		UG_CATCH_THROW("TurbulentViscosityData: trial space missing, Reference Object: "
				<<roid<<", Trial Space: CROUZEIX_RAVIART, refDim="<<refDim);

		number kinViscValues[max_number_of_ips];
		(*m_imKinViscosity)(kinViscValues,
				vGlobIP,
				time, si,
				elem,
				vCornerCoords,
				vLocIP,
				nip,
				u,
				vJT);
		for (size_t ip=0;ip < nip;ip++){
			// UG_LOG("turbVis(" << ip << ")=" << vValue[ip] << "+" << kinViscValues[ip] << "\n");
			vValue[ip] += kinViscValues[ip];
		}
	}

	static const size_t max_number_of_ips = 20;

	void update();

  };

template <typename TGridFunction>
class CRDynamicTurbViscData
: public StdTurbulentViscosityData<number, TGridFunction::dim,
  CRDynamicTurbViscData<TGridFunction>,TGridFunction >, virtual public INewtonUpdate
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

	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// side iterator
	typedef typename TGridFunction::template traits<side_type>::const_iterator SideIterator;

	/// attachment accessor types
	typedef MathMatrix<dim,dim> dimMat;
	typedef Attachment<dimMat> ATensor;

	typedef MathVector<dim> vecDim;
	typedef Attachment<vecDim> AMathVectorDim;

	typedef PeriodicAttachmentAccessor<side_type,ANumber > aSideNumber;
	typedef PeriodicAttachmentAccessor<side_type,ATensor > aSideTensor;
	typedef PeriodicAttachmentAccessor<side_type,AMathVectorDim > aSideDimVector;

	  public:
	/**
	 * This method sets the kinematic viscosity value. Kinematic viscosity is added to turbulent viscosity in evaluation routine.
	 *
	 * \param[in]	data		kinematic Viscosity
	 */
	///	\{
	void set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > user){
		m_imKinViscosity = user;
	}
	void set_kinematic_viscosity(number val){
		m_viscosityNumber = val;
		set_kinematic_viscosity(CreateSmartPtr(new ConstUserNumber<dim>(val)));
	}
#ifdef UG_FOR_LUA
	void set_kinematic_viscosity(const char* fctName){
		set_kinematic_viscosity(LuaUserDataFactory<number, dim>::create(fctName));
	}
#endif
	///	\}

	  private:
	///	Data import for kinematic viscosity
	SmartPtr<CplUserData<number,dim> > m_imKinViscosity;

	number m_viscosityNumber;

	static const number m_small = 1e-8;

	  private:
	// grid function
	SmartPtr<TGridFunction> m_u;

	grid_type* m_grid;

	//  turbulent viscosity attachment
	aSideNumber m_acTurbulentViscosity;
	ANumber m_aTurbulentViscosity;

	//  turbulent model parameter attachment
	aSideNumber m_acTurbulentC;
	ANumber m_aTurbulentC;

	//  new turbulent model parameter attachment (needed in time filtering)
	aSideNumber m_acTurbulentCNew;
	ANumber m_aTurbulentCNew;

	//  volume attachment
	aSideNumber m_acVolume;
	ANumber m_aVolume;

	//  coarser grid volume attachment
	aSideNumber m_acVolumeHat;
	ANumber m_aVolumeHat;

	//  filtered u attachment
	aSideDimVector m_acUHat;
	AMathVectorDim m_aUHat;

	//	deformation tensor attachment
	aSideTensor m_acDeformation;
	ATensor m_aDeformation;

	//	coarser grid deformation tensor attachment
	aSideTensor m_acDeformationHat;
	ATensor m_aDeformationHat;

	//	Leonard tensor attachment
	aSideTensor m_acLij;
	ATensor m_aLij;

	//	Mij tensor attachment
	aSideTensor m_acMij;
	ATensor m_aMij;

	//	approximation space for level and surface grid
	SmartPtr<ApproximationSpace<domain_type> > m_spApproxSpace;

	//  periodic boundary manager
	PeriodicBoundaryManager* m_pbm;

	  public:
	void init(SmartPtr<ApproximationSpace<domain_type> > approxSpace,SmartPtr<TGridFunction> spGridFct,bool spaceFilter,number timeFilterEps){
		m_u = spGridFct;
		this->m_uInfo = m_u;
		m_spApproxSpace = approxSpace;
		domain_type& domain = *m_u->domain().get();
		grid_type& grid = *domain.grid();
		m_grid = &grid;
		m_pbm = m_grid->periodic_boundary_manager();
		// attachments
		grid.template attach_to<side_type>(m_aTurbulentViscosity);
		grid.template attach_to<side_type>(m_aTurbulentC);
		grid.template attach_to<side_type>(m_aTurbulentCNew);
		grid.template attach_to<side_type>(m_aVolume);
		grid.template attach_to<side_type>(m_aVolumeHat);
		grid.template attach_to<side_type>(m_aUHat);
		grid.template attach_to<side_type>(m_aDeformation);
		grid.template attach_to<side_type>(m_aDeformationHat);
		grid.template attach_to<side_type>(m_aLij);
		grid.template attach_to<side_type>(m_aMij);
		// accessors
		m_acTurbulentViscosity.access(grid,m_aTurbulentViscosity);
		m_acTurbulentC.access(grid,m_aTurbulentC);
		m_acTurbulentCNew.access(grid,m_aTurbulentCNew);
		m_acVolume.access(grid,m_aVolume);
		m_acVolumeHat.access(grid,m_aVolumeHat);
		m_acUHat.access(grid,m_aUHat);
		m_acDeformation.access(grid,m_aDeformation);
		m_acDeformationHat.access(grid,m_aDeformationHat);
		m_acLij.access(grid,m_aLij);
		m_acMij.access(grid,m_aMij);
		m_spaceFilter=spaceFilter;
		m_timeFilterEps=timeFilterEps;
		if (timeFilterEps==1) m_timeFilter=false;
	}
	/// constructor
	CRDynamicTurbViscData(SmartPtr<ApproximationSpace<domain_type> > approxSpace,SmartPtr<TGridFunction> spGridFct){
		// use default settings for filtering of model constant c, use space filter, no time filter (timeFilter eps=1)
		init(approxSpace,spGridFct,true,1);
	}

	CRDynamicTurbViscData(SmartPtr<ApproximationSpace<domain_type> > approxSpace,SmartPtr<TGridFunction> spGridFct,bool spaceFilter,number timeFilter){
		init(approxSpace,spGridFct,spaceFilter,timeFilter);
	}

	CRDynamicTurbViscData(SmartPtr<ApproximationSpace<domain_type> > approxSpace,SmartPtr<TGridFunction> spGridFct,bool spaceFilter,bool timeFilter){
		init(approxSpace,spGridFct,spaceFilter,1);
		set_time_filter(timeFilter);
	}

	virtual ~CRDynamicTurbViscData() {
		domain_type& domain = *m_u->domain().get();
		grid_type& grid = *domain.grid();
		grid.template detach_from<side_type>(m_aTurbulentViscosity);
		grid.template detach_from<side_type>(m_aTurbulentC);
		grid.template detach_from<side_type>(m_aVolume);
		grid.template detach_from<side_type>(m_aVolumeHat);
		grid.template detach_from<side_type>(m_aUHat);
		grid.template detach_from<side_type>(m_aDeformation);
		grid.template detach_from<side_type>(m_aDeformationHat);
		grid.template detach_from<side_type>(m_aLij);
		grid.template detach_from<side_type>(m_aMij);
	};

	template <int refDim>
	inline void evaluate(number vValue[],
	                     const MathVector<dim> vGlobIP[],
	                     number time, int si,
	                     GeometricObject* elem,
	                     const MathVector<dim> vCornerCoords[],
	                     const MathVector<refDim> vLocIP[],
	                     const size_t nip,
	                     LocalVector* u,
	                     const MathMatrix<refDim, dim>* vJT = NULL) const
	{
		//	reference object id
		ReferenceObjectID roid = elem->reference_object_id();

		typename grid_type::template traits<side_type>::secure_container sides;

		UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Unsupported element type");

		m_grid->associated_elements_sorted(sides, static_cast<elem_type*>(elem) );

		//	get trial space
		try{
			const LocalShapeFunctionSet<refDim>& rTrialSpace =
					LocalFiniteElementProvider::get<refDim>(roid, LFEID(LFEID::CROUZEIX_RAVIART, refDim, 1));

			//	memory for shapes
			std::vector<number> vShape;

			//	loop ips
			for(size_t ip = 0; ip < nip; ++ip)
			{
				//	evaluate at shapes at ip
				rTrialSpace.shapes(vShape, vLocIP[ip]);

				// 	compute solution at integration point
				vValue[ip] = 0.0;
				for(size_t sh = 0; sh < vShape.size(); ++sh)
				{
					number valSH = m_acTurbulentViscosity[sides[sh]];
					vValue[ip] += valSH * vShape[sh];
				}
			}

		}
		UG_CATCH_THROW("TurbulentViscosityData: trial space missing, Reference Object: "
				<<roid<<", Trial Space: CROUZEIX_RAVIART, refDim="<<refDim);

		number kinViscValues[max_number_of_ips];
		(*m_imKinViscosity)(kinViscValues,
				vGlobIP,
				time, si,
				elem,
				vCornerCoords,
				vLocIP,
				nip,
				u,
				vJT);
		for (size_t ip=0;ip < nip;ip++){
			//  UG_LOG("turbVis(" << ip << ")=" << vValue[ip] << "+" << kinViscValues[ip] << "\n");
			vValue[ip] += kinViscValues[ip];
		}
	}

	static const size_t max_number_of_ips = 20;

	bool m_spaceFilter;
	bool m_timeFilter;
	number m_timeFilterEps;

	void update();

	void set_space_filter(bool b){
		m_spaceFilter=b;
	}
	void set_time_filter(bool b){
		m_timeFilter=b;
		m_timeFilterEps=0.001;
	}
	void set_time_filter_eps(number eps){
		if (eps!=1)
			m_timeFilter=true;
		else
			m_timeFilter=false;
		m_timeFilterEps=eps;
	}
  };


} // namespace NavierStokes
} // end namespace ug

// include implementation
#include "turbulent_viscosity_fvcr_impl.h"

#endif /* __H__UG__NAVIER_STOKES_TURBULENT_VISCOSITY_DATA__ */
