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
#include "lib_disc/local_finite_element/local_shape_function_set.h"
#include "lib_disc/spatial_disc/user_data/user_data.h"

namespace ug{
namespace NavierStokes{

template <typename TData, int dim, typename TImpl,typename TGridFunction>
class StdTurbulentViscosityData
	: 	public UserData<TData,dim>
{
	public:
		////////////////
		// one value
		////////////////
		virtual void operator() (TData& value,
		                         const MathVector<dim>& globIP,
		                         number time, int si) const
		{
			UG_THROW("StdTurbulentViscosityData: Need element.");
		}

		virtual void operator() (TData& value,
		                         const MathVector<dim>& globIP,
		                         number time, int si,
		                         LocalVector& u,
		                         GeometricObject* elem,
		                         const MathVector<dim> vCornerCoords[],
		                         const MathVector<1>& locIP) const
		{
			getImpl().template evaluate<1>(&value,&globIP,time,si,u,elem,vCornerCoords,&locIP, 1, NULL);
		}

		virtual void operator() (TData& value,
		                         const MathVector<dim>& globIP,
		                         number time, int si,
		                         LocalVector& u,
		                         GeometricObject* elem,
		                         const MathVector<dim> vCornerCoords[],
		                         const MathVector<2>& locIP) const
		{
			getImpl().template evaluate<2>(&value,&globIP,time,si,u,elem,vCornerCoords,&locIP, 1, NULL);
		}

		virtual void operator() (TData& value,
		                         const MathVector<dim>& globIP,
		                         number time, int si,
		                         LocalVector& u,
		                         GeometricObject* elem,
		                         const MathVector<dim> vCornerCoords[],
		                         const MathVector<3>& locIP) const
		{
			getImpl().template evaluate<3>(&value,&globIP,time,si,u,elem,vCornerCoords,&locIP, 1, NULL);
		}

		////////////////
		// vector of values
		////////////////

		virtual void operator() (TData vValue[],
		                         const MathVector<dim> vGlobIP[],
		                         number time, int si, const size_t nip) const
		{
			UG_THROW("StdTurbulentViscosityData: Need element.");
		}


		virtual void operator()(TData vValue[],
		                        const MathVector<dim> vGlobIP[],
		                        number time, int si,
		                        LocalVector& u,
		                        GeometricObject* elem,
		                        const MathVector<dim> vCornerCoords[],
		                        const MathVector<1> vLocIP[],
		                        const size_t nip,
		                        const MathMatrix<1, dim>* vJT = NULL) const
		{
			getImpl().template evaluate<1>(vValue,vGlobIP,time,si,u,elem,
			                               vCornerCoords,vLocIP,nip, vJT);
		}

		virtual void operator()(TData vValue[],
		                        const MathVector<dim> vGlobIP[],
		                        number time, int si,
		                        LocalVector& u,
		                        GeometricObject* elem,
		                        const MathVector<dim> vCornerCoords[],
		                        const MathVector<2> vLocIP[],
		                        const size_t nip,
		                        const MathMatrix<2, dim>* vJT = NULL) const
		{
			getImpl().template evaluate<2>(vValue,vGlobIP,time,si,u,elem,
			                               vCornerCoords,vLocIP,nip, vJT);
		}

		virtual void operator()(TData vValue[],
		                        const MathVector<dim> vGlobIP[],
		                        number time, int si,
		                        LocalVector& u,
		                        GeometricObject* elem,
		                        const MathVector<dim> vCornerCoords[],
		                        const MathVector<3> vLocIP[],
		                        const size_t nip,
		                        const MathMatrix<3, dim>* vJT = NULL) const
		{
			getImpl().template evaluate<3>(vValue,vGlobIP,time,si,u,elem,
			                               vCornerCoords,vLocIP,nip, vJT);
		}

		virtual void compute(LocalVector* u, GeometricObject* elem, bool bDeriv = false)
		{
			UG_THROW("Not implemented.");
		}

		virtual bool update(const TGridFunction& u) = 0;

	protected:
	///	access to implementation
		TImpl& getImpl() {return static_cast<TImpl&>(*this);}

	///	const access to implementation
		const TImpl& getImpl() const {return static_cast<const TImpl&>(*this);}
};


template <typename TGridFunction>
class CRSmagorinskyTurbViscData
	: public StdTurbulentViscosityData<number, TGridFunction::dim,
	  	  	  CRSmagorinskyTurbViscData<TGridFunction>,TGridFunction >
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

		typedef typename Grid::AttachmentAccessor<side_type,Attachment<number > > aNumber;

	private:
		
		grid_type* m_grid;

	//  turbulent viscosity attachment
		aNumber m_aTurbulentViscosity;
		
		bool m_init;
		
	//  Smagorinsky model parameter, typical values [0.01 0.1]
		number m_c;

		void initTVAttachment(TGridFunction u){

			//	get domain of grid function
			domain_type& domain = *u.domain().get();

			//	get grid type of domain
			typedef typename domain_type::grid_type grid_type;

			//	get grid of domain
			grid_type& grid = *domain.grid();

			m_grid = grid;
						
			grid.template attach_to<side_type>(m_aTurbulentViscosity);
			
		};

	public:
	/// constructor
		CRSmagorinskyTurbViscData(number c = 0.05) : m_init(false) {
			m_c = c;
		}
		
		virtual ~CRSmagorinskyTurbViscData() {};

		void set_model_parameter(number c){
			m_c = c;
		}

		void reset(){ m_init=false; }
		
		template <int refDim>
		inline void evaluate(number vValue[],
		                     const MathVector<dim> vGlobIP[],
		                     number time, int si,
		                     LocalVector& u,
		                     GeometricObject* elem,
		                     const MathVector<dim> vCornerCoords[],
		                     const MathVector<refDim> vLocIP[],
		                     const size_t nip,
		                     const MathMatrix<refDim, dim>* vJT = NULL) const
		{
		//	reference object id
			ReferenceObjectID roid = elem->reference_object_id();

			typename grid_type::template traits<side_type>::secure_container sides;

			UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Only elements of type elem_type are currently supported");

			m_grid->associated_elements_sorted(sides, static_cast<elem_type*>(elem) );

		//	get trial space
			try{
			const LocalShapeFunctionSet<refDim>& rTrialSpace =
					LocalShapeFunctionSetProvider::get<refDim>(roid, LFEID(LFEID::CROUZEIX_RAVIART, 1));

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
					const number valSH = m_aTurbulentViscosity[sides[sh]];
					 vValue[ip] += valSH * vShape[sh];
				}
			}

			}catch(UGError_LocalShapeFunctionSetNotRegistered& ex){
				UG_THROW("TurbulentViscosityData: "<< ex.get_msg()<<", Reference Object: "
				         <<roid<<", Trial Space: CROUZEIX_RAVIART, refDim="<<refDim);
			}
		}
		
		bool update(const TGridFunction& u);
};

template <typename TGridFunction>
class CRDynamicTurbViscData
	: public StdTurbulentViscosityData<number, TGridFunction::dim,
	  	  	  CRDynamicTurbViscData<TGridFunction>,TGridFunction >
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

		typedef typename Grid::AttachmentAccessor<side_type,Attachment<number > > aNumber;

	private:
		
		grid_type* m_grid;

	//  turbulent viscosity attachment
		aNumber m_aTurbulentViscosity;
		
		bool m_init;
		
		void initTVAttachment(TGridFunction u){

			//	get domain of grid function
			domain_type& domain = *u.domain().get();

			//	get grid type of domain
			typedef typename domain_type::grid_type grid_type;

			//	get grid of domain
			grid_type& grid = *domain.grid();

			m_grid = grid;
						
			grid.template attach_to<side_type>(m_aTurbulentViscosity);
			
		};

	public:
	/// constructor
		CRDynamicTurbViscData(){ m_init = false; }
		
		virtual ~CRDynamicTurbViscData() {};

		void reset(){ m_init=false; }
		
		template <int refDim>
		inline void evaluate(number vValue[],
		                     const MathVector<dim> vGlobIP[],
		                     number time, int si,
		                     LocalVector& u,
		                     GeometricObject* elem,
		                     const MathVector<dim> vCornerCoords[],
		                     const MathVector<refDim> vLocIP[],
		                     const size_t nip,
		                     const MathMatrix<refDim, dim>* vJT = NULL) const
		{
		//	reference object id
			ReferenceObjectID roid = elem->reference_object_id();

			typename grid_type::template traits<side_type>::secure_container sides;

			UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Only elements of type elem_type are currently supported");

			m_grid->associated_elements_sorted(sides, static_cast<elem_type*>(elem) );

		//	get trial space
			try{
			const LocalShapeFunctionSet<refDim>& rTrialSpace =
					LocalShapeFunctionSetProvider::get<refDim>(roid, LFEID(LFEID::CROUZEIX_RAVIART, 1));

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
					const number valSH = m_aTurbulentViscosity[sides[sh]];
					 vValue[ip] += valSH * vShape[sh];
				}
			}

			}catch(UGError_LocalShapeFunctionSetNotRegistered& ex){
				UG_THROW("TurbulentViscosityData: "<< ex.get_msg()<<", Reference Object: "
				         <<roid<<", Trial Space: CROUZEIX_RAVIART, refDim="<<refDim);
			}
		}
		
		bool update(const TGridFunction& u);
};



} // namespace NavierStokes
} // end namespace ug

// include implementation
#include "turbulent_viscosity_data_impl.h"

#endif /* __H__UG__NAVIER_STOKES_TURBULENT_VISCOSITY_DATA__ */
