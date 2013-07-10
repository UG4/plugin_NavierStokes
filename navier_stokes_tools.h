/*
 * navier_stokes_tools.h
 *
 *  Created on: 13.12.2012
 *      Author: Christian Wehner
 */

#ifndef __H__UG__LIB_DISC__NAVIER_STOKES_TOOLS__
#define __H__UG__LIB_DISC__NAVIER_STOKES_TOOLS__

#include <vector>
#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_grid/lg_base.h"
#include "lib_grid/tools/periodic_boundary_manager.h"
#include "common/profiler/profiler.h"

namespace ug{

// compute vorticity
// for velocity field (u,v,w) the vorticity is \partial_x v - \partial_y u
template <typename TGridFunction>
void vorticityFVCR(TGridFunction& vort,TGridFunction& u)
{
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;

	///	grid type
	typedef typename domain_type::grid_type grid_type;

	///	world dimension
	static const int dim = domain_type::dim;

	// get grid
	grid_type& grid = *u.domain()->grid();

	// position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;

	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object elem_type;

    /// side type
	typedef typename elem_type::side side_type;

	//  volume attachment
	typedef PeriodicAttachmentAccessor<side_type,ANumber > aSideNumber;
	aSideNumber m_acVolume;
	ANumber m_aVolume;

	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// side iterator
	typedef typename TGridFunction::template traits<side_type>::const_iterator SideIterator;

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	position_accessor_type& posAcc = u.domain()->position_accessor();

	DimCRFVGeometry<dim> geo;

	grid.template attach_to<side_type>(m_aVolume);
	m_acVolume.access(grid,m_aVolume);

	SetAttachmentValues(m_acVolume,grid.template begin<side_type>(), grid.template end<side_type>(), 0);

	if (vort.local_finite_element_id(0) != LFEID(LFEID::CROUZEIX_RAVIART, dim, 1)){
				UG_THROW("Component " << 0 << " in approximation space of parameter 1 must be of Crouzeix-Raviart type.");
	}

	for (int d=0;d<dim;d++){
		if (u.local_finite_element_id(d) != LFEID(LFEID::CROUZEIX_RAVIART, dim, 1)){
			UG_THROW("Component " << d << " in approximation space of parameter 2 must be of Crouzeix-Raviart type.");
		}
	}

	PeriodicBoundaryManager* pbm = grid.periodic_boundary_manager();

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	//	create Multiindex
	std::vector<MultiIndex<2> > multInd;

	for(int si = 0; si < u.num_subsets(); ++si)
	{
		if (si>0) continue;
		//	get iterators
		ElemIterator iter = u.template begin<elem_type>(si);
		ElemIterator iterEnd = u.template end<elem_type>(si);

		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			elem_type* elem = *iter;
			//	get vertices and extract corner coordinates
			const size_t numVertices = elem->num_vertices();
			for(size_t i = 0; i < numVertices; ++i){
				vVrt[i] = elem->vertex(i);
				coCoord[i] = posAcc[vVrt[i]];
				// UG_LOG("co_coord(" << i<< "+1,:)=" << coCoord[i] << "\n");
			}
			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), u.domain()->subset_handler().get());

			size_t nofsides = geo.num_scv();

			typename grid_type::template traits<side_type>::secure_container sides;

			UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Only elements of type elem_type are currently supported");

			grid.associated_elements_sorted(sides, static_cast<elem_type*>(elem) );

			static const size_t MaxNumSidesOfElem = 18;

			typedef MathVector<dim> MVD;
			std::vector<MVD> uValue(MaxNumSidesOfElem);

			for (size_t s=0;s < nofsides;s++)
			{
				for (int d=0;d<dim;d++){
					//	get indices of function fct on vertex
					u.multi_indices(sides[s], d, multInd);
					//	read value of index from vector
					uValue[s][d]=DoFRef(u,multInd[0]);
				}
			}

			for (size_t s=0;s < nofsides;s++)
			{
				number localvort = 0;
				const typename DimCRFVGeometry<dim>::SCV& scv = geo.scv(s);

				for(size_t sh = 0 ; sh < nofsides; ++sh){
					localvort += uValue[sh][1] * scv.global_grad(sh)[0] - uValue[sh][0] * scv.global_grad(sh)[1];
				};

				number vol = scv.volume();

				localvort*=vol;

				vort.multi_indices(sides[s], 0, multInd);
				DoFRef(vort,multInd[0])+=localvort;
				m_acVolume[sides[s]] += vol;
			}
		}
		// average vorticity
		SideIterator sideIter = vort.template begin<side_type>(si);
		SideIterator sideIterEnd = vort.template end<side_type>(si);
		number maxvort = 0;
		for(  ;sideIter !=sideIterEnd; sideIter++)
		{
			//	get Elem
			side_type* elem = *sideIter;
			// if periodic slave continue
			if (pbm && pbm->is_slave(elem)) continue;
			vort.multi_indices(elem, 0, multInd);
			DoFRef(vort,multInd[0])/=m_acVolume[elem];
			if (DoFRef(vort,multInd[0])*DoFRef(vort,multInd[0])>maxvort*maxvort){
				maxvort = DoFRef(vort,multInd[0]);
			}
		}
	}
	grid.template detach_from<side_type>(m_aVolume);
}

template <typename TGridFunction>
void vorticityFV1(TGridFunction& vort,TGridFunction& u)
{
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;

	///	grid type
	typedef typename domain_type::grid_type grid_type;

	///	world dimension
	static const int dim = domain_type::dim;

	// get grid
	grid_type& grid = *u.domain()->grid();

	// position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;

	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object elem_type;

	//  volume attachment
	typedef PeriodicAttachmentAccessor<VertexBase,ANumber > aSideNumber;
	aSideNumber m_acVolume;
	ANumber m_aVolume;

	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// side iterator
	typedef typename TGridFunction::template traits<VertexBase>::const_iterator VertexBaseConstIterator;

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	position_accessor_type& posAcc = u.domain()->position_accessor();

	DimFV1Geometry<dim> geo;

	grid.template attach_to<VertexBase>(m_aVolume);
	m_acVolume.access(grid,m_aVolume);

	SetAttachmentValues(m_acVolume,grid.template begin<VertexBase>(), grid.template end<VertexBase>(), 0);

	if (vort.local_finite_element_id(0) != LFEID(LFEID::LAGRANGE, dim, 1)){
				UG_THROW("Component " << 0 << " in approximation space of parameter 1 must be of Crouzeix-Raviart type.");
	}

	for (int d=0;d<dim;d++){
		if (u.local_finite_element_id(d) != LFEID(LFEID::LAGRANGE, dim, 1)){
			UG_THROW("Component " << d << " in approximation space of parameter 2 must be of Crouzeix-Raviart type.");
		}
	}

	PeriodicBoundaryManager* pbm = grid.periodic_boundary_manager();

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	//	create Multiindex
	std::vector<MultiIndex<2> > multInd;

	for(int si = 0; si < u.num_subsets(); ++si)
	{
		if (si>0) continue;
		//	get iterators
		ElemIterator iter = u.template begin<elem_type>(si);
		ElemIterator iterEnd = u.template end<elem_type>(si);

		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			elem_type* elem = *iter;
			//	get vertices and extract corner coordinates
			const size_t numVertices = elem->num_vertices();
			for(size_t i = 0; i < numVertices; ++i){
				vVrt[i] = elem->vertex(i);
				coCoord[i] = posAcc[vVrt[i]];
				// UG_LOG("co_coord(" << i<< "+1,:)=" << coCoord[i] << "\n");
			}
			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), u.domain()->subset_handler().get());

			static const size_t MaxNumSidesOfElem = 18;

			typedef MathVector<dim> MVD;
			std::vector<MVD> uValue(MaxNumSidesOfElem);

			for (size_t co=0;co < numVertices;co++)
			{
				for (int d=0;d<dim;d++){
					//	get indices of function fct on vertex
					u.multi_indices(elem->vertex(co), d, multInd);
					//	read value of index from vector
					uValue[co][d]=DoFRef(u,multInd[0]);
				}
			}

			for (size_t co=0;co < numVertices;co++)
			{
				number localvort = 0;
				const typename DimFV1Geometry<dim>::SCV& scv = geo.scv(co);

				for(size_t sh = 0 ; sh < scv.num_sh(); ++sh){
					localvort += uValue[sh][1] * scv.global_grad(sh)[0] - uValue[sh][0] * scv.global_grad(sh)[1];
				};

				number vol = scv.volume();

				localvort*=vol;

				vort.multi_indices(elem->vertex(co), 0, multInd);
				DoFRef(vort,multInd[0])+=localvort;
				m_acVolume[elem->vertex(co)] += vol;
			}
		}
		// average vorticity
		VertexBaseConstIterator vertexIter = vort.template begin<VertexBase>(si);
		VertexBaseConstIterator vertexIterEnd = vort.template end<VertexBase>(si);
		number maxvort = 0;
		for(  ;vertexIter !=vertexIterEnd; vertexIter++)
		{
			//	get Elem
			VertexBase* vrt = *vertexIter;
			// if periodic slave continue
			if (pbm && pbm->is_slave(vrt)) continue;
			vort.multi_indices(vrt, 0, multInd);
			DoFRef(vort,multInd[0])/=m_acVolume[vrt];
			if (DoFRef(vort,multInd[0])*DoFRef(vort,multInd[0])>maxvort*maxvort){
				maxvort = DoFRef(vort,multInd[0]);
			}
		}
	}
	grid.template detach_from<VertexBase>(m_aVolume);
}

template <typename TGridFunction>
void vorticity(TGridFunction& vort,TGridFunction& u){
	if (vort.local_finite_element_id(0) == LFEID(LFEID::LAGRANGE, TGridFunction::dim, 1)){
		vorticityFV1<TGridFunction>(vort,u);
	} else
		if (vort.local_finite_element_id(0) == LFEID(LFEID::CROUZEIX_RAVIART, TGridFunction::dim, 1)){
			vorticityFVCR<TGridFunction>(vort,u);
		} else {
			UG_LOG("Function type " << vort.local_finite_element_id(0) << " not supported in vorticity computation.");
		}
}

// array0 = array1
void copyGhiaNumbers(number array0[17],number array1[17]){
	for (size_t i=0;i<17;i++) array0[i]=array1[i];
}

void copyErturkNumbers(number array0[23],number array1[23]){
	for (size_t i=0;i<23;i++) array0[i]=array1[i];
}

template <typename TGridFunction>
void drivenCavityEvaluationErturk(TGridFunction& u,size_t Re){
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;

	///	grid type
	typedef typename domain_type::grid_type grid_type;

	///	world dimension
	static const int dim = domain_type::dim;

	// position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;

	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object elem_type;

	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	position_accessor_type& posAcc = u.domain()->position_accessor();

	bool foundRe = false;

	// data from Erturk, Corke, Gökcöl paper
	number yco[23]={   1.00000,
   0.99000,
   0.98000,
   0.97000,
   0.96000,
   0.95000,
   0.94000,
   0.93000,
   0.92000,
   0.91000,
   0.90000,
   0.50000,
   0.20000,
   0.18000,
   0.16000,
   0.14000,
   0.12000,
   0.10000,
   0.08000,
   0.06000,
   0.04000,
   0.02000,
   0.00000
	};
   number xco[23]={1.00000,
   0.98500,
   0.97000,
   0.95500,
   0.94000,
   0.92500,
   0.91000,
   0.89500,
   0.88000,
   0.86500,
   0.85000,
   0.50000,
   0.15000,
   0.13500,
   0.12000,
   0.10500,
   0.09000,
   0.07500,
   0.06000,
   0.04500,
   0.03000,
   0.01500,
   0.00000};
	number xLineReferenceValue[23];
	number yLineReferenceValue[23];
	number yLineReferenceValue1000[23]={   0.00000,
  -0.09730,
  -0.21730,
  -0.34000,
  -0.44170,
  -0.50520,
  -0.52630,
  -0.51320,
  -0.48030,
  -0.44070,
  -0.40280,
   0.02580,
   0.37560,
   0.37050,
   0.36050,
   0.34600,
   0.32730,
   0.30410,
   0.27460,
   0.23490,
   0.17920,
   0.10190,
   0.00000};
   number yLineReferenceValue2500[23]={
		   0.0000,-0.1675,-0.3725,-0.5192,-0.5603,-0.5268,-0.4741,-0.4321,-0.4042,-0.3843,-0.3671,
		0.0160,0.3918,0.4078,0.4187,0.4217,0.4142,0.3950,0.3649,0.3238,0.2633,0.1607,0.0000};
	number yLineReferenceValue5000[23]={
		0.0000,-0.2441,-0.5019,-0.5700,-0.5139,-0.4595,-0.4318,-0.4147,-0.3982,-0.3806,-0.3624,
		0.0117,0.3699,0.3878,0.4070,0.4260,0.4403,0.4426,0.4258,0.3868,0.3263,0.2160,0.0000};
	number yLineReferenceValue7500[23]={
		0.0000,-0.2991,-0.5550,-0.5434,-0.4748,-0.4443,-0.4283,-0.4118,-0.3938,-0.3755,-0.3574,
		0.0099,0.3616,0.3779,0.3950,0.4137,0.4337,0.4495,0.4494,0.4210,0.3608,0.2509,0.0000};
	number yLineReferenceValue10000[23]={
		0.0000,-0.3419,-0.5712,-0.5124,-0.4592,-0.4411,-0.4256,-0.4078,-0.3895,-0.3715,-0.3538,
		0.0088,0.3562,0.3722,0.3885,0.4056,0.4247,0.4449,0.4566,0.4409,0.3844,0.2756,0.0000};
	number yLineReferenceValue12500[23]={
		0.0000,-0.3762,-0.5694,-0.4899,-0.4534,-0.4388,-0.4221,-0.4040,-0.3859,-0.3682,-0.3508,
		0.0080,0.3519,0.3678,0.3840,0.4004,0.4180,0.4383,0.4563,0.4522,0.4018,0.2940,0.0000};
	number yLineReferenceValue15000[23]={
		0.0000,-0.4041,-0.5593,-0.4754,-0.4505,-0.4361,-0.4186,-0.4005,-0.3828,-0.3654,-0.3481,
		0.0074,0.3483,0.3641,0.3801,0.3964,0.4132,0.4323,0.4529,0.4580,0.4152,0.3083,0.0000};
	number yLineReferenceValue17500[23]={
		0.0000,-0.4269,-0.5460,-0.4664,-0.4482,-0.4331,-0.4153,-0.3975,-0.3800,-0.3627,-0.3457,
		0.0069,0.3452,0.3608,0.3767,0.3929,0.4093,0.4273,0.4484,0.4602,0.4254,0.3197,0.0000};
	number yLineReferenceValue20000[23]={
		0.0000,-0.4457,-0.5321,-0.4605,-0.4459,-0.4300,-0.4122,-0.3946,-0.3774,-0.3603,-0.3434,
		0.0065,0.3423,0.3579,0.3736,0.3897,0.4060,0.4232,0.4438,0.4601,0.4332,0.3290,0.0000};
	number yLineReferenceValue21000[23]={
		0.0000,-0.4522,-0.5266,-0.4588,-0.4449,-0.4287,-0.4110,-0.3936,-0.3764,-0.3593,-0.3425,
		0.0063,0.3413,0.3567,0.3725,0.3885,0.4048,0.4218,0.4420,0.4596,0.4357,0.3323,0.0000};
	number xLineReferenceValue1000[23]={
	   1.00000,
   0.84860,
   0.70650,
   0.59170,
   0.51020,
   0.45820,
   0.42760,
   0.41010,
   0.39930,
   0.39130,
   0.38380,
  -0.06200,
  -0.37560,
  -0.38690,
  -0.38540,
  -0.36900,
  -0.33810,
  -0.29600,
  -0.24720,
  -0.19510,
  -0.13920,
  -0.07570,
   0.00000
   };
   number xLineReferenceValue2500[23]={
		1.0000,0.7704,0.5924,0.4971,0.4607,0.4506,0.4470,0.4424,0.4353,0.4256,0.4141,
		-0.0403,-0.3228,-0.3439,-0.3688,-0.3965,-0.4200,-0.4250,-0.3979,-0.3372,-0.2547,-0.1517,
		0.0000};
	number xLineReferenceValue5000[23]={
		1.0000,0.6866,0.5159,0.4749,0.4739,0.4738,0.4683,0.4582,0.4452,0.4307,0.4155,
		-0.0319,-0.3100,-0.3285,-0.3467,-0.3652,-0.3876,-0.4168,-0.4419,-0.4272,-0.3480,-0.2223,
		0.0000};
	number xLineReferenceValue7500[23]={
		1.0000,0.6300,0.4907,0.4817,0.4860,0.4824,0.4723,0.4585,0.4431,0.4275,0.4123,
		-0.0287,-0.3038,-0.3222,-0.3406,-0.3587,-0.3766,-0.3978,-0.4284,-0.4491,-0.3980,-0.2633,
		0.0000};
	number xLineReferenceValue10000[23]={
		1.0000,0.5891,0.4837,0.4891,0.4917,0.4843,0.4711,0.4556,0.4398,0.4243,0.4095,
		-0.0268,-0.2998,-0.3179,-0.3361,-0.3543,-0.3721,-0.3899,-0.4142,-0.4469,-0.4259,-0.2907,
		0.0000};
	number xLineReferenceValue12500[23]={
		1.0000,0.5587,0.4833,0.4941,0.4937,0.4833,0.4684,0.4523,0.4366,0.4216,0.4070,
		-0.0256,-0.2967,-0.3146,-0.3326,-0.3506,-0.3685,-0.3859,-0.4054,-0.4380,-0.4407,-0.3113,
		0.0000};
	number xLineReferenceValue15000[23]={
		1.0000,0.5358,0.4850,0.4969,0.4937,0.4811,0.4653,0.4492,0.4338,0.4190,0.4047,
		-0.0247,-0.2942,-0.3119,-0.3297,-0.3474,-0.3652,-0.3827,-0.4001,-0.4286,-0.4474,-0.3278,
		0.0000};
	number xLineReferenceValue17500[23]={
		1.0000,0.5183,0.4871,0.4982,0.4925,0.4784,0.4622,0.4463,0.4312,0.4166,0.4024,
		-0.0240,-0.2920,-0.3096,-0.3271,-0.3446,-0.3622,-0.3797,-0.3965,-0.4206,-0.4490,-0.3412,
		0.0000};
	number xLineReferenceValue20000[23]={
		1.0000,0.5048,0.4889,0.4985,0.4906,0.4754,0.4592,0.4436,0.4287,0.4142,0.4001,
		-0.0234,-0.2899,-0.3074,-0.3248,-0.3422,-0.3595,-0.3769,-0.3936,-0.4143,-0.4475,-0.3523,
		0.0000};
	number xLineReferenceValue21000[23]={
		1.0000,0.5003,0.4895,0.4983,0.4897,0.4742,0.4580,0.4425,0.4277,0.4132,0.3992,
		-0.0232,-0.2892,-0.3066,-0.3239,-0.3412,-0.3585,-0.3758,-0.3925,-0.4121,-0.4463,-0.3562,
		0.0000
	};
	switch (Re){
	case 1000:
		foundRe=true;
		copyErturkNumbers(xLineReferenceValue,xLineReferenceValue1000);
		copyErturkNumbers(yLineReferenceValue,yLineReferenceValue1000);
		break;
	case 2500:
		foundRe=true;
		copyErturkNumbers(xLineReferenceValue,xLineReferenceValue2500);
		copyErturkNumbers(yLineReferenceValue,yLineReferenceValue2500);
		break;
	case 5000:
		foundRe=true;
		copyErturkNumbers(xLineReferenceValue,xLineReferenceValue5000);
		copyErturkNumbers(yLineReferenceValue,yLineReferenceValue5000);
		break;
	case 7500:
		foundRe=true;
		copyErturkNumbers(xLineReferenceValue,xLineReferenceValue7500);
		copyErturkNumbers(yLineReferenceValue,yLineReferenceValue7500);
		break;
	case 10000:
		foundRe=true;
		copyErturkNumbers(xLineReferenceValue,xLineReferenceValue10000);
		copyErturkNumbers(yLineReferenceValue,yLineReferenceValue10000);
		break;
	case 12500:
		foundRe=true;
		copyErturkNumbers(xLineReferenceValue,xLineReferenceValue12500);
		copyErturkNumbers(yLineReferenceValue,yLineReferenceValue12500);
		break;
	case 15000:
		foundRe=true;
		copyErturkNumbers(xLineReferenceValue,xLineReferenceValue15000);
		copyErturkNumbers(yLineReferenceValue,yLineReferenceValue15000);
		break;
	case 17500:
		foundRe=true;
		copyErturkNumbers(xLineReferenceValue,xLineReferenceValue17500);
		copyErturkNumbers(yLineReferenceValue,yLineReferenceValue17500);
		break;
	case 20000:
		foundRe=true;
		copyErturkNumbers(xLineReferenceValue,xLineReferenceValue20000);
		copyErturkNumbers(yLineReferenceValue,yLineReferenceValue20000);
		break;
	case 21000:
		foundRe=true;
		copyErturkNumbers(xLineReferenceValue,xLineReferenceValue21000);
		copyErturkNumbers(yLineReferenceValue,yLineReferenceValue21000);
		break;
	}

	if (foundRe==false){
		UG_LOG("Reynolds number " << Re << " results not in Erturk data set.\n\n");
		return;
	}

	MathVector<dim> xLineEvalPos;
	MathVector<dim> yLineEvalPos;
	xLineEvalPos[0]=0.5;
	yLineEvalPos[1]=0.5;

	number xLineValue[23];
	number yLineValue[23];

	static const number unhandled = 1537454.54345435;

	for (size_t i=0;i<23;i++){
		xLineValue[i]=unhandled;
		yLineValue[i]=unhandled;
	}

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	//	create Multiindex
	std::vector<MultiIndex<2> > multInd;

	for(int si = 0; si < u.num_subsets(); ++si)
	{
		if (si>0) continue;
		//	get iterators
		ElemIterator iter = u.template begin<elem_type>(si);
		ElemIterator iterEnd = u.template end<elem_type>(si);

		static const number bigNumber = 1e+8;

		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			elem_type* elem = *iter;
			//	get vertices and extract corner coordinates
			const size_t numVertices = elem->num_vertices();
			number xmin=bigNumber,xmax=-bigNumber,ymin=bigNumber,ymax=-bigNumber;
			//UG_LOG("-----------------------------------------------------------------\n");
			// find element bounds
			for(size_t i = 0; i < numVertices; ++i){
				vVrt[i] = elem->vertex(i);
				coCoord[i] = posAcc[vVrt[i]];
				//UG_LOG("co(" << i+1 << ",:)=[" << coCoord[i][0] << "," << coCoord[i][1] << "];" << "\n");
				if (coCoord[i][0]<xmin) xmin=coCoord[i][0];
				if (coCoord[i][0]>xmax) xmax=coCoord[i][0];
				if (coCoord[i][1]<ymin) ymin=coCoord[i][1];
				if (coCoord[i][1]>ymax) ymax=coCoord[i][1];
			};
			//UG_LOG("co(" << numVertices+1 << ",:)=[" << coCoord[0][0] << "," << coCoord[0][1] << "];" << "\n");

			//UG_LOG("bound=[" << xmin << " " << xmax << " " << ymin << " " << ymax << "]\n");

			bool checkx=true;
			bool checky=true;

			if (xmin>0.5) checkx=false;
			else if (xmax<0.5) checkx=false;

			if (ymin>0.5) checky=false;
			else if (ymax<0.5) checky=false;

			if ((checkx==false)&&(checky==false)) continue;

			number interpolation;

			if (checkx==true){
				for (size_t i=0;i<23;i++)
				{
					if ((yco[i]>=ymin)&&(yco[i]<=ymax)){
						xLineEvalPos[1]=yco[i];
						//UG_LOG("p=[" << xLineEvalPos[0] << "," << xLineEvalPos[1] << "]\n");
						if ( ContainsPoint(elem, xLineEvalPos, posAcc) )
						{
							//UG_LOG("***************************************************\n");
							//UG_LOG(xLineEvalPos << "\n");
							if (yco[i]==0){
								//UG_LOG("##################################################################\n");
								//UG_LOG("##################################################################\n");
								//UG_LOG("##################################################################\n");
								//UG_LOG("##################################################################\n");
							}
							//	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
							ReferenceObjectID roid = (ReferenceObjectID) elem->reference_object_id();
							//	get the reference mapping for the element using global corners
							DimReferenceMapping<dim, dim>& mapping
								= ReferenceMappingProvider::get<dim, dim>(roid, coCoord);
							MathVector<dim> localPos;

							//	compute global integration points
							mapping.global_to_local(localPos, xLineEvalPos);

							//	local finite element id
							const LFEID m_id = u.local_finite_element_id(0);

							//	get trial space
							const LocalShapeFunctionSet<dim>& rTrialSpace =
								LocalFiniteElementProvider::get<dim>(roid, m_id);

							//	number of dofs on element
							const size_t num_sh = rTrialSpace.num_sh();

							//	get multiindices of element
							std::vector<MultiIndex<2> > ind;  // 	aux. index array
							u.multi_indices(elem, 0, ind);
							// 	compute approximated solution at integration point
							interpolation = 0.0;
							for(size_t sh = 0; sh < num_sh; ++sh)
							{
								//	get value at shape point (e.g. corner for P1 fct)
								const number valSH = DoFRef(u, ind[sh]);

								//	add shape fct at ip * value at shape
								interpolation += valSH * rTrialSpace.shape(sh, localPos);
							}
							// if there is a previous value, then average
							if (xLineValue[i]!=unhandled) xLineValue[i] = 0.5*(xLineValue[i] + interpolation);
							else xLineValue[i]=interpolation;
						}
					}
				}
			}

			if (checky==true){
				for (size_t i=0;i<23;i++)
				{
					if ((xco[i]>=xmin)&&(xco[i]<=xmax)){
						yLineEvalPos[0]=xco[i];
						if ( ContainsPoint(elem, yLineEvalPos, posAcc) )
						{
							//	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
							ReferenceObjectID roid = (ReferenceObjectID) elem->reference_object_id();
							//	get the reference mapping for the element using global corners
							DimReferenceMapping<dim, dim>& mapping
								= ReferenceMappingProvider::get<dim, dim>(roid, coCoord);
							MathVector<dim> localPos;

							//	compute global integration points
							mapping.global_to_local(localPos, yLineEvalPos);

							//	local finite element id
							const LFEID m_id = u.local_finite_element_id(0);

							//	get trial space
							const LocalShapeFunctionSet<dim>& rTrialSpace =
								LocalFiniteElementProvider::get<dim>(roid, m_id);

							//	number of dofs on element
							const size_t num_sh = rTrialSpace.num_sh();

							//	get multiindices of element
							std::vector<MultiIndex<2> > ind;  // 	aux. index array
							u.multi_indices(elem, 1, ind);
							// 	compute approximated solution at integration point
							interpolation = 0.0;
							for(size_t sh = 0; sh < num_sh; ++sh)
							{
								//	get value at shape point (e.g. corner for P1 fct)
								const number valSH = DoFRef(u, ind[sh]);

								//	add shape fct at ip * value at shape
								interpolation += valSH * rTrialSpace.shape(sh, localPos);
							}
							// if there is a previous value, then average
							if (yLineValue[i]!=unhandled) yLineValue[i] = 0.5*(yLineValue[i] + interpolation);
							else yLineValue[i]=interpolation;
						}
					}
				}
			}
		};// for(  ;iter !=iterEnd; ++iter)
	};// for(int si = 0; si < u.num_subsets(); ++si)
	UG_LOG("Comparison with Erturk data\n");
	UG_LOG("\nData evaluation for Re=" << Re << ":\n\n");
	UG_LOG("u values on line through x=0.5:" << "\n\n");
	number maxdiff = 0;
	number diffsum = 0;
	for (size_t i=0;i<23;i++){
		number localdiff=std::abs(xLineReferenceValue[i]-xLineValue[i]);
		UG_LOG("y1(" << i+1 << ") = " << yco[i] << "; u1(" << i+1 << ") = " << xLineValue[i] << "; u_ert(" << i+1 << ") = " << xLineReferenceValue[i] << "; uerror_ert(" << i+1 << ") = " << localdiff << ";\n");
		if (localdiff>maxdiff) maxdiff = localdiff;
		diffsum += localdiff;
	}
	UG_LOG("u max difference: " << maxdiff << "\nu average difference: " << (number)diffsum/23.0);
	UG_LOG("\n\nv values on line through y=0.5:" << "\n\n");
	maxdiff = 0;
	diffsum = 0;
	for (size_t i=0;i<23;i++){
		number localdiff=std::abs(yLineReferenceValue[i]-yLineValue[i]);
		UG_LOG("x1(" << i+1 << ") = " << xco[i] << "; v1(" << i+1 << ") = " << yLineValue[i] << "; v_ert(" << i+1 << ") = " << yLineReferenceValue[i] << "; verror_ert(" << i+1 << ") = " << std::abs(yLineReferenceValue[i]-yLineValue[i]) << ";\n");
		if (localdiff>maxdiff) maxdiff = localdiff;
		diffsum += localdiff;
	}
	UG_LOG("v max difference: " << maxdiff << "\nv average difference: " << (number)diffsum/23.0);
	UG_LOG("\n\n");
}


template <typename TGridFunction>
void drivenCavityEvaluation(TGridFunction& u,size_t Re){
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;

	///	grid type
	typedef typename domain_type::grid_type grid_type;

	///	world dimension
	static const int dim = domain_type::dim;

	// position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;

	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object elem_type;

	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	position_accessor_type& posAcc = u.domain()->position_accessor();

	bool foundRe = false;

	// data from Ghia paper, see also cast3m implementation
	number yco[17]={0.0000,0.0547,0.0625,0.0703,0.1016,0.1719,0.2813,0.4531,0.5000,0.6172,0.7344,0.8516,0.9531,0.9609,0.9688,0.9766,1.0000};
	number xco[17]={0.0000,0.0625,0.0703,0.0781,0.0938,0.1563,0.2266,0.2344,0.5000,0.8047,0.8594,0.9063,0.9453,0.9531,0.9609,0.9688,1.0000};
	number xLineReferenceValue[17];
	number yLineReferenceValue[17];
	number xLineReferenceValue100[17]={0.,-0.03717,-0.04192,-0.04775,-0.06434,-0.10150,-0.15662,-0.2109,-0.20581,-0.13641,0.00332,0.23151,0.68717,
						0.73722,0.78871,0.84123,1.};
	number yLineReferenceValue100[17]={0.00000,0.09233,0.10091,0.10890,0.12317,0.16077,0.17507
						,0.17527
						,0.05454,-0.24533,-0.22445,-0.16914,-0.10313,-0.08864
						,-0.07391,-0.05906,0.00000};
	number xLineReferenceValue400[17]={0.,-0.08186,-0.09266,-0.10338,-0.14612,-0.24299,-0.32726,-0.17119,-0.11477,0.02135,0.16256,0.29093,0.55892,
					0.61756,0.68439,0.75837,1.};
				number yLineReferenceValue400[17]={0.00000,0.18360,0.19713,0.20920,0.22965,0.28124,0.30203
						,0.30174
						,0.05186,-0.38598,-0.44993,-0.3827,-0.22847,-0.19254
						,-0.15663,-0.12146,0.00000};
			number xLineReferenceValue1000[17]={0.,-0.18109,-0.20196,-0.22220,-0.29730,-0.38289,-0.27805,-0.10648,-0.06080,0.05702,0.18719,0.33304,0.46604,
					0.51117,0.57492,0.65928,1.};
			number yLineReferenceValue1000[17]={0.00000,0.27485,0.29012,0.30353,0.32627,0.37095,0.33075
					,0.32235
					,0.02526,-0.31966,-0.42665,-0.51550,-0.39188,-0.33714
					,-0.27669,-0.21388,0.00000};
			number xLineReferenceValue3200[17]={0.00000,-0.32407,-0.35344,-0.37827,-0.41933,-0.34323
					,-0.24427,-0.086636
					,-0.04272,0.07156,0.19791,0.34682,0.46101,0.46547,0.48296
					,0.53236,1.00000};
			number yLineReferenceValue3200[17]={0.00000,0.39560,0.40917,0.41906,0.42768,0.37119,0.29030
					,0.28188
					,0.00999,-0.31184,-0.37401,-0.44307,-0.54053,-0.52357
					,-0.47425,-0.39017,0.00000};
		number xLineReferenceValue5000[17]={0.00000,-0.41165,-0.42901,-0.43643,-0.40435,-0.33050
				,-0.22855,-0.07404
				,-0.03039,0.08183,0.20087,0.33556,0.46036,0.45992,0.46120
				,0.48223,1.00000};
		number yLineReferenceValue5000[17]={0.00000,0.42447,0.43329,0.43648,0.42951,0.35368,0.28066
				,0.27280
				,0.00945,-0.30018,-0.36214,-0.41442,-0.52876,-0.55408
				,-0.55069,-0.49774,0.00000};
		number xLineReferenceValue7500[17]={0.00000,-0.43154,-0.43590,-0.43025,-0.38324,-0.32393
				,-0.23176,-0.07503
				,-0.03800,0.08342,0.20591,0.34228,0.47167,0.47323,0.47048
				,0.47244,1.00000};
		number yLineReferenceValue7500[17]={0.00000,0.43979,0.44030,0.43564,0.41824,0.35060,0.28117
				,0.27348
				,0.00824,-0.30448,-0.36213,-0.41050,-0.48590,-0.52347
				,-0.55216,-0.53858,0.00000};
		number xLineReferenceValue10000[17]={0.00000,-0.42735,-0.42537,-0.41657,-0.38000,-0.32709
				,-0.23186,-0.07540
				,-0.03111,0.08344,0.20673,0.34635,0.47804,0.48070,0.47783
				,0.47221,1.00000};
		number yLineReferenceValue10000[17]={0.00000,0.43983,0.43733,0.43124,0.41487,0.35070,0.28003
				,0.27224
				,0.00831,-0.30719,-0.36737,-0.41496,-0.45863,-0.49099
				,-0.52987,-0.54302,0.00000};
		// Botella reference data for Re=1000
		number xLineReferenceValue1000Botella[17]={0.0000000 , -0.1812881 , -0.2023300 , -0.2228955 , -0.3004561 , -0.3885691 , -0.2803696 , -0.1081999 , -0.0620561 , 
		0.0570178 , 0.1886747 , 0.3372212 , 0.4723329 , 0.5169277 , 0.5808359 , 0.6644227 , 1.0000000
		};
		number yLineReferenceValue1000Botella[17]={0.0000000 , 0.2807056 , 0.2962703 , 0.3099097 , 0.3330442 , 0.3769189 , 0.3339924 , 0.3253592 , 0.0257995 ,
		 -0.3202137 , -0.4264545 , -0.5264392 , -0.4103754 , -0.3553213 , -0.2936869 , -0.2279225 , 0.0000000};
	switch (Re){
	case 100:
		foundRe=true;
		copyGhiaNumbers(xLineReferenceValue,xLineReferenceValue100);
		copyGhiaNumbers(yLineReferenceValue,yLineReferenceValue100);
		break;
	case 400:
		foundRe=true;
		copyGhiaNumbers(xLineReferenceValue,xLineReferenceValue400);
		copyGhiaNumbers(yLineReferenceValue,yLineReferenceValue400);
		break;
	case 1000:
		foundRe=true;
		copyGhiaNumbers(xLineReferenceValue,xLineReferenceValue1000);
		copyGhiaNumbers(yLineReferenceValue,yLineReferenceValue1000);
		break;
	case 3200:
		foundRe=true;
		copyGhiaNumbers(xLineReferenceValue,xLineReferenceValue3200);
		copyGhiaNumbers(yLineReferenceValue,yLineReferenceValue3200);
		break;
	case 5000:
		foundRe=true;
		copyGhiaNumbers(xLineReferenceValue,xLineReferenceValue5000);
		copyGhiaNumbers(yLineReferenceValue,yLineReferenceValue5000);
		break;
	case 7500:
		foundRe=true;
		copyGhiaNumbers(xLineReferenceValue,xLineReferenceValue7500);
		copyGhiaNumbers(yLineReferenceValue,yLineReferenceValue7500);
		break;
	case 10000:
		foundRe=true;
		copyGhiaNumbers(xLineReferenceValue,xLineReferenceValue10000);
		copyGhiaNumbers(yLineReferenceValue,yLineReferenceValue10000);
		break;
	}

	if (foundRe==false){
		UG_LOG("Reynolds number " << Re << " results not in Ghia data set.\n\n");
		return;
	}

	MathVector<dim> xLineEvalPos;
	MathVector<dim> yLineEvalPos;
	xLineEvalPos[0]=0.5;
	yLineEvalPos[1]=0.5;

	number xLineValue[17];
	number yLineValue[17];

	static const number unhandled = 1537454.54345435;

	for (size_t i=0;i<17;i++){
		xLineValue[i]=unhandled;
		yLineValue[i]=unhandled;
	}

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	//	create Multiindex
	std::vector<MultiIndex<2> > multInd;

	for(int si = 0; si < u.num_subsets(); ++si)
	{
		if (si>0) continue;
		//	get iterators
		ElemIterator iter = u.template begin<elem_type>(si);
		ElemIterator iterEnd = u.template end<elem_type>(si);

		static const number bigNumber = 1e+8;

		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			elem_type* elem = *iter;
			//	get vertices and extract corner coordinates
			const size_t numVertices = elem->num_vertices();
			number xmin=bigNumber,xmax=-bigNumber,ymin=bigNumber,ymax=-bigNumber;
			//UG_LOG("-----------------------------------------------------------------\n");
			// find element bounds
			for(size_t i = 0; i < numVertices; ++i){
				vVrt[i] = elem->vertex(i);
				coCoord[i] = posAcc[vVrt[i]];
				//UG_LOG("co(" << i+1 << ",:)=[" << coCoord[i][0] << "," << coCoord[i][1] << "];" << "\n");
				if (coCoord[i][0]<xmin) xmin=coCoord[i][0];
				if (coCoord[i][0]>xmax) xmax=coCoord[i][0];
				if (coCoord[i][1]<ymin) ymin=coCoord[i][1];
				if (coCoord[i][1]>ymax) ymax=coCoord[i][1];
			};
			//UG_LOG("co(" << numVertices+1 << ",:)=[" << coCoord[0][0] << "," << coCoord[0][1] << "];" << "\n");

			//UG_LOG("bound=[" << xmin << " " << xmax << " " << ymin << " " << ymax << "]\n");

			bool checkx=true;
			bool checky=true;

			if (xmin>0.5) checkx=false;
			else if (xmax<0.5) checkx=false;

			if (ymin>0.5) checky=false;
			else if (ymax<0.5) checky=false;

			if ((checkx==false)&&(checky==false)) continue;

			number interpolation;

			if (checkx==true){
				for (size_t i=0;i<17;i++)
				{
					if ((yco[i]>=ymin)&&(yco[i]<=ymax)){
						xLineEvalPos[1]=yco[i];
						//UG_LOG("p=[" << xLineEvalPos[0] << "," << xLineEvalPos[1] << "]\n");
						if ( ContainsPoint(elem, xLineEvalPos, posAcc) )
						{
							//UG_LOG("***************************************************\n");
							//UG_LOG(xLineEvalPos << "\n");
							if (yco[i]==0){
								//UG_LOG("##################################################################\n");
								//UG_LOG("##################################################################\n");
								//UG_LOG("##################################################################\n");
								//UG_LOG("##################################################################\n");
							}
							//	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
							ReferenceObjectID roid = (ReferenceObjectID) elem->reference_object_id();
							//	get the reference mapping for the element using global corners
							DimReferenceMapping<dim, dim>& mapping
								= ReferenceMappingProvider::get<dim, dim>(roid, coCoord);
							MathVector<dim> localPos;

							//	compute global integration points
							mapping.global_to_local(localPos, xLineEvalPos);

							//	local finite element id
							const LFEID m_id = u.local_finite_element_id(0);

							//	get trial space
							const LocalShapeFunctionSet<dim>& rTrialSpace =
								LocalFiniteElementProvider::get<dim>(roid, m_id);

							//	number of dofs on element
							const size_t num_sh = rTrialSpace.num_sh();

							//	get multiindices of element
							std::vector<MultiIndex<2> > ind;  // 	aux. index array
							u.multi_indices(elem, 0, ind);
							// 	compute approximated solution at integration point
							interpolation = 0.0;
							for(size_t sh = 0; sh < num_sh; ++sh)
							{
								//	get value at shape point (e.g. corner for P1 fct)
								const number valSH = DoFRef(u, ind[sh]);

								//	add shape fct at ip * value at shape
								interpolation += valSH * rTrialSpace.shape(sh, localPos);
							}
							// if there is a previous value, then average
							if (xLineValue[i]!=unhandled) xLineValue[i] = 0.5*(xLineValue[i] + interpolation);
							else xLineValue[i]=interpolation;
						}
					}
				}
			}

			if (checky==true){
				for (size_t i=0;i<17;i++)
				{
					if ((xco[i]>=xmin)&&(xco[i]<=xmax)){
						yLineEvalPos[0]=xco[i];
						if ( ContainsPoint(elem, yLineEvalPos, posAcc) )
						{
							//	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
							ReferenceObjectID roid = (ReferenceObjectID) elem->reference_object_id();
							//	get the reference mapping for the element using global corners
							DimReferenceMapping<dim, dim>& mapping
								= ReferenceMappingProvider::get<dim, dim>(roid, coCoord);
							MathVector<dim> localPos;

							//	compute global integration points
							mapping.global_to_local(localPos, yLineEvalPos);

							//	local finite element id
							const LFEID m_id = u.local_finite_element_id(0);

							//	get trial space
							const LocalShapeFunctionSet<dim>& rTrialSpace =
								LocalFiniteElementProvider::get<dim>(roid, m_id);

							//	number of dofs on element
							const size_t num_sh = rTrialSpace.num_sh();

							//	get multiindices of element
							std::vector<MultiIndex<2> > ind;  // 	aux. index array
							u.multi_indices(elem, 1, ind);
							// 	compute approximated solution at integration point
							interpolation = 0.0;
							for(size_t sh = 0; sh < num_sh; ++sh)
							{
								//	get value at shape point (e.g. corner for P1 fct)
								const number valSH = DoFRef(u, ind[sh]);

								//	add shape fct at ip * value at shape
								interpolation += valSH * rTrialSpace.shape(sh, localPos);
							}
							// if there is a previous value, then average
							if (yLineValue[i]!=unhandled) yLineValue[i] = 0.5*(yLineValue[i] + interpolation);
							else yLineValue[i]=interpolation;
						}
					}
				}
			}
		};// for(  ;iter !=iterEnd; ++iter)
	};// for(int si = 0; si < u.num_subsets(); ++si)
	UG_LOG("\nData evaluation for Re=" << Re << ":\n\n");
	UG_LOG("u values on line through x=0.5:" << "\n\n");
	number maxdiff = 0;
	number diffsum = 0;
	for (size_t i=0;i<17;i++){
		number localdiff=std::abs(xLineReferenceValue[i]-xLineValue[i]);
		UG_LOG("y(" << i+1 << ") = " << yco[i] << "; u(" << i+1 << ") = " << xLineValue[i] << "; u_ghia(" << i+1 << ") = " << xLineReferenceValue[i] << "; uerror_ghia(" << i+1 << ") = " << localdiff << ";\n");
		if (localdiff>maxdiff) maxdiff = localdiff;
		diffsum += localdiff;
	}
	UG_LOG("u max difference: " << maxdiff << "\naverage difference: " << (number)diffsum/17.0);
	UG_LOG("\n\nv values on line through y=0.5:" << "\n\n");
	maxdiff = 0;
	diffsum = 0;
	for (size_t i=0;i<17;i++){
		number localdiff=std::abs(yLineReferenceValue[i]-yLineValue[i]);
		UG_LOG("x(" << i+1 << ") = " << xco[i] << "; v(" << i+1 << ") = " << yLineValue[i] << "; v_ghia(" << i+1 << ") = " << yLineReferenceValue[i] << "; verror_ghia(" << i+1 << ") = " << std::abs(yLineReferenceValue[i]-yLineValue[i]) << ";\n");
		if (localdiff>maxdiff) maxdiff = localdiff;
		diffsum += localdiff;
	}
	UG_LOG("v max difference: " << maxdiff << "\naverage difference: " << (number)diffsum/17.0);
	UG_LOG("\n\n");
	// if Re == 1000 also data from Botella/Peyret paper is available
	if (Re==1000){
		copyGhiaNumbers(xLineReferenceValue,xLineReferenceValue1000Botella);
		copyGhiaNumbers(yLineReferenceValue,yLineReferenceValue1000Botella);
		UG_LOG("Comparison with Botella/Peyret data:\n");
		UG_LOG("\nData evaluation for Re=" << Re << ":\n\n");
		UG_LOG("u values on line through x=0.5:" << "\n\n");
		number maxdiff = 0;
		number diffsum = 0;
		for (size_t i=0;i<17;i++){
			number localdiff=std::abs(xLineReferenceValue[i]-xLineValue[i]);
			UG_LOG("y(" << i+1 << ") = " << yco[i] << "; u(" << i+1 << ") = " << xLineValue[i] << "; u_bot(" << i+1 << ") = " << xLineReferenceValue[i] << "; uerror_bot(" << i+1 << ") = " << localdiff << ";\n");
			if (localdiff>maxdiff) maxdiff = localdiff;
			diffsum += localdiff;
		}
		UG_LOG("u max difference: " << maxdiff << "\nv average difference: " << (number)diffsum/17.0);
		UG_LOG("\n\nv values on line through y=0.5:" << "\n\n");
		maxdiff = 0;
		diffsum = 0;
		for (size_t i=0;i<17;i++){
			number localdiff=std::abs(yLineReferenceValue[i]-yLineValue[i]);
			UG_LOG("x(" << i+1 << ") = " << xco[i] << "; v(" << i+1 << ") = " << yLineValue[i] << "; v_bot(" << i+1 << ") = " << yLineReferenceValue[i] << "; verror_bot(" << i+1 << ") = " << std::abs(yLineReferenceValue[i]-yLineValue[i]) << ";\n");
			if (localdiff>maxdiff) maxdiff = localdiff;
			diffsum += localdiff;
		}
		UG_LOG("v max difference: " << maxdiff << "\nu average difference: " << (number)diffsum/17.0);
		UG_LOG("\n\n");
	}
	// check for Erturk data
	drivenCavityEvaluationErturk(u,Re);
}

/*
// Computation of maximum CFL-number over elements
// General CFL number is given as |u|*deltaT/h where u is velocity, deltaT time step length
// and h mesh width.
// Here following approximation for element local CFL-number is used:
// In element compute interpolated velocity v_bary in barycenter
// go over all "edges" between Crouzeix-Raviart nodes c_i and c_j (middle of edge/face)
// approximate velocity between nodes is 1/|c_i-c_j|*|(c_i-c_j)*v_bary|
// step length is |c_i-c_j| Therefore overall local cfl number is
// 1/|c_i-c_j|^2*|(c_i-c_j)*v_bary|
// Function computes maximum over all elements.
 */
template<typename TGridFunction>
void cflNumber(TGridFunction& u,number deltaT){
	/// domain type
	typedef typename TGridFunction::domain_type domain_type;

	///	world dimension
	static const int dim = domain_type::dim;

	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object elem_type;

	/// side type
	typedef typename elem_type::side side_type;

	///	grid type
	typedef typename domain_type::grid_type grid_type;

	// cfl number
	number maxCfl=0;
	//	get domain
	domain_type& domain = *u.domain().get();

	DimCRFVGeometry<dim> geo;

	//	create Multiindex
	std::vector<MultiIndex<2> > multInd;
	
	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();

	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		//	get iterators
		ElemIterator iter = u.template begin<elem_type>(si);
		ElemIterator iterEnd = u.template end<elem_type>(si);

		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			elem_type* elem = *iter;
			//	get vertices and extract corner coordinates
			const size_t numVertices = elem->num_vertices();
			MathVector<dim> bary,localBary;
			bary = 0;

			for(size_t i = 0; i < numVertices; ++i){
				vVrt[i] = elem->vertex(i);
				coCoord[i] = posAcc[vVrt[i]];
				bary += coCoord[i];
			};
			bary /= numVertices;

			//	reference object id
			ReferenceObjectID roid = elem->reference_object_id();

			//	get trial space
			const LocalShapeFunctionSet<dim>& rTrialSpace =
			LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::CROUZEIX_RAVIART, dim, 1));

			//	get Reference Mapping
			DimReferenceMapping<dim, dim>& map = ReferenceMappingProvider::get<dim, dim>(roid, coCoord);

			map.global_to_local(localBary,bary);

			typename grid_type::template traits<side_type>::secure_container sides;

			UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Only elements of type elem_type are currently supported");

			domain.grid()->associated_elements_sorted(sides, static_cast<elem_type*>(elem) );

			//	memory for shapes
			std::vector<number> vShape;

			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

			rTrialSpace.shapes(vShape, localBary);

			size_t nofsides = geo.num_scv();

			MathVector<dim> baryV;
			baryV = 0;
			for (size_t s=0;s<nofsides;s++){
				MathVector<dim> localbaryV;
				for (int d=0;d<dim;d++){
					u.multi_indices(sides[s], d, multInd);
					localbaryV[d]=DoFRef(u,multInd[0]);
				}
				localbaryV *= vShape[s];
				baryV += localbaryV;
			}
			for (size_t i=0;i<nofsides;i++){
				const typename DimCRFVGeometry<dim>::SCV& scvi = geo.scv(i);
				MathVector<dim> iCoord = scvi.global_ip();
				for (size_t j=i+1;j<nofsides;j++){
					const typename DimCRFVGeometry<dim>::SCV& scvj = geo.scv(j);
					MathVector<dim> jCoord = scvj.global_ip();
					MathVector<dim> subVec;
					VecSubtract(subVec,iCoord,jCoord);
					number localCfl=deltaT*(number)1.0/VecTwoNormSq(subVec)*std::abs(subVec*baryV);
					if (localCfl>maxCfl) maxCfl=localCfl;
				}
			}
		}
	}
	UG_LOG("Max CFL number is " << maxCfl << "\n");
}

// Compute 1/|\Omega| \int_{\Omega} \hat{u_i} \hat{u_i} dx , where \Omega is the computational domain.
// Elementwise approximation of integral by interpolation to barycenter.
template<typename TGridFunction>
void kineticEnergy(TGridFunction& u){
	// total kinetic energy
	number totalE=0;
	// total volume
	number totalVol=0;

	/// domain type
	typedef typename TGridFunction::domain_type domain_type;

	///	world dimension
	static const int dim = domain_type::dim;

	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object elem_type;

	/// side type
	typedef typename elem_type::side side_type;

	///	grid type
	typedef typename domain_type::grid_type grid_type;

	//	get domain of grid function
	domain_type& domain = *u.domain().get();
	DimCRFVGeometry<dim> geo;

	 std::vector<MultiIndex<2> > multInd;

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();

	// assemble deformation tensor fluxes
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		//	get iterators
		ElemIterator iter = u.template begin<elem_type>(si);
		ElemIterator iterEnd = u.template end<elem_type>(si);

		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			elem_type* elem = *iter;
			//	get vertices and extract corner coordinates
			const size_t numVertices = elem->num_vertices();
			MathVector<dim> bary,localBary;
			bary = 0;

			for(size_t i = 0; i < numVertices; ++i){
				vVrt[i] = elem->vertex(i);
				coCoord[i] = posAcc[vVrt[i]];
				bary += coCoord[i];
			};
			bary /= numVertices;

			//	reference object id
			ReferenceObjectID roid = elem->reference_object_id();

			//	get trial space
			const LocalShapeFunctionSet<dim>& rTrialSpace =
			LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::CROUZEIX_RAVIART, dim, 1));

			//	get Reference Mapping
			DimReferenceMapping<dim, dim>& map = ReferenceMappingProvider::get<dim, dim>(roid, coCoord);

			map.global_to_local(localBary,bary);

			typename grid_type::template traits<side_type>::secure_container sides;

			UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Only elements of type elem_type are currently supported");

			domain.grid()->associated_elements_sorted(sides, static_cast<elem_type*>(elem) );

			//	memory for shapes
			std::vector<number> vShape;

			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

			rTrialSpace.shapes(vShape, localBary);

			size_t nofsides = geo.num_scv();

			MathVector<dim> value;
			value = 0;
			number elementVolume = 0;
			for (size_t s=0;s<nofsides;s++){
				const typename DimCRFVGeometry<dim>::SCV& scv = geo.scv(s);
				MathVector<dim> localValue;
				for (int d=0;d<dim;d++){
					u.multi_indices(sides[s], d, multInd);
					localValue[d]=DoFRef(u,multInd[0]);
				}
				localValue *= vShape[s];
				value += localValue;
				elementVolume += scv.volume();
			}
			for (int d=0;d<dim;d++){
				totalE += elementVolume*value[d]*value[d];
			}

			totalVol+=elementVolume;
		}
	}
	// average
	totalE/=(number)totalVol;
	UG_LOG("Total kinetic energy in domain is " << totalE << "\n");
}


} // end namespace ug

#endif /* __H__UG__LIB_DISC__NAVIER_STOKES_TOOLS__ */
