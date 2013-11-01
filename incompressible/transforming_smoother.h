/*
 * transforming_smoother.h
 *
 *  Created on: 01.10.2013
 *      Author: andreasvogel
 */

#ifndef __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__TRANSFORMING_SMOOTHER__
#define __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__TRANSFORMING_SMOOTHER__

#include "lib_algebra/operator/interface/linear_iterator.h"
#include "lib_algebra/operator/interface/matrix_operator.h"
#include "lib_algebra/operator/debug_writer.h"


namespace ug{


template <typename TDomain, typename TAlgebra>
class AssembledTransformingSmoother :
	public ILinearIterator<typename TAlgebra::vector_type>,
	public DebugWritingObject<TAlgebra>
{
	public:
	///	Domain
		typedef TDomain domain_type;

	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Vector type
		typedef typename algebra_type::vector_type vector_type;

	///	Matrix type
		typedef typename algebra_type::matrix_type matrix_type;

		using DebugWritingObject<TAlgebra>::write_debug;
		using ILinearIterator<typename TAlgebra::vector_type>::damping;

	public:
	/// constructor setting approximation space
		AssembledTransformingSmoother(SmartPtr<IAssemble<TAlgebra> > spRightTrafoAss,
		                              SmartPtr<IAssemble<TAlgebra> > spTrafoSystemAss,
		                              SmartPtr<ILinearIterator<vector_type> > spSmoother)
			: m_spRightTrafoAss(spRightTrafoAss),
			  m_spRightTrafoMat(new MatrixOperator<matrix_type, vector_type>),
			  m_spTrafoSystemAss(spTrafoSystemAss),
			  m_spTrafoSystemMat(new MatrixOperator<matrix_type, vector_type>),
			  m_spSmoother(spSmoother)
		{}

	///	name
		virtual const char* name() const {return "Assembled Transform Smoother";}

	/// Prepare for Operator J(u) and linearization point u (current solution)
		virtual bool init(SmartPtr<ILinearOperator<vector_type> > J, const vector_type& u);

	///	Prepare for Linear Operartor L
		virtual bool init(SmartPtr<ILinearOperator<vector_type> > L);

	///	Compute new correction c = B*d
		virtual bool apply(vector_type& c, const vector_type& d);

	///	Compute new correction c = B*d and return new defect d := d - A*c
		virtual bool apply_update_defect(vector_type& c, vector_type& d);

	///	Clone
		SmartPtr<ILinearIterator<vector_type> > clone();

	protected:
	///	matrix for original system
		SmartPtr<ILinearOperator<vector_type> > m_spOriginalSystemMat;

	///	assembling for right-transformation
		SmartPtr<IAssemble<TAlgebra> > m_spRightTrafoAss;

	///	matrix operator for right-transformation
		SmartPtr<MatrixOperator<matrix_type, vector_type> > m_spRightTrafoMat;

	///	assembling for transformed system
		SmartPtr<IAssemble<TAlgebra> > m_spTrafoSystemAss;

	///	matrix operator for transformed system
		SmartPtr<MatrixOperator<matrix_type, vector_type> > m_spTrafoSystemMat;

	///	smoother used on transformed system
		SmartPtr<ILinearIterator<vector_type> > m_spSmoother;

};

template <typename TDomain, typename TAlgebra>
bool AssembledTransformingSmoother<TDomain, TAlgebra>::
apply(vector_type &c, const vector_type& d)
{
//	temporary vector for defect
	vector_type dTmp; dTmp.resize(d.size());

//	copy defect
	dTmp = d;

//	work on copy
	return apply_update_defect(c, dTmp);
}

template <typename TDomain, typename TAlgebra>
bool AssembledTransformingSmoother<TDomain, TAlgebra>::
apply_update_defect(vector_type &c, vector_type& d)
{
	SmartPtr<vector_type> y = c.clone_without_values();

//	apply smoother of transformed system
	if(!m_spSmoother->apply(*y, d))
	{
		UG_LOG("AssembledTransformingSmoother: Smoother applied incorrectly.\n");
		return false;
	}

//	apply right-transformation
	m_spRightTrafoMat->apply(c, *y);

#ifdef UG_PARALLEL
	c.change_storage_type(PST_CONSISTENT);
#endif

	const number damp = this->damping()->damping();
	c *= damp;

//	compute updated defect
	m_spOriginalSystemMat->apply_sub(d, c);

	return true;
}

template <typename TDomain, typename TAlgebra>
bool AssembledTransformingSmoother<TDomain, TAlgebra>::
init(SmartPtr<ILinearOperator<vector_type> > J, const vector_type& u)
{
	m_spOriginalSystemMat = J;

	GridLevel gridLevel;

	const GridFunction<TDomain, TAlgebra>* pSol =
		dynamic_cast<const GridFunction<TDomain, TAlgebra>*>(&u);
	if(pSol){
		gridLevel = pSol->dof_distribution()->grid_level();
	}

	try{
		m_spRightTrafoAss->assemble_jacobian(*m_spRightTrafoMat, u, gridLevel);
	}
	UG_CATCH_THROW("AssembledTransformingSmoother: "
					" Cannot assemble right-transformation matrix.");

	write_debug(*m_spRightTrafoMat, "RightTrafo");

	try{
		m_spTrafoSystemAss->assemble_jacobian(*m_spTrafoSystemMat, u, gridLevel);
	}
	UG_CATCH_THROW("AssembledTransformingSmoother: "
					" Cannot assemble transformed system matrix.");

	write_debug(*m_spTrafoSystemMat, "TrafoSystem");

	if(!m_spSmoother->init(m_spTrafoSystemMat, u))
		UG_THROW("AssembledTransformingSmoother: "
				" Cannot init smoother for transformed system matrix.");

	return true;
}

template <typename TDomain, typename TAlgebra>
bool AssembledTransformingSmoother<TDomain, TAlgebra>::
init(SmartPtr<ILinearOperator<vector_type> > L)
{
	UG_THROW("Not implemented.")
}

template <typename TDomain, typename TAlgebra>
SmartPtr<ILinearIterator<typename TAlgebra::vector_type> >
AssembledTransformingSmoother<TDomain, TAlgebra>::
clone()
{
	SmartPtr<AssembledTransformingSmoother<TDomain, TAlgebra> > clone(
		new AssembledTransformingSmoother<TDomain, TAlgebra>(m_spRightTrafoAss, m_spTrafoSystemAss, m_spSmoother));

	return clone;
}


} // end namespace ug

#endif /* __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__TRANSFORMING_SMOOTHER__ */
