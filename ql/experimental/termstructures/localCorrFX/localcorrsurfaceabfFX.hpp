/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C)  2017 Cord Harms

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file localcorrsurface.hpp
    \brief Local Correlation surface derived.
	J. Guyon, 2013, A new Class of local correlation models
*/

#ifndef quantlib_localcorrsurfaceabffx_hpp
#define quantlib_localcorrsurfaceabffx_hpp

#include <ql/experimental/termstructures/localcorrsurfaceabf.hpp>

#include <ql\math\matrixutilities\pseudosqrt.hpp>
#include <ql/experimental/templatemodels/auxilliaries/svdT.hpp>
#include <ql/math/matrixutilities/SymmetricSchurDecomposition.hpp>
#include <ql/experimental/templatemodels/auxilliaries/choleskyfactorisationT.hpp>

namespace QuantLib {

    //! Local Correlation surface derived 
    /*! For details about this implementation refer to
        
        \bug this class is untested, probably unreliable.
    */
    class LocalCorrSurfaceABFFX : public LocalCorrSurfaceABF {
      public:
        LocalCorrSurfaceABFFX(const std::vector<boost::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>>& processes,
							  const boost::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>&			  processToCal);
		LocalCorrSurfaceABFFX(const std::vector<boost::shared_ptr<QuantLib::HestonSLVProcess>>& processes,
			const boost::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>&			  processToCal,
			const RealStochasticProcess::MatA											  correlations);

		//@}
		//! \name Visitability
		//@{
		virtual void accept(AcyclicVisitor&);
        //@}
		virtual QuantLib::Real localA(Time t, const RealStochasticProcess::VecA& assets,
			bool extrapolate = false) const = 0;
		virtual QuantLib::Real localB(Time t, const RealStochasticProcess::VecA& assets,
			bool extrapolate = false) const = 0;

		Matrix getLocalCorrelationSurface(Time t, std::vector<Real> assetGrid1, std::vector<Real> assetGrid2);
		virtual QuantLib::Real localFStrike(Time t, const RealStochasticProcess::VecA& X0);

      protected:
		  virtual Real localCorrImplTeq0(Time t, const RealStochasticProcess::VecA& X0, bool extrapolate = false);
		  virtual QuantLib::Real checkLambdaValue(QuantLib::Real lambda);
	  private:
		  RealStochasticProcess::MatA getPureHestonImpliedCorrelationMatrix();
    };

	namespace
	{

		// Matrix infinity norm. See Golub and van Loan (2.3.10) or
		// <http://en.wikipedia.org/wiki/Matrix_norm>
		Real normInf(const Matrix& M) {
			Size rows = M.rows();
			Size cols = M.columns();
			Real norm = 0.0;
			for (Size i = 0; i<rows; ++i) {
				Real colSum = 0.0;
				for (Size j = 0; j<cols; ++j)
					colSum += std::fabs(M[i][j]);
				norm = std::max(norm, colSum);
			}
			return norm;
		}

		//J. Higham, Computating the nearest correlation matrix - a problem from finance
		void projectSymmetricToCorrelation(RealStochasticProcess::MatA& Y)
		{
			//Set W:=I and compute Algorithm 3.3

			RealStochasticProcess::MatA S = RealStochasticProcess::MatA(Y);
			Matrix R = Matrix(Y.size(), Y.size());
			Matrix X = Matrix(Y.size(), Y.size());
			Matrix Xp = Matrix(Y.size(), Y.size());
			Matrix SV = Matrix(Y.size(), Y.size(), 0.0);

			int cc = 0;
			Real tol = 10e-8;
			Real err = QL_MAX_REAL;

			while (err > tol)
			{
				for (size_t i = 0; i < Y.size(); i++)
				{
					for (size_t j = 0; j < Y.size(); j++)
					{
						if (cc == 0) S[i][j] = 0; //initialize
						R[i][j] = Y[i][j] - S[i][j]; //Dykstra's correction
						Xp[i][j] = X[i][j];
					}
				}

				//projectino onto S
				//SVD tmp = SVD(R);

				SymmetricSchurDecomposition dec = SymmetricSchurDecomposition(R);

				for (size_t i = 0; i < Y.size(); i++)
				{
					SV[i][i] = std::max(dec.eigenvalues()[i], 0.0);
				}

				X = dec.eigenvectors()*SV*transpose(dec.eigenvectors());

				for (size_t i = 0; i < Y.size(); i++)
				{
					for (size_t j = 0; j < Y.size(); j++)
					{
						QL_ASSERT(abs(X[i][j] - X[j][i]) < 10e-8, "X not symmetric.");
						//X[i][j] = X[j][i]; //we need QL_EPSILON accuracy
						Y[i][j] = i == j ? 1 : X[i][j]; //projektion onto U
						S[i][j] = X[i][j] - R[i][j];
					}
				}
				cc++;
				QL_ASSERT(cc < 100, "Correlation matrix projection does not converge.");
				err = normInf(X - Xp) / normInf(X);
			}
			//check
			TemplateAuxilliaries::performCholesky(RealStochasticProcess::MatA(Y), Y.size(), true);
		}

	}


}

#endif
