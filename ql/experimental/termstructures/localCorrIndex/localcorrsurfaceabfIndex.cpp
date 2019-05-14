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

#include <vector>
#include <ql\errors.hpp>
#include <ql/experimental/templatemodels/auxilliaries/choleskyfactorisationT.hpp>
#include <ql\experimental\termstructures\localCorrIndex\localcorrsurfaceabfIndex.hpp>
#include <ql\experimental\termstructures\Helper\ParticleMethodUtils.hpp>
#include <ql\math\matrixutilities\pseudosqrt.hpp>
#include <ql/experimental/templatemodels/auxilliaries/svdT.hpp>
#include <ql/math/matrixutilities/SymmetricSchurDecomposition.hpp>
#include <ql/experimental/termstructures/localCorrFX/localcorrsurfaceabfFX.hpp>

namespace QuantLib {

    LocalCorrSurfaceABFIndex::LocalCorrSurfaceABFIndex(
		const std::vector<boost::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>>& processes,
		const boost::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>&  		    processToCal,
		const RealStochasticProcess::MatA&												corr0,
		const RealStochasticProcess::MatA&												corr1,
		const RealStochasticProcess::VecA&												indexWeights,
		bool																			possibleNegativeIndex,
		double																		    processToCalBlackVolShift)
    : LocalCorrSurfaceABF(processes, processToCal) {
		corr0_ = RealStochasticProcess::MatA(corr0);
		corr1_ = RealStochasticProcess::MatA(corr1);
		indexWeights_ = RealStochasticProcess::VecA(indexWeights);
		possibleNegativeIndex_ = possibleNegativeIndex;
		processToCalBlackVolShift_ = processToCalBlackVolShift;

		QL_ASSERT(processes.size() == indexWeights.size(), "processes and indexWeights do not fit.");
		QL_ASSERT(processes.size() == corr0.size(), "processes and corr0 do not fit.");
		QL_ASSERT(processes.size() == corr0[0].size(), "processes and corr0 do not fit.");
		QL_ASSERT(processes.size() == corr1.size(), "processes and corr1 do not fit.");
		QL_ASSERT(processes.size() == corr1[0].size(), "processes and corr1 do not fit.");
    }

	LocalCorrSurfaceABFIndex::LocalCorrSurfaceABFIndex(
		const std::vector<boost::shared_ptr<QuantLib::HestonSLVProcess>>&				processes,
		const boost::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>&  		    processToCal,
		const RealStochasticProcess::MatA&												corr0,
		const RealStochasticProcess::MatA&												corr1,
		const RealStochasticProcess::VecA&												 indexWeights,
		bool																			possibleNegativeIndex,
		double																		    processToCalBlackVolShift)
		: LocalCorrSurfaceABF(processes, processToCal) {


		checkHestonConsistence(corr0);
		checkHestonConsistence(corr1);

		corr0_ = RealStochasticProcess::MatA(corr0);
		corr1_ = RealStochasticProcess::MatA(corr1);
		indexWeights_ = RealStochasticProcess::VecA(indexWeights);
		possibleNegativeIndex_ = possibleNegativeIndex;
		processToCalBlackVolShift_ = processToCalBlackVolShift;

		QL_ASSERT(processes.size() == indexWeights.size(), "processes and indexWeights do not fit.");
		QL_ASSERT(2*processes.size() == corr0.size(), "processes and corr0 do not fit.");
		QL_ASSERT(2*processes.size() == corr0[0].size(), "processes and corr0 do not fit for Heston.");
		QL_ASSERT(2*processes.size() == corr1.size(), "processes and corr1 do not fit for Heston.");
		QL_ASSERT(2*processes.size() == corr1[0].size(), "processes and corr1 do not fit.");

		projectSymmetricToCorrelation(corr0_);
		projectSymmetricToCorrelation(corr1_);
		
	}


	QuantLib::Real LocalCorrSurfaceABFIndex::localFStrike(Time t, const RealStochasticProcess::VecA& X0) {
		double weightedSum = 0;
		for (size_t i = 0; i < indexWeights_.size(); i++)
		{
			weightedSum += processes_[i]->s0()->value() * std::exp(X0[i]) * indexWeights_[i];
		}
		return weightedSum;
	}
	
	void LocalCorrSurfaceABFIndex::accept(AcyclicVisitor& v) {
		Visitor<LocalCorrSurfaceABFIndex>* v1 =
			dynamic_cast<Visitor<LocalCorrSurfaceABFIndex>*>(&v);
		if (v1 != 0)
			v1->visit(*this);
		else
			LocalCorrSurfaceABF::accept(v);
	}

	Real LocalCorrSurfaceABFIndex::localCorrImplTeq0(Time t, const RealStochasticProcess::VecA& X0, bool extrapolate) {
		
		//smiled surface will through an error, therefore assume one minute ahead
		t = 1.0 / (365 * 24 * 60);

		RealStochasticProcess::VecA s0 = RealStochasticProcess::VecA(processes_.size());
		RealStochasticProcess::VecA vol = RealStochasticProcess::VecA(processes_.size());

		for (size_t i = 0; i < s0.size(); i++)
		{
			s0[i] = processes_[i]->s0()->value() * std::exp(X0[i]);
			vol[i] = processes_[i]->leverageFct()->localVol(t, s0[i], extrapolate);

			if (!processes_[i]->isLocalVolProcess()) vol[i] *= (X0[i+ processes_.size()] <= 0 ? 0.001 : std::sqrt(X0[i+ processes_.size()]));
			if (vol[i] != vol[i]) QL_FAIL(std::string("leverage function of asset ") + std::to_string(i) + std::string(" does have non-real values"));
		}

		Real cov1 = getIndexCovariance(CTSIndexCovarianceType::CORR1, s0, vol);
		Real cov0 = getIndexCovariance(CTSIndexCovarianceType::CORR0, s0, vol);

		Real indexValue = localFStrike(t,X0);
		if(possibleNegativeIndex_){ 
			Real volIndex = ParticleMethodUtils::getLocalVolFromPriceFormula(processToCal_,indexValue,t,processToCalBlackVolShift_);
			return (volIndex - cov0) / (cov1 - cov0);
		}
		else {
			Real volIndex = processToCal_->localVolatility()->localVol(t, indexValue, extrapolate);
			return (indexValue*indexValue*volIndex*volIndex - cov0) / (cov1 - cov0);
		}
	}
	

	Matrix LocalCorrSurfaceABFIndex::getLocalCorrelationSurface2dim(Time t,
		std::vector<Real> assetGrid1, std::vector<Real> assetGrid2) {

		QL_ASSERT(corr0_.size() == 2, "function only works for index with two assets");

		Matrix result(assetGrid1.size(), assetGrid2.size());
		std::vector<std::vector<Real>> corrM;
		if(processes_[0]->isLocalVolProcess())
			corrM = std::vector<std::vector<Real>>(2, std::vector<Real>(2));
		else
			corrM = std::vector<std::vector<Real>>(4, std::vector<Real>(4));

		std::vector<Real> x0(2);

		for (size_t i = 0; i < assetGrid1.size(); i++)
		{
			for (size_t j = 0; j < assetGrid2.size(); j++)
			{
				x0[0] = log(assetGrid1[i] / processes_[0]->s0()->value());
				x0[1] = log(assetGrid2[j] / processes_[1]->s0()->value());

				localCorr(corrM, t, x0, true);
				result[i][j] = corrM[0][1];
			}
		}
		return result;
	}

	QuantLib::Real LocalCorrSurfaceABFIndex::checkLambdaValue(QuantLib::Real lambda) {
	
		if (lambda != lambda)
			QL_FAIL("lambda is erroneous.");

		if (lambda > 1) return 1;
		if (lambda < 0) return 0;
		return lambda;
	
	}

	Real LocalCorrSurfaceABFIndex::getIndexCovariance(CTSIndexCovarianceType type, RealStochasticProcess::VecA& s0, RealStochasticProcess::VecA& vol) {
		switch (type)
		{
		case QuantLib::LocalCorrSurfaceABFIndex::CORR0:
			return getIndexCovariance(corr0_,s0,vol);
			break;
		case QuantLib::LocalCorrSurfaceABFIndex::CORR1:
			return getIndexCovariance(corr1_, s0, vol);
			break;
		default:
			QL_FAIL("given CTSIndexCovarianceType not considered.");
			break;
		}
	}


	Real LocalCorrSurfaceABFIndex::getIndexCovariance(RealStochasticProcess::MatA& corrMatrix, RealStochasticProcess::VecA& s0, RealStochasticProcess::VecA& vol) {
		Real cov = 0;

		for (size_t i = 0; i < corrMatrix.size(); i++)
		{
			for (size_t j = 0; j < corrMatrix.size(); j++)
			{
				cov += corrMatrix[i][j] * indexWeights_[i] * indexWeights_[j] * vol[i] * vol[j] * s0[i] * s0[j];
			}
		}
		
		return cov;
	}

	void LocalCorrSurfaceABFIndex::checkHestonConsistence(const RealStochasticProcess::MatA& corrMatrix) {
		for (size_t i = 0; i < 2 * processes_.size(); i++)
		{
			for (size_t j = 0; j < 2 * processes_.size(); j++)
			{
				if (i == j) {
					QL_ASSERT(corrMatrix[i][j] == 1, "Error: Correlation matrix has ones on diagonals.");
				}
				else if (i == j + processes_.size() || i + processes_.size() == j) {
					int assetIndex = std::min(i, j);
					QL_ASSERT(abs(corrMatrix[i][j] - processes_[assetIndex]->rho())<QL_EPSILON, std::string( "correlation of asset " ) + std::to_string(assetIndex) + std::string(" does not fit to rho of heston process."));
				}
			}
		}
	}

}

