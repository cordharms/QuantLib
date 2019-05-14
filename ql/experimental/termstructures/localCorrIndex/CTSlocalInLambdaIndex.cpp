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

#include <ql\experimental\termstructures\localCorrIndex\CTSlocalInLambdaIndex.hpp>


namespace QuantLib {

	CTSlocalInLambdaIndex::CTSlocalInLambdaIndex(
		const std::vector<boost::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>>& processes,
		const boost::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>&			    processToCal,
		const RealStochasticProcess::MatA&												corr0,
		const RealStochasticProcess::MatA&												corr1,
		const RealStochasticProcess::VecA&												indexWeights,
		bool																			possibleNegativeIndex,
		double																			processToCalBlackVolShift)
    : LocalCorrSurfaceABFIndex(processes,processToCal,corr0,corr1,indexWeights,possibleNegativeIndex,processToCalBlackVolShift){
		//initializeF();
		//setInterpolation<Linear>();
	}

	CTSlocalInLambdaIndex::CTSlocalInLambdaIndex(
		const std::vector<boost::shared_ptr<QuantLib::HestonSLVProcess>>&				processes,
		const boost::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>&			    processToCal,
		const RealStochasticProcess::MatA&												corr0,
		const RealStochasticProcess::MatA&												corr1,
		const RealStochasticProcess::VecA&												indexWeights,
		bool																			possibleNegativeIndex,
		double																			processToCalBlackVolShift)
		: LocalCorrSurfaceABFIndex(processes, processToCal, corr0, corr1, indexWeights, possibleNegativeIndex, processToCalBlackVolShift) {
		//initializeF();
		//setInterpolation<Linear>();
	}

	QuantLib::Real CTSlocalInLambdaIndex::localA(Time t, const RealStochasticProcess::VecA& assets,
		bool extrapolate) const {
		return 0;
	}

	QuantLib::Real CTSlocalInLambdaIndex::localB(Time t, const RealStochasticProcess::VecA& assets,
		bool extrapolate) const {
		return 1;
	}
	void CTSlocalInLambdaIndex::accept(AcyclicVisitor& v) {
		Visitor<CTSlocalInLambdaIndex>* v1 =
			dynamic_cast<Visitor<CTSlocalInLambdaIndex>*>(&v);
		if (v1 != 0)
			v1->visit(*this);
		else
			LocalCorrSurfaceABFIndex::accept(v);
	}


	//void CTSlocalInCrossCorrelationFX::initializeF() {
	//	calibratorLocalCorr_->calibrateFX(strikes_, times_, surfaceF_, processes_, processToCal_);
	//	interpolatorStrikesF_.resize(times_.size());
	//	valuesSecInt_.resize(times_.size());
	//}
}

