/*
 * Copyright (C) 2000, 2001
 * Ferdinando Ametrano, Luigi Ballabio, Adolfo Benin, Marco Marchioro
 * 
 * This file is part of QuantLib.
 * QuantLib is a C++ open source library for financial quantitative
 * analysts and developers --- http://quantlib.sourceforge.net/
 *
 * QuantLib is free software and you are allowed to use, copy, modify, merge,
 * publish, distribute, and/or sell copies of it under the conditions stated 
 * in the QuantLib License.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the license for more details.
 *
 * You should have received a copy of the license along with this file;
 * if not, contact ferdinando@ametrano.net
 *
 * QuantLib license is also available at 
 * http://quantlib.sourceforge.net/LICENSE.TXT
*/
/*! \file MontecarloPricers.i
	
	$Source$
	$Name$
	$Log$
	Revision 1.7  2001/02/26 17:05:30  lballabio
	Ultimate Array interface and typemap for SWIG

	Revision 1.6  2001/02/13 10:06:50  marmar
	Ambiguous variable name underlyingGrowthRate changed in
	unambiguos dividendYield
	
	Revision 1.5  2001/02/07 10:15:57  marmar
	Interface for Himalaya-type option pricer
	
	Revision 1.4  2001/02/05 16:55:12  marmar
	McAsianPricer replaced by AveragePriceAsian and AverageStrikeAsian
	
	Revision 1.3  2001/02/02 10:58:39  marmar
	MonteCarloTools.i added
	
	Revision 1.2  2001/01/15 13:43:20  lballabio
	Using PyArray typemap
	
	Revision 1.1  2001/01/04 17:31:23  marmar
	Alpha version of the Monte Carlo tools.
	
*/

#ifndef quantlib_Montecarlo_Pricers_i
#define quantlib_Montecarlo_Pricers_i

%module MontecarloPricers

%{
#include "quantlib.h"
%}

%include QLArray.i
%include Vectors.i
%include Matrix.i

#if defined(SWIGPYTHON)

%{
#include "standardpathgenerator.h"
using QuantLib::MonteCarlo::StandardPathGenerator;
%}

class StandardPathGenerator{
    public:
	StandardPathGenerator(int dimension, long seed=0);
	Array next() const; // Note that currently Path and Array are equivalent
	double weight() const;
};

%{
#include "mceuropeanpricer.h"
using QuantLib::Pricers::McEuropeanPricer;
%}

class McEuropeanPricer{
    public:
	McEuropeanPricer(OptionType type, double underlying, double strike, 
		 Rate dividendYield,   Rate riskFreeRate, double residualTime, 
		 double volatility,	int timesteps, int confnumber, long seed);
	double value() const;
	double errorEstimate() const;
};

%{
#include "geometricasianoption.h"
using QuantLib::Pricers::GeometricAsianOption;
%}

class GeometricAsianOption {
  public:
	GeometricAsianOption(OptionType type, double underlying, double strike, 
		Rate dividendYield, Rate exerciseRate,
		double residualTime, double volatility);
	double value() const;
};

%{
#include "averagepriceasian.h"
using QuantLib::Pricers::AveragePriceAsian;
%}

class AveragePriceAsian{
    public:
	AveragePriceAsian(OptionType type, double underlying, double strike, 
		 Rate dividendYield,   Rate riskFreeRate, double residualTime, 
		 double volatility,	int timesteps, int confnumber, long seed);
	double value() const;
	double errorEstimate() const;
};
%{
#include "averagestrikeasian.h"
using QuantLib::Pricers::AverageStrikeAsian;
%}

class AverageStrikeAsian{
    public:
	AverageStrikeAsian(OptionType type, double underlying, double strike, 
		 Rate dividendYield,   Rate riskFreeRate, double residualTime, 
		 double volatility,	int timesteps, int confnumber, long seed);
	double value() const;
	double errorEstimate() const;
};


%{
#include "plainbasketoption.h"
using QuantLib::Pricers::PlainBasketOption;
%}

class PlainBasketOption{
    public:
    PlainBasketOption(const Array &underlying, 
        const Array &dividendYield, const Matrix &covariance, 
        Rate riskFreeRate,  double residualTime, 
        int timesteps, long samples, long seed = 0);
	double value() const;
	double errorEstimate() const;
};

%{
#include "himalaya.h"
using QuantLib::Pricers::Himalaya;
%}

class Himalaya{
    public:
    Himalaya(const Array& underlying, const Array& dividendYield, 
        const Matrix &covariance, Rate riskFreeRate, double strike, 
        const DoubleVector &timeDelays, long samples, long seed=0);
	double value() const;
	double errorEstimate() const;
};

#elif defined(SWIGJAVA)
// export relevant functions
#endif

#endif
