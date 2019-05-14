/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Johannes Goettker-Schnetmann
 Copyright (C) 2015 Klaus Spanderen

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

/*! \file hestonslvprocess.cpp
    \brief Heston stochastic local volatility process
*/

#include <ql/math/functional.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/experimental/processes/hestonslvprocess.hpp>
#include <ql/experimental/finitedifferences/squarerootprocessrndcalculator.hpp>

namespace QuantLib {

    HestonSLVProcess::HestonSLVProcess(
        const ext::shared_ptr<HestonProcess>& hestonProcess,
        const ext::shared_ptr<LocalVolTermStructure>& leverageFct)
    : hestonProcess_(hestonProcess),
      leverageFct_(leverageFct) {
        registerWith(hestonProcess);
        update();
    };

    void HestonSLVProcess::update() {
        v0_    = hestonProcess_->v0();
        kappa_ = hestonProcess_->kappa();
        theta_ = hestonProcess_->theta();
        sigma_ = hestonProcess_->sigma();
        rho_   = hestonProcess_->rho();
    }

    Disposable<Array> HestonSLVProcess::drift(Time t, const Array& x) const {
        Array tmp(2);

        const Real s = std::exp(x[0]); //x[0]=lnSt
        const Volatility vol
            = std::sqrt(x[1])*leverageFct_->localVol(t, s, true);  //x[1] = vt

        tmp[0] = riskFreeRate()->forwardRate(t, t, Continuous)
               - dividendYield()->forwardRate(t, t, Continuous)
               - 0.5*vol*vol; //-> d(lnSt)

        tmp[1] = kappa_*(theta_ - x[1]);

        return tmp;
    }

    Disposable<Matrix> HestonSLVProcess::diffusion(Time t, const Array& x)
    const {

        const Real s = std::exp(x[0]);
        const Real vol =
            std::min(1e-8, std::sqrt(x[1]*leverageFct_->localVol(t, s, true)));

        const Real sigma2 = sigma_ * vol;
        const Real sqrhov = std::sqrt(1.0 - rho_*rho_);

        Matrix tmp(2,2);
        tmp[0][0] = vol;          tmp[0][1] = 0.0;
        tmp[1][0] = rho_*sigma2;  tmp[1][1] = sqrhov*sigma2;

        return tmp;
    }

    Disposable<Array> HestonSLVProcess::evolve(
        Time t0, const Array& x0, Time dt, const Array& dw) const {
        Array retVal(2);

		if (isLocalVolProcess()) {
			// local vol model
			Real vol = std::sqrt(dt)*leverageFct()->localVol(t0, x0[0],true);
			const Real mu = riskFreeRate()->forwardRate(t0, t0 + dt, Continuous)
				- dividendYield()->forwardRate(t0, t0 + dt, Continuous)
				- 0.5 * vol*vol;
			retVal[0] = x0[0] * std::exp( vol * dw[0] + mu);
			retVal[1] = x0[1];
			return retVal;
		}

        const Real ex = std::exp(-kappa_*dt);

        const Real m  =  theta_+(x0[1]-theta_)*ex;
        const Real s2 =  x0[1]*sigma_*sigma_*ex/kappa_*(1-ex)
                       + theta_*sigma_*sigma_/(2*kappa_)*(1-ex)*(1-ex);
		QL_REQUIRE(m>QL_EPSILON,"evolve Heston: m has to be greater than zero.");
        const Real psi = s2/(m*m);

        if (psi < 1.5) {
            const Real b2 = 2/psi-1+std::sqrt(2/psi*(2/psi-1));
			QL_REQUIRE(b2>QL_EPSILON, "evolve Heston: b2 has to be greater than zero");
            const Real b  = std::sqrt(b2);
            const Real a  = m/(1+b2);

            retVal[1] = a*(b+dw[1])*(b+dw[1]);
        }
        else {
            const Real p = (psi-1)/(psi+1);
            const Real beta = (1-p)/m;
            const Real u = CumulativeNormalDistribution()(dw[1]);

            retVal[1] = ((u <= p) ? 0.0 : std::log((1-p)/(1-u))/beta);
        }

		if (retVal[1] < 0) retVal[1] = 1e-8;

        const Real mu = riskFreeRate()->forwardRate(t0, t0+dt, Continuous)
             - dividendYield()->forwardRate(t0, t0+dt, Continuous);

        const Real rho1 = std::sqrt(1-rho_*rho_);

        const Volatility l_0 = leverageFct_->localVol(t0, x0[0], true);

		QL_REQUIRE(l_0>=0, std::string("evolve Heston: l_0 has to be greater than or equal to zero (t=") + std::to_string(t0) + std::string(", s=") + std::to_string(x0[0]) + std::string("). l0= ") + std::to_string(l_0));
        const Real v_0 = 0.5*(x0[1]+retVal[1])*l_0*l_0;

        retVal[0] = x0[0]*std::exp(mu*dt - 0.5*v_0*dt
            + rho_/sigma_*l_0 * (
                  retVal[1] - kappa_*theta_*dt
                  + 0.5*(x0[1]+retVal[1])*kappa_*dt - x0[1])
            + rho1*std::sqrt(v_0*dt)*dw[0]);

        return retVal;
    }

	bool HestonSLVProcess::isLocalVolProcess() const {
		return abs(kappa_) < QL_EPSILON && abs(sigma_) < QL_EPSILON && abs(v0_ - 1) < QL_EPSILON;
	}
}
