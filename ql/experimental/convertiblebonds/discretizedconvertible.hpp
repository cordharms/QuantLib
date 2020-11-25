/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2005, 2006 Theo Boafo
 Copyright (C) 2006, 2007 StatPro Italia srl

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

/*! \file discretizedconvertible.hpp
    \brief discretized convertible
*/

#ifndef quantlib_discretized_convertible_hpp
#define quantlib_discretized_convertible_hpp

#include <ql/discretizedasset.hpp>
#include <ql/experimental/convertiblebonds/convertiblebond.hpp>
#include <ql/processes/blackscholesprocess.hpp>

namespace QuantLib {

    class DiscretizedConvertible : public DiscretizedAsset {
      public:
        DiscretizedConvertible(
             const ConvertibleBond::option::arguments&,
             const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
             const TimeGrid& grid = TimeGrid());

        void reset(Size size);

        const Array& conversionProbability() const {
            return conversionProbability_;
        }
        Array& conversionProbability() { return conversionProbability_; }

        const Array& spreadAdjustedRate() const { return spreadAdjustedRate_; }
        Array& spreadAdjustedRate() { return spreadAdjustedRate_; }

        const Array& dividendValues() const { return dividendValues_; }
        Array& dividendValues() { return dividendValues_; }

        std::vector<Time> mandatoryTimes() const {
            std::vector<Time> result;
            std::copy(stoppingTimes_.begin(), stoppingTimes_.end(),
                      std::back_inserter(result));
            std::copy(callabilityTimes_.begin(), callabilityTimes_.end(),
                      std::back_inserter(result));
            std::copy(couponTimes_.begin(), couponTimes_.end(),
                      std::back_inserter(result));
            return result;
        }

      protected:
        void postAdjustValuesImpl();
        Array conversionProbability_, spreadAdjustedRate_, dividendValues_;
		Disposable<Array> adjustedGrid() const;
		ConvertibleBond::option::arguments arguments_;
		std::vector<Time> stoppingTimes_;
		std::vector<Time> callabilityTimes_;
		std::vector<Time> couponTimes_;
		ext::shared_ptr<GeneralizedBlackScholesProcess> process_;

      private:
		void addCoupon(Size);
        void applyConvertibility();
        void applyCallability(Size, bool convertible);
        std::vector<Time> dividendTimes_;
    };

	class DiscretizedCoCo : public DiscretizedConvertible {
	public:
		DiscretizedCoCo(
			const ContingentConvertible::optionCoCo::arguments&,
			const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
			const TimeGrid& grid = TimeGrid());
		void reset(Size size);
	protected:
		ContingentConvertible::optionCoCo::arguments arguments_;
		void postAdjustValuesImpl();
	private:
		void addCoupon(Size);
		void computeDefaultedCoupon();
		void applyConvertibility();
		void applyCallability(Size, bool convertible);
		Real defaultedCoupon_;
	};

}


#endif

