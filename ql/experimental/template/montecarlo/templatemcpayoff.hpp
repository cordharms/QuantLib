/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015, Sebastian Schlenkrich

*/

/*! \file templatemcpayoff.hpp
    \brief generic payoff interface for MC simulation
	
*/


#ifndef quantlib_templatemcpayoff_hpp
#define quantlib_templatemcpayoff_hpp


#include <ql/experimental/template/montecarlo/templatemcsimulation.hpp>



namespace QuantLib {

	// Base class for template payoffs
    template <class DateType, class PassiveType, class ActiveType>
	class TemplateMCPayoff {
	protected:
		typedef TemplateMCSimulation<DateType, PassiveType, ActiveType>        SimulationType;
		typedef typename TemplateMCSimulation<DateType, PassiveType, ActiveType>::Path  PathType;

	    DateType observationTime_;
	public:
	    TemplateMCPayoff( const DateType observationTime ) : observationTime_(observationTime) { }
		// inspectors
		inline DateType observationTime() { return observationTime_; }
		// generic payoff(observationTime, p) needs to be implemented by derived classes
        inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) = 0;
		// discounted payoff for NPV valuation
		inline virtual ActiveType discountedAt(const boost::shared_ptr<PathType>& p) { return at(p) / p->numeraire(observationTime_); }

		// generic pricer
		class Pricer {
		protected:
			std::vector< boost::shared_ptr<TemplateMCPayoff> > payoffs_;
			boost::shared_ptr<SimulationType>                  simulation_;
		public:

			inline static ActiveType NPV(const std::vector< boost::shared_ptr<TemplateMCPayoff> >& payoffs, 
				                         const boost::shared_ptr<SimulationType>&                  simulation) {
				ActiveType npv = 0.0;
				for (size_t k=0; k<payoffs.size(); ++k) {
					ActiveType npvk = 0;
					for (size_t n=0; n<simulation->nPaths(); ++n) {
						npvk += payoffs[k]->discountedAt(simulation->path(n));
					}
					npv += npvk;
				}
				return npv / simulation->nPaths();
			}


			Pricer (const std::vector< boost::shared_ptr<TemplateMCPayoff> >& payoffs, 
				    const boost::shared_ptr<SimulationType>&                  simulation)
					: payoffs_(payoffs), simulation_(simulation) { }

			inline ActiveType NPV() {
				/*
				ActiveType npv = 0.0;
				for (size_t k=0; k<payoffs_.size(); ++k) {
					ActiveType npvk = 0;
					for (size_t n=0; n<simulation_->nPaths(); ++n) {
						npvk += payoffs_[k]->discountedAt(simulation_->path(n));
					}
					npv += npvk;
				}
				return npv / simulation_->nPaths();
				*/
				return NPV(payoffs_,simulation_);
			}
		};


		// particular payoffs

		// simple cash payment
		class Cash : public TemplateMCPayoff {
		protected:
			DateType payTime_;
		public:
			Cash( DateType obsTime, DateType payTime ) : TemplateMCPayoff(obsTime), payTime_(payTime) { }
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				//if (payTime_<=observationTime()) return (ActiveType)1.0;
				return p->zeroBond(observationTime(),payTime_);
			}
		};

		// call or put exercised at observation time and settled at pay time
		class VanillaOption : public TemplateMCPayoff {
		protected:
			DateType    payTime_;
			PassiveType callOrPut_;
			PassiveType strike_;
		public:
			VanillaOption( DateType obsTime, DateType payTime, PassiveType strike, PassiveType callOrPut ) : TemplateMCPayoff(obsTime), payTime_(payTime), strike_(strike), callOrPut_(callOrPut) { }
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				if (payTime_<observationTime()) return (ActiveType)0.0;
				ActiveType DF = p->zeroBond(observationTime(),payTime_);
				ActiveType S  = p->asset(observationTime());
				ActiveType V  = callOrPut_ * DF * (S - strike_);
				return (V>0.0) ? (V) : ((ActiveType)0.0);
			}
		};


		// annuity
		class Annuity : public TemplateMCPayoff {
		protected:
			std::vector<DateType>    payTimes_;
			std::vector<PassiveType> payWeights_;  // these are typically year fractions
		public:
			Annuity( DateType                        obsTime, 
				     const std::vector<DateType>&    payTimes,
				     const std::vector<PassiveType>& payWeights)
				: TemplateMCPayoff(obsTime), payTimes_(payTimes), payWeights_(payWeights) { }
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				ActiveType res = 0.0;
				size_t N = (payTimes_.size()<payWeights_.size()) ? (payTimes_.size()) : (payWeights_.size());
				for (size_t k=0; k<N; ++k) {
					if (payTimes_[k]>observationTime()) {
						res += payWeights_[k] * p->zeroBond(observationTime(),payTimes_[k]);
					}
				}
				return res;
			}
		};

		// prototypical physically settled European swaption
		class ModelSwaption : public TemplateMCPayoff {
		protected:
			std::vector<DateType>    times_;        // T_0, .., T_N
			std::vector<PassiveType> payWeights_;   // tau_0, .., tau_N-1
			PassiveType              strikeRate_;   // option strike
			PassiveType              payOrRec_;     // call (+1) or put (-1) option on swap rate
			bool                     isConsistent_; // check consistency of input
		public:
			ModelSwaption( DateType                        obsTime, 
				           const std::vector<DateType>&    times,
				           const std::vector<PassiveType>& payWeights,
					       PassiveType                     strikeRate,
					       PassiveType                     payOrRec      )
				: TemplateMCPayoff(obsTime), times_(times), payWeights_(payWeights), strikeRate_(strikeRate), payOrRec_(payOrRec) { 
				isConsistent_ = true;
				if (times_.size()<2) isConsistent_ = false;
				for (size_t k=0; k<times_.size(); ++k) if (times_[k]<observationTime()) isConsistent_ = false;
				// default weights
				if (payWeights_.size()!=times_.size()-1) {
					payWeights_.resize(times_.size()-1);
					for (size_t k=0; k<payWeights_.size(); ++k) payWeights_[k] = times_[k+1] - times_[k];
				}
				// finished
			}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				ActiveType res = 0.0;
				if (!isConsistent_) return res;
				// annuity...
				for (size_t k=0; k<payWeights_.size(); ++k) res += payWeights_[k] * p->zeroBond(observationTime(),times_[k+1]);
				// floatleg - fixedleg...
				res = p->zeroBond(observationTime(),times_[0]) - p->zeroBond(observationTime(),times_[times_.size()-1]) - strikeRate_*res;
				// payer or receiver swap...
				res *= payOrRec_;
				// exercise option...
				res = (res>0) ? res : 0.0;
				return res;
			}
		};

		// prototypical physically settled European swaption
		class GeneralSwaption : public TemplateMCPayoff {
		protected:
			std::vector<DateType>    floatTimes_;     // T_1, .., T_M
			std::vector<PassiveType> floatWeights_;   // u_1, .., u_M
			std::vector<DateType>    fixedTimes_;     // T_1, .., T_N
			std::vector<PassiveType> fixedWeights_;   // w_1, .., w_N
			PassiveType              strikeRate_;   // option strike
			PassiveType              payOrRec_;     // call (+1) or put (-1) option on swap rate
		public:
			GeneralSwaption( DateType                        obsTime, 
				             const std::vector<DateType>&    floatTimes,
				             const std::vector<PassiveType>& floatWeights,
				             const std::vector<DateType>&    fixedTimes,
				             const std::vector<PassiveType>& fixedWeights,
					         PassiveType                     strikeRate,
					         PassiveType                     payOrRec      )
				: TemplateMCPayoff(obsTime),  floatTimes_(floatTimes), floatWeights_(floatWeights),
				  fixedTimes_(fixedTimes), fixedWeights_(fixedWeights), strikeRate_(strikeRate), payOrRec_(payOrRec) { 
			    // check consistency of swap
			    // float leg
			    QL_REQUIRE(floatWeights.size()>0,"TemplateQGSwaptionModel: empty float weights.");
			    QL_REQUIRE(floatTimes.size()==floatWeights.size(),"TemplateQGSwaptionModel: float sizes mismatch.");
			    QL_REQUIRE(floatTimes[0]>0,"TemplateQGSwaptionModel: future float times required");
			    for (size_t k=1; k<floatTimes.size(); ++k) QL_REQUIRE(floatTimes[k]>floatTimes[k-1],"TemplateQGSwaptionModel: ascending float times required");
			    // fixed leg
			    QL_REQUIRE(fixedWeights.size()>0,"TemplateQGSwaptionModel: empty fixed weights.");
			    QL_REQUIRE(fixedTimes.size()==fixedWeights.size(),"TemplateQGSwaptionModel: fixed sizes mismatch.");
			    QL_REQUIRE(fixedTimes[0]>0,"TemplateQGSwaptionModel: future fixed times required");
			    for (size_t k=1; k<fixedTimes.size(); ++k) QL_REQUIRE(fixedTimes[k]>fixedTimes[k-1],"TemplateQGSwaptionModel: ascending fixed times required");
				// finished
			}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				ActiveType floatleg = 0.0;
				ActiveType annuity  = 0.0;
				// float leg
				for (size_t k=0; k<floatTimes_.size(); ++k) floatleg += floatWeights_[k] * p->zeroBond(observationTime(),floatTimes_[k]);
				// annuity
				for (size_t k=0; k<fixedTimes_.size(); ++k) annuity  += fixedWeights_[k] * p->zeroBond(observationTime(),fixedTimes_[k]);
				// floatleg - fixedleg...
				ActiveType res = floatleg - strikeRate_*annuity;
				// payer or receiver swap...
				res *= payOrRec_;
				// exercise option...
				res = (res>0) ? res : 0.0;
				return res;
			}
		};

		// undiscounted correlation between prototypical physically settled European swaption
		class ModelCorrelation : public TemplateMCPayoff {
		protected:
			std::vector<DateType> times_;
			DateType T1_, T2_;
			ActiveType swapRate(const boost::shared_ptr<PathType>& p, const DateType t, const DateType TN ) {
				ActiveType num = p->zeroBond(t,t) - p->zeroBond(t,TN);
				ActiveType den = 0.0;
				for (ActiveType Ti = t; Ti<TN; Ti+=1.0) {
					ActiveType T = (Ti+1.0>TN) ? (TN) : (Ti+1.0);
					den += (T-Ti) * p->zeroBond(t,T);
				}
				return num / den;
			}
		public:
			ModelCorrelation( const std::vector<DateType>&    times,   // observation times
				              const DateType                  T1,      // swap term one
							  const DateType                  T2 )     // swap term two
							  : TemplateMCPayoff(0.0), times_(times), T1_(T1), T2_(T2) {
				QL_REQUIRE(times_.size()>1,"ModelCorrelation: At least two observation times required.");
			}
		    // payoff should NOT be discounted
		    inline virtual ActiveType discountedAt(const boost::shared_ptr<PathType>& p) { return at(p); }
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				std::vector<ActiveType> dS1(times_.size()-1), dS2(times_.size()-1);
				ActiveType EdS1 = 0.0, EdS2 = 0.0; 
				for (size_t i=1; i<times_.size(); ++i) {
					dS1[i-1] =  swapRate(p,times_[i],times_[i]+T1_) - swapRate(p,times_[i-1],times_[i-1]+T1_);
					dS2[i-1] =  swapRate(p,times_[i],times_[i]+T2_) - swapRate(p,times_[i-1],times_[i-1]+T2_);
					EdS1 += dS1[i-1];
					EdS2 += dS2[i-1];
				}
				EdS1 /= dS1.size();
				EdS2 /= dS2.size();
				ActiveType Var1=0.0, Var2=0.0, Cov=0.0;
				for (size_t i=0; i<times_.size()-1; ++i) {
					Var1 += (dS1[i] - EdS1)*(dS1[i] - EdS1);
					Var2 += (dS2[i] - EdS2)*(dS2[i] - EdS2);
					Cov  += (dS1[i] - EdS1)*(dS2[i] - EdS2);
				}
				return Cov / sqrt(Var1*Var2);
			}
		};

	};



}

#endif  /* ifndef quantlib_templatemcpayoff_hpp */
