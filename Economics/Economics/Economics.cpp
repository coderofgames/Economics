// Economics.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <cmath>
#include <vector>
#include <iostream>

using namespace std;

#include "finance\all_cc_progs\fin_recipes.h"
#include "finance\all_cc_progs\normdist.h"


#ifdef _HAVE_NEWMAT_
#include "fin_recipes_newmat.h"
#include "examples_mean_variance_cxx_newmat.cc"
#include "examples_implicit_finite_diff_using_newmat.cc"
#endif 
#ifdef _HAVE_ITPP_
#include "fin_recipes_itpp.h"
#include "examples_mean_variance_cxx_itpp.cc"
#include "examples_implicit_finite_diff_using_itpp.cc"
#endif 

#ifdef _HAVE_GSL_
#include "fin_recipes_gsl.h"
#endif

//====================================================================================================================================
// ALTERNATE FORMULAS
//====================================================================================================================================
void test_merton_jump_diff_call(){
	double S = 100; double K = 100; double r = 0.05;
	double sigma = 0.3;
	double time_to_maturity = 1;
	double lambda = 0.5;
	double kappa = 0.5;
	double delta = 0.5;
	cout << " Merton Jump diffusion call = "
		<< option_price_call_merton_jump_diffusion(S, K, r, sigma, time_to_maturity, lambda, kappa, delta)
		<< endl;
};

void test_heston(){
#ifdef _HAVE_GSL_
	double S = 100;
	double K = 100;
	double r = 0.01;
	double v = 0.01;
	double tau = 0.5;
	double rho = 0;
	double kappa = 2;
	double lambda = 0.0;
	double theta = 0.01;
	double sigma = 0.01;
	cout << "heston call price " << heston_call_option_price(S, K, r, v, tau, rho, kappa, lambda, theta, sigma) << endl;
#endif
};

void alternative_formulas_examples(){
	cout << "-----------------------------" << endl;
	cout << "Alternative formulas " << endl;
	cout << "-----------------------------" << endl;
	test_merton_jump_diff_call();
	test_heston();
};


//====================================================================================================================================
// APPROXIMATIONS
//====================================================================================================================================
void tst_johnson_approximation_am_put(){
	double r = 0.125;    double S = 1.1;    double X = 1;
	double sigma = 0.5;    double time = 1;
	cout << " American put price using Johnson approximation = "
		<< option_price_american_put_approximated_johnson(S, X, r, sigma, time)
		<< endl;
};

void test_baw_approximation_call(){
	double S = 100;   double X = 100;     double sigma = 0.20;
	double r = 0.08;  double b = -0.04;   double time = 0.25;
	cout << " Call price using Barone-Adesi Whaley approximation = "
		<< option_price_american_call_approximated_baw(S, X, r, b, sigma, time) << endl;
};

void approximations_examples(){
	cout << "------------------------------------" << endl;
	cout << "Approximations chapter " << endl;
	cout << "------------------------------------" << endl;
	tst_johnson_approximation_am_put();
	test_baw_approximation_call();
};

//====================================================================================================================================
// AVERAGE AND LOOPBACK OPTIONS
//====================================================================================================================================
void test_bermudan_option(){
	double S = 80;         double K = 100;          double r = 0.20;
	double time = 1.0;   double sigma = 0.25;
	int steps = 500;
	double q = 0.0;
	vector<double> potential_exercise_times;  potential_exercise_times.push_back(0.25);
	potential_exercise_times.push_back(0.5);  potential_exercise_times.push_back(0.75);
	cout << " Bermudan put price = "
		<< option_price_put_bermudan_binomial(S, K, r, q, sigma, time, potential_exercise_times, steps)
		<< endl;
};

void test_analytical_geometric_average(){
	double S = 100;  double K = 100; double q = 0;
	double r = 0.06; double sigma = 0.25; double time = 1.0;
	cout << " Analytical geometric average = "
		<< option_price_asian_geometric_average_price_call(S, K, r, q, sigma, time)
		<< endl;
};

void test_exotics_lookback(){
	double S = 100; double Smin = S; double q = 0; double r = 0.06;
	double sigma = 0.346; double time = 1.0;
	cout << " Lookback call price = "
		<< option_price_european_lookback_call(S, Smin, r, q, sigma, time) << endl;
};

//====================================================================================================================================
// EXAMPLES BINOMIAL
//====================================================================================================================================

void test_bin_eur_call_ud(){
	double S = 100.0;    double K = 100.0;    double r = 0.025;
	double u = 1.05;     double d = 1 / u;
	cout << " one period european call =  "
		<< option_price_call_european_binomial_single_period(S, K, r, u, d) << endl;
	int no_periods = 2;
	cout << " two period european call =  "
		<< option_price_call_european_binomial_multi_period_given_ud(S, K, r, u, d, no_periods) << endl;

};

void binomial_examples(){
	cout << "----------------------------" << endl;
	cout << "Binomial Chapter " << endl;
	cout << "----------------------------" << endl;
	test_bin_eur_call_ud();
};

//====================================================================================================================================
// EXAMPLES BINOMIAL APPROXIMATIONS
//====================================================================================================================================
void test_binomial_approximations_option_pricing(){
	double S = 100.0;    double K = 100.0;
	double r = 0.1;      double sigma = 0.25;
	double time = 1.0;
	int no_steps = 100;
	cout << " european call = "
		<< option_price_call_european_binomial(S, K, r, sigma, time, no_steps)
		<< endl;
	cout << " american call = "
		<< option_price_call_american_binomial(S, K, r, sigma, time, no_steps)
		<< endl;
};

void test_binomial_approximations_option_price_partials(){
	double S = 100.0;    double K = 100.0;
	double r = 0.1;      double sigma = 0.25;
	double time = 1.0;     int no_steps = 100;

	double delta, gamma, theta, vega, rho;
	option_price_partials_american_call_binomial(S, K, r, sigma, time, no_steps,
		delta, gamma, theta, vega, rho);
	cout << " Call price partials " << endl;
	cout << "  delta = " << delta << endl;
	cout << "  gamma = " << gamma << endl;
	cout << "  theta = " << theta << endl;
	cout << "  vega  = " << vega << endl;
	cout << "  rho   = " << rho << endl;
};

void test_binomial_approximations_option_price_dividends(){
	double S = 100.0;    double K = 100.0;
	double r = 0.10;      double sigma = 0.25;
	double time = 1.0;
	int no_steps = 100;
	double d = 0.02;
	cout << " call price with continuous dividend payout = "
		<< option_price_call_american_binomial(S, K, r, d, sigma, time, no_steps) << endl;
	vector<double> dividend_times;      vector<double> dividend_yields;
	dividend_times.push_back(0.25);     dividend_yields.push_back(0.025);
	dividend_times.push_back(0.75);     dividend_yields.push_back(0.025);
	cout << " call price with proportial dividend yields at discrete dates = "
		<< option_price_call_american_proportional_dividends_binomial(S, K, r, sigma, time, no_steps,
		dividend_times, dividend_yields)
		<< endl;

	vector<double> dividend_amounts;
	dividend_amounts.push_back(2.5);
	dividend_amounts.push_back(2.5);
	cout << " call price with proportial dividend amounts at discrete dates = "
		<< option_price_call_american_discrete_dividends_binomial(S, K, r, sigma, time, no_steps,
		dividend_times, dividend_amounts)
		<< endl;
};

void test_binomial_approximations_futures_options(){
	double F = 50.0;    double K = 45.0;
	double r = 0.08;    double sigma = 0.2;
	double time = 0.5;
	int no_steps = 100;
	cout << " european futures call option = "
		<< futures_option_price_call_american_binomial(F, K, r, sigma, time, no_steps) << endl;
};

void test_binomial_approximations_currency_options(){
	double S = 50.0;      double K = 52.0;
	double r = 0.08;      double rf = 0.05;
	double sigma = 0.2;   double time = 0.5;

	int no_steps = 100;
	cout << " european currency option call = "
		<< currency_option_price_call_american_binomial(S, K, r, rf, sigma, time, no_steps) << endl;
};

void binomial_approximations_examples(){
	cout << "-------------------------------------" << endl;
	cout << "Binomial Approximations examples" << endl;
	cout << "-------------------------------------" << endl;
	test_binomial_approximations_option_pricing();
	test_binomial_approximations_option_price_partials();
	test_binomial_approximations_option_price_dividends();
	test_binomial_approximations_futures_options();
	test_binomial_approximations_currency_options();
};
//====================================================================================================================================
// EXMAPLES BINOMIAL TERM STRUCTURE MODELS
//====================================================================================================================================
void test_rendleman_bartter_zero_coupon_call() {
	double K = 950; double S = 0.15; double M = 0.05; double interest = 0.10;
	double option_maturity = 4; double bond_maturity = 5; double bond_maturity_payment = 1000;
	int no_steps = 100;
	cout << " Rendleman Bartter price of option on zero coupon bond: ";
	double c = bond_option_price_call_zero_american_rendleman_bartter(K, option_maturity, S, M,
		interest, bond_maturity,
		bond_maturity_payment, no_steps);
	cout << " c = " << c << endl;
};


void binomial_term_structure_models_examples(){
	cout << "---------------------------------------" << endl;
	cout << "Binomial term structure examples " << endl;
	cout << "---------------------------------------" << endl;
	test_rendleman_bartter_zero_coupon_call();
};

//====================================================================================================================================
// EXMAPLES BLACK SHOLES
//====================================================================================================================================
void test_option_price_call_black_scholes(){
	double S = 50; double K = 50; double r = 0.10;
	double sigma = 0.30; double time = 0.50;
	cout << " Black Scholes call price = ";
	cout << option_price_call_black_scholes(S, K, r, sigma, time) << endl;
};
void test_black_scholes_partials_call(){
	cout << " Black Scholes call partial derivatives " << endl;
	double S = 50; double K = 50; double r = 0.10;
	double sigma = 0.30; double time = 0.50;
	double Delta, Gamma, Theta, Vega, Rho;
	option_price_partials_call_black_scholes(S, K, r, sigma, time,
		Delta, Gamma, Theta, Vega, Rho);
	cout << "  Delta = " << Delta << endl;
	cout << "  Gamma = " << Gamma << endl;
	cout << "  Theta = " << Theta << endl;
	cout << "  Vega  = " << Vega << endl;
	cout << "  Rho   = " << Rho << endl;
};

void test_black_scholes_implied_volatility(){
	double S = 50; double K = 50; double r = 0.10; double time = 0.50;
	double C = 2.5;
	cout << " Black Scholes implied volatility using Newton search = ";
	cout << option_price_implied_volatility_call_black_scholes_newton(S, K, r, time, C) << endl;
	cout << " Black Scholes implied volatility using bisections = ";
	cout << option_price_implied_volatility_call_black_scholes_bisections(S, K, r, time, C) << endl;
};

void black_scholes_examples(){
	cout << "----------------------------------------" << endl;
	cout << "Examples in Black Scholes chapter " << endl;
	cout << "----------------------------------------" << endl;
	test_option_price_call_black_scholes();
	test_black_scholes_partials_call();
	test_black_scholes_implied_volatility();
};
//====================================================================================================================================
// EXAMPLES BLACK SHOLES EXTENSIONS
//====================================================================================================================================
void test_black_scholes_with_dividends(){
	double S = 100.0;    double K = 100.0;
	double r = 0.1;      double sigma = 0.25;
	double time = 1.0;
	double dividend_yield = 0.05;
	vector<double> dividend_times;     vector<double> dividend_amounts;
	dividend_times.push_back(0.25);    dividend_amounts.push_back(2.5);
	dividend_times.push_back(0.75);    dividend_amounts.push_back(2.5);
	cout << " european stock call option with contininous dividend = "
		<< option_price_european_call_payout(S, K, r, dividend_yield, sigma, time) << endl;
	cout << " european stock call option with discrete dividend =  "
		<< option_price_european_call_dividends(S, K, r, sigma, time, dividend_times, dividend_amounts) << endl;
};

void test_rgw_price_am_call_div(){
	double S = 100.0;    double K = 100.0;
	double r = 0.1;      double sigma = 0.25;
	double tau = 1.0;    double tau1 = 0.5;
	double D1 = 10.0;
	cout << " american call price with one dividend = "
		<< option_price_american_call_one_dividend(S, K, r, sigma, tau, D1, tau1) << endl;
};

void test_futures_option_price_black(){
	double F = 50.0;    double K = 45.0;
	double r = 0.08;    double sigma = 0.2;
	double time = 0.5;
	cout << " european futures call option = "
		<< futures_option_price_put_european_black(F, K, r, sigma, time) << endl;
};

void test_currency_option_european_call(){
	double S = 50.0;      double K = 52.0;
	double r = 0.08;      double rf = 0.05;
	double sigma = 0.2;   double time = 0.5;
	cout << " european currency call option = "
		<< currency_option_price_call_european(S, K, r, rf, sigma, time) << endl;
};

void test_option_price_perpetual_american_call(){
	double S = 50.0;    double K = 40.0;
	double r = 0.05;    double q = 0.02;
	double sigma = 0.05;
	double price = option_price_american_perpetual_call(S, K, r, q, sigma);
	cout << " perpetual call price = " << price << endl;
};

void black_scholes_extensions_examples(){
	cout << "-------------------------------" << endl;
	cout << "Black Scholes extensions" << endl;
	cout << "-------------------------------" << endl;
	test_black_scholes_with_dividends();
	test_futures_option_price_black();
	test_currency_option_european_call();
	test_rgw_price_am_call_div();
	test_option_price_perpetual_american_call();
};

//====================================================================================================================================
// EXAMPLES BOND OPTIONS
//====================================================================================================================================
void test_bond_option_gbm_pricing(){
	double B = 100;
	double K = 100;
	double r = 0.05;
	double sigma = 0.1;
	double time = 1;
	cout << " zero coupon put option price = "
		<< bond_option_price_put_zero_black_scholes(B, K, r, sigma, time) << endl;

	vector<double> coupon_times; coupon_times.push_back(0.5);
	vector<double> coupons; coupons.push_back(1);
	cout << " coupon bond put option price = "
		<< bond_option_price_put_coupon_bond_black_scholes(B, K, r, sigma, time, coupon_times, coupons);
	cout << endl;

	int steps = 100;
	cout << " zero coupon american put option price, binomial = "
		<< bond_option_price_put_american_binomial(B, K, r, sigma, time, steps) << endl;
};

void bond_options_examples(){
	cout << "---------------------------------------" << endl;
	cout << "Bond option pricing, simple (GBM) case " << endl;
	cout << "---------------------------------------" << endl;
	test_bond_option_gbm_pricing();
};

//====================================================================================================================================
// EXAMPLES BOND PRICING FLAT TERM STRUCTURE
//====================================================================================================================================
void example_bond_pricing_flat_term_structure(){
	vector<double> cflows; cflows.push_back(10); cflows.push_back(10); cflows.push_back(110);
	vector<double> times;  times.push_back(1);   times.push_back(2);   times.push_back(3);
	double r = 0.09;
	cout << "Discrete discounting" << endl;
	double B = bonds_price_discrete(times, cflows, r);
	cout << " bonds price    = " << B << endl;
	cout << " bond yield to maturity = " << bonds_yield_to_maturity_discrete(times, cflows, B) << endl;
	cout << " bond duration  = " << bonds_duration_discrete(times, cflows, r) << endl;
	cout << " bond duration modified = " << bonds_duration_modified_discrete(times, cflows, B) << endl;
	cout << " bond convexity =" << bonds_convexity_discrete(times, cflows, r) << endl;
	cout << " new bond price = " << bonds_price_discrete(times, cflows, 0.1) << endl;
	cout << "Continous discounting" << endl;
	B = bonds_price(times, cflows, r);
	cout << " bonds price    = " << B << endl;
	cout << " bond yield to maturity = " << bonds_yield_to_maturity(times, cflows, B) << endl;
	cout << " bond duration  = " << bonds_duration(times, cflows, r) << endl;
	cout << " bond convexity =" << bonds_convexity(times, cflows, r) << endl;
	cout << " new bond price = " << bonds_price(times, cflows, 0.08) << endl;
}

void bond_pricing_flat_term_structure_examples(){
	cout << "----------------------------" << endl;
	cout << "Bond Pricing with a flat term structure Chapter " << endl;
	cout << "----------------------------" << endl;
	example_bond_pricing_flat_term_structure();
};
//====================================================================================================================================
// EXAMPLES CREDIT RISK
//====================================================================================================================================
void test_credit_risk(){
	cout << " Credit Risk Calculation " << endl;
	double V = 100; double F = 90; double r = 0.05; double T = 1; double sigma = 0.25;
	double p = option_price_put_black_scholes(V, F, r, sigma, T);
	cout << " Debt value = " << exp(-r*T)*F - p << endl;
};

void credit_risk_examples(){
	cout << "------------------------------" << endl;
	cout << "Credit Risk Examples " << endl;
	cout << "------------------------------" << endl;
	test_credit_risk();
};
//====================================================================================================================================
// EXAMPLES FINITE DIFFERENCES
//====================================================================================================================================

void test_explicit_finite_differences(){
	double S = 50.0;
	double K = 50.0;
	double r = 0.1;
	double sigma = 0.4;
	double time = 0.4167;
	int no_S_steps = 20;
	int no_t_steps = 11;
	cout << " explicit finite differences, european put price = ";
	cout << option_price_put_european_finite_diff_explicit(S, K, r, sigma, time, no_S_steps, no_t_steps)
		<< endl;
	cout << " explicit finite differences, american put price = ";
	cout << option_price_put_american_finite_diff_explicit(S, K, r, sigma, time, no_S_steps, no_t_steps)
		<< endl;
};

void finite_differences_examples(){
	cout << "----------------------------" << endl;
	cout << "Finite Differences examples " << endl;
	cout << "----------------------------" << endl;
	test_explicit_finite_differences();
};

//====================================================================================================================================
// EXAMPLES FORWARDS FUTURES
//====================================================================================================================================
void test_futures_price(){
	double S = 100;    double r = 0.10;    double time = 0.5;
	cout << " futures price = " << futures_price(S, r, time) << endl;
};

void forwards_futures_examples(){
	cout << "---------------------------" << endl;
	cout << "Futures/Forwards pricing   " << endl;
	cout << "---------------------------" << endl;
	test_futures_price();
};

//====================================================================================================================================
// EXAMPLES GENERIC BINOMIAL 
//====================================================================================================================================
void example_binomial_generic_standard_put_and_calls(){

	double S = 100.0;
	double K = 100.0;
	double r = 0.1;
	double sigma = 0.25;
	double time_to_maturity = 1.0;
	int steps = 100;
	cout << " american call price = "
		<< option_price_generic_binomial(S, K, payoff_call, r, sigma, time_to_maturity, steps)
		<< endl;
	cout << " american put price = "
		<< option_price_generic_binomial(S, K, payoff_put, r, sigma, time_to_maturity, steps)
		<< endl;
};
void example_binomial_generic_binary(){
	double S = 100.0;
	double K = 120.0;
	double r = 0.1;
	double sigma = 0.25;
	double time_to_maturity = 1.0;
	int steps = 100;
	cout << " binary option price = "
		<< option_price_generic_binomial(S, K, payoff_binary_call, r, sigma, time_to_maturity, steps)
		<< endl;
};

void examples_generic_binomial(){
	cout << "----------------------------------------" << endl;
	cout << "Generic binomial pricing                " << endl;
	cout << "----------------------------------------" << endl;
	example_binomial_generic_standard_put_and_calls();
	example_binomial_generic_binary();
};

//====================================================================================================================================
// EXAMPLES IMPLICIT FINITE DIFFERENCE USING ITPP
//====================================================================================================================================
#ifdef _HAVE_ITPP_
void test_implicit_finite_differences_using_itpp(){
	double S = 50.0;
	double K = 50.0;
	double r = 0.1;
	double sigma = 0.4;
	double time = 0.5;
	int no_S_steps = 200;
	int no_t_steps = 200;
	cout << " black scholes put price = " << option_price_put_black_scholes(S, K, r, sigma, time) << endl;
	cout << " implicit American put price = ";
	cout << option_price_put_american_finite_diff_implicit_itpp(S, K, r, sigma, time, no_S_steps, no_t_steps) << endl;
};

void examples_finite_diffs_using_itpp(){
	cout << "----------------------------" << endl;
	cout << "Finite Differences examples using it++ " << endl;
	cout << "----------------------------" << endl;
	test_implicit_finite_differences_using_itpp();
};
#endif
//====================================================================================================================================
// EXAMPLES IMPLICIT FINITE DIFFERENCE USING NEWMAT
//====================================================================================================================================
#ifdef _HAVE_NEWMAT_

void test_implicit_finite_differences_using_newmat(){
	double S = 50.0;
	double K = 50.0;
	double r = 0.1;
	double sigma = 0.4;
	double time = 0.5;
	int no_S_steps = 200;
	int no_t_steps = 200;
	cout << " black scholes put price = " << option_price_put_black_scholes(S, K, r, sigma, time) << endl;
	cout << " implicit Euro put price = ";
	cout << option_price_put_european_finite_diff_implicit(S, K, r, sigma, time, no_S_steps, no_t_steps) << endl;
	cout << " implicit American put price = ";
	cout << option_price_put_american_finite_diff_implicit(S, K, r, sigma, time, no_S_steps, no_t_steps) << endl;
};

void examples_finite_diffs_using_newmat(){
	cout << "----------------------------" << endl;
	cout << "Finite Differences examples using newmat " << endl;
	cout << "----------------------------" << endl;
	test_implicit_finite_differences_using_newmat();
};
#endif
//====================================================================================================================================
// EXAMPLES INTREST RATE TREES GBM
//====================================================================================================================================
void example_interest_rate_trees_gbm_build_tree(){
	vector< vector<double> > tree = interest_rate_trees_gbm_build(0.1, 1.02, 0.99, 3);
	cout << " Interest rate tree: " << endl;
	cout << " Time 0: " << tree[0][0] << endl;
	cout << " Time 1: " << tree[1][0] << "  " << tree[1][1] << endl;
	cout << " Time 2: " << tree[2][0] << "  " << tree[2][1] << "  " << tree[2][2] << endl;
};

void example_interest_rate_trees_gbm_price_bond(){
	double r0 = 0.1;
	double u = 1.02;  double d = 0.99;
	int n = 3;
	double q = 0.5;
	vector< vector<double> > tree = interest_rate_trees_gbm_build(r0, u, d, n);
	vector<double> cashflows;
	cashflows.push_back(0); cashflows.push_back(10); cashflows.push_back(10); cashflows.push_back(110);
	cout << "Bond price B = " << interest_rate_trees_gbm_value_of_cashflows(cashflows, tree, q) << endl;
};

void example_interest_rate_trees_gbm_price_callable_bond(){
	double r0 = 0.06;
	double u = 1.2;  double d = 0.9;
	int n = 10;
	double q = 0.5;
	vector< vector<double> > tree = interest_rate_trees_gbm_build(r0, u, d, n);
	vector<double> cashflows;
	cashflows.push_back(0);
	for (int t = 1; t <= 9; ++t){ cashflows.push_back(6); };
	cashflows.push_back(106);
	cout << "Straight bond price  = " << interest_rate_trees_gbm_value_of_cashflows(cashflows, tree, q) << endl;
	int first_call_time = 6;
	double call_price = 106;
	cout << "Callable bond price = "
		<< interest_rate_trees_gbm_value_of_callable_bond(cashflows, tree, q, first_call_time, call_price) << endl;
};

void examples_interest_rate_trees_gbm(){
	cout << "-----------------------------" << endl;
	cout << "Interest rate trees " << endl;
	cout << "-----------------------------" << endl;
	example_interest_rate_trees_gbm_build_tree();
	example_interest_rate_trees_gbm_price_bond();
	example_interest_rate_trees_gbm_price_callable_bond();

};
//====================================================================================================================================
// EXAMPLES INTREST RATE TREES HO LEE
//====================================================================================================================================

void test_ho_lee_pricing_european_call(){
	double delta = 0.98;
	double pi = 0.5;
	double r = 0.1;
	term_structure_class* initial = new term_structure_class_flat(r);
	vector<double> times; times.push_back(5.0);
	vector<double> cflows; cflows.push_back(100);
	double K = 80;
	double time_to_maturity = 3;
	cout << " Flat term structure " << endl;
	cout << " c= " << price_european_call_option_on_bond_using_ho_lee(initial, delta, pi, times, cflows, K, time_to_maturity);
	cout << endl;
	delete (initial);
	double beta0 = 0.09;  double beta1 = 0.01; double beta2 = 0.01; double lambda = 5.0;
	initial = new term_structure_class_nelson_siegel(beta0, beta1, beta2, lambda);
	cout << " Nelson Siegel term structure " << endl;
	cout << " c= " << price_european_call_option_on_bond_using_ho_lee(initial, delta, pi, times, cflows, K, time_to_maturity);
	cout << endl;
	delete (initial);
};

void examples_interest_rate_trees_ho_lee(){
	cout << "------------------------------------" << endl;
	cout << "---- interest rate trees ho lee ----" << endl;
	cout << "------------------------------------" << endl;
	test_ho_lee_pricing_european_call();
};
//====================================================================================================================================
// EXAMPLES MEAN VARIANCE CXX ITPP
//====================================================================================================================================

#ifdef _HAVE_ITPP_

#include "fin_recipes_itpp.h"

void test_mean_variance_calculations_itpp(){
	cout << " testing basic mean variance calculations using IT++ " << endl;
	vec e = " 0.05 0.1 0.075";
	mat V = "0.9  -0.2 0.0; -0.2 1.0 -0.3; 0.0 -0.3 0.6";
	vec w = "0.3333 0.3333 0.333";
	cout << " mean " << mv_calculate_mean(e, w) << endl;
	cout << " variance " << mv_calculate_variance(V, w) << endl;
	cout << " stdev " << mv_calculate_st_dev(V, w) << endl;
};

void test_mean_variance_calculations_portfolio_itpp(){
	cout << " Example mean variance optimal portfolis using itpp " << endl;
	vec e = " 0.05 0.1 0.075";
	mat V = "0.9  -0.2 0.0; -0.2 1.0 -0.3; 0.0 -0.3 0.6";
	double r = 0.075;
	mat w = mv_calculate_portfolio_given_mean_unconstrained(e, V, r);
	cout << " suggested portfolio " << endl;
	cout << w << endl;
};


#endif

void examples_mean_variance_itpp() {
	cout << "-------------------------------" << endl;
	cout << "Mean Variance Examples using itpp " << endl;
	cout << "-------------------------------" << endl;
#ifdef _HAVE_ITPP_
	test_mean_variance_calculations_itpp();
	test_mean_variance_calculations_portfolio_itpp();
#endif
};
//====================================================================================================================================
// EXAMAPLES MEAN VARIANCE CXX NEWMAT
//====================================================================================================================================

#ifdef _HAVE_NEWMAT_
#include "newmat.h"
#include "fin_recipes_newmat.h"

void test_mean_variance_calculations_newmat(){
	cout << "Simple example of mean variance calculations, using newmat " << endl;
	Matrix e(2, 1); e(1, 1) = 0.05; e(2, 1) = 0.1;
	Matrix V(2, 2);
	V(1, 1) = 1.0; V(2, 1) = 0.0;
	V(1, 2) = 0.0; V(2, 2) = 1.0;
	Matrix w(2, 1);
	w(1, 1) = 0.5;
	w(2, 1) = 0.5;
	cout << " mean " << mv_calculate_mean(e, w) << endl;
	cout << " variance " << mv_calculate_variance(V, w) << endl;
	cout << " stdev " << mv_calculate_st_dev(V, w) << endl;
};

void test_mean_variance_portfolio_calculation_newmat(){
	cout << "Calculate mean variance portfolio using newmat: " << endl;
	Matrix e(2, 1);
	e(1, 1) = 0.05; e(2, 1) = 0.1;
	Matrix V(2, 2);
	V(1, 1) = 1.0; V(2, 1) = 0.0;
	V(1, 2) = 0.0; V(2, 2) = 1.0;
	double r = 0.075;
	Matrix w = mv_calculate_portfolio_given_mean_unconstrained(e, V, r);
	cout << " suggested portfolio:  ";
	cout << " w1 = " << w(1, 1) << " w2 = " << w(2, 1) << endl;
};
#endif

void examples_mean_variance_newmat() {
	cout << "-------------------------------" << endl;
	cout << "Mean Variance Examples using newmat " << endl;
	cout << "-------------------------------" << endl;
#ifdef _HAVE_NEWMAT_
	test_mean_variance_calculations_newmat();
	test_mean_variance_portfolio_calculation_newmat();
#endif
};

//====================================================================================================================================
// EXAMPLES PRESENT VALUE
//====================================================================================================================================

void test_present_value(){
	vector<double> cflows; cflows.push_back(-100.0); cflows.push_back(10.0); cflows.push_back(110.0);
	vector<double> times; times.push_back(0.0); times.push_back(1); times.push_back(2);
	double r = 0.05;
	cout << " present value, 5\% discretely compounded interest = ";
	cout << cash_flow_pv_discrete(times, cflows, r) << endl;
	cout << " internal rate of return, discrete compounding = ";
	cout << cash_flow_irr_discrete(times, cflows) << endl;
	cout << " present value, 5\% continously compounded interest = ";
	cout << cash_flow_pv(times, cflows, r) << endl;
	cout << " internal rate of return, continous compounding = ";
	cout << cash_flow_irr(times, cflows) << endl;
};

void present_value_examples(){
	cout << "----------------------------" << endl;
	cout << "Present Value Chapter " << endl;
	cout << "----------------------------" << endl;
	test_present_value();
};


//====================================================================================================================================
// EXAMPLES SIMULATION
//====================================================================================================================================
void test_simulation_pricing() {
	double S = 100.0;  double K = 100.0; double r = 0.1; double sigma = 0.25;
	double time = 1.0; int no_sims = 5000;
	cout << " call:  black scholes price = " << option_price_call_black_scholes(S, K, r, sigma, time) << endl;
	cout << "        simulated price     = "
		<< option_price_call_european_simulated(S, K, r, sigma, time, no_sims) << endl;
	cout << " put:  black scholes price = " << option_price_put_black_scholes(S, K, r, sigma, time) << endl;
	cout << "       simulated price     = "
		<< option_price_put_european_simulated(S, K, r, sigma, time, no_sims) << endl;
};

void test_simulation_pricing_delta(){
	double S = 100.0;  double K = 100.0; double r = 0.1; double sigma = 0.25;
	double time = 1.0; int no_sims = 5000;
	cout << " call: bs delta  = " << option_price_delta_call_black_scholes(S, K, r, sigma, time)
		<< "       sim delta = " << option_price_delta_call_european_simulated(S, K, r, sigma, time, no_sims)
		<< endl;
	cout << " put: bs delta  = " << option_price_delta_put_black_scholes(S, K, r, sigma, time)
		<< "      sim delta = " << option_price_delta_put_european_simulated(S, K, r, sigma, time, no_sims)
		<< endl;
};


void test_simulation_bs_case_using_generic_routine(){
	double S = 100; double X = 100; double r = 0.1;
	double sigma = 0.25; double time = 1.0;  int no_sims = 50000;
	cout << "Black Scholes call option price = " << option_price_call_black_scholes(S, X, r, sigma, time) << endl;
	cout << "Simulated call option price     = "
		<< derivative_price_simulate_european_option_generic(S, X, r, sigma, time, payoff_call, no_sims)
		<< endl;
	cout << "Black Scholes put option price  = " << option_price_put_black_scholes(S, X, r, sigma, time) << endl;
	cout << "Simulated put option price      = "
		<< derivative_price_simulate_european_option_generic(S, X, r, sigma, time, payoff_put, no_sims)
		<< endl;
};
void test_simulation_bs_case_using_generic_routine_improving_efficiency(){
	double S = 100; double K = 100; double r = 0.1;
	double sigma = 0.25; double time = 1; int no_sims = 50000;
	cout << "Black Scholes call option price = "
		<< option_price_call_black_scholes(S, K, r, sigma, time) << endl;
	cout << "Simulated call option price     = "
		<< derivative_price_simulate_european_option_generic(S, K, r, sigma, time, payoff_call, no_sims)
		<< endl;
	cout << "Simulated call option price, CV = "
		<< derivative_price_simulate_european_option_generic_with_control_variate(S, K, r, sigma, time,
		payoff_call, no_sims)
		<< endl;
	cout << "Simulated call option price, AV = "
		<< derivative_price_simulate_european_option_generic_with_antithetic_variate(S, K, r, sigma, time,
		payoff_call, no_sims)
		<< endl;
	cout << "Black Scholes put option price  = " << option_price_put_black_scholes(S, K, r, sigma, time) << endl;
	cout << "Simulated put option price      = "
		<< derivative_price_simulate_european_option_generic(S, K, r, sigma, time, payoff_put, no_sims) << endl;
	cout << "Simulated put option price,  CV = "
		<< derivative_price_simulate_european_option_generic_with_control_variate(S, K, r, sigma, time,
		payoff_put, no_sims)
		<< endl;
	cout << "Simulated put option price,  AV = "
		<< derivative_price_simulate_european_option_generic_with_antithetic_variate(S, K, r, sigma, time,
		payoff_put, no_sims)
		<< endl;
};

void test_simulation_binary_options(){
	double S = 100.0; double K = 100.0; double r = 0.1; double sigma = 0.25;
	double time = 1.0;  int no_sims = 5000;
	cout << " cash or nothing, Q=1: "
		<< derivative_price_simulate_european_option_generic(S, K, r, sigma, time,
		payoff_cash_or_nothing_call,
		no_sims)
		<< endl;
	cout << " control_variate "
		<< derivative_price_simulate_european_option_generic_with_control_variate(S, K, r, sigma, time,
		payoff_cash_or_nothing_call,
		no_sims)
		<< endl;
	cout << " antithetic_variate "
		<< derivative_price_simulate_european_option_generic_with_antithetic_variate(S, K, r, sigma, time,
		payoff_cash_or_nothing_call,
		no_sims)
		<< endl;
	cout << " asset or nothing: "
		<< derivative_price_simulate_european_option_generic(S, K, r, sigma, time,
		payoff_asset_or_nothing_call,
		no_sims)
		<< endl;
	cout << " control_variate "
		<< derivative_price_simulate_european_option_generic_with_control_variate(S, K, r, sigma, time,
		payoff_asset_or_nothing_call,
		no_sims)
		<< endl;
	cout << " antithetic_variate "
		<< derivative_price_simulate_european_option_generic_with_antithetic_variate(S, K, r, sigma, time,
		payoff_asset_or_nothing_call,
		no_sims)
		<< endl;
};

void simulation_examples(){
	cout << "--------------------------" << endl;
	cout << " Simulation examples " << endl;
	cout << "--------------------------" << endl;
	test_simulation_pricing();
	test_simulation_pricing_delta();
	test_simulation_bs_case_using_generic_routine();
	test_simulation_bs_case_using_generic_routine_improving_efficiency();
	test_simulation_binary_options();
};

//====================================================================================================================================
// EXAMPLES TERM STRUCTURE
//====================================================================================================================================
void test_termstru_transforms(){
	double t1 = 1;  double r_t1 = 0.05; double d_t1 = term_structure_discount_factor_from_yield(r_t1, t1);
	cout << " a " << t1 << " period spot rate of " << r_t1
		<< " corresponds to a discount factor of " << d_t1 << endl;
	double t2 = 2;  double d_t2 = 0.9;
	double r_t2 = term_structure_yield_from_discount_factor(d_t2, t2);
	cout << " a " << t2 << " period discount factor of " << d_t2
		<< " corresponds to a spot rate of " << r_t2 << endl;
	cout << " the forward rate between " << t1 << " and " << t2
		<< " is " << term_structure_forward_rate_from_discount_factors(d_t1, d_t2, t2 - t1)
		<< " using discount factors " << endl;
	cout << "  and is " << term_structure_forward_rate_from_yields(r_t1, r_t2, t1, t2)
		<< " using yields " << endl;
};

void test_term_structure_class_flat(){
	cout << "flat term structure class " << endl;
	term_structure_class_flat ts(0.05);
	double t1 = 1;
	cout << " discount factor t1 = " << t1 << ":" << ts.d(t1) << endl;
	double t2 = 2;
	cout << " discount factor t2 = " << t2 << ":" << ts.d(t2) << endl;
	cout << " spot rate t = " << t1 << ":" << ts.r(t1) << endl;
	cout << " spot rate t = " << t2 << ":" << ts.r(t2) << endl;
	cout << " forward rate from t1= " << t1 << " to t2= " << t2 << ":" << ts.f(t1, t2) << endl;
};

void test_termstru_interpolated(){
	vector<double> times;
	vector<double> yields;
	times.push_back(0.1);  times.push_back(0.5);  times.push_back(1);
	yields.push_back(0.1); yields.push_back(0.2); yields.push_back(0.3);
	times.push_back(5);    times.push_back(10);
	yields.push_back(0.4); yields.push_back(0.5);
	cout << " testing interpolated term structure" << endl;
	cout << " yields at times: " << endl;
	cout << " t=.1 " << term_structure_yield_linearly_interpolated(0.1, times, yields) << endl;
	cout << " t=0.5 " << term_structure_yield_linearly_interpolated(0.5, times, yields) << endl;
	cout << " t=1 " << term_structure_yield_linearly_interpolated(1, times, yields) << endl;
	cout << " t=3 " << term_structure_yield_linearly_interpolated(3, times, yields) << endl;
	cout << " t=5 " << term_structure_yield_linearly_interpolated(5, times, yields) << endl;
	cout << " t=10 " << term_structure_yield_linearly_interpolated(10, times, yields) << endl;
};


void test_term_structure_class_interpolated(){
	vector<double> times;     times.push_back(0.1);
	vector<double> spotrates; spotrates.push_back(0.05);
	times.push_back(1);       times.push_back(5);
	spotrates.push_back(0.07); spotrates.push_back(0.08);
	term_structure_class_interpolated ts(times, spotrates);
	double t1 = 1;
	cout << "discount factor t1 = " << t1 << ":" << ts.d(t1) << endl;
	double t2 = 2;
	cout << "discount factor t2 = " << t2 << ":" << ts.d(t2) << endl;
	cout << "spot rate t = " << t1 << ":" << ts.r(t1) << endl;
	cout << "spot rate t = " << t2 << ":" << ts.r(t2) << endl;
	cout << "forward rate from t1= " << t1 << " to t2= " << t2 << ":" << ts.f(t1, t2) << endl;
};


void test_term_structure_class_bond_calculations(){
	vector <double> times;     times.push_back(1);       times.push_back(2);
	vector <double> cashflows; cashflows.push_back(10);  cashflows.push_back(110);
	term_structure_class_flat tsflat(0.1);
	cout << " price = " << bonds_price(times, cashflows, tsflat) << endl;
	cout << " duration = " << bonds_duration(times, cashflows, tsflat) << endl;
	cout << " convexity = " << bonds_convexity(times, cashflows, tsflat) << endl;
};


void term_structure_examples(){
	cout << "--------------------------------------------" << endl;
	cout << " Term structure of interest rates examples  " << endl;
	cout << "--------------------------------------------" << endl;
	test_termstru_transforms();
	test_term_structure_class_flat();
	test_termstru_interpolated();
	test_term_structure_class_interpolated();
	test_term_structure_class_bond_calculations();
};
//====================================================================================================================================
// EXAMPLES TERM STRUCTURE DERIVATIVES
//====================================================================================================================================
void test_vasicek_option_pricing(){
	double a = 0.1;  double b = 0.1; double sigma = 0.02; double r = 0.05; double X = 0.9;
	cout << " Vasicek call option price "
		<< bond_option_price_call_zero_vasicek(X, r, 1, 5, a, b, sigma) << endl;
};

void term_structure_derivatives_examples(){
	cout << "--------------------------------------" << endl;
	cout << " term structure derivatives examples  " << endl;
	cout << "--------------------------------------" << endl;
	test_vasicek_option_pricing();
}
//====================================================================================================================================
// EXAMPLES TERM STRUCTURE MODEL
//====================================================================================================================================
void test_term_structure_nelson_siegel(){
	double beta0 = 0.09;  double beta1 = 0.01; double beta2 = 0.01; double lambda = 5.0;
	double t = 1.0;
	cout << "Example calculations using the Nelson Siegel term structure approximation" << endl;
	cout << " direct calculation, yield = "
		<< term_structure_yield_nelson_siegel(t, beta0, beta1, beta2, lambda) << endl;

	term_structure_class_nelson_siegel ns(beta0, beta1, beta2, lambda);
	cout << " using a term structure class" << endl;
	cout << " yield (t=1) = " << ns.r(t) << endl;
	cout << " discount factor (t=1) = " << ns.d(t) << endl;
	cout << " forward rate (t1=1, t2=2) = " << ns.f(1, 2) << endl;
};

void test_term_structure_cir(){
	cout << "Example calculations using the Cox Ingersoll Ross term structure model " << endl;
	double r = 0.05; double kappa = 0.01; double sigma = 0.1; double theta = 0.08; double lambda = 0.0;
	cout << " direct calculation, discount factor (t=1): "
		<< term_structure_discount_factor_cir(1, r, kappa, lambda, theta, sigma) << endl;
	cout << " using a class " << endl;
	term_structure_class_cir cir(r, kappa, lambda, theta, sigma);
	cout << " yield (t=1) = " << cir.r(1) << endl;
	cout << " discount factor (t=1) = " << cir.d(1) << endl;
	cout << " forward (t1=1, t2=2) = " << cir.f(1, 2) << endl;
};

void test_term_structure_cubic_spline(){
	cout << "Example term structure calculations using a cubic spline " << endl;
	double b = 0.1;  double c = 0.1; double d = -0.1;
	vector<double> f;     f.push_back(0.01);  f.push_back(0.01);  f.push_back(-0.01);
	vector<double> knots; knots.push_back(2); knots.push_back(7); knots.push_back(12);
	cout << " direct calculation, discount factor (t=1) "
		<< term_structure_discount_factor_cubic_spline(1, b, c, d, f, knots) << endl;
	cout << " Using a term structure class " << endl;
	term_structure_class_cubic_spline cs(b, c, d, f, knots);
	cout << " yield (t=1) = " << cs.r(1) << endl;
	cout << " discount factor (t=1) = " << cs.d(1) << endl;
	cout << " forward (t1=1, t2=2) = " << cs.f(1, 2) << endl;
};

void test_term_structure_vasicek() {
	cout << "Example term structure calculation using the Vasicek term structure model" << endl;
	double r = 0.05;  double a = -0.1; double b = 0.1; double sigma = 0.1;
	cout << " direct calculation, discount factor (t=1): "
		<< term_structure_discount_factor_vasicek(1, r, a, b, sigma) << endl;
	term_structure_class_vasicek vc(r, a, b, sigma);
	cout << " using a term structure class " << endl;
	cout << " yield (t=1) = " << vc.r(1) << endl;
	cout << " discount factor (t=1) = " << vc.d(1) << endl;
	cout << " forward rate (t1=1, t2=2) = " << vc.f(1, 2) << endl;
}

void test_term_structure_svensson(){
	double beta0 = 0.01;  double beta1 = 0.01; double beta2 = 0.01; double beta3 = 0.01;
	double tau1 = 5.0;     double tau2 = 5.0;
	double t = 1.0;
	cout << "Example calculations using the Svensson term structure approximation" << endl;
	cout << " direct calculation, yield = "
		<< term_structure_yield_svensson(t, beta0, beta1, beta2, beta3, tau1, tau2)
		<< endl;

	term_structure_class_svensson s(beta0, beta1, beta2, beta3, tau1, tau2);
	cout << " using a term structure class" << endl;
	cout << " yield (t=1) = " << s.r(t) << endl;
	cout << " discount factor (t=1) = " << s.d(t) << endl;
	cout << " forward rate (t1=1, t2=2) = " << s.f(1, 2) << endl;
};

void term_structure_model_examples(){
	cout << "----------------------------" << endl;
	cout << "Further term structure models " << endl;
	cout << "----------------------------" << endl;
	test_term_structure_nelson_siegel();
	test_term_structure_svensson();
	test_term_structure_cir();
	test_term_structure_cubic_spline();
	test_term_structure_vasicek();
};
//====================================================================================================================================
// EXAMPLES TRINOMIAL
//====================================================================================================================================
void test_trinomial_trees(){
	double S = 100.0; double K = 100.0;
	double r = 0.1;   double sigma = 0.25;
	double time = 1.0;  double q = 0;
	int no_steps = 100;
	cout << " european put using binomial = "
		<< option_price_put_european_binomial(S, K, r, sigma, time, no_steps)
		<< endl;
	cout << " american put using binomial = "
		<< option_price_put_american_binomial(S, K, r, sigma, time, no_steps)
		<< endl;
	cout << " american put using trinomial = "
		<< option_price_put_american_trinomial(S, K, r, q, sigma, time, no_steps)
		<< endl;
};

void examples_trinomial_trees(){
	cout << "----------------------------" << endl;
	cout << "-- trinomial trees ---------" << endl;
	cout << "----------------------------" << endl;
	test_trinomial_trees();
};

//====================================================================================================================================
// WARRENTS
//====================================================================================================================================
void test_warrant_price_adjusted_black_scholes(){
	double S = 48;  double K = 40; double r = 0.08;  double sigma = 0.30;
	double time = 0.5;   double m = 1000;   double n = 10000;
	double w = warrant_price_adjusted_black_scholes(S, K, r, sigma, time, m, n);
	cout << " warrant price = " << w << endl;
};

void warrant_examples(){
	cout << "----------------------------------------" << endl;
	cout << "Warrant pricing chapter" << endl;
	cout << "----------------------------------------" << endl;
	test_warrant_price_adjusted_black_scholes();
};

//====================================================================================================================================
// EXAMPLE FUTURES PRICE
//====================================================================================================================================


//====================================================================================================================================
// NORMAL DISTRIBUTION
//====================================================================================================================================
void test_cumulative_normal() {
	cout << " N(0) = " << N(0) << endl;
	cout << " N(0,0,0) = " << N(0, 0, 0) << endl;
};

void test_random_normal(){
	cout << " 5 random uniform numbers between 0 and 1: " << endl << "    ";
	for (int i = 0; i<5; ++i){ cout << " " << random_uniform_0_1(); }; cout << endl;
	cout << " 5 random normal(0,1) numbers: " << endl << "    ";
	for (int i = 0; i<5; ++i){ cout << " " << random_normal(); }; cout << endl;
};

void normal_distribution_examples(){
	cout << "----------------------------------" << endl;
	cout << " Normal distribution calculations " << endl;
	cout << "----------------------------------" << endl;
	test_cumulative_normal();
	test_random_normal();
};





int main(int argc, char** argv)

{

	present_value_examples();
	bond_pricing_flat_term_structure_examples();
	term_structure_examples();
	forwards_futures_examples();
	binomial_examples();
	black_scholes_examples();
	black_scholes_extensions_examples();
	binomial_approximations_examples();
	warrant_examples();
	simulation_examples();
	finite_differences_examples();
	approximations_examples();
	alternative_formulas_examples();
//	average_and_lookback_options_examples();
	examples_generic_binomial();
	examples_trinomial_trees();
	bond_options_examples();
	credit_risk_examples();
	term_structure_model_examples();
	binomial_term_structure_models_examples();
	examples_interest_rate_trees_gbm();
	examples_interest_rate_trees_ho_lee();
	term_structure_derivatives_examples();
	normal_distribution_examples();

#ifdef _HAVE_NEWMAT_
	examples_mean_variance_newmat();
	examples_finite_diffs_using_newmat();
#endif
#ifdef _HAVE_ITPP_
	examples_mean_variance_itpp();
	examples_finite_diffs_using_itpp();
#endif
	return 0;
}

