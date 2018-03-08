#include "barrier_option.hpp"


double barrier_option::d_calculate(double x, double y){
	return (log(x / y) + (r + sigma * sigma / 2.0) * tau) / (sigma * sqrt(tau));
}
double barrier_option::d_calculate_minus(double x, double y){
	return (log(x / y) + (r - sigma * sigma / 2.0) * tau) / (sigma * sqrt(tau));
}



double barrier_option::barrier_down(){
  double current_value = St;
  double minimum_price = St;
  double h = tau / number_iterations;
  for (int i = 0; i < number_iterations; i++){
    current_value = stock_price_single(current_value, t + (i + 1) * h, t + i * h);
    if (minimum_price > current_value) minimum_price = current_value;
    if (minimum_price < barrier) return 0;
 }
 if (current_value > K) return current_value;
   return 0.0;
}
double barrier_option::barrier_up(){
  double current_value = St;
  double maximum_price = St;
  double h = tau / number_iterations;
  for (int i = 0; i < number_iterations; i++){
    current_value = stock_price_single(current_value, t + (i + 1) * h, t + i * h);
    if (maximum_price < current_value) maximum_price = current_value;
    if (maximum_price > barrier) return 0;
  }
  if (current_value > K) return current_value;
  return 0.0;
}
double barrier_option::barrier_option_price(){
	double result = 0.0;
	double tem;
	for (int i = 0; i < number_iterations; i++){
		if(barrier<St)tem = barrier_down();
		else tem = barrier_up();
		if (tem!= 0)result += exp(-r*(tau))*(tem - K);
	}
	return result / number_iterations;
}



double barrier_option::delta_lr_step(double x){
	if (x == 0 )return 0.0;
	double mu = (r - sigma*sigma / 2) / sigma;
	double f = (normal_pdf((log(x / St) - sigma*mu*(tau)) / (sigma*sqrt(tau))) - exp(2 * mu*log(barrier / St) / sigma)*normal_pdf((log(x / St) - 2 * log(barrier / St) - sigma*mu*(tau)) / (sigma*sqrt(tau)))) / sqrt(tau);
	double df = exp(-0.5*pow((log(x / St) / sigma - mu*(tau)), 2) / (tau))*(log(x / St) / sigma - mu*(tau)) / (sqrt(2 * M_PI*pow((tau), 3))*sigma*St);
	df += sqrt(2 / M_PI) *mu*barrier*pow(barrier / St, (2 * mu / sigma - 1))*exp(-pow((-2 * log(barrier / St) / sigma + log(x / St) / sigma - mu*(tau)), 2) / (2 * (tau))) / (sigma*St*St*sqrt(tau));
	df += pow(barrier / St, (2 * mu / sigma))*(-2 * log(barrier / St) / sigma + log(x / St) / sigma - mu*(tau))*exp(-pow((-2 * log(barrier / St) / sigma + log(x / St) / sigma - mu*(tau)), 2) / (2 * (tau))) / (sqrt(2 * M_PI*pow((tau), 3))*sigma*St);
	double result = exp(-r*(tau))*(x - K)*df / f;
	return result;
}
double barrier_option::delta_lr(){
	double result = 0.0;
	if (barrier<St){
		for (int i = 0; i < number_iterations; i++){
			result += delta_lr_step(barrier_down());
		}
	}
	else{
		for (int i = 0; i < number_iterations; i++){
			result += delta_lr_step(barrier_up());
		}
	}
	return result / number_iterations;
}

double barrier_option::barrier_down_theoric(){
	double mu = r - sigma*sigma / 2.0;
	return payoff_theoretic(St, K) - pow(St / barrier, -2.0 * mu / (sigma*sigma))*payoff_theoretic(barrier*barrier / St, K);
}
double barrier_option::barrier_up_theoric(){
	if (barrier < K)return 0;
	double mu = r - sigma*sigma / 2.0;
	return payoff_theoretic(St, K) - payoff_theoretic(St, barrier) - (barrier - K)*exp(-r*(tau))*normal_cdf(d_calculate_minus(St, barrier)) - pow(barrier / St, 2.0 * mu / (sigma*sigma))*(payoff_theoretic(barrier*barrier / St, K) - payoff_theoretic(barrier*barrier / St, barrier) - (barrier - K)*exp(-r*(tau))*normal_cdf(d_calculate_minus(barrier, St)));
}
double barrier_option::price_theoretic(){
	if (barrier < St)return barrier_down_theoric();
	else return barrier_up_theoric();
}
double dnorm(double x){
	return -exp(-0.5*x*x)*x / sqrt(2.0*M_PI);
}
double barrier_option::delta_theoretic(){
	double mu = r - sigma*sigma / 2.0;
	double result = 0.0;
	if (barrier < St){
		result = normal_cdf((log(St / K) + (r + sigma*sigma / 2.0)*(tau)) / (sigma*sqrt(tau)));
		result -= pow((barrier / St), (r / (sigma*sigma) - 1))*(-barrier*barrier* normal_cdf((log(barrier*barrier / (St*K)) + mu*(tau)) / (sigma*sqrt(tau)) + sigma*sqrt(tau)) / (St*St) - 2.0 * mu*payoff_theoretic(barrier*barrier / St, K) / (St*sigma*sigma));

	}
	else{
		result = delta_theoretic_call(St, K) - delta_theoretic_call(St, barrier);
		result -= (barrier - K)*exp(-r*(tau))*normal_pdf(d_calculate_minus(St,barrier))/(sigma*St*sqrt(tau));
		result += (2.0*mu*pow(barrier / St, 2.0*mu / pow(sigma, 2)) / (pow(sigma, 2)*St))*(payoff_theoretic(barrier*barrier / St, K) - payoff_theoretic(barrier*barrier / St, barrier) - (barrier - K)*exp(-r*(tau))*normal_cdf(d_calculate_minus(barrier, St)));
		result -= pow(barrier / St, 2.0*mu / pow(sigma, 2))*(-barrier*barrier / (St*St)*delta_theoretic_call(barrier*barrier / St, K) + barrier*barrier / (St*St)*delta_theoretic_call(barrier*barrier / St, barrier) + (barrier - K)*exp(-r*(tau))*normal_pdf(d_calculate_minus(barrier, St)) / (sigma*St*sqrt(tau)));
	}
	return result;
}
double barrier_option::gamma_lr_step(double x){
	if (x == 0)return 0.0;
	double mu = (r - sigma*sigma / 2.0) / sigma;
	double f = (normal_pdf((log(x / St) - sigma*mu*(tau)) / (sigma*sqrt(tau))) - exp(2 * mu*log(barrier / St) / sigma)*normal_pdf((log(x / St) - 2 * log(barrier / St) - sigma*mu*(tau)) / (sigma*sqrt(tau)))) / sqrt(tau);
	double df = exp(-pow((log(x / St) / sigma - mu*(tau)), 2) / (2.0*(tau)))*pow(log(x / St) / sigma - mu*tau, 2) / (sigma*sigma*St*St*pow(tau, 2.5)*sqrt(2*M_PI));
	df -= exp(-pow(log(x / St) / sigma - mu*(tau), 2) / (2.0*(tau))) / (sqrt(2.0 * M_PI*pow((tau), 3))*sigma*sigma*St*St);
	df -= exp(-pow(log(x / St) / sigma - mu*(tau), 2) / (2.0*(tau)))*(log(x / St) / sigma - mu*tau) / (sqrt(2.0 * M_PI*pow((tau), 3))*sigma*St*St);
	df -= sqrt(2.0 / M_PI)*mu*barrier*barrier*(2.0 * mu / sigma - 1.0)*pow(barrier / St, 2.0 * mu / sigma - 2.0)*exp(-pow(-2.0*log(barrier / St) / sigma + log(x / St) / sigma - mu*(tau), 2) / (2.0*(tau))) / (sigma*pow(St, 4)*sqrt(tau));
	df -= sqrt(2.0 / M_PI)*mu*barrier*pow(barrier / St, 2.0*mu / sigma - 1.0)*(-2.0*log(barrier / St) / sigma + log(x / St) / sigma - mu*(tau))*exp(-pow(-2.0*log(barrier / St) / sigma + log(x / St) / sigma - mu*(tau), 2) / (2.0*(tau))) / (sigma*sigma*pow(St, 3)*pow(T - t, 1.5));
	df -= 2.0*sqrt(2.0 / M_PI)*mu*barrier*pow(barrier / St, 2 * mu / sigma - 1.0)*exp(-pow(-2.0*log(barrier / St) / sigma + log(x / St) / sigma - mu*(tau), 2) / (2.0*tau)) / (sigma*pow(St, 3)*sqrt(tau));
	df -= sqrt(2.0 / M_PI)*mu*barrier*pow(barrier / St, 2.0*mu / sigma - 1.0)*(-2.0*log(barrier / St) / sigma + log(x / St) / sigma - mu*(tau))*exp(-pow(-2.0*log(barrier / St) / sigma + log(x / St) / sigma - mu*(tau), 2) / (2.0*(tau))) / (sigma*sigma*pow(T - t, 1.5)*pow(St, 3));
	df -= pow(barrier / St, 2 * mu / sigma)*pow(-2.0*log(barrier / St) / sigma + log(x / St) / sigma - mu*(tau), 2)*exp(-pow(-2.0*log(barrier / St) / sigma + log(x / St) / sigma - mu*(tau), 2) / (2.0*(tau))) / (sqrt(2 * M_PI*pow(T - t, 5))*sigma*sigma*St*St);
	df += pow(barrier / St, 2 * mu / sigma)*exp(-pow(-2.0*log(barrier / St) / sigma + log(x / St) / sigma - mu*(tau), 2) / (2.0*(tau))) / (sqrt(2 * M_PI*pow(T - t, 3))*sigma*sigma*St*St);
	df -= pow(barrier / St, 2 * mu / sigma)*(-2.0*log(barrier / St) / sigma + log(x / St) / sigma - mu*(tau))*exp(-pow(-2.0*log(barrier / St) / sigma + log(x / St) / sigma - mu*(tau), 2) / (2.0*(tau))) / (sqrt(2 * M_PI*pow(T - t, 3))*sigma*St*St);
	return exp(-r*(tau))*(x - K)*df / f;
}
double barrier_option::gamma_lr(){
	double result = 0.0;
	if (barrier < St){
		for (int i = 0; i < number_iterations; i++){
			result += gamma_lr_step(barrier_down());
		}
	}
	else{
		for (int i = 0; i < number_iterations; i++){
			result += gamma_lr_step(barrier_up());
		}
	}

	return result / number_iterations;
}
double barrier_option::gamma_theoretic(){
	if (barrier < St){
		double mu = r - sigma*sigma / 2.0;
		double d = (log(St / K) + (r + sigma*sigma / 2)*(tau)) / (sigma*sqrt(tau));
		double d_calculate_minus = (log(barrier*barrier / (St*K)) + mu*(tau)) / (sigma*sqrt(tau)) + sigma*sqrt(tau);
		double delta = -barrier*barrier / (St*St)*normal_cdf(d_calculate_minus);
		double gamma = (barrier*barrier / pow(St, 3))*(2.0*normal_cdf(d_calculate_minus) + normal_pdf(d_calculate_minus) / (sigma*sqrt(tau)));
		double result = normal_pdf(d) / (St*sigma*sqrt(tau));
		result -= pow(barrier / St, 2.0*mu / (sigma*sigma))*(2.0*mu*(2.0*mu + sigma*sigma)*payoff_theoretic(barrier*barrier / St, K) / (St*St*pow(sigma, 4)) + gamma - 4.0*mu / (St*sigma*sigma)*delta);
		return result;
	}
	else {
		double mu = r - sigma*sigma / 2.0;
		double result = gamma_theoretic_call(St, K) - gamma_theoretic_call(St, barrier) - (barrier - K)*exp(-r*(tau))*dnorm(d_calculate_minus(St, barrier)) / (sigma*sigma*St*St*(tau)) + (barrier - K)*exp(-r*(tau))*normal_pdf(d_calculate_minus(St, barrier)) / (sigma*St*St*sqrt(tau));
		result -= (4.0*mu*(mu + 0.5*sigma*sigma)*pow(barrier / St, 2.0*mu / sigma*sigma) / (pow(sigma, 4)*St*St))*(payoff_theoretic(barrier*barrier / St, K) - payoff_theoretic(barrier*barrier / St, barrier) - (barrier - K)*exp(-r*(tau))*normal_cdf(d_calculate_minus(barrier, St)));
		result += (4.0*mu*pow(barrier / St, 2.0*mu / pow(sigma, 2)) / (pow(sigma, 2)*St))*(-barrier*barrier / (St*St)*delta_theoretic_call(barrier*barrier / St, K) + barrier*barrier / (St*St)*delta_theoretic_call(barrier*barrier / St, barrier) + (barrier - K)*exp(-r*(tau))*normal_pdf(d_calculate_minus(barrier, St)) / (sigma*St*sqrt(tau)));
		result -= pow(barrier / St, 2.0*mu / pow(sigma, 2))*(pow(barrier / St, 4)*gamma_theoretic_call(barrier*barrier / St, K) +2.0*barrier*barrier/pow(St,3)*delta_theoretic_call(barrier*barrier/St,K)- pow(barrier / St, 4)*gamma_theoretic_call(barrier*barrier / St, barrier)-2.0*barrier*barrier/pow(St,3)*delta_theoretic_call(barrier*barrier/St,barrier) - (barrier - K)*exp(-r*(tau))*dnorm(d_calculate_minus(barrier, St)) / (sigma*sigma*St*St*(tau)) - (barrier - K)*exp(-r*(tau))*normal_pdf(d_calculate_minus(barrier, St)) / (sigma*St*St*sqrt(tau)));
		return result;
	}
}

double barrier_option::vega_lr_step(double x){
	if (x == 0 )return 0.0;
	double mu = (r - sigma*sigma / 2) / sigma;
	double f = (normal_pdf((log(x / St) - sigma*mu*(tau)) / (sigma*sqrt(tau))) - exp(2 * mu*log(barrier / St) / sigma)*normal_pdf((log(x / St) - 2 * log(barrier / St) - sigma*mu*(tau)) / (sigma*sqrt(tau)))) / sqrt(tau);
	f /= (sigma*x);
	double df = -exp(-pow(log(x / St) / sigma - (tau)*(r / sigma - sigma / 2.0), 2) / (2.0*(tau)))*((tau)*(r / pow(sigma, 2) + 0.5) - log(x / St) / pow(sigma, 2))*(log(x / St) / sigma - (tau)*(r / sigma - sigma / 2.0)) / (sqrt(2.0*M_PI*pow(T - t, 3))*sigma*x);
	df -= exp(-pow(log(x / St) / sigma - (tau)*(r / sigma - sigma / 2.0), 2) / (2.0*(tau))) / (sqrt(2.0*M_PI*(tau))*pow(sigma, 2)*x);
	df += (4.0*r*log(barrier / St) / pow(sigma, 3) + (2.0*log(barrier / St) / pow(sigma, 2) + (tau)*(r / pow(sigma, 2) + 0.5) - log(x / St) / pow(sigma, 2))*(-2.0*log(barrier / St) / sigma + (tau)*(-r / sigma + sigma / 2.0) + log(x / St) / sigma) / (tau))*exp(2.0*(r / pow(sigma, 2) - 0.5)*log(barrier / St) - pow(-2.0*log(barrier / St) / sigma + (tau)*(-r / sigma + sigma / 2.0) + log(x / St) / sigma, 2) / (2.0*(tau))) / (sqrt(2.0*M_PI*(tau))*sigma*x);
	df += exp(2.0*(r / pow(sigma, 2) - 0.5)*log(barrier / St) - pow(-2.0*log(barrier / St) / sigma + (tau)*(-r / sigma + sigma / 2.0) + log(x / St) / sigma, 2) / (2.0*(tau))) / (sqrt(2.0*M_PI*(tau))*pow(sigma, 2)*x);
	return exp(-r*(tau))*(x - K)*df / f;
}


double barrier_option::vega_lr(){
	double result = 0.0;
	if (barrier < St){
		for (int i = 0; i < number_iterations; i++){
			result += vega_lr_step(barrier_down());
		}
	}
	else{
		for (int i = 0; i < number_iterations; i++){
			result += vega_lr_step(barrier_up());
		}
	}
	return result / number_iterations;
}
double barrier_option::vega_theoretic(){
	if (barrier < St){
		double mu = r - sigma*sigma / 2.0;
		double d = (log(St / K) + (r + sigma*sigma / 2)*(tau)) / (sigma*sqrt(tau));
		double d_calculate_minus = (log(barrier*barrier / (St*K)) + mu*(tau)) / (sigma*sqrt(tau));
		double vega = sqrt(tau)*K*exp(-r*(tau))*normal_pdf(d_calculate_minus);
		double result = St*normal_pdf(d)*sqrt(tau);
		result -= pow(barrier / St, 2.0*mu / (sigma*sigma))*(vega - 4.0*r*payoff_theoretic(barrier*barrier / St, K)*log(barrier / St) / pow(sigma, 3));
		return result;
	}
	else{
		double mu = r - sigma*sigma / 2.0;
		double result = vega_theoretic_call(St, K) - vega_theoretic_call(St, barrier) + (barrier - K)*exp(-r*(tau))*normal_pdf(d_calculate_minus(St, barrier))*((tau)*(sigma*sigma + mu) + log(St / barrier)) / (sigma*sigma*sqrt(tau));
		result += 4.0*r*log(barrier / St)*pow(barrier / St, 2 * r / pow(sigma, 2) - 1) / pow(sigma, 3)*(payoff_theoretic(barrier*barrier / St, K) - payoff_theoretic(barrier*barrier / St, barrier) - (barrier - K)*exp(-r*(tau))*normal_cdf(d_calculate_minus(barrier, St)));
		result -= pow(barrier / St, 2.0*mu / pow(sigma, 2))*(sqrt(tau)*exp(-r*(tau))*(K*normal_pdf(d_calculate_minus(barrier*barrier / St, K)) - barrier*normal_pdf(d_calculate_minus(barrier*barrier / St, barrier))) + (barrier - K)*exp(-r*(tau))*normal_pdf(d_calculate_minus(barrier, St))*((tau)*(sigma*sigma + mu) + log(barrier / St)) / (sigma*sigma*sqrt(tau)));
		return result;
	}
	return 0.0;
}



double barrier_option::price(){
  return barrier_option_price();
}

double barrier_option::delta(std::string method){
    if (method=="th") {return delta_theoretic();}
    return delta_lr();
}
double barrier_option::delta(){return delta("lr");}

double barrier_option::gamma(){return gamma("lr");}

double barrier_option::gamma(std::string method){
    if (method=="th") {return gamma_theoretic();}
    return gamma_lr();
}

double barrier_option::vega(){return vega("lr");}

double barrier_option::vega(std::string method){
    if (method=="th") {return vega_theoretic();}
    return vega_lr();
}
