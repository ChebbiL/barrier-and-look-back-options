#include "barrier_option.hpp"

double barrier_option::barrier_down(){
  double current_value = St;
  double minimum_price = St;
  double h = tau / number_iterations;
  for (int i = 0; i < number_iterations; i++){
    current_value = stock_price_single(t + (i + 1) * h, t + i * h);
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
    current_value = stock_price_single(t + (i + 1) * h, t + i * h);
    if (maximum_price < current_value) maximum_price = current_value;
    if (maximum_price > barrier) return 0;
  }
  if (current_value > K) return current_value;
  return 0.0;
}

double barrier_option::d_calculate(double x, double y){
	return (log(x / y) + (r + sigma * sigma / 2.0) * tau) / (sigma * sqrt(tau));
}
double barrier_option::d_calculate_minus(double x, double y){
	return (log(x / y) + (r - sigma * sigma / 2.0) * tau) / (sigma * sqrt(tau));
}


double barrier_option::payoff_theoretic(double initial_stock_price, double strike){
  double d_reference = (log(initial_stock_price / strike) + (r + sigma*sigma / 2)*tau) / (sigma*sqrt(tau));
	return normal_cdf(d_reference) * initial_stock_price - normal_cdf(d_reference - sigma*sqrt(tau)) * strike * discount;
}

double barrier_option::barrier_down_theoric(){
	double mu = r - sigma * sigma / 2.0;
	return payoff_theoretic(St, K) - pow(St / barrier, -2.0 * mu / (sigma*sigma))*payoff_theoretic(barrier * barrier / St, K);
}
double barrier_option::barrier_up_theoric(){
	if (barrier < K) return 0;
	double mu = r - sigma * sigma / 2.0;
	return payoff_theoretic(St, K) - payoff_theoretic(St, barrier) - (barrier - K)*exp(-r*(tau))*normal_cdf(d_calculate_minus(St, barrier)) - pow(barrier / St, 2.0 * mu / (sigma*sigma))*(payoff_theoretic(barrier*barrier / St, K) - payoff_theoretic(barrier*barrier / St, barrier) - (barrier - K)*exp(-r*(tau))*normal_cdf(d_calculate_minus(barrier, St)));
}


double barrier_option::delta_lr_step(double current_value){
	if (current_value == 0) return .0;
	double mu = (r - sigma*sigma / 2) / sigma;

	double f = (normal_pdf((log(current_value / St) - sigma*mu*(tau)) / (sigma*sqrt(tau))) - exp(2 * mu*log(barrier / St) / sigma)*normal_pdf((log(current_value / St) - 2 * log(barrier / St) - sigma*mu*(tau)) / (sigma*sqrt(tau)))) / sqrt(tau);

	double df = exp(-0.5*pow((log(current_value / St) / sigma - mu*(tau)), 2) / (tau))*(log(current_value / St) / sigma - mu*(tau)) / (sqrt(2 * M_PI*pow((tau), 3))*sigma*St);

	df += sqrt(2 / M_PI) *mu*barrier*pow(barrier / St, (2 * mu / sigma - 1))*exp(-pow((-2 * log(barrier / St) / sigma + log(current_value / St) / sigma - mu*(tau)), 2) / (2 * (tau))) / (sigma*St*St*sqrt(tau));

	df += pow(barrier / St, (2 * mu / sigma))*(-2 * log(barrier / St) / sigma + log(current_value / St) / sigma - mu*(tau))*exp(-pow((-2 * log(barrier / St) / sigma + log(current_value / St) / sigma - mu*(tau)), 2) / (2 * (tau))) / (sqrt(2 * M_PI*pow((tau), 3)) * sigma * St);

	return exp(-r*(tau))*(current_value - K)*df / f;
}

double barrier_option::gamma_lr_step(double current_value){
	if (current_value == 0)return 0.0;
	double mu = (r - sigma*sigma / 2.0) / sigma;
	double f = (normal_pdf((log(current_value / St) - sigma*mu*(tau)) / (sigma*sqrt(tau))) - exp(2 * mu*log(barrier / St) / sigma)*normal_pdf((log(current_value / St) - 2 * log(barrier / St) - sigma*mu*(tau)) / (sigma*sqrt(tau)))) / sqrt(tau);
	double df = exp(-pow((log(current_value / St) / sigma - mu*(tau)), 2) / (2.0*(tau)))*pow(log(current_value / St) / sigma - mu*(tau), 2) / (sigma*sigma*St*St*pow((tau), 2.5)*sqrt(2*M_PI));
	df -= exp(-pow(log(current_value / St) / sigma - mu*(tau), 2) / (2.0*(tau))) / (sqrt(2.0 * M_PI*pow((tau), 3))*sigma*sigma*St*St);
	df -= exp(-pow(log(current_value / St) / sigma - mu*(tau), 2) / (2.0*(tau)))*(log(current_value / St) / sigma - mu*(tau)) / (sqrt(2.0 * M_PI*pow((tau), 3))*sigma*St*St);
	df -= sqrt(2.0 / M_PI)*mu*barrier*barrier*(2.0 * mu / sigma - 1.0)*pow(barrier / St, 2.0 * mu / sigma - 2.0)*exp(-pow(-2.0*log(barrier / St) / sigma + log(current_value / St) / sigma - mu*(tau), 2) / (2.0*(tau))) / (sigma*pow(St, 4)*sqrt(tau));
	df -= sqrt(2.0 / M_PI)*mu*barrier*pow(barrier / St, 2.0*mu / sigma - 1.0)*(-2.0*log(barrier / St) / sigma + log(current_value / St) / sigma - mu*(tau))*exp(-pow(-2.0*log(barrier / St) / sigma + log(current_value / St) / sigma - mu*(tau), 2) / (2.0*(tau))) / (sigma*sigma*pow(St, 3)*pow(tau, 1.5));
	df -= 2.0*sqrt(2.0 / M_PI)*mu*barrier*pow(barrier / St, 2 * mu / sigma - 1.0)*exp(-pow(-2.0*log(barrier / St) / sigma + log(current_value / St) / sigma - mu*(tau), 2) / (2.0*(tau))) / (sigma*pow(St, 3)*sqrt(tau));
	df -= sqrt(2.0 / M_PI)*mu*barrier*pow(barrier / St, 2.0*mu / sigma - 1.0)*(-2.0*log(barrier / St) / sigma + log(current_value / St) / sigma - mu*(tau))*exp(-pow(-2.0*log(barrier / St) / sigma + log(current_value / St) / sigma - mu*(tau), 2) / (2.0*(tau))) / (sigma*sigma*pow(tau, 1.5)*pow(St, 3));
	df -= pow(barrier / St, 2 * mu / sigma)*pow(-2.0*log(barrier / St) / sigma + log(current_value / St) / sigma - mu*(tau), 2)*exp(-pow(-2.0*log(barrier / St) / sigma + log(current_value / St) / sigma - mu*(tau), 2) / (2.0*(tau))) / (sqrt(2 * M_PI*pow(tau, 5))*sigma*sigma*St*St);
	df += pow(barrier / St, 2 * mu / sigma)*exp(-pow(-2.0*log(barrier / St) / sigma + log(current_value / St) / sigma - mu*(tau), 2) / (2.0*(tau))) / (sqrt(2 * M_PI*pow(tau, 3))*sigma*sigma*St*St);
	df -= pow(barrier / St, 2 * mu / sigma)*(-2.0*log(barrier / St) / sigma + log(current_value / St) / sigma - mu*(tau))*exp(-pow(-2.0*log(barrier / St) / sigma + log(current_value / St) / sigma - mu*(tau), 2) / (2.0*(tau))) / (sqrt(2 * M_PI*pow(tau, 3))*sigma*St*St);
	return exp(-r*(tau))*(current_value - K)*df / f;
}








double barrier_option::barrier_option_price(){
  double result = .0;
  double current_value;
  for (int i = 0; i < number_iterations; i++){
    barrier < St ? current_value = barrier_down() : current_value = barrier_down();
    if (current_value!=0) result += discount * (current_value - K);
  }
  return result / number_iterations;
}

double barrier_option::payoff_theoretic() {
  double result;
  (barrier < St) ? result = barrier_down_theoric() : result = barrier_up_theoric();
  return result;
}

double barrier_option::delta_lr(){
  double result = .0;
  if (barrier < St){
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

double barrier_option::delta_theoric(){
	double mu = r - sigma*sigma / 2.0;
	double res = 0.0;
	if (barrier < St){
		res = normal_cdf((log(St / K) + (r + sigma*sigma / 2.0)*(tau)) / (sigma*sqrt(tau)));
		res -= pow((barrier / St), (r / (sigma*sigma) - 1))*(-barrier*barrier* normal_cdf((log(barrier*barrier / (St*K)) + mu*(tau)) / (sigma*sqrt(tau)) + sigma*sqrt(tau)) / (St*St) - 2.0 * mu*stock_price_single(barrier*barrier / St, K) / (St*sigma*sigma));

	}
	else{
		res = normal_cdf(d_calculate(St, K)) - normal_cdf(d_calculate_minus(St, barrier));
		res -= (barrier - K)*exp(-r*(tau))*normal_pdf(d_calculate(St,barrier))/(sigma*St*sqrt(T-t));
		res += (2.0*mu*pow(barrier / St, 2.0*mu / pow(sigma, 2)) / (pow(sigma, 2)*St))*(stock_price_single(barrier*barrier / St, K) - stock_price_single(barrier*barrier / St, barrier) - (barrier - K)*exp(-r*(tau))*normal_cdf(d_calculate_minus(barrier, St)));
		res -= pow(barrier / St, 2.0*mu / pow(sigma, 2))*(-barrier*barrier / (St*St)*normal_cdf(d_calculate(barrier*barrier / St, K)) + barrier * barrier / (St*St)*normal_cdf(d_calculate(barrier * barrier / St, barrier)) + (barrier - K)*exp(-r*(tau))*normal_pdf(d_calculate_minus(barrier, St)) / (sigma*St*sqrt(tau)));
	}
	return res;
}



double barrier_option::gamma_lr(){
	double result = .0;
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

double barrier_option::gamma_theoretic_dnorm(double x){
  return -exp(-0.5*x*x)*x / sqrt(2.0*M_PI);
}


double barrier_option::gamma_theoretic(){
	if (barrier < St){
		double mu = r - sigma*sigma / 2.0;
		double d = (log(St / K) + (r + sigma*sigma / 2)*tau) / (sigma*sqrt(tau));
		double d_calculate = (log(barrier*barrier / (St*K)) + mu*tau) / (sigma*sqrt(tau)) + sigma*sqrt(tau);
		double delta = -barrier*barrier / (St*St)*normal_cdf(d_calculate);
		double gamma = (barrier*barrier / pow(St, 3))*(2.0*normal_cdf(d_calculate) + normal_pdf(d_calculate) / (sigma*sqrt(tau)));
		double result = normal_pdf(d) / (St*sigma*sqrt(tau));
		result -= pow(barrier / St, 2.0*mu / (sigma*sigma))*(2.0*mu*(2.0*mu + sigma*sigma)*stock_price_single(barrier*barrier / St, K) / (St*St*pow(sigma, 4)) + gamma - 4.0*mu / (St*sigma*sigma)*delta);
		return result;
	}
	else {
		double mu = r - sigma*sigma / 2.0;
		double result = gamma_theoretic_call(St, K) - gamma_theoretic_call(St, barrier) - (barrier - K)*discount*gamma_theoretic_dnorm(d_calculate_minus(St, barrier)) / (sigma*sigma*St*St*tau) + (barrier - K)*discount*normal_pdf(d_calculate_minus(St, barrier)) / (sigma*St*St*sqrt(tau));
		result -= (4.0*mu*(mu + 0.5*sigma*sigma)*pow(barrier / St, 2.0*mu / sigma*sigma) / (pow(sigma, 4)*St*St))*(stock_price_single(barrier*barrier / St, K) - stock_price_single(barrier*barrier / St, barrier) - (barrier - K)*discount*normal_cdf(d_calculate_minus(barrier, St)));
		result += (4.0*mu*pow(barrier / St, 2.0*mu / pow(sigma, 2)) / (pow(sigma, 2)*St))*(-barrier*barrier / (St*St)*delta_theoretic_call(barrier*barrier / St, K) + barrier*barrier / (St*St)*delta_theoretic_call(barrier*barrier / St, barrier) + (barrier - K)*discount*normal_pdf(d_calculate_minus(barrier, St)) / (sigma*St*sqrt(tau)));
		result -= pow(barrier / St, 2.0*mu / pow(sigma, 2))*(pow(barrier / St, 4)*gamma_theoretic_call(barrier*barrier / St, K) +2.0*barrier*barrier/pow(St,3)*delta_theoretic_call(barrier*barrier/St,K)- pow(barrier / St, 4)*gamma_theoretic_call(barrier*barrier / St, barrier)-2.0*barrier*barrier/pow(St,3)*delta_theoretic_call(barrier*barrier/St,barrier) - (barrier - K)*discount*gamma_theoretic_dnorm(d_calculate_minus(barrier, St)) / (sigma*sigma*St*St*tau) - (barrier - K)*discount*normal_pdf(d_calculate_minus(barrier, St)) / (sigma*St*St*sqrt(tau)));
		return result;
	}
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



double barrier_option::vega_lr_step(double current_value){
	if (current_value == 0 ) return 0.0;
	double mu = (r - sigma*sigma / 2) / sigma;
	double f = (normal_pdf((log(current_value / St) - sigma*mu*(tau)) / (sigma*sqrt(tau))) - exp(2 * mu*log(barrier / St) / sigma)*normal_pdf((log(current_value / St) - 2 * log(barrier / St) - sigma*mu*(tau)) / (sigma*sqrt(tau)))) / sqrt(tau);
	f /= (sigma*current_value);
	double df = -exp(-pow(log(current_value / St) / sigma - (tau)*(r / sigma - sigma / 2.0), 2) / (2.0*(tau)))*((tau)*(r / pow(sigma, 2) + 0.5) - log(current_value / St) / pow(sigma, 2))*(log(current_value / St) / sigma - (tau)*(r / sigma - sigma / 2.0)) / (sqrt(2.0*M_PI*pow(tau, 3))*sigma*current_value);
	df -= exp(-pow(log(current_value / St) / sigma - (tau)*(r / sigma - sigma / 2.0), 2) / (2.0*(tau))) / (sqrt(2.0*M_PI*(tau))*pow(sigma, 2)*current_value);
	df += (4.0*r*log(barrier / St) / pow(sigma, 3) + (2.0*log(barrier / St) / pow(sigma, 2) + (tau)*(r / pow(sigma, 2) + 0.5) - log(current_value / St) / pow(sigma, 2))*(-2.0*log(barrier / St) / sigma + (tau)*(-r / sigma + sigma / 2.0) + log(current_value / St) / sigma) / (tau))*exp(2.0*(r / pow(sigma, 2) - 0.5)*log(barrier / St) - pow(-2.0*log(barrier / St) / sigma + (tau)*(-r / sigma + sigma / 2.0) + log(current_value / St) / sigma, 2) / (2.0*(tau))) / (sqrt(2.0*M_PI*(tau))*sigma*current_value);
	df += exp(2.0*(r / pow(sigma, 2) - 0.5)*log(barrier / St) - pow(-2.0*log(barrier / St) / sigma + (tau)*(-r / sigma + sigma / 2.0) + log(current_value / St) / sigma, 2) / (2.0*(tau))) / (sqrt(2.0*M_PI*(tau))*pow(sigma, 2)*current_value);
	return exp(-r*(tau))*(current_value - K)*df / f;
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
		double d_calculate = (log(barrier*barrier / (St*K)) + mu*(tau)) / (sigma*sqrt(tau));
		double vega = sqrt(tau)*K*exp(-r*(tau))*normal_pdf(d_calculate);
		double result = St*normal_pdf(d)*sqrt(tau);
		result -= pow(barrier / St, 2.0*mu / (sigma*sigma))*(vega - 4.0*r*payoff_theoretic(barrier*barrier / St, K)*log(barrier / St) / pow(sigma, 3));
		return result;
	}
	else{
		double mu = r - sigma*sigma / 2.0;
		double result = vega_theoretic_call(St, K) - vega_theoretic_call(St, barrier) + (barrier - K)*exp(-r*(tau))*normal_pdf(d_calculate(St, barrier))*((tau)*(sigma*sigma + mu) + log(St / barrier)) / (sigma*sigma*sqrt(tau));
		result += 4.0*r*log(barrier / St)*pow(barrier / St, 2 * r / pow(sigma, 2) - 1) / pow(sigma, 3)*(payoff_theoretic(barrier*barrier / St, K) - payoff_theoretic(barrier*barrier / St, barrier) - (barrier - K)*exp(-r*(tau))*normal_cdf(d_calculate(barrier, St)));
		result -= pow(barrier / St, 2.0*mu / pow(sigma, 2))*(sqrt(tau)*exp(-r*(tau))*(K*normal_pdf(d_calculate(barrier*barrier / St, K)) - barrier*normal_pdf(d_calculate(barrier*barrier / St, barrier))) + (barrier - K)*exp(-r*(tau))*normal_pdf(d_calculate(barrier, St))*((tau)*(sigma*sigma + mu) + log(barrier / St)) / (sigma*sigma*sqrt(tau)));
		return result;
	}
	return 0.0;
}
