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
	return payoff_theoretic(St, K) - payoff_theoretic(St, barrier) - (barrier - K)*exp(-r*(tau))*normal_cdf(d_calculate(St, barrier)) - pow(barrier / St, 2.0 * mu / (sigma*sigma))*(payoff_theoretic(barrier*barrier / St, K) - payoff_theoretic(barrier*barrier / St, barrier) - (barrier - K)*exp(-r*(tau))*normal_cdf(d_calculate(barrier, St)));
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
		res = normal_cdf(d_calculate(St, K)) - normal_cdf(d_calculate(St, barrier));
		res -= (barrier - K)*exp(-r*(tau))*normal_pdf(stock_price_single(St,barrier))/(sigma*St*sqrt(T-t));
		res += (2.0*mu*pow(barrier / St, 2.0*mu / pow(sigma, 2)) / (pow(sigma, 2)*St))*(stock_price_single(barrier*barrier / St, K) - stock_price_single(barrier*barrier / St, barrier) - (barrier - K)*exp(-r*(tau))*normal_cdf(stock_price_single(barrier, St)));
		res -= pow(barrier / St, 2.0*mu / pow(sigma, 2))*(-barrier*barrier / (St*St)*normal_cdf(d_calculate(barrier*barrier / St, K)) + barrier * barrier / (St*St)*normal_cdf(d_calculate(barrier * barrier / St, barrier)) + (barrier - K)*exp(-r*(tau))*normal_pdf(stock_price_single(barrier, St)) / (sigma*St*sqrt(tau)));
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
    //if (method=="th") {return gamma_theoretic();}
    return gamma_lr();
}
