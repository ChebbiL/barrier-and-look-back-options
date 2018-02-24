# User Manual

## Features

* Compute the European call option price
* Compute the European call option greeks:
  * delta
  * gamma
  * vega


## Usage

Before you use any of the functions of the package, you must create a `call` object to interact with. This can be done using the following method. You will need the following parameters:
* `S_t0` the value of the underlying stock at time t0. Default is `100`.
* `strike` the value of the strike of the European option. Default is `100`.
* `interest_rate` the decimal value of the interest rate. For an interest rate of 15%, enter `0.15`.  Default is `0.05`.
* `volatility` the decimal value of the volatility of the underlying stock. For an volatility of 15%, enter `0.15`.  Default is `0.40`.
* `time_maturity_T` the time (double) of maturity of the call. Default value is `1
.0`.
* `initial_time_t0` the initial time (double) at which `S_t0` was recorded. Default value is `0.0`.
```
call my_european_call_option(S_t0, strike, interest_rate, volatility, time_maturity_T, time_initial_t);
```

You can now access all the functionalities of the program by accessing the functions as follows `my_european_call_option.the_function_i_use()`.

### Call Option price
In order to compute the call option price, you can use the following method. By default (with any string input or no input at all), the method used is _likelihood ratios_.
```
double my_european_call_option.delta();
```
You can alsospecify the method you want to use.
- To use the _pathwise estimates_ method, enter the argument `"pw"`.
- To use the _likelihood ratios_ method, enter the argument `"lr"` or alternatively do not enter any argument. 
```
double my_european_call_option.delta("pw");
double my_european_call_option.delta("lr");
```
