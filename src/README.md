# User Manual

This package allows you to perform various computations on European call options and on barrier options in the `C++` language.

## Features

* Compute the European call option price
* Compute the European call option greeks:
  * delta
  * gamma
  * vega
* Provide statistical informations

## Installation

In order to compile and run the project, you can use the following command in the `src` folder.

```
make
./run
```

### Initialising the European call option

Before you use any of the functions of the package, you must create a `call` object to interact with. This can be done using the following method. You will need the following parameters:
* `S_t0 (double)` the value of the underlying stock at time t0. Default is `100`.
* `strike (double)` the value of the strike of the European option. Default is `100`.
* `interest_rate (double)` the decimal value of the interest rate. For an interest rate of 15%, enter `0.15`.  Default is `0.05`.
* `volatility (double)` the decimal value of the volatility of the underlying stock. For an volatility of 15%, enter `0.15`.  Default is `0.40`.
* `time_maturity_T (double)` the time of maturity of the call. Default value is `1
.0`.
* `initial_time_t0 (double)` the initial time at which `S_t0` was recorded. Default value is `0.0`.
* `number_iterations_approximation (int)` the number of iterations for the approximation methods. The default value is 10,000. Note than the larger this number is, the slower but more accurate computations are. Industry standards are 100,000.
```
european_option call(S_t0, strike, interest_rate, volatility, time_maturity_T, time_initial_t, number_iterations_approximation);
```

### Initialising the Barrier option

---

## Usage

### European call option

You can access most the functionalities of the program by accessing the functions as follows `my_european_call_option.the_function_i_use()`.

#### Price

The price is computed using the Monte Carlo method with the number of iterations specified when initialising the `european_call_option` class. Recall that the value is 10,000 iterations by default.
```
double call.price();
```
You can however specify any other number of iterations by specifying an `int` number in the parameters as follows.
```
double call.price(100000);
```

#### Delta
In order to compute the call option delta, you can use the following method. By default (with any string input or no input at all), the method used is _likelihood ratios_.
```
double call.delta();
```
You can alsospecify the method you want to use.
- To use the _pathwise estimates_ method, enter the argument `"pw"`.
- To use the _likelihood ratios_ method, enter the argument `"lr"` or alternatively do not enter any argument.
```
double call.delta("pw");
double call.delta("lr");
```

#### Gamma
In order to compute the call option gamma, you can use the following method. By default (with any string input or no input at all), the method used is _pathwise estimates - likelihood ratios_.
```
double call.gamma();
```
You can alsospecify the method you want to use.
- To use the _likelihood ratios - pathwise estimates_ method, enter the argument `"lrpw"`.
- To use the _likelihood ratios - likelihood ratios_ method, enter the argument `"lrlr"`.
- To use the _pathwise estimates - likelihood ratios_ method, enter the argument `"pwlr"` or alternatively do not enter any argument.
```
double call.gamma("lrpw");
double call.gamma("lrlr");
double call.gamma("pwlr");
```

### Barrier option
