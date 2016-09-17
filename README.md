# Codes

A repository for the paper "Granularity, Network Asymmetry, and Aggregate Volatility". Codes reproduce full analysis in the paper, but with fake data. Real data can be accessed through [http://www.statcan.gc.ca/eng/cder/data](CDER) at Statistics Canada.

## To use

The general strategy: generate fake plant characteristics, use those to calculate region characteristics, then generate demand shares that are consistent with those characteristics. Generate trade costs and other relevant parameters. From there, you can start the analysis. 

Given observed demand characteristics, and the elasticity of substitution \sigma, calculate \Lambda and \Gamma. Then see how closely idiosyncratic shocks to demand and productivity can match observed covariance matrix and aggregate volatility.

### Data

* ```R``` number of regions
* ```N``` number of plants
* ```\beta``` vector of value-added shares, ```Nx1```
* ```s``` vector of output, ```Nx1```
* ```z``` vector of productivities, ```Nx1```
* ```ir``` data.frame containing a plant's id and associated region, ```Nx2```, or ```NxR``` in sparse matrix form
* ```Er``` matrix of extensive region-plant demand characteristics, ```RxN```
* ```En``` matrix of extensive plant-plant demand characteristics, ```NxN```

### To solve

* ```\Lambda``` matrix of region-plant demand characteristics, ```RxN```
* ```\Gamma``` matrix of plant-plant demand characteristics, ```NxN```


## Organization:

* ```main.R``` sets up the environment. It loads the functions in the ```./R``` subdirectory.
* ```helpers.R``` has useful functions to convert data frames to sparse matrices and back.
* ```initialize_functions.R``` includes the functions that create fake data. Feel free to change the parameters to get more or less sparse functions, more or less plants, more or less variance in plant parameters, and so on.
* ```benchmark.R``` given certain parameters, calculate observed region-plant demand shares ```A``` and plant-plant demand shares ```G``` that are consistent with the parameters.
* ```solve_lambda_gamma.R``` given the data, calculate unobserved ```\Lambda``` and ```\Gamma``` that are consistent with the data.
* ```???``` Unnamed function to calibrate the elasticity of substitution
* ```???``` Unnamed function to calibrate the variances of idiosyncratic shocks that match the observed covariance matrix and aggregate volatility.




