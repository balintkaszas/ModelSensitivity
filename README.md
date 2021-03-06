# Model Sensitivity


## License

This software is made public for research use only. It may be modified and redistributed under the terms of the GNU General Public License.



## Motivation

Calculation of bounds on prediction errors in dynamical systems subject to modeling uncertainties. 
The theory is discussed in the manuscript B. K., G. Haller, __Universal Upper Estimate for Prediction Errors under Moderate Model Uncertainty__ (submitted) [1].

Given a known dynamical system in the form 
<p align="center"><img src="svgs/8ae26bb45da4a1da385d8831c338005d.svg?invert_in_darkmode" align=middle width=393.49880625pt height=16.438356pt/></p>

we can quantify the prediction errors in the presence of modeling errors of the form

<p align="center"><img src="svgs/ae9ecea711fefb62a659a2baba56db9f.svg?invert_in_darkmode" align=middle width=170.88747059999997pt height=16.438356pt/></p>

where <img src="svgs/749e5d74b13d610763eb7779886568ef.svg?invert_in_darkmode" align=middle width=26.66961824999999pt height=24.65753399999998pt/> represents white noise. These bounds are expressed as functions of the invariants of the Cauchy-Green strain tensor of the idealized model:
<p align="center"><img src="svgs/1459bdf98910ccdf3ce84c76ea32948e.svg?invert_in_darkmode" align=middle width=224.35356899999996pt height=23.5253469pt/></p>

Its leading eigenvalue is denoted by <img src="svgs/b68b7fd0730749f66010e7a278f16778.svg?invert_in_darkmode" align=middle width=47.64853334999999pt height=26.085962100000025pt/>.

If <img src="svgs/bcaccb82f389300ad559e9be322690d3.svg?invert_in_darkmode" align=middle width=15.94753544999999pt height=26.76175259999998pt/> is the idealized solution of <img src="svgs/a8cf3097890381cfa1311a7a0838c48f.svg?invert_in_darkmode" align=middle width=82.15744514999999pt height=24.65753399999998pt/> and <img src="svgs/84614bcbba656f441da3dfdd6c082d31.svg?invert_in_darkmode" align=middle width=15.60562409999999pt height=21.839370299999988pt/> is the real solution, the leading order error is bounded by

<p align="center"><img src="svgs/b0e9dffda3866d1b9edeaf4fc6ee3a62.svg?invert_in_darkmode" align=middle width=611.17729035pt height=45.361174649999995pt/></p>

The quantities <img src="svgs/ea26fdf44c4f592467343f74a20a06fa.svg?invert_in_darkmode" align=middle width=26.80375499999999pt height=22.465723500000017pt/> and <img src="svgs/855ac9885de30be8c9d1359561c892b3.svg?invert_in_darkmode" align=middle width=26.80375499999999pt height=22.465723500000017pt/> denote the maximal value of the modeling errors 
<p align="center"><img src="svgs/a2bee8faec81269c3230c77d7966165a.svg?invert_in_darkmode" align=middle width=165.80475779999998pt height=16.438356pt/></p>

<p align="center"><img src="svgs/702ea200c58f41d105b95ee5e2e2f2ab.svg?invert_in_darkmode" align=middle width=230.1787488pt height=18.7598829pt/></p>

Then, the system's sensitivity with respect to modeling errors can be characterized by the following scalar over the time interval <img src="svgs/1036fd308dc681d0e5922e143269d0ce.svg?invert_in_darkmode" align=middle width=35.68498559999999pt height=24.65753399999998pt/>, at the point <img src="svgs/e714a3139958da04b41e3e607a544455.svg?invert_in_darkmode" align=middle width=15.94753544999999pt height=14.15524440000002pt/>, 
<p align="center"><img src="svgs/f6182ea3ade85b84cc667b3a0e098fd7.svg?invert_in_darkmode" align=middle width=379.6676433pt height=45.361174649999995pt/></p>

where <img src="svgs/be923a44af90ca78e993ae8bed6e4a6e.svg?invert_in_darkmode" align=middle width=92.43919904999998pt height=26.76175259999998pt/>.

## Installation

To install and run the examples
 
- Clone the repository `git clone https://github.com/balintkaszas/ModelSensitivity.git`
- In MATLAB, run the script `addPath.m`

## Implementation 

The software implements the calculation of the scalar field <img src="svgs/820a653cdb14a665cb4a3749bd54dc1b.svg?invert_in_darkmode" align=middle width=65.95915754999999pt height=28.493173500000005pt/>. 
The idealized <img src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg?invert_in_darkmode" align=middle width=9.86687624999999pt height=14.15524440000002pt/>-dimensional model must be provided in the form of a function handle. For example, the damped-driven Duffing oscillator (see below) may be specified as
```MATLAB
function dy = f0(t, x)
    dy(1) = x(2);
    dy(2) = x(1) - x(1)^3 - 0.15*x(2) + 0.3*cos(t);
end
```

The Jacobian of the system must also be given. For example, 
```MATLAB
function dyGrad = grad_f0(t, x)
    dyGrad = [0, 1; 1 - 2*x(1)^2, -0.15];
end
```
The invariants of the Cauchy-Green strain tensor may be computed from 

- finite differencing 
- using the equation of variations.

### Finite differencing 

For finite differencing, an auxiliary grid is used. If the difference between gridpoints is <img src="svgs/ca5b783e4830121f2344b3259c731083.svg?invert_in_darkmode" align=middle width=8.98692299999999pt height=27.342470100000007pt/>, then around each gridpoint, an auxiliary grid is created with size <img src="svgs/29532cac2a859e208b1278a3be5f274a.svg?invert_in_darkmode" align=middle width=72.91950105pt height=27.342470100000007pt/>. 
In this case, the Jacobian is not used for the calculation.

### Equation of variations

If the equation of variations is used, then the matrix-differential equation 

<p align="center"><img src="svgs/a48e7f7bd0fb445298d995cfdb3e3d7b.svg?invert_in_darkmode" align=middle width=144.06387435pt height=20.8872774pt/></p>
is solved. The solution obeying the initial condition
<p align="center"><img src="svgs/f8cacfe7193a583adc457440f4b03390.svg?invert_in_darkmode" align=middle width=74.95433055pt height=16.438356pt/></p>

is the flowmap-gradient. 


## Example

As an example, let us assume that we wish to calculate MS for the damped-driven Duffing equation given by

<p align="center"><img src="svgs/0481a492d99d0bf1214e8fabeb9ff04a.svg?invert_in_darkmode" align=middle width=44.52802695pt height=14.174856299999998pt/></p>
<p align="center"><img src="svgs/d41b284a4837dc7b27c922126e1b7cf6.svg?invert_in_darkmode" align=middle width=183.90170865pt height=17.399144399999997pt/></p>
with <img src="svgs/3dd19a025b3d253c43db1187e62c7953.svg?invert_in_darkmode" align=middle width=59.06955449999999pt height=22.831056599999986pt/> and <img src="svgs/eab9b48fe3d2ac3deb50c9c260a9fa41.svg?invert_in_darkmode" align=middle width=55.25107169999998pt height=22.465723500000017pt/>. 

Also assume that the derivative is available as a function handle in the file `d_duffing.m` and the Jacobian is contained in `d_duffing_grad.m`

We first create a DynSystem object, by specifying the function handle, the dimensions of the phase space and the value of the model uncertainty: <img src="svgs/3b468fc90e1e51b9213b87f119e64ec2.svg?invert_in_darkmode" align=middle width=26.80375499999999pt height=22.465723500000017pt/> and <img src="svgs/e0ba81a24a00ccda69bbff19e04c095a.svg?invert_in_darkmode" align=middle width=26.80375499999999pt height=22.465723500000017pt/>. For simplicity, both deltas are chosen to be 1 here.

```
duffing = DynSystem(@(t,x) d_duffing(t,x), 2, [1,1], @(t,x) d_duffing_grad(t,x));
```

Next, we set up the computational domain by creating a Grid object. We specify a grid of 250 by 250, over the domain <img src="svgs/5451d22ed4f28a2f76eb9081de995e8c.svg?invert_in_darkmode" align=middle width=220.78394085pt height=24.65753399999998pt/>.

```
resolution = [250, 250]; 
domain = [-1.5, 1.5; -1.5, 1.5];
init = Grid(2, [1,2], resolution, domain, 1e-3);
timeSpan = [0, pi];
```

Then we can call the wrapper to calculate MS by providing the necessary arguments. The fourth argument enables parallelization, the fifth specifies the desired accuracy for the integration. The last argument is the method for the calculation of the flow-map gradient, this can be either 'finitedifference' or 'eov'. 

```
ms = modelSensitivity(duffing, init, timeSpan, true, 1e-7, 'finitedifference');
```

The resulting MS field is 

<img src="test/duffing_ms.png" alt="" width="500"/>


## Validation 

An example, in which the MS field is computable analytically is used for validation. The computation is detailed [here](test/validation.ipynb).

## References

[1] B. Kaszás, G. Haller, Universal Upper Estimate for Prediction Errors under Moderate Model Uncertainty, [arXiv:2007.07330 [nlin.CD]](https://arxiv.org/abs/2007.07330) (2020)
