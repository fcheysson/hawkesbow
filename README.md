
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hawkesbow

<!-- badges: start -->

<!-- badges: end -->

In this package, we implement a spectral approach to the parametric
estimation of Hawkes processes from binned observations.

## Installation

You can install the released version of hawkesbow from
[GitHub](https://github.com/fcheysson/hawkesbow) with:

``` r
devtools::install_github("fcheysson/hawkesbow")
```

## The Hawkes process

Hawkes processes form a family of models for point processes for which
the occurrence of an event temporarily increases the probability of
future events occurring (Hawkes 1971). Formally, a Hawkes process
![N](https://latex.codecogs.com/png.latex?N "N") on ![\\mathbb
R](https://latex.codecogs.com/png.latex?%5Cmathbb%20R "\\mathbb R") can
be defined from its conditional intensity
![\\lambda(\\cdot)](https://latex.codecogs.com/png.latex?%5Clambda%28%5Ccdot%29
"\\lambda(\\cdot)") :  where the variables
![T\_i](https://latex.codecogs.com/png.latex?T_i "T_i") denote the
arrival times of the process, the *immigration intensity*
![\\eta](https://latex.codecogs.com/png.latex?%5Ceta "\\eta") is a
positive constant, and the *reproduction function* ![h: \\mathbb
R\_{\\ge 0} \\rightarrow \\mathbb
R\_{\\ge 0}](https://latex.codecogs.com/png.latex?h%3A%20%5Cmathbb%20R_%7B%5Cge%200%7D%20%5Crightarrow%20%5Cmathbb%20R_%7B%5Cge%200%7D
"h: \\mathbb R_{\\ge 0} \\rightarrow \\mathbb R_{\\ge 0}") is a
measurable function. Alternatively, the Hawkes process can be
constructed as a poissonian cluster process (Hawkes and Oakes 1974). The
process consists of a flow of *immigrants*, the cluster centres,
arriving according to a homogeneous Poisson process of intensity
![\\eta](https://latex.codecogs.com/png.latex?%5Ceta "\\eta"). Then, an
immigrant arriving at time
![T\_i](https://latex.codecogs.com/png.latex?T_i "T_i") generates
*children* according to an inhomogeneous Poisson process of intensity
![h(\\cdot -
T\_i)](https://latex.codecogs.com/png.latex?h%28%5Ccdot%20-%20T_i%29
"h(\\cdot - T_i)"). These in turn independently generate children
according to the same process, and so on *ad infinitum*. The processes
consisting of an immigrant and all its descendants are therefore
branching processes, and are independent of each other. Finally, the
process resulting from the superposition of all these branching
processes is a Hawkes process of conditional intensity
![\\lambda(\\cdot)](https://latex.codecogs.com/png.latex?%5Clambda%28%5Ccdot%29
"\\lambda(\\cdot)") (see the following figure of a realisation of an
exponential Hawkes process with ![\\\\eta
= 1](https://latex.codecogs.com/png.latex?%5C%5Ceta%20%3D%201
"\\\\eta = 1"), ![h(t) =
e^{-2t}](https://latex.codecogs.com/png.latex?h%28t%29%20%3D%20e%5E%7B-2t%7D
"h(t) = e^{-2t}")).

![Realisation of an exponential Hawkes process](man/figures/hawkes.png)

## Usage

## Reproduction kernels

Fourier transform:   
![\\widetilde {h^\\ast}(\\omega) = \\int\_\\mathbb R \\exp(-i\\omega t)
h(t) \\mathrm
dt](https://latex.codecogs.com/png.latex?%5Cwidetilde%20%7Bh%5E%5Cast%7D%28%5Comega%29%20%3D%20%5Cint_%5Cmathbb%20R%20%5Cexp%28-i%5Comega%20t%29%20h%28t%29%20%5Cmathrm%20dt
"\\widetilde {h^\\ast}(\\omega) = \\int_\\mathbb R \\exp(-i\\omega t) h(t) \\mathrm dt")  

### The exponential kernel

  
![h^\\ast(t) = \\beta e^{-\\beta t} 1\_{\\{t
\> 0\\}}](https://latex.codecogs.com/png.latex?h%5E%5Cast%28t%29%20%3D%20%5Cbeta%20e%5E%7B-%5Cbeta%20t%7D%201_%7B%5C%7Bt%20%3E%200%5C%7D%7D
"h^\\ast(t) = \\beta e^{-\\beta t} 1_{\\{t \> 0\\}}")  

Fourier transform   
![\\widetilde {h^\\ast}(\\omega) = \\beta \\frac{1}{\\beta + i\\omega} =
\\beta \\frac{\\beta - i \\omega}{\\beta^2 +
\\omega^2}](https://latex.codecogs.com/png.latex?%5Cwidetilde%20%7Bh%5E%5Cast%7D%28%5Comega%29%20%3D%20%5Cbeta%20%5Cfrac%7B1%7D%7B%5Cbeta%20%2B%20i%5Comega%7D%20%3D%20%5Cbeta%20%5Cfrac%7B%5Cbeta%20-%20i%20%5Comega%7D%7B%5Cbeta%5E2%20%2B%20%5Comega%5E2%7D
"\\widetilde {h^\\ast}(\\omega) = \\beta \\frac{1}{\\beta + i\\omega} = \\beta \\frac{\\beta - i \\omega}{\\beta^2 + \\omega^2}")  

The exponential kernel can be specified with the string “Exponential”

Both mle and whittle methods are available for exponential kernels.

### The Gaussian kernel

  
![h^\\ast(t) = \\frac{1}{\\sigma
\\sqrt(2\\pi)}\\exp\\left(-\\frac{(t-\\nu)^2}{2\\sigma^2}\\right)](https://latex.codecogs.com/png.latex?h%5E%5Cast%28t%29%20%3D%20%5Cfrac%7B1%7D%7B%5Csigma%20%5Csqrt%282%5Cpi%29%7D%5Cexp%5Cleft%28-%5Cfrac%7B%28t-%5Cnu%29%5E2%7D%7B2%5Csigma%5E2%7D%5Cright%29
"h^\\ast(t) = \\frac{1}{\\sigma \\sqrt(2\\pi)}\\exp\\left(-\\frac{(t-\\nu)^2}{2\\sigma^2}\\right)")  

Fourier transform   
![\\widetilde {h^\\ast}(\\omega) =
\\exp\\left(-\\frac{\\sigma^2\\omega^2}{2}-i\\nu\\omega\\right)](https://latex.codecogs.com/png.latex?%5Cwidetilde%20%7Bh%5E%5Cast%7D%28%5Comega%29%20%3D%20%5Cexp%5Cleft%28-%5Cfrac%7B%5Csigma%5E2%5Comega%5E2%7D%7B2%7D-i%5Cnu%5Comega%5Cright%29
"\\widetilde {h^\\ast}(\\omega) = \\exp\\left(-\\frac{\\sigma^2\\omega^2}{2}-i\\nu\\omega\\right)")  

The gaussian kernel can be specified with the string “Gaussian”

Only whittle method is available for gaussian kernels.

### The Pareto kernels

  
![h^\\ast(t) = \\theta a^\\theta t^{-\\theta - 1} 1\_{\\{t \>
a\\}}](https://latex.codecogs.com/png.latex?h%5E%5Cast%28t%29%20%3D%20%5Ctheta%20a%5E%5Ctheta%20t%5E%7B-%5Ctheta%20-%201%7D%201_%7B%5C%7Bt%20%3E%20a%5C%7D%7D
"h^\\ast(t) = \\theta a^\\theta t^{-\\theta - 1} 1_{\\{t \> a\\}}")  

Fourier transform   
![\\widetilde {h^\\ast}(\\omega) = \\theta a^\\theta E\_{\\theta + 1}
(i\\omega)
](https://latex.codecogs.com/png.latex?%5Cwidetilde%20%7Bh%5E%5Cast%7D%28%5Comega%29%20%3D%20%5Ctheta%20a%5E%5Ctheta%20E_%7B%5Ctheta%20%2B%201%7D%20%28i%5Comega%29%20
"\\widetilde {h^\\ast}(\\omega) = \\theta a^\\theta E_{\\theta + 1} (i\\omega) ")  
where   
![E\_{\\theta+1}(i\\omega) = \\int\_a^{+\\infty} \\frac{\\exp(-i\\omega
t)}{t^{\\theta+1}} \\mathrm
dt](https://latex.codecogs.com/png.latex?E_%7B%5Ctheta%2B1%7D%28i%5Comega%29%20%3D%20%5Cint_a%5E%7B%2B%5Cinfty%7D%20%5Cfrac%7B%5Cexp%28-i%5Comega%20t%29%7D%7Bt%5E%7B%5Ctheta%2B1%7D%7D%20%5Cmathrm%20dt
"E_{\\theta+1}(i\\omega) = \\int_a^{+\\infty} \\frac{\\exp(-i\\omega t)}{t^{\\theta+1}} \\mathrm dt")  

Only Pareto kernels with ![\\theta
= 1](https://latex.codecogs.com/png.latex?%5Ctheta%20%3D%201
"\\theta = 1"), ![2](https://latex.codecogs.com/png.latex?2 "2"), and
![3](https://latex.codecogs.com/png.latex?3 "3") have been implemented
and can specified with the strngs “Pareto1”, “Pareto2” and “Pareto3”
respectively.

Only whittle method is available for Pareto kernels.

### Power law kernel

  
![h^\\ast(t) = \\theta a^\\theta (t+a)^{-\\theta-1} 1\_{\\{\\theta \> 0
\\}}](https://latex.codecogs.com/png.latex?h%5E%5Cast%28t%29%20%3D%20%5Ctheta%20a%5E%5Ctheta%20%28t%2Ba%29%5E%7B-%5Ctheta-1%7D%201_%7B%5C%7B%5Ctheta%20%3E%200%20%5C%7D%7D
"h^\\ast(t) = \\theta a^\\theta (t+a)^{-\\theta-1} 1_{\\{\\theta \> 0 \\}}")  

Fourier transform   
![\\widetilde{h\_\\theta^\\ast}(\\omega) = \\theta a^\\theta
\\exp(i\\omega a)E\_{\\theta + 1} (i\\omega)
](https://latex.codecogs.com/png.latex?%5Cwidetilde%7Bh_%5Ctheta%5E%5Cast%7D%28%5Comega%29%20%3D%20%5Ctheta%20a%5E%5Ctheta%20%5Cexp%28i%5Comega%20a%29E_%7B%5Ctheta%20%2B%201%7D%20%28i%5Comega%29%20
"\\widetilde{h_\\theta^\\ast}(\\omega) = \\theta a^\\theta \\exp(i\\omega a)E_{\\theta + 1} (i\\omega) ")  

The power law kernel can be specified with the string “PowerLaw”

Both mle and whittle methods are available for power law kernels.

<div id="refs" class="references">

<div id="ref-Hawkes1971a">

Hawkes, Alan G. 1971. “Spectra of Some Self-Exciting and Mutually
Exciting Point Processes.” *Biometrika* 58 (1): 83–90.
<https://doi.org/10.2307/2334319>.

</div>

<div id="ref-Hawkes1974">

Hawkes, Alan G., and David Oakes. 1974. “A cluster process
representation of a self-exciting process.” *J. Appl. Probab.* 11 (03):
493–503. <https://doi.org/10.2307/3212693>.

</div>

</div>
