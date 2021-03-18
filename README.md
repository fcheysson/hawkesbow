
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hawkesbow

<!-- badges: start -->

<!-- badges: end -->

In this package, we implement a spectral approach to the parametric
estimation of Hawkes processes from binned observations.

## Installation

You can install the released version of hawkesbow from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("hawkesbow")
```

## The Hawkes process

Hawkes processes form a family of models for point processes for which
the occurrence of an event temporarily increases the probability of
future events occurring (Hawkes 1971). Formally, a Hawkes process \(N\)
on \(\mathbb R\) can be defined from its conditional intensity
\(\lambda(\cdot)\) :  where the variables \(T_i\) denote the arrival
times of the process, the *immigration intensity* \(\eta\) is a positive
constant, and the *reproduction function*
\(h: \mathbb R_{\ge 0} \rightarrow \mathbb R_{\ge 0}\) is a measurable
function. Alternatively, the Hawkes process can be constructed as a
poissonian cluster process (Hawkes and Oakes 1974). The process consists
of a flow of *immigrants*, the cluster centres, arriving according to a
homogeneous Poisson process of intensity \(\eta\). Then, an immigrant
arriving at time \(T_i\) generates *children* according to an
inhomogeneous Poisson process of intensity \(h(\cdot - T_i)\). These in
turn independently generate children according to the same process, and
so on *ad infinitum*. The processes consisting of an immigrant and all
its descendants are therefore branching processes, and are independent
of each other. Finally, the process resulting from the superposition of
all these branching processes is a Hawkes process of conditional
intensity \((\lambda(\cdot)\) (Figure ).

``` r
knitr::include_graphics("man/figures/hawkes.eps")
```

<div class="figure" style="text-align: center">

<img src="man/figures/hawkes.eps" alt="Realization of an exponential Hawkes process, with $\eta = 1$, $h(t) = e^{-2t}$. The crosses represent the arrival times of the process. (Top) Representation of the conditional intensity: each event increases the probability of occurrence of future events, according to the reproduction function $h$. (Bottom) Representation of branching: each immigrant (black squares, of generation 0) can generate children (red dots, of generation 1), which can in turn generate children, and so on." width="80%" />

<p class="caption">

Realization of an exponential Hawkes process, with \(\eta = 1\),
\(h(t) = e^{-2t}\). The crosses represent the arrival times of the
process. (Top) Representation of the conditional intensity: each event
increases the probability of occurrence of future events, according to
the reproduction function \(h\). (Bottom) Representation of branching:
each immigrant (black squares, of generation 0) can generate children
(red dots, of generation 1), which can in turn generate children, and so
on.

</p>

</div>

## Reproduction kernels

Fourier transform:
\[\widetilde {h^\ast}(\omega) = \int_\mathbb R \exp(-i\omega t) h(t) \mathrm dt\]

### The exponential kernel

\[h^\ast(t) = \beta e^{-\beta t} 1_{\{t > 0\}}\]

Fourier transform
\[\widetilde {h^\ast}(\omega) = \beta \frac{1}{\beta + i\omega} = \beta \frac{\beta - i \omega}{\beta^2 + \omega^2}\]

The exponential kernel can be specified with the string “Exponential”

Both mle and whittle methods are available for exponential kernels.

### The Gaussian kernel

\[h^\ast(t) = \frac{1}{\sigma \sqrt(2\pi)}\exp\left(-\frac{(t-\nu)^2}{2\sigma^2}\right)\]

Fourier transform
\[\widetilde {h^\ast}(\omega) = \exp\left(-\frac{\sigma^2\omega^2}{2}-i\nu\omega\right)\]

The gaussian kernel can be specified with the string “Gaussian”

Only whittle method is available for gaussian kernels.

### The Pareto kernels

\[h_\theta^\ast(t) = \theta a^\theta t^{-\theta - 1} 1_{\{t > a\}}\]

Fourier transform
\[\widetilde{h_\theta^\ast}(\omega) = \theta a^\theta E_{\theta + 1} (i\omega) \]
where
\[E_{\theta+1}(i\omega) = \int_a^{+\infty} \frac{\exp(-i\omega t)}{t^{\theta+1}} \mathrm dt\]

Only Pareto kernels with \(\theta = 1\), \(2\), and \(3\) have been
implemented and can specified with the strngs “Pareto1”, “Pareto2” and
“Pareto3” respectively.

Only whittle method is available for Pareto kernels.

### Power law kernel

\[h^\ast(t) = \theta a^\theta (t+a)^{-\theta-1} 1_{\{\theta > 0 \}}\]

Fourier transform
\[\widetilde{h_\theta^\ast}(\omega) = \theta a^\theta \exp(i\omega a)E_{\theta + 1} (i\omega) \]

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
