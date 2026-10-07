# cS: Contaminated Survival Models

`cS` fits the **contaminated survival (cS) model**, a two-component mixture for
time-to-event data in which both components share the same reference
distribution and differ only in speed. A proportion `alpha` of the
observations belongs to a **delayed** component whose hazard is multiplied by
`eta` (between 0 and 1), so these observations take longer to experience the
event.

The model was developed to identify delayed birth registrations from two
administrative dates: the date of birth and the date of registration. It
tells you

1. **how many** registrations follow a slower regime (`alpha`),
2. **how much slower** that regime is (`eta`), and
3. **which** records are likely to belong to it (a posterior probability for
   every birth).

Concomitant variables can be added to see **which groups** are more likely to
be delayed.

Two reference distributions are implemented: the exponential (giving the
contaminated exponential, cE) and the Weibull (giving the contaminated
Weibull, cW).

## Installation

```r
# install.packages("devtools")
library(devtools)
install_github("arnootto/cS")
library(cS)
```

## The model

The survival function of the cS model is

$$
S_C(t) = (1-\alpha)\,S_R(t) + \alpha\,S_R(t)^{\eta},
$$

where $S_R(t)$ is the reference survival function, $\alpha \in (0,1)$ is the
proportion of observations in the delayed component, and $\eta \in (0,1)$ is
the hazard multiplier of that component. For the contaminated Weibull, every
quantile of the delayed component (the median, for example) is
$\eta^{-1/k}$ times the corresponding quantile of the reference component.
The functions return this as `time.ratio`.

With concomitant variables $\mathbf{x}_i$, the proportion becomes

$$
\alpha_i = \frac{\exp(\mathbf{x}_i'\boldsymbol\gamma)}{1+\exp(\mathbf{x}_i'\boldsymbol\gamma)}.
$$

All models are fitted by maximum likelihood with the EM algorithm.

## Main functions

| Function | What it does |
|---|---|
| `fit.e()`, `fit.w()` | Fit the exponential and Weibull reference models |
| `ml.cexp()`, `ml.cweib()` | Fit the cE and cW models |
| `ml.cexp.cv()`, `ml.cweib.cv()` | Fit the cE and cW models with concomitant variables |
| `lrt()` | Likelihood ratio test between two nested concomitant models |
| `births` | January 2024 Ecuadorian birth registration data |

## Example: Ecuadorian birth registrations

The `births` data contain the 17,144 live births that occurred in Ecuador in
January 2024 (source: INEC), with the date of birth, the date of registration,
the lag in days, and four indicators (`home`, `rural`, `indigenous`,
`out_province`). The code below reproduces the results in the paper.

### Step 1: Prepare the survival time

The registry records whole days, and some births are registered on the day of
birth (`lag = 0`). We use the midpoint of the day, `lag + 0.5`.

```r
library(cS)
data(births)
births$t <- births$lag + 0.5
```

### Step 2: Is there a delayed regime?

Fit the two reference models and their contaminated versions, and compare
them with BIC (lower is better).

```r
fe  <- fit.e(births$t)
fw  <- fit.w(births$t)
fce <- ml.cexp(births$t)
fcw <- ml.cweib(births$t)

c(E = fe$BIC, W = fw$BIC, cE = fce$BIC, cW = fcw$BIC)
#>        E        W       cE       cW
#> 152687.7 150858.6 149143.3 148943.4
```

Both contaminated models fit far better than their reference models, and the
cW fits best.

### Step 3: How large and how slow is the delayed regime?

```r
fcw$results
#>       estimates    se
#> sigma    14.661 0.227
#> k         1.134 0.009
#> alpha     0.270 0.008
#> eta       0.140 0.004

fcw$time.ratio
#> 5.67
```

About 27% of births belong to the delayed regime, and every quantile of the
delayed component is about 5.7 times the corresponding quantile of the
routine component. The routine median is
$\hat\sigma(\ln 2)^{1/\hat k} \approx 10.6$ days, and the delayed median is
about $5.67 \times 10.6 \approx 60$ days.

### Step 4: Which births are delayed?

`fcw$v` holds the posterior probability that each birth belongs to the
delayed component. Classify a birth as delayed if this probability exceeds
0.5.

```r
births$prob.delayed <- fcw$v
births$delayed <- ifelse(fcw$v > 0.5, 1, 0)
table(births$delayed)
#>     0     1
#> 13811  3333
```

### Step 5: Which groups are more likely to be delayed?

Let the proportion `alpha` depend on the four indicators, and test whether
they improve the fit.

```r
f0 <- ml.cweib.cv(t ~ 1, ~ 1, data = births)
f1 <- ml.cweib.cv(t ~ 1, ~ home + rural + indigenous + out_province,
                  data = births)

lrt(f0, f1)
#>    loglike0  loglike1     stat df             p
#> 1 -74452.22 -74198.71 507.0136  4 2.037444e-108

f1$results
exp(f1$gamma[-1])     # odds ratios for belonging to the delayed component
#>         home        rural   indigenous out_province
#>         6.27         1.35         2.13         8.15
```

A birth registered outside its province of birth has about eight times the
odds of belonging to the delayed component, and a home birth about six times.
`f1$alpha` gives the fitted probability of delayed registration for every
birth.

## Using your own data

All the model needs is two dates per birth. If your data frame `mydata` has
columns `birth_date` and `registration_date`:

```r
mydata$birth_date        <- as.Date(mydata$birth_date)
mydata$registration_date <- as.Date(mydata$registration_date)

mydata$lag <- as.numeric(mydata$registration_date - mydata$birth_date)
mydata <- mydata[mydata$lag >= 0, ]
mydata$t <- mydata$lag + 0.5

fit <- ml.cweib(mydata$t)
fit$results
fit$time.ratio
```

Then follow Steps 2–5 above, replacing `births` with `mydata`.

## Return values

`ml.cexp()` and `ml.cweib()` return a list with, among others:

| Component | Description |
|---|---|
| `results` | Estimates and standard errors |
| `alpha`, `eta` | Proportion delayed and hazard multiplier |
| `lambda` (cE) or `sigma`, `k` (cW) | Parameters of the reference distribution |
| `time.ratio` | How many times slower the delayed component is |
| `v` | Posterior probability that each observation is delayed |
| `loglike`, `AIC`, `BIC` | Log-likelihood and information criteria |
| `AIC.ref`, `BIC.ref` | Information criteria of the reference model |

`ml.cexp.cv()` and `ml.cweib.cv()` also return `gamma` (the coefficients for
`alpha`, on the logit scale), `alpha` (the fitted proportion for every
observation) and `npar` (the number of parameters, used by `lrt()`).

## Other functions

| Function | Description |
|---|---|
| `dcexp()`, `scexp()`, `hcexp()` | Density, survival and hazard of the cE distribution |
| `dcweib()`, `scweib()`, `hcweib()` | Density, survival and hazard of the cW distribution |
| `h0weib()` | Hazard of the reference Weibull distribution |
| `pi.cexp()` | Proportion of delayed observations still at risk at time *t* (cE) |
| `rcexp()`, `rcweib()` | Simulate from the cE and cW distributions |
| `rcexp.cv()`, `rcweib.cv()` | Simulate from the models with concomitant variables |

For more details, type `?function_name` in R (e.g., `?ml.cweib`).

## Citation

Otto, A. (2026). Who registers late? Identifying delayed birth registrations
from two administrative dates.
