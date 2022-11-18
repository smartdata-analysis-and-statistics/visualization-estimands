
# Visualization as a diagnostic tool for estimands: A proof-of-concept and case study with pairwise matching of multiple treatments

Propensity score matching is among the most popular confounding
adjustment methods for comparative effectiveness research in
non-randomized settings. Although propensity score methods were
originally introduced in the context of comparisons of two treatments,
they are increasingly being used for multiple comparisons of more than
two treatments. When pairwise matching is applied to estimate treatment
effects for multiple pairs of treatments, the underlying target
population changes across the comparisons. This is attributable to the
differences in covariate distributions across treatment groups.
Consequently, this results in different treatment effect interpretations
for the different underlying populations across the multiple pairwise
comparisons. The interpretation of treatment effect estimates in
relation to the matching (i.e., the target estimand) is rarely clarified
and consequently might lead to erroneous conclusions about real-world
effectiveness of different treatments. Based on empirical research, we
illustrate that multiple pairwise matching for the investigation of
comparative effectiveness of more than two treatments can lead to
targeting different estimands. We use a proof-of-concept simulation
study and a case study in multiple sclerosis. We propose visualization
tools to illustrate the problem and clarify the connection between
estimand, target population and pairwise matching to avoid
misinterpretations and treatment decision-making errors in clinical
practice.

All R code from this project is located in this repository. To run the
simulation, use `runSimulation.r`. For the visualization functions, use
`/R/visualisation.r`.

The code contains two scenarios:

- Scenario 1: Absence of treatment effect heterogeneity

- Scenario 2: Presence of treatment effect heterogeneity (used in the
  appendix of the manuscript)

Depending on the selected scenario, the values for $\beta_5$ and
$\beta_6$ were modified.

------------------------------------------------------------------------

### Data-generating mechanism

We generated data for a super-population of $N=100,000$ individuals
treated with treatment A, B, or C. We generated two continuous variables
$x_1$ and $x_2$ for each patient as baseline covariates, labeled as age
(standardized) and disease severity (z-score). We generated a potential
outcome, labeled as a test score for which higher values are better,
under each of the three treatments. The test scores were generated from
a linear model with age, disease severity, and the treatment indicators
as main effects and an interaction term between the treatment indicators
and age to introduce treatment effect heterogeneity. The outcome $y$ is
continuous, and generated as follows:

$y_i \sim N(\mu_i, \sigma^2)$

with

$\mu_i = \beta_0 + \beta_1 x_{1i} + \beta_2 x_{2i} + \beta_3 I_{Bi} + \beta_4 I_{Ci} + (\beta_5 I_{Bi} + \beta_6 I_{Ci}) x_{1i}$

Using the following regression coefficients:

- Baseline outcome risk ($\beta_0$): **0**
- Confounder effect of $x_1$ ($\beta_1$): **1**
- Confounder effect of $x_2$ ($\beta_2$): **0.25**
- Treatment effect of B versus A ($\beta_3$): **-0.5**
- Treatment effect of C versus A ($\beta_4$): **-0.25**
- Standard deviation of the residual error ($\sigma$): **1**

Two scenarios were evaluated in the simulations by modifying the values
for $\beta_5$ and $\beta_6$.

**Scenario 1**: Absence of treatment effect heterogeneity (HTE)

- Effect modification between B and x1 ($\beta_5$): **0**
- Effect modification between C and x1 ($\beta_6$): **0**

**Scenario 2**: Presence of HTE

- Effect modification between B and x1 ($\beta_5$): **0.25**
- Effect modification between C and x1 ($\beta_6$): **0.125**

We assigned the treatment received as a function of the two baseline
covariates $x_1$ and $x_2$. Thus, the observable dataset consisted of
the two baseline covariates, the treatment received and one outcome –
the one generated under the treatment that the patient received.

### Implementation of pairwise matching

We applied 1:1 PS nearest-neighbor matching with a caliper and with
replacement for each pairwise treatment comparison. We estimated the
propensity score with a logistic regression model as a function of the
two baseline covariates separately for each treatment comparison to
mimic how pairwise matching would be repetitively applied in the context
of three treatments being compared. We used a caliper of 0.05 standard
deviation of the logit of the estimated PS. We assessed covariates
balance pre- and post-matching with absolute standardized mean
differences.

For each treatment comparison, we generated created two matched samples,
separately targeting eachither of the ATTs. For example, for the
comparison of treatment A vs B, we first match individuals treated with
B to individuals treated with A, thus targeting the ATT-A. A first
matched sample results from this matching. Second, starting from the
original sample again, we now match individuals treated with A to
individuals treated with B, targeting the ATT-B, which results in a
second matched sample. For each treatment comparison and in each matched
sample, we estimate the treatment effect (e.g., difference in means at a
certain timepoint) using a linear regression model. We note that, in an
application of pairwise matching to real data, researchers would
typically generate only one matched sample, targeting one estimand or
the other.

### Evaluation and visualization

For the evaluation and visualization of covariate overlap, we compared
the target populations before and after matching for all pairwise
treatment comparisons with the bivariate ellipses. We also compared the
estimated treatment effects in the matched sample.
