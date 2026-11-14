# Venus flytrap demographic model

`flytrap_ipm.py` fits a compact, size-structured integral projection model (IPM)
to `Datasets/demoVFT.csv`. It requires only Python's standard library.

```bash
python3 flytrap_ipm.py --output flytrap_results.json
python3 graph_flytrap.py flytrap_results.json --output flytrap_growth_vs_fire.svg
```

The JSON report contains the fitted coefficients, scaling means and standard
deviations, observation counts, the site-year burn ages, and projected annual
population multipliers for several deterministic fire-return intervals. A
multiplier below 1 means decline; above 1 means increase.
The graph shows those multipliers by fire-return interval, marks the highest
tested value, and includes a horizontal reference at 1 (population stability).
It is an SVG file that opens in a browser or can be inserted into a document.

## What λ means

**λ (lambda) is the annual population multiplier.** If λ = 1, the projected
number of plants is stable; λ = 1.05 means roughly 5% growth per year; λ = 0.95
means roughly 5% decline per year. It describes the whole population, not the
change in an individual plant's size. In this implementation, each simulated
year has a multiplier `N_next / N_current`. The reported λ for a fire interval
is the geometric mean of those multipliers over the final
`min(--years, 5 × interval)` simulated years:

```text
λ = exp(mean(log(N_next / N_current)))
```

The projection normalizes the plant-size distribution after each year, but
retains each year's multiplier before normalizing. This avoids numerical
underflow without changing λ. The main paper instead reports a transient
growth rate over its historical or future climate period,
`(N_final / N_initial)^(1 / number_of_years)`, so its published λ is not
directly comparable to this fixed-climate, late-cycle estimate.

## Script inputs and outputs

`flytrap_ipm.py` accepts these command-line inputs:

| Input | Default | Meaning |
| --- | --- | --- |
| `--data PATH` | `Datasets/demoVFT.csv` | Plant-level observational CSV. |
| `--burn-events PATH` | None | Optional CSV of actual fire dates or census years to check against the dataset's burn ages. |
| `--intervals N [N ...]` | `1 2 3 5 8 10 12 15 20` | Positive integer years between simulated fires. |
| `--years N` | `100` | Positive integer number of projected years per interval. |
| `--output PATH` | Standard output | Write the JSON report to this path. |

The plant CSV needs `site` and `year`, plus these columns (missing values may
be written as `NA`):

| Column | Meaning and unit | Used for |
| --- | --- | --- |
| `size_prevTS` | Previous census plant size: natural log of leaf count × longest leaf length in mm | Predictor in all fitted rates; size bins in the projection. |
| `size` | Current census log plant size, same measure | Growth response for survivors. |
| `libsur_prevTS` | Survival from previous census, 0 or 1 | Survival response; identifies survivors for other fits. |
| `binrep` | Fruited this year, 0 or 1 | Fruiting response among survivors. |
| `fr` | Fruit count this year | Fruit-number response among fruiting survivors. |
| `burn` | Integer years since most recent fire; 1 is the first census after a burn | Fire predictor; checked for consistency within each site-year. |
| `annual.mean.temp.ibutton` | Annual mean local soil temperature, °C | Climate predictor. |
| `annual.precip.metdata` | Annual precipitation, mm | Climate predictor. |

`ID`, `region`, and the other climate columns are present in this dataset but
are not used by this compact refit. Rows missing a required predictor are
excluded from that vital-rate fit. The fixed projection climate is the mean
across the unique observed site-year temperature and precipitation records,
so sites with more plant rows do not receive extra climate weight.

The JSON report contains:

| Key | Meaning |
| --- | --- |
| `projected_lambda_by_fire_return_interval` | Map from interval in years to projected annual population multiplier; this is what `graph_flytrap.py` plots. |
| `projection_settings`, `lambda_definition` | Simulation length, size-bin count, scenario assumptions, and precise λ calculation. |
| `fixed_climate` | Temperature and precipitation used in every simulated year. |
| `burn_site_years`, `burn_source` | Observed or checked burn ages used for fitting, and their source. |
| `features` | Ordered predictor names used for the coefficient arrays. |
| `scaling_mean_sd` | Mean and standard deviation used to standardize each continuous predictor. |
| `coefficients` | Fitted coefficient arrays for the four vital rates, plus the scalar pooled `growth_variance`. |
| `sample_sizes` | Number of plant records used to fit each vital rate. |
| `recruitment_per_fruit` | Fixed value, 0.01818, taken from Appendix S1. |
| `model`, `warning` | Model identity and its main limitation. |

For direct Python use, `Features(rows)` builds the predictor encoding,
`fit_vital_rates(rows, features)` returns `(model, sample_sizes)`, and
`project(model, features, climate, interval, years=100, mesh_points=100)`
returns one λ as a `float`. `climate` is a two-element sequence containing
annual local temperature in °C and annual precipitation in mm. `model` is the
dictionary returned by `fit_vital_rates`; `interval` is a positive integer
number of years between fires. `project` does not write a file.

## What is fitted

For each plant transition, the model uses previous log plant size, years since
fire, a squared years-since-fire term, annual local temperature, annual
precipitation, and a site indicator. Continuous inputs are centered and scaled
using the fitting data. A small ridge penalty stabilizes the coefficients.

It fits:

1. Annual survival with a logistic regression.
2. Next-year log size among survivors with a linear regression and a pooled
   residual variance.
3. Fruiting probability among survivors with a logistic regression.
4. The log of fruits per unit original-scale plant size among fruiting survivors
   with a linear regression.

The fertility part of the IPM multiplies survival, fruiting probability,
predicted fruit number, and the paper's constant **0.01818 seedlings per
fruit**. Seedlings enter with the paper's log-size mean 3.38 and variance 0.33.
The size domain uses 100 bins. Survival-growth and fertility are integrated
over these bins, and the final five burn cycles give the reported annual
population multiplier. Each scenario uses the same mean observed climate and
repeats a fire on a fixed interval. For example, interval 5 uses years since
fire 1, 2, 3, 4, 5, then 1 again.

## Integrating burn history

`demoVFT.csv` already contains `burn`, meaning **years since fire** for each
site-year. This is the most direct burn input for the observed flytrap plants.
The program verifies that every plant at a given site and year has the same
burn age. To check it against a separate event record, supply a CSV with
`site,year` or `site,burn_date` columns:

Illustrative event record (the exact date is only an example):

```csv
site,burn_date
coastal.1,2018-04-10
```

Dates must be ISO `YYYY-MM-DD`. A fire on or before the June 15 census is
assigned to that year's transition; one after June 15 is assigned to the next
census year. The first census following the fire has `burn=1`. Only site names
matching `demoVFT.csv` are accepted. The program checks any overlapping
site-years against the existing `burn` values and reports disagreements.

```bash
python3 flytrap_ipm.py --burn-events my_flytrap_burn_events.csv \
  --intervals 1 3 5 10 12 20 --output flytrap_results.json
```

`Datasets/burnData.csv` cannot be directly joined to the flytrap observations:
its sites are `S1`–`S7` from a separate study, while `demoVFT.csv` uses six
named monitoring sites. It also records days since fire at one survey, rather
than a fire history through the flytrap censuses. A defensible join would
require a documented site mapping and matching dates; the program deliberately
rejects unknown sites.

## Scope and limits

This is a **refit**, not a numerical reproduction of the article's published
forecasts. The paper selected different climate variables and interactions for
each vital rate using AICc, modeled growth variance as its own function, used
multiple climate scenarios and climate models, and propagated coefficient
uncertainty. Its Table S3 coefficients use scaled predictors, but the scaling
constants and all projection inputs are not supplied in the two PDFs. The
compact fit here uses a common predictor set and a constant growth variance
so that it remains estimable from the 18 observed site-year climate/burn
combinations. It should be used to explore mechanisms and compare burn
scenarios under the observed mean climate, not to reproduce the published
10-year and 12-year optima or make a management forecast.

Sources: `Documents/SensitivityCoefficientsDataSource.pdf`, pages 5–8, and
Appendix S1 Table S3 in the supplied attachment.
