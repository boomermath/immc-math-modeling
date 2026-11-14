"""Reproducible, dependency-free approximation of the Louthan et al. flytrap IPM.

Fits vital rates to Datasets/demoVFT.csv, then projects a size-structured
population under a deterministic fire-return interval. This is a refit with a
compact predictor set, not a reproduction of the paper's AICc-selected fits.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "Datasets" / "demoVFT.csv"
RECRUITMENT_PER_FRUIT = 0.01818  # Appendix S1, Table S3
SEEDLING_MEAN = 3.38             # main paper, log-size scale
SEEDLING_VARIANCE = 0.33
CLIMATE = ("annual.mean.temp.ibutton", "annual.precip.metdata")


def number(value):
    if value is None or value in ("", "NA"):
        return None
    return float(value)


def read_rows(path):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def burn_history(rows):
    """Validate that years-since-fire is constant within each site-year."""
    by_group = defaultdict(set)
    for row in rows:
        value = number(row["burn"])
        if value is not None:
            by_group[row["site"], int(row["year"])].add(int(value))
    bad = {key: values for key, values in by_group.items() if len(values) != 1}
    if bad:
        raise ValueError(f"Conflicting burn ages within site-year: {bad}")
    return {key: next(iter(values)) for key, values in by_group.items()}


def read_burn_events(path, sites):
    """Read site,year or site,burn_date; calendar year is the census burn year."""
    from datetime import date

    events = defaultdict(set)
    for row in read_rows(path):
        site = row.get("site", "").strip()
        if site not in sites:
            raise ValueError(f"Unknown burn-event site {site!r}; expected {sorted(sites)}")
        if row.get("year", "").strip():
            year = int(row["year"])
        elif row.get("burn_date", "").strip():
            day = date.fromisoformat(row["burn_date"].strip())
            # The main study uses a June 15 census. A burn after it affects
            # the next annual transition, so first appears as burn=1 then.
            year = day.year if (day.month, day.day) <= (6, 15) else day.year + 1
        else:
            raise ValueError("Burn events need a year or ISO burn_date column")
        events[site].add(year)
    return events


def integrate_burn_events(rows, path):
    sites = {row["site"] for row in rows}
    observed = burn_history(rows)
    events = read_burn_events(path, sites)
    result = {}
    for (site, year), existing in observed.items():
        prior = [event for event in events[site] if event <= year]
        if prior:
            inferred = year - max(prior) + 1
            if inferred != existing:
                raise ValueError(
                    f"Burn-history conflict at {site}, {year}: dataset={existing}, "
                    f"events imply {inferred}. Check the June 15 census boundary."
                )
            result[site, year] = inferred
        else:
            result[site, year] = existing
    return result


def solve(matrix, vector):
    """Gaussian elimination with partial pivoting for the small fit matrices."""
    a = [row[:] + [value] for row, value in zip(matrix, vector)]
    n = len(a)
    for col in range(n):
        pivot = max(range(col, n), key=lambda i: abs(a[i][col]))
        if abs(a[pivot][col]) < 1e-12:
            raise ValueError("Model matrix is singular")
        a[col], a[pivot] = a[pivot], a[col]
        scale = a[col][col]
        for j in range(col, n + 1):
            a[col][j] /= scale
        for i in range(n):
            if i == col:
                continue
            factor = a[i][col]
            for j in range(col, n + 1):
                a[i][j] -= factor * a[col][j]
    return [a[i][n] for i in range(n)]


def fit(x, y, logistic=False, penalty=2.0):
    p = len(x[0])
    beta = [0.0] * p
    if logistic:
        avg = min(0.999, max(0.001, sum(y) / len(y)))
        beta[0] = math.log(avg / (1 - avg))
    for _ in range(60 if logistic else 1):
        hessian = [[0.0] * p for _ in range(p)]
        gradient = [0.0] * p
        for row, target in zip(x, y):
            eta = sum(a * b for a, b in zip(row, beta))
            if logistic:
                prob = sigmoid(eta)
                weight = max(1e-8, prob * (1 - prob))
                error = target - prob
            else:
                weight = 1.0
                error = target - eta
            for i in range(p):
                gradient[i] += row[i] * error
                for j in range(p):
                    hessian[i][j] += row[i] * row[j] * weight
        for i in range(1, p):
            hessian[i][i] += penalty
            gradient[i] -= penalty * beta[i]
        delta = solve(hessian, gradient)
        beta = [value + step for value, step in zip(beta, delta)]
        if max(abs(step) for step in delta) < 1e-7:
            break
    return beta


def sigmoid(value):
    if value >= 0:
        return 1 / (1 + math.exp(-value))
    exp_value = math.exp(value)
    return exp_value / (1 + exp_value)


class Features:
    def __init__(self, rows):
        self.sites = sorted({row["site"] for row in rows})
        values = [number(row["size_prevTS"]) for row in rows]
        values = [value for value in values if value is not None]
        self.size_bounds = (min(values), max(values))
        self.scaling = {}
        for key in ("size_prevTS", "burn", *CLIMATE):
            data = [number(row[key]) for row in rows]
            data = [value for value in data if value is not None]
            mean = sum(data) / len(data)
            sd = math.sqrt(sum((value - mean) ** 2 for value in data) / len(data))
            self.scaling[key] = (mean, sd or 1.0)
        self.names = ["intercept", "size", "burn", "burn_squared", "annual_temp", "annual_precip"]
        self.names += ["site:" + site for site in self.sites[1:]]

    def vector(self, row):
        z = number(row["size_prevTS"])
        fire = number(row["burn"])
        temp = number(row[CLIMATE[0]])
        precip = number(row[CLIMATE[1]])
        if None in (z, fire, temp, precip):
            return None
        scale = lambda key, value: (value - self.scaling[key][0]) / self.scaling[key][1]
        f = scale("burn", fire)
        return [1.0, scale("size_prevTS", z), f, f * f,
                scale(CLIMATE[0], temp), scale(CLIMATE[1], precip)] + [
                    float(row["site"] == site) for site in self.sites[1:]
                ]


def fit_vital_rates(rows, features):
    """Fit four vital rates from CSV rows.

    Returns (model, counts). ``model`` contains coefficient arrays in
    ``features.names`` order for survival, growth, fruiting and fruit_count,
    plus a scalar growth_variance. ``counts`` gives each fit's sample size.
    """
    samples = {name: ([], []) for name in ("survival", "growth", "fruiting", "fruit_count")}
    for row in rows:
        x = features.vector(row)
        if x is None:
            continue
        survival = number(row["libsur_prevTS"])
        size = number(row["size"])
        fruiting = number(row["binrep"])
        fruits = number(row["fr"])
        if survival is not None:
            samples["survival"][0].append(x)
            samples["survival"][1].append(survival)
        if survival == 1 and size is not None:
            samples["growth"][0].append(x)
            samples["growth"][1].append(size)
        if fruiting is not None and survival == 1:
            samples["fruiting"][0].append(x)
            samples["fruiting"][1].append(fruiting)
        if fruits is not None and fruits > 0 and survival == 1:
            samples["fruit_count"][0].append(x)
            # Paper's response: log(fruits / size), with size on original scale.
            samples["fruit_count"][1].append(math.log(fruits) - number(row["size_prevTS"]))
    model = {}
    for name, (x, y) in samples.items():
        if not x:
            raise ValueError(f"No usable observations for {name}")
        model[name] = fit(x, y, logistic=name in ("survival", "fruiting"))
    x, y = samples["growth"]
    residuals = [target - sum(a * b for a, b in zip(row, model["growth"]))
                 for row, target in zip(x, y)]
    model["growth_variance"] = max(0.05, sum(r * r for r in residuals) / len(residuals))
    return model, {name: len(y) for name, (_, y) in samples.items()}


def normal_pdf(value, mean, variance):
    return math.exp(-0.5 * (value - mean) ** 2 / variance) / math.sqrt(2 * math.pi * variance)


def project(model, features, climate, interval, years=100, mesh_points=100):
    """Return annual population multiplier for one deterministic burn interval.

    ``model`` and ``features`` come from fit_vital_rates and Features.
    ``climate`` is [annual local temperature (C), annual precipitation (mm)].
    ``interval`` is years between burns; ``years`` is simulated years;
    ``mesh_points`` is the number of log-size bins. The return value is the
    geometric mean of yearly N_next/N_current over the last min(years,
    5*interval) years. A value of 1 is a stable population.
    """
    lo, hi = features.size_bounds
    width = hi - lo
    lo -= 0.1 * width
    hi += 0.3 * width
    step = (hi - lo) / mesh_points
    mesh = [lo + (i + 0.5) * step for i in range(mesh_points)]
    recruit = [normal_pdf(z, SEEDLING_MEAN, SEEDLING_VARIANCE) for z in mesh]
    recruit_sum = sum(recruit)
    recruit = [value / recruit_sum for value in recruit]
    population = recruit[:]
    kernels = {}
    for age in range(1, interval + 1):
        columns = []
        for z in mesh:
            vital = []
            for site in features.sites:
                row = {"site": site, "size_prevTS": str(z), "burn": str(age),
                       CLIMATE[0]: str(climate[0]), CLIMATE[1]: str(climate[1])}
                x = features.vector(row)
                predict = lambda key: sum(a * b for a, b in zip(x, model[key]))
                survival = sigmoid(predict("survival"))
                fruiting = sigmoid(predict("fruiting"))
                mean_size = predict("growth")
                fruits = math.exp(min(20, predict("fruit_count") + z))
                vital.append((survival, fruiting, mean_size, fruits))
            growth = [0.0] * mesh_points
            reproduction = 0.0
            for survival, fruiting, mean_size, fruits in vital:
                density = [normal_pdf(next_z, mean_size, model["growth_variance"])
                           for next_z in mesh]
                total = sum(density)
                if total < 1e-300:
                    density = [float(i == min(range(mesh_points), key=lambda j: abs(mesh[j] - mean_size)))
                               for i in range(mesh_points)]
                    total = 1.0
                for i in range(mesh_points):
                    growth[i] += survival * density[i] / total / len(vital)
                # Fruit production is conditional on being alive at the next census.
                reproduction += survival * fruiting * fruits * RECRUITMENT_PER_FRUIT / len(vital)
            columns.append([growth[i] + reproduction * recruit[i] for i in range(mesh_points)])
        kernels[age] = columns
    growth_factors = []
    for year in range(years):
        age = year % interval + 1
        columns = kernels[age]
        next_pop = [0.0] * mesh_points
        for j, count in enumerate(population):
            if count:
                for i in range(mesh_points):
                    next_pop[i] += columns[j][i] * count
        old_total = sum(population)
        new_total = sum(next_pop)
        if new_total <= 0:
            return 0.0
        growth_factors.append(new_total / old_total)
        population = [value / new_total for value in next_pop]
    # Geometric mean over the final 5 full burn cycles, after transient dynamics.
    tail = growth_factors[-min(years, 5 * interval):]
    return math.exp(sum(math.log(value) for value in tail) / len(tail))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DATA)
    parser.add_argument("--burn-events", type=Path, help="CSV with site,year or site,burn_date")
    parser.add_argument("--intervals", type=int, nargs="+", default=[1, 2, 3, 5, 8, 10, 12, 15, 20])
    parser.add_argument("--years", type=int, default=100)
    parser.add_argument("--output", type=Path, help="Write JSON report here")
    args = parser.parse_args()
    if args.years < 1 or any(i < 1 for i in args.intervals):
        parser.error("Years and fire-return intervals must be positive")
    rows = read_rows(args.data)
    history = integrate_burn_events(rows, args.burn_events) if args.burn_events else burn_history(rows)
    for row in rows:
        key = row["site"], int(row["year"])
        if key in history:
            row["burn"] = str(history[key])
    features = Features(rows)
    model, counts = fit_vital_rates(rows, features)
    site_year = {(row["site"], row["year"]): row for row in rows if all(number(row[k]) is not None for k in CLIMATE)}
    climate = [sum(number(row[key]) for row in site_year.values()) / len(site_year)
               for key in CLIMATE]
    estimates = {str(i): project(model, features, climate, i, args.years) for i in args.intervals}
    report = {
        "model": "compact observational-data refit of Louthan et al. IPM",
        "warning": "Not the authors' AICc-selected model or their historical/future climate projection.",
        "burn_source": str(args.burn_events) if args.burn_events else "demoVFT.csv burn column",
        "burn_site_years": {f"{site}:{year}": age for (site, year), age in sorted(history.items())},
        "features": features.names,
        "scaling_mean_sd": features.scaling,
        "coefficients": model,
        "sample_sizes": counts,
        "fixed_climate": dict(zip(CLIMATE, climate)),
        "recruitment_per_fruit": RECRUITMENT_PER_FRUIT,
        "projection_settings": {"years": args.years, "mesh_points": 100,
                                "climate": "mean across observed site-years",
                                "fire": "deterministic interval; years since fire reset to 1"},
        "lambda_definition": "Geometric mean of yearly population multipliers over the final min(years, 5*interval) simulated years",
        "projected_lambda_by_fire_return_interval": estimates,
    }
    output = json.dumps(report, indent=2)
    if args.output:
        args.output.write_text(output + "\n", encoding="utf-8")
        print(f"Wrote {args.output}")
    else:
        print(output)


if __name__ == "__main__":
    main()
