# Simulation Code: Multivariate Longitudinal Finite Mixture Models

Code for the simulation study in:

> *Robustness of multivariate longitudinal finite mixture models to covariance misspecification*  
> (manuscript under review)

Three models are compared. The general model is the **Multivariate Covariance Pattern Growth Mixture Model (MCPGMM)**; GBMTM and MLCGA are nested special cases differing only in residual covariance structure.

| Abbreviation | Full name | Covariance structure |
|:---|:---|:---|
| **MCPGMM** | Multivariate Covariance Pattern Growth Mixture Model | Time-varying residual variance + Toeplitz within-outcome serial correlation |
| **GBMTM** | Multivariate Group-Based Trajectory Model | Constant residual variance, no serial correlation |
| **MLCGA** | Multivariate Latent Class Growth Analysis | Time-varying residual variance, no serial correlation |

The study evaluates class enumeration accuracy and classification performance across 64 crossed design conditions (sample size × class separation × between-outcome correlation × covariance structure), with 200 replications per condition.

---

## Repository structure

```
simulation_code/
├── README.md
├── requirements.R                    Install and report all required R packages
├── run_models.R                      Standalone script to generate .inp files and run Mplus
│
├── data_generation/                  Data-generating and Mplus execution scripts
│   ├── GBMTM_DGM.R
│   ├── MLCGA_DGM.R
│   └── MCPGMM_DGM.R
│
├── mplus_templates/
│   ├── enum/                         MplusAutomation templates: class enumeration (K = 1..6)
│   │   ├── GBMTM_enum.txt
│   │   ├── MLCGA_enum.txt
│   │   └── MCPGMM_enum.txt
│   └── mc/                           MplusAutomation templates: MC classification (K = 3 fixed)
│       ├── GBMTM_mc.txt
│       ├── MLCGA_mc.txt
│       └── MCPGMM_mc.txt
│
├── ovl/                              Overlap Coefficient (OVL) utility
│   ├── ovl_trajectories.R
│   └── tests/
│       └── run_ovl_tests.R
│
├── abc/                              Area Between Curves (ABC) label-switching utility
│   ├── abc_trajectories.R
│   └── tests/
│       ├── test_abc_trajectories.R   testthat unit tests
│       ├── run_tests.R               Standalone test runner
│       └── run_demo.R                Interactive worked example
│
└── application/                      Empirical application: fitted models for the real-data example
    ├── GBMTM_7class.inp              7-class GBMTM
    ├── MLCGA_7class.inp              7-class MLCGA
    └── MCPGMM_7class.inp             7-class MCPGMM
```

---

## Simulation design

### Base parameters

| Parameter | Value |
|:---|:---|
| Replications per condition | 200 |
| Sample sizes (*N*) | 500, 1000 |
| Repeated measurements (*m*) | 5 |
| True number of classes (*K*) | 3 |
| Class proportions | 25% / 50% / 25% |
| Trajectory type | Quadratic polynomial |
| Random effects | None (fixed-effects mixture) |

### Data-generating covariance structures

| Feature | GBMTM | MLCGA | MCPGMM |
|:---|:---|:---|:---|
| Within-outcome serial correlation | None | None | Toeplitz (lag-1: 0.45, lag-2: 0.40, lag-3: 0.35, lag-4: 0.30) |
| Residual variance | Constant over time | Time-varying | Time-varying |
| Between-outcome correlation | Contemporaneous (ρ = 0.45) | Contemporaneous (ρ = 0.45) | Contemporaneous (ρ = 0.45) |

### Crossed design factors

- **Sample size**: N = 500 or N = 1000
- **Class separation**: high or low (varied via fixed-effect coefficients)
- **Between-outcome correlation**: high (ρ = 0.45) or low (ρ = 0.10)
- **Covariance (mis)specification**: within-outcome serial correlation and time-dependence of residual variance

---

## Requirements

### R

R 4.5.2 (developed and tested). Versions ≥ 4.1 will likely work but are untested.

Run once before using any scripts:

```r
source("requirements.R")
```

This installs all required packages and reports installed versions.

| Package | Version | Purpose |
|:---|:---|:---|
| `MplusAutomation` | 1.1.1 | Generate `.inp` files, run Mplus, extract output |
| `MASS` | 7.3.64 | Multivariate normal simulation (`mvrnorm`) |
| `stringr` | 1.5.1 | String helpers for parsing Mplus output titles |
| `plyr` | 1.8.9 | Collating results across replications (`ldply`) |
| `testthat` | 3.2.3 | Unit tests for `abc/` and `ovl/` only |

The `ovl/` and `abc/` modules have no package dependencies beyond base R and can be sourced independently.

For fully reproducible environments with all transitive dependencies pinned:

```r
install.packages("renv")
renv::restore()
```

### Mplus

[Mplus](https://www.statmodel.com/) v8.7 or later (commercial licence required).

The DGM scripts detect the Mplus executable automatically: on Windows they search the default install location (`C:/Program Files/Mplus/Mplus.exe`); on Mac/Linux they expect `mplus` on `PATH`.

---

## Usage

### Step 1: Generate data and run enumeration models

1. Open the relevant DGM script (`data_generation/GBMTM_DGM.R`, `MLCGA_DGM.R`, or `MCPGMM_DGM.R`).
2. Set `wd` to your working directory. This must contain a `Templates/` subfolder.  
   **Windows:** use a short path with no spaces (e.g. `D:/sim/GBMTM`). Mplus wraps DATA file paths longer than approximately 80 characters across lines, silently dropping the space at the break and causing a file-not-found error.
3. Copy the appropriate enumeration template from `mplus_templates/enum/` into `Templates/`.
4. Replace all `PATH_TO_DATA_DIR` placeholders in the template with your `wd` path.
5. Source the script. It will simulate 200 datasets, write `.dat` files, generate and run all Mplus `.inp` files, then extract fit statistics (AIC, BIC, entropy, etc.) and save them as `df.Rda`.

### Step 2: Run MC classification models

1. Copy the appropriate MC template from `mplus_templates/mc/` into a `Templates/` subfolder.
2. Replace `PATH_TO_DATA_DIR` and `PATH_TO_MC_DIR` placeholders.
3. Source `run_models.R` with `wd` set to your MC directory.

### Adapting to a different condition

All manipulable parameters are in the `# CONFIGURATION` block at the top of each DGM script.

| Parameter | Variable | Notes |
|:---|:---|:---|
| Sample size | `N` | Subjects per dataset |
| Class proportions | `pr` | Numeric vector, length K, sums to 1 |
| Trajectory shapes (Y1) | `bmeanmat` | 3 × K matrix of fixed-effect coefficients |
| Trajectory shapes (Y2) | `bmeanmat2` | 3 × K matrix of fixed-effect coefficients |
| Residual variance | `vary1`, `vary2` or `vary1_tv`, `vary2_tv` | Scalar (GBMTM) or time-point vector (MLCGA, MCPGMM) |
| Within-outcome correlation | `rho_within` | Toeplitz structure (MCPGMM only) |
| Between-outcome correlation | `rho_between` | Scalar contemporaneous correlation |

---

## Mplus template syntax

Templates use [MplusAutomation](https://michaelhallquist.github.io/MplusAutomation/) template syntax. The `[[init]]` block defines iterators over replications and class solutions; `[[classes = n]]` blocks insert model syntax conditionally for a specific K.

| Setting | Enum templates | MC templates | Notes |
|:---|:---|:---|:---|
| `PATH_TO_DATA_DIR` | Replace | Replace | Directory containing `Mplus_1.dat` ... `Mplus_n.dat` |
| `PATH_TO_MC_DIR` | n/a | Replace | Output path for MC `.out` files and saved posterior probabilities |
| `STARTS` | 1000 250 | 5000 1250 | MC stage uses more starts to reduce local-maxima risk |
| `PROCESSORS` | 4 | 16 | Adjust to available CPU cores |

---

## Overlap Coefficient (OVL)

`ovl/ovl_trajectories.R` is a self-contained utility for computing pairwise Overlap Coefficients between latent classes in a longitudinal mixture model. It is used in the paper to characterise class separation across simulation conditions (Supplementary Materials, Section S.11).

The OVL is the integral of the minimum of two class-specific normal densities at each measurement occasion, averaged across occasions. OVL = 0 indicates perfectly separated classes; OVL = 1 indicates complete overlap. When multiple outcomes are provided, matrices are returned per outcome and optionally averaged. It is used in the paper to characterise class separation across simulation conditions and to assist practitioners in identifying over-extracted solutions.

### Quick start

```r
source("ovl/ovl_trajectories.R")

# Time points and design matrix matching the DGM (quadratic, 5 occasions)
t_vec <- seq(0, 7, length.out = 5)   # {0, 1.75, 3.5, 5.25, 7}
X     <- cbind(1, t_vec, t_vec^2)

# Class-specific trajectory means (evaluated at observed time points)
means <- list(
  Class_1 = X %*% c(0,  1.752, -0.140),
  Class_2 = X %*% c(0,  0.000,  0.000),
  Class_3 = X %*% c(0, -1.752,  0.140)
)

# Shared residual covariance (scalar variance, independent over time)
covs <- replicate(3, 3.125 * diag(5), simplify = FALSE)

ovl <- compute_ovl(means, covs)
print_ovl(ovl)
```

Run `Rscript ovl/tests/run_ovl_tests.R` to execute the test suite (19 tests).

**Reference:** Nowakowska, E., Koronacki, J., & Lipovetsky, S. (2014). Tractable measure of component overlap for Gaussian mixture models. *arXiv:1407.7172 [math.ST]*. https://arxiv.org/abs/1407.7172

---

## Area Between Curves (ABC): label switching

`abc/abc_trajectories.R` is a self-contained utility for resolving label switching in latent class mixture models via an Area Between Curves criterion.

### The problem

Class labels in finite mixture models are statistically unidentified: the class labelled "1" in one replication may correspond to population class "3" in another. Before summarising results across replications (e.g. computing parameter recovery bias or classification accuracy), estimated labels must be aligned with the population labels used to generate the data.

### The algorithm

For a single replication with K classes and P outcomes:

1. Enumerate all K! relabellings of the K estimated classes.
2. For each relabelling, compute the mean Area Between Curves between the estimated and population polynomial trajectories across all classes and outcomes.
3. Select the relabelling that minimises this mean ABC.
4. Return the optimal permutation and (optionally) a relabelled vector of class memberships.

The ABC is computed **analytically**, not via numerical quadrature. For polynomial trajectories of degree d, the difference polynomial has degree d with antiderivative F(t) = ∑ᵢ Δbᵢ · tⁱ⁺¹/(i+1). Sign-change points within the integration domain are found via `polyroot()` (base R). This yields machine-precision results (ABC = 0 exactly when estimated equals population) and runs efficiently inside large simulation loops. Population and estimated trajectories may have different polynomial degrees; the shorter coefficient vector is zero-padded to the degree of the longer.

### Usage

```r
source("abc/abc_trajectories.R")

# Population coefficient matrices from the DGM (K × d: rows = classes, cols = coefficients)
pop_y1 <- t(bmeanmat)    # transpose from d × K to K × d
pop_y2 <- t(bmeanmat2)

for (i in seq_len(n_replications)) {

  # Estimated trajectory coefficients for this replication (K × d matrices)
  est_y1 <- ...
  est_y2 <- ...

  # Estimated class memberships (integer vector, values in 1:K)
  est_classes <- ...

  abc_res <- compute_abc(
    pop_betas         = list(Y1 = pop_y1, Y2 = pop_y2),
    est_betas         = list(Y1 = est_y1, Y2 = est_y2),
    t_min             = 0,
    t_max             = 7,
    class_assignments = est_classes
  )

  aligned_classes    <- abc_res$relabeled_assignments
  est_y1_aligned     <- est_y1[abc_res$optimal_permutation, ]
  est_y2_aligned     <- est_y2[abc_res$optimal_permutation, ]
}
```

### Diagnostics and visualisation

```r
# Formatted summary: optimal permutation, grand mean ABC, per-class and per-outcome breakdown,
# and the full K!-length vector of mean ABC values for every permutation tried
print_abc(abc_res)

# P × K panel grid: population trajectory (solid), optimally-relabelled estimated
# trajectory (dashed), shaded ABC region per panel
plot_abc(abc_res,
         pop_betas = list(Y1 = pop_y1, Y2 = pop_y2),
         est_betas = list(Y1 = est_y1, Y2 = est_y2))
```

Run `Rscript abc/tests/run_tests.R` to execute the test suite (21 tests, including cubic and mixed-degree cases).

---

## Empirical application

`application/` contains the Mplus `.inp` files for the real-data example reported in the paper. All three models are fitted to the same dataset (`all.csv`) with 7 latent classes, 2 outcomes (`dm`, `fp`), and 4 measurement occasions (0, 4, 8, 12).

| File | Model | Residual variance | Serial correlation |
|:---|:---|:---|:---|
| `GBMTM_7class.inp` | GBMTM | Constant | None |
| `MLCGA_7class.inp` | MLCGA | Time-varying | None |
| `MCPGMM_7class.inp` | MCPGMM | Time-varying | Toeplitz (lags 1-3) |

All three models use quadratic polynomial trajectories, no random effects, and a class-invariant contemporaneous cross-outcome correlation. The MCPGMM additionally estimates lag-1, lag-2, and lag-3 within-outcome correlations via `MODEL CONSTRAINT`, with covariances parameterised as `rho(lag) * SD(t1) * SD(t2)`.

To run: update the `FILE IS` path in each `.inp` to point to your copy of `all.csv`, then run via Mplus or `MplusAutomation::runModels()`.

---

## Citation

If you use this code, please cite the associated paper:

> *Robustness of multivariate longitudinal finite mixture models to covariance misspecification* (manuscript under review)

This section will be updated with full citation details upon publication.

```bibtex
@article{,
  title   = {Robustness of multivariate longitudinal finite mixture models
             to covariance misspecification},
  journal = {},
  year    = {},
  doi     = {}
}
```
