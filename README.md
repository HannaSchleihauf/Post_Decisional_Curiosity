# Chimpanzees and Human Children Explore Unchosen Options

Data, analysis code, and power simulations for:

> Schleihauf, H., Haux, L., Ghavamizadeh, N., Fischer, J., Engelmann, J. M.\*, & Herrmann, E.\* (2026). Chimpanzees and human children explore unchosen options. *Proceedings of the Royal Society B*. DOI: *to be added*
>
> \* shared last authorship

Supplementary materials and example videos: **OSF – <https://osf.io/nu67c/>**

---

## Overview

Post-decisional information search, meaning looking for information about an option you didn't choose after the choice has been made, is a behavioural marker of **epistemic curiosity**. When this search goes up after a person **had control over the decision** or **received a bad outcome**, it also indicates **counterfactual curiosity**: wanting to know what would have happened had you chosen differently.

In two preregistered studies we compared semi-free-ranging chimpanzees (*Pan troglodytes*) with 2- to 5-year-old children using a nearly identical, non-verbal paradigm:

| Study | Sample | Design | Question |
|---|---|---|---|
| **1A** | 16 chimpanzees (14–37 yrs) | Within-subject; choice / no-choice / control; 54 trials each | Do participants peek into the unchosen box? Does peeking depend on decisional control and on the value of the reward received? |
| **1B** | 144 children (2.05–5.98 yrs) | Between-subject; choice / no-choice / control; 4 trials each | Same as 1A, plus age |
| **2A** | 9 chimpanzees (subset of 1A) | Forced choice between peeking at a box that *had been* available vs. one that *never was*; 18 trials each | Is curiosity directed specifically at the forgone option? |
| **2B** | 48 children (4.04–5.99 yrs) | Same as 2A; 4 trials each | Same as 2A |

**Main findings.** Both species searched after making a decision, and much less when box contents were already known (control condition). This is evidence of epistemic curiosity. In children, search increased with age. Only children searched more after receiving no reward. In Study 2, both species peeked above chance into the box that had been available to choose. Children therefore show the full profile of counterfactual curiosity (sensitive to both control and outcome); chimpanzees show sensitivity to control only.

**Preregistrations (OSF):** Study 1 – [chimpanzees](https://osf.io/5hxcj/), [children](https://osf.io/kpae4/) · Study 2 – [chimpanzees](https://osf.io/3ngqt/), [children](https://osf.io/4zpu3/)

---

## Repository structure

```
.
├── Post_Decisional_Curiosity.Rproj   # open this first; sets the working directory
├── Plot_1A.pdf                       # Figure 2 (Study 1A results)
├── Plot_1B.pdf                       # Figure 3 (Study 1B results)
│
├── chimpanzees/
│   ├── study1/                       # Study 1A
│   │   ├── data/       data_study1A_chimpanzees.csv
│   │   ├── scripts/    Analysis_Study_1A.R
│   │   │               Study_1A_Power_Simulation_with_control.Rmd
│   │   │               Study_1A_Power_Simulation_test_conditions.Rmd
│   │   ├── functions/  helper functions (see below)
│   │   └── R_objects/  analysis_1.RData, power-simulation results
│   └── study2/                       # Study 2A
│       ├── data/       data_study2A_chimpanzees.csv
│       ├── scripts/    Analysis_2A.R
│       │               Study_2A_Power_Simulation.Rmd
│       │               Combine_figures.R
│       ├── functions/
│       └── R_objects/  analysis_2.RData, power-simulation results
│
└── children/
    ├── study1/                       # Study 1B
    │   ├── data/       data_study1B_children.csv
    │   ├── materials/  Instructions.docx, Counterbalancing_Sheet_Study1_Kids.xlsx
    │   ├── scripts/    Analysis_Study_1B.R
    │   │               Study_1B_Power_Simulation_with_control.Rmd
    │   │               Study_1B_Power_Simulation_test_conditions.Rmd
    │   ├── functions/
    │   └── R_objects/  analysis_1B.RData, power-simulation results
    └── study2/                       # Study 2B
        ├── data/       data_study2B_children.csv
        ├── scripts/    Analysis_2B.R
        │               Study_2B_Power_Simulation.Rmd
        ├── functions/
        └── R_objects/  Analysis_2.RData, power-simulation results
```

### Helper functions (`*/functions/`)

Based on GLMM utility functions by Roger Mundry; used across all analyses:

| File | Purpose |
|---|---|
| `diagnostic_fcns.r` | Model diagnostics (overdispersion, random-effect distributions, collinearity) |
| `glmm_stability.r` | Model stability: refits the model leaving out one level of the random effect at a time |
| `drop1_para.r` | Likelihood-ratio tests for individual fixed effects (parallelisable `drop1`) |
| `boot_glmm.r` | Parametric bootstrap of GLMM estimates and 95% CIs (used for the figures) |

---

## Reproducing the analyses

### Requirements

- **R ≥ 4.3** (analyses in the paper were run in R 4.3.3)
- Packages: `lme4`, `tidyverse`, `tidyselect`, `readxl`, `parallel`, `optimx`, `emmeans`, `car`, `irr`, `effects`, `gghalves`, `ggthemes`, `cowplot`, `ggpubr`
- Additional packages for the power simulations: `data.table`, `kyotil`

The analysis scripts install any missing core packages automatically.

### Steps

1. **Open `Post_Decisional_Curiosity.Rproj`** in RStudio. All analysis scripts use paths relative to the project root (e.g. `./chimpanzees/study1/data/...`). If you're not using RStudio, run `setwd()` to the repository root.
2. Run the analysis script for the study you're interested in:

   | Study | Script | Reported in |
   |---|---|---|
   | 1A | `chimpanzees/study1/scripts/Analysis_Study_1A.R` | Results Study 1A, Figure 2 |
   | 1B | `children/study1/scripts/Analysis_Study_1B.R` | Results Study 1B, Figure 3 |
   | 2A | `chimpanzees/study2/scripts/Analysis_2A.R` | Results Study 2A, Figure 4b |
   | 2B | `children/study2/scripts/Analysis_2B.R` | Results Study 2B, Figure 4c |

3. **To save time:** the bootstrapped confidence intervals (1,000 resamples) can take a long time to compute. Each script has a commented-out line near the top that loads the saved workspace (`R_objects/*.RData`) containing all fitted models and bootstraps. Uncomment it to explore the results without rerunning everything.

### What each analysis script does

1. Loads and cleans the data, and excludes familiarisation trials
2. Computes inter-rater reliability (percentage agreement and Cohen's κ on the 25% of trials double-coded)
3. Reports sample demographics and descriptive statistics
4. Fits the preregistered binomial GLMMs (`lme4::glmer`, `bobyqa` optimiser):
   - **Study 1:** peeking (yes/no) ~ condition × reward value (× age for children) + sex + trial, with random intercepts for individual and random slopes where identifiable; full–null model comparison followed by likelihood-ratio tests via `drop1`
   - **Study 2:** intercept-only model with z-transformed trial, testing whether peeking into the previously available box differs from chance (0.5)
5. Checks model assumptions and stability
6. Bootstraps CIs and produces the figures
7. Runs robustness analyses excluding the control condition (Study 1)
8. Saves the workspace to `R_objects/`

### Power simulations

The `*_Power_Simulation*.Rmd` files contain the a priori simulation-based power analyses reported in the preregistrations and the SI. They use paths relative to their own `scripts/` folder, so knit them from there. Their results are stored in the corresponding `R_objects/` folders.

---

## Data

All data are stored as one row per trial. Note that the delimiters differ:

| File | Delimiter | Read with |
|---|---|---|
| `data_study1A_chimpanzees.csv` | `;` | `read.csv2()` |
| `data_study1B_children.csv` | `;` | `read.csv2()` |
| `data_study2A_chimpanzees.csv` | `,` | `read.csv()` |
| `data_study2B_children.csv` | `;` | `read.csv2()` |

### Key variables

**General conventions**
- Columns ending in **`.reli`** hold the second coder's reliability coding of the preceding variable (filled in for the randomly selected 25% of trials).
- Reward values are coded as **`high`** (chimpanzees: three banana pieces; children: a small toy or sticker), **`low`** (chimpanzees: one carrot piece; children: a plain dot sticker), or **`none`**.

**Study 1 (1A and 1B)**

| Variable | Description |
|---|---|
| `id` | Participant identifier |
| `condition` | `choice`, `no-choice` (children: `no.choice`), or `control` |
| `trial` | Trial number |
| `left`, `right` | Reward value in the left and right box |
| `side.chosen` | Box the participant pulled |
| `food.received` / `reward.received` | Value of the reward obtained |
| `food.remaining` / `reward.remaining` | Value of the reward in the unchosen box |
| **`searched`** | **Outcome: did the participant peek into the unchosen box? (`yes`/`no`)** |
| `seconds.until.searched` | Latency to peek (s) |
| `value.remaining.seen` | Whether the remaining reward was better, the same, or worse |
| `turned.away`, `time.from.last.gaze...` | (1A only) Whether, and how quickly, the chimpanzee turned away after peeking; used as an indicator that they did not expect to get the remaining reward |
| `age`, `gender` | Age in years; sex/gender |
| `condition.trial` | (1B only) Trial label; `Fam_*` marks familiarisation trials, which are excluded from the analyses |

**Study 2 (2A and 2B)**

| Variable | Description |
|---|---|
| `left`, `middle`, `right` | Reward values; the middle box always contained a high-value reward |
| `side.available` | (2A) Which outer box had been available to choose |
| **`searched.in.available`** (2A) / **`searched.in.targeted.box`** (2B) | **Outcome: did the participant peek into the box that had been available? (`yes`/`no`)** |
| `trial.repeated`, `count.of.repeated.trials` | (2B) Repeated trials (participant chose an outer box or did not peek) |

**Children's demographics** (1B, 2B): `ethnic.info`, `education.parent1/2`, and `family.income` are parent-reported. `NA` or "prefer not to say" means no information was provided.

---

## Ethics

Chimpanzee research was approved by the relevant institutional ethics committee, the Chimpanzee Sanctuary & Wildlife Conservation Trust, and national wildlife and science authorities. It followed the PASA Primate Veterinary Healthcare Manual and ASAB guidelines. Chimpanzees were never deprived of food or water and participated voluntarily. Children's participation was approved by the institutional ethics committee, with written parental consent and verbal assent from the children. Child data are de-identified.

## Funding

H.S. was supported by the European Union's Horizon 2020 programme under a Marie Skłodowska-Curie grant. L.H. was supported by the Deutsche Forschungsgemeinschaft (DFG; HE 2768/11-1).

## Contact

Hanna Schleihauf – Department of Developmental Psychology, Utrecht University – <h.schleihauf@uu.nl>

## License

*To be added* (e.g. code under MIT, data under CC BY 4.0).
