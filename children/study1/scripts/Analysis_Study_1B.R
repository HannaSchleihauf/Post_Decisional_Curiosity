# Analysis Study 1 - Children

# Load packages -------------------------------------------------------------
# Define required packages
required_packages <- c(
  "lme4","readxl","tidyverse","tidyselect","parallel","optimx","emmeans","car",
  "irr"
)
# Check and install missing packages
installed_packages <- rownames(installed.packages())

for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}
# Source custom functions
source("./children/study1/functions/diagnostic_fcns.r")
source("./children/study1/functions/glmm_stability.r")
source("./children/study1/functions/drop1_para.r")
source("./children/study1/functions/boot_glmm.r")
options(scipen = 9999)

# Load data -------------------------------------------------------------

# To load all the objects, including model results and bootstraps, you can run the follwing line:
# load("./children/study1/R_objects/analysis_1B.RData")

xdata <-
  read.csv2("./children/study1/data/data_study1B_children.csv",
           header = TRUE, na = c("NA", "")
  )
# remove spaces at the end of the data
xdata <- as.data.frame(
  lapply(xdata, function(col) {
    if (is.character(col)) {
      sub("\\s+$", "", col)
    } else {
      col
    }
  }),
  stringsAsFactors = FALSE
)

# Data Wrangling --------------------------------------------------------
# Ignore familiarization trials
xdata <-
  subset(xdata, xdata$condition.trial == "choice" |
           xdata$condition.trial == "no.choice" |
           xdata$condition.trial == "control")
xdata$condition <-
  droplevels(as.factor(xdata$condition))

# Data Wrangling
xdata$searched <-
  as.factor(xdata$searched)
ftable(searched ~
         condition + reward.received, xdata)
xdata$condition <-
  factor(xdata$condition,
         levels = c("control", "choice", "no.choice"),
         labels = c("control", "choice", "no.choice")
  )
xdata$value.remaining.seen <-
  as.factor(xdata$value.remaining.seen)
xdata$value.remaining.seen <-
  relevel(xdata$value.remaining.seen,
          ref = "worse"
  )
xdata$trial.per.child <-
  as.numeric(ave(xdata$id,
                 list(xdata$id, xdata$id),
                 FUN = seq_along
  ))

# Participant information ------------------------------------------
table(xdata$age.group)
tapply(xdata$age, xdata$age.group, mean)
tapply(xdata$age, xdata$age.group, min)
tapply(xdata$age, xdata$age.group, max)
mean(xdata$age)
sd(xdata$age)

# Interrater Reliability ------------------------------------------------
interrater <-
  xdata %>%
  filter(!is.na(searched.reli))
# percentage of reliability coding 25%
nrow(interrater) / nrow(xdata)

# number of dyads for which the reliability coding was not the same
sum(interrater$searched.reli !=
      interrater$searched)
(sum(interrater$searched.reli ==
       interrater$searched)) / nrow(interrater)

# Cohen's Kappa
kappa_results <-
  kappa2(cbind(
    as.factor(interrater$searched.reli),
    interrater$searched
  ), "unweighted")
print(kappa_results)

# Descriptives
ftable(searched ~ condition + reward.received, xdata)
round(prop.table(ftable(searched ~ condition, xdata), margin = 1) * 100, 2)
round(prop.table(ftable(searched ~ reward.received, xdata), margin = 1) * 100, 2)

# Prepare data for model fitting------------------------------------
xx.fe.re <- fe.re.tab(
  fe.model = "searched ~ age*condition*reward.received +
                   gender + trial.per.child",
  re = "(1|id)", other.vars = c("age.group"),
  data = xdata
)
xx.fe.re$summary
t.data <- xx.fe.re$data
str(t.data)
# We do not have enough variation levels of reward received that have more than 2 observation per subject.
# Thus, we cannot include that random slope into the model.

# Means and SDs
t.data$searched_num <- as.numeric(t.data$searched) - 1

per_subject <- t.data %>%
  group_by(id, condition) %>%
  summarise(prop_searched = mean(searched_num, na.rm = TRUE),
            n_trials = n(),
            .groups = "drop")
per_subject_summary <- per_subject %>%
  group_by(condition) %>%
  summarise(mean_prop = mean(prop_searched, na.rm = TRUE),
            sd_prop   = sd(prop_searched, na.rm = TRUE),
            n_subjects = n(),
            .groups = "drop")

per_subject <- t.data %>%
  group_by(id, reward.received) %>%
  summarise(prop_searched = mean(searched_num, na.rm = TRUE),
            n_trials = n(),
            .groups = "drop")
per_subject_summary <- per_subject %>%
  group_by(reward.received) %>%
  summarise(mean_prop = mean(prop_searched, na.rm = TRUE),
            sd_prop   = sd(prop_searched, na.rm = TRUE),
            n_subjects = n(),
            .groups = "drop")

per_subject <- t.data %>%
  group_by(id, condition, reward.received) %>%
  summarise(prop_searched = mean(searched_num, na.rm = TRUE),
            n_trials = n(),
            .groups = "drop")

per_subject_summary <- per_subject %>%
  group_by(condition, reward.received) %>%
  summarise(mean_prop  = mean(prop_searched, na.rm = TRUE),
            sd_prop    = sd(prop_searched, na.rm = TRUE),
            n_subjects = n(),
            .groups = "drop")

## Center dummy variables$
## (necessary for the random effects in the model or plotting)
t.data$reward.received.low.code <-
  t.data$reward.received.low - mean(t.data$reward.received.low)
t.data$reward.received.none.code <-
  t.data$reward.received.none - mean(t.data$reward.received.none)
t.data$condition.no.choice.code <-
  t.data$condition.no.choice - mean(t.data$condition.no.choice)
t.data$condition.choice.code <-
  t.data$condition.choice - mean(t.data$condition.choice)
t.data$gender.code <-
  t.data$gender.male - mean(t.data$gender.male)
t.data$z.age <- scale(t.data$age)
t.data$z.trial <- scale(t.data$trial.per.child)

# Fitting models ----------------------------------------------------

contr <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000000))

full <-
  glmer(
    searched ~
      (reward.received + condition + z.age)^2  + gender + z.trial +
      (1 + z.trial | id), #(reward.received.low.code + reward.received.none.code)
    data = t.data, control = contr,
    family = binomial(link = "logit")
  )

red <-
  glmer(
    searched ~
      reward.received + (condition + z.age)^2 + gender + z.trial +
      (1 + z.trial | id), #(reward.received.low.code + reward.received.none.code)
    data = t.data, control = contr,
    family = binomial(link = "logit")
  )

main <-
  glmer(
    searched ~
      reward.received + condition + z.age + gender + z.trial +
      (1 + z.trial | id), #(reward.received.low.code + reward.received.none.code)
    data = t.data, control = contr,
    family = binomial(link = "logit")
  )

null <-
  glmer(
    searched ~
      1 +
      (1 + z.trial | id), #(reward.received.low.code + reward.received.none.code)
    data = t.data, control = contr,
    family = binomial(link = "logit")
  )

summary(full)$varcor
round(summary(full)$coefficients, 3)

tapply(
  (as.numeric(t.data$searched) - 1),
  list(t.data$condition, t.data$reward.received), mean
)

round(summary(full)$coefficients, 3)

## Test assumptions ----
overdisp.test(full) # no overdispersion
vif(main) # no colliniarity
ranef.diagn.plot(full)

## Checking model stability
m.stab.b <-
  glmm.model.stab(model.res = full, contr = contr, use = c("id"))
m.stab.b$detailed$warnings
as.data.frame(round(m.stab.b$summary[, -1], 3))
m.stab.plot(round(m.stab.b$summary[, -1], 3))
# Model stablility only with models that did not give a warning message
m.stab.b$detailed <- m.stab.b$detailed %>%
  filter(warnings == "none")
cbind(
  apply(m.stab.b$detailed, 2, min),
  apply(m.stab.b$detailed, 2, max)
)

# Model Comparisons -------------------------------------------------------

anova(full, null, test = "LRT")

## Reduced model comparisons
tests.full <- drop1p(
  model.res = full,
  contr = contr
)
round(tests.full$drop1.res, 3)
tests.red <- drop1p(
  model.res = red,
  contr = contr
)
round(tests.red$drop1.res, 3)
tests.main <- drop1p(
  model.res = main,
  contr = contr
)
round(tests.main$drop1.res, 3)

## First peek at effects
library("effects")
plot(effect("condition:z.age", red))
plot(effect("z.trial", red))
plot(effect("gender", red))
plot(effect(c("reward.received"), red))

## Pairwise comparisons
emm1 <- emmeans(main, ~condition)
summary(emm1, type = "response")
emmeans(main, pairwise ~ condition)
summary(pairs(emm1), type = "response")

emm2 <- emmeans(main, ~reward.received)
summary(emm2, type = "response")
emmeans(red, pairwise ~reward.received)
summary(pairs(emm2), type = "response")

# Bootstraps --------------------------------------------------------------
## Bootstraps of full model
boot.full <- boot.glmm.pred(
  model.res = full, excl.warnings = T,
  nboots = 1000, para = F, level = 0.95,
  use = c("condition", "reward.received", "z.age")
)
round(boot.full$ci.estimates, 3)
as.data.frame(round(boot.full$ci.estimates, 3))
m.stab.plot(round(boot.full$ci.estimates, 3))
boot.full$ci.predicted

## Bootstraps of plot.model
model.condition.age <-
  glmer(
    searched ~ condition * z.age +
      (reward.received.low.code + reward.received.none.code) + gender.male + z.trial +
      (1 + z.trial | id), #(reward.received.low.code + reward.received.none.code)
    data = t.data, control = contr,
    family = binomial(link = "logit")
  )

boot.condition.age <- boot.glmm.pred(
  model.res = model.condition.age,
  excl.warnings = T,
  nboots = 1000, para = F, level = 0.95,
  use = c("condition", "z.age")
)
round(boot.condition.age$ci.estimates, 3)
as.data.frame(round(boot.condition.age$ci.estimates, 3))
m.stab.plot(round(boot.condition.age$ci.estimates, 3))
head(boot.condition.age$ci.predicted)

model.reward <-
  glmer(
    searched ~ reward.received +
      z.age * (condition.no.choice.code + condition.choice.code) + gender.male + z.trial +
      (1 + z.trial | id),
    data = t.data, control = contr,
    family = binomial(link = "logit")
  )

summary(model.reward)$coefficients

boot.reward <-
  boot.glmm.pred(
    model.res = model.reward, excl.warnings = T,
    nboots = 1000, para = F,
    level = 0.95, use = c("reward.received")
  )
round(boot.reward$ci.estimates, 3)
as.data.frame(round(boot.reward$ci.estimates, 3))
m.stab.plot(round(boot.reward$ci.estimates, 3))
head(boot.reward$ci.predicted)

# Plotting Reward model-----------------------------------------------------
library(gghalves)
library(ggthemes)
library(cowplot)

xdata.agg <- xdata %>%
  mutate(searched.numeric = as.numeric(searched) - 1) %>%
  group_by(id, reward.received) %>%
  summarise(mean.resp = mean(searched.numeric, na.rm = T)) %>%
  ungroup()

xdata.agg$reward.received2 <-
  jitter(as.numeric(as.factor(xdata.agg$reward.received)), amount = 0.13)
xdata.agg$mean.resp2 <-
  jitter(xdata.agg$mean.resp, amount = 0.04)

exp1_plot_reward <-
  ggplot(data = xdata.agg, aes(x = factor(reward.received,
                                          levels = c("high", "low", "none")
  ), y = mean.resp2)) +

  # geom_line(aes(x = reward.received2, y = mean.resp2, group = id),
  #           color = "gray", lty = 1, alpha = .3) +

  geom_point(
    data = xdata.agg %>% filter(reward.received == "high"),
    aes(x = reward.received2), color = "chartreuse4", size = 1.5,
    alpha = .6
  ) +
  geom_point(
    data = xdata.agg %>% filter(reward.received == "low"),
    aes(x = reward.received2), color = "gold3", size = 1.5,
    alpha = .6
  ) +
  geom_point(
    data = xdata.agg %>% filter(reward.received == "none"),
    aes(x = reward.received2), color = "red3", size = 1.5,
    alpha = .4
  ) +
  geom_violin(
    data = xdata.agg %>% filter(reward.received == "high"),
    aes(x = reward.received, y = mean.resp),
    position = position_nudge(x = 0.00),
    fill = "chartreuse4", alpha = .2
  ) +
  geom_violin(
    data = xdata.agg %>% filter(reward.received == "low"),
    aes(x = reward.received, y = mean.resp),
    position = position_nudge(x = 0.00),
    fill = "gold3", alpha = .2
  ) +
  geom_violin(
    data = xdata.agg %>% filter(reward.received == "none"),
    aes(x = reward.received, y = mean.resp),
    position = position_nudge(x = 0.00),
    fill = "red3", alpha = .2
  ) +
  geom_errorbar(
    data = boot.reward$ci.predicted %>%
      filter(reward.received == "high"),
    aes(
      x = as.numeric(reward.received) + 0.25, y = fitted,
      ymin = lower.cl, ymax = upper.cl
    ), color = "chartreuse4",
    width = 0.1, size = 1
  ) +
  geom_errorbar(
    data = boot.reward$ci.predicted %>%
      filter(reward.received == "low"),
    aes(
      x = as.numeric(reward.received) + 0.25, y = fitted,
      ymin = lower.cl, ymax = upper.cl
    ), color = "gold3",
    width = 0.1, size = 1
  ) +
  geom_errorbar(
    data = boot.reward$ci.predicted %>%
      filter(reward.received == "none"),
    aes(
      x = as.numeric(reward.received) + 0.25, y = fitted,
      ymin = lower.cl, ymax = upper.cl
    ), color = "red3",
    width = 0.1, size = 1
  ) +
  geom_point(
    data = boot.reward$ci.predicted %>%
      filter(reward.received == "high"),
    aes(x = as.numeric(reward.received) + 0.25, y = fitted),
    color = "darkgreen", size = 2.5
  ) +
  geom_point(
    data = boot.reward$ci.predicted %>%
      filter(reward.received == "low"),
    aes(x = as.numeric(reward.received) + 0.25, y = fitted),
    color = "gold4", size = 2.5
  ) +
  geom_point(
    data = boot.reward$ci.predicted %>%
      filter(reward.received == "none"),
    aes(x = as.numeric(reward.received) + 0.25, y = fitted),
    color = "darkred", size = 2.5
  ) +
  scale_x_discrete(
    limits =
      c(
        "high",
        "low",
        "none"
      ),
    name =
      "Reward received",
    labels =
      c(
        "High-value",
        "Low-value",
        "No reward"
      )
  ) +
  scale_y_continuous(
    name =
      "Post-decisional information-search",
    labels = scales::percent
  ) +
  labs(
    title =
      "Effect of received reward value"
  ) +
  theme_classic() +
  theme(
    axis.ticks.x =
      element_blank(),
    axis.title =
      element_text(size = 15),
    axis.text =
      element_text(size = 13),
    strip.text.x =
      element_text(size = 13),
    plot.margin =
      unit(c(1, 1, 1, 1), "cm"),
    axis.title.y.left =
      element_text(vjust = 3),
    plot.title =
      element_text(color = "black", size = 15, face = "bold"),
    axis.title.x =
      element_text(
        margin =
          margin(t = 10, r = 0, b = 0, l = 0)
      )
  )
exp1_plot_reward

# Plotting Condition model---------------------------------------------------
xdata$z.age <- scale(xdata$age)

xdata.agg.condition <- xdata %>%
  mutate(searched.numeric = as.numeric(searched) - 1) %>%
  group_by(id, z.age, condition) %>%
  summarise(mean.resp = mean(searched.numeric, na.rm = T)) %>%
  ungroup()
xdata.agg.age.group.condition$age.group <- as.numeric(xdata.agg.age.group.condition$age.group) + 0.5

xdata.agg.condition$mean.resp2 <-
  jitter(xdata.agg.condition$mean.resp, amount = 0.04)

xdata.agg.age.group.condition <- xdata %>%
  mutate(searched.numeric = as.numeric(searched) - 1) %>%
  group_by(age.group, condition) %>%
  summarise(mean.resp = mean(searched.numeric, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    # if age.group is a factor like "2","3",..., go through as.character first
    age.group.mid = as.numeric(as.character(age.group)) + 0.5,
    # put it on the SAME z-scale as z.age
    z.age.group   = (age.group.mid - mean(xdata$age)) / sd(xdata$age)
  )

exp1_plot_condition <-
  ggplot() +
  geom_point(
    data = xdata.agg.condition, shape = 21,
    aes(x = z.age, y = mean.resp2, fill = condition, color = condition),
    size = 1.7, alpha = .4
  ) +
  # geom_point(
  #   data = xdata.agg.age.group.condition, shape = 23,
  #   aes(x = z.age.group, y = mean.resp, fill = condition, color = condition),
  #   size = 3, alpha = 1
  # ) +
  geom_ribbon(
    data = boot.condition.age$ci.predicted,
    aes(x = z.age, ymin = lower.cl, ymax = upper.cl,
        group = condition, fill = condition, color = condition),
    alpha = 0.3) +
  scale_fill_manual(values = c("dodgerblue", "darkorange", "darkorange")) +
  geom_line(
    data = boot.condition.age$ci.predicted,
    aes(
      x = z.age,
      y = fitted,
      color = condition), size = 0.8) +
  scale_color_manual(values = c("dodgerblue4", "darkorange4", "darkorange4")) +
  #facet_wrap(~condition, strip.position = "bottom") +
  facet_wrap(
    ~condition, strip.position = "bottom",
    labeller = as_labeller(c(
      control = "Control",
      choice = "Choice",
      no.choice = "No choice"
    ))
  ) +
  scale_x_continuous(
    name = "Condition * Age",
    breaks = c(
      (2 - mean(xdata$age)) / sd(xdata$age),
      (3 - mean(xdata$age)) / sd(xdata$age),
      (4 - mean(xdata$age)) / sd(xdata$age),
      (5 - mean(xdata$age)) / sd(xdata$age),
      (6 - mean(xdata$age)) / sd(xdata$age)
    ),
    labels = c(2, 3, 4, 5, 6)
  ) +
  scale_y_continuous(
    name = "Post-decisional information-search",
    labels = scales::percent
  ) +
  labs(title = "Effect of the condition*age interaction") +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_text(
      size = 13,
      margin = margin(t = 0, r = 0, b = 8, l = 0)
    ),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.title.y.left = element_text(vjust = 3),
    plot.title = element_text(color = "black", size = 15, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))
  )
exp1_plot_condition

# Combine plots
library(ggpubr)
# theme_set(theme_pubr())
figure <- ggarrange(exp1_plot_condition, exp1_plot_reward,
                    labels = c("(a)", "(b)"),
                    hjust = -1, # hjust = 0 for left alignment
                    vjust = 2.5, # hjust = 0 for left alignment
                    ncol = 2, nrow = 1
)
figure

# Post-hoc Analysis without Control Condition ----------------------------------

## Subset the data
without_contr <- subset(t.data, t.data$condition == "choice" |
                          t.data$condition == "no.choice")

## Fit the same models as above for both test conditions only
full.no.cont <-
  glmer(
    searched ~
      (reward.received + condition + z.age)^2 + gender + z.trial +
      (1 + z.trial | id),
    data = without_contr, control = contr,
    family = binomial(link = "logit")
  )

main.no.cont <-
  glmer(
    searched ~
      reward.received + condition + z.age + gender + z.trial +
      (1 + z.trial | id),
    data = without_contr, control = contr,
    family = binomial(link = "logit")
  )

null.no.cont <-
  glmer(
    searched ~
      1 +
      (1 + z.trial | id),
    data = without_contr, control = contr,
    family = binomial(link = "logit")
  )

round(summary(full.no.cont)$coefficients, 3)

## Test assumptions ----
overdisp.test(full.no.cont) # no overdispersion
vif(main.no.cont) # no colliniarity
ranef.diagn.plot(full.no.cont)

## Checking model stability
m.stab.b <-
  glmm.model.stab(model.res = full.no.cont, contr = contr, use = c("id"))
m.stab.b$detailed$warnings
as.data.frame(round(m.stab.b$summary[, -1], 3))
m.stab.plot(round(m.stab.b$summary[, -1], 3))
# Model stablility only with models that did not give a warning message
m.stab.b$detailed <- m.stab.b$detailed %>%
  filter(warnings == "none")
cbind(
  apply(m.stab.b$detailed, 2, min),
  apply(m.stab.b$detailed, 2, max))


# Model Comparisons -------------------------------------------------------

## Full-null models comparison
round(anova(full.no.cont, null.no.cont, test = "Chisq"), 3)

## Reduced model comparisons
tests.full <- drop1p(
  model.res = full.no.cont,
  contr = contr
)
round(tests.full$drop1.res, 3)

tests.main <- drop1p(
  model.res = main.no.cont,
  contr = contr
)
round(tests.main$drop1.res, 3)

## First peek at effects
library("effects")
plot(effect("condition:z.age", red))
plot(effect(c("reward.received"), main.no.cont))
plot(effect(c("z.trial"), main.no.cont))

## Bootstraps of full model
boot.full <- boot.glmm.pred(
  model.res = full.no.cont, excl.warnings = T,
  nboots = 1000, para = F, level = 0.95,
  use = c("condition", "reward.received", "z.age")
)
round(boot.full$ci.estimates, 3)
as.data.frame(round(boot.full$ci.estimates, 3))
m.stab.plot(round(boot.full$ci.estimates, 3))
boot.full$ci.predicted

save.image("./children/study1/R_objects/analysis_1B.RData")
