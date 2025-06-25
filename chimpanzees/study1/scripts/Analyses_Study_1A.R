# Analysis Study 1 - Chimpanzees

# Load packages -------------------------------------------------------------
# Define required packages
required_packages <- c(
  "lme4", "readxl", "tidyverse", "tidyselect", "parallel", "optimx", "emmeans", "car",
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
source("./study1/functions/diagnostic_fcns.r")
source("./study1/functions/glmm_stability.r")
source("./study1/functions/drop1_para.r")
source("./study1/functions/boot_glmm.r")

# Load data -------------------------------------------------------------

# To load all the objects, including model results and bootstraps, you can run the follwing line:
# load("./study1/R_objects/analysis_1.RData")

xdata <-
  read.csv2("./study1/data/data_study1A_chimpanzees.csv",
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
xdata$searched <-
  as.factor(xdata$searched)
xdata$condition <-
  factor(xdata$condition,
    levels = c("control", "choice", "no-choice")
  )
xdata$value.remaining.seen <-
  as.factor(xdata$value.remaining.seen)
xdata$value.remaining.seen <-
  relevel(xdata$value.remaining.seen, ref = "worse")

xdata$turned.away.within.2 <- NA
xdata$turned.away.within.2[xdata$searched == "yes" &
  xdata$time.from.last.gaze.to.back.of.head.to.camera >= 2] <- "yes"
xdata$turned.away.within.2[xdata$searched == "yes" &
  xdata$turned.away == "no"] <- "no"
xdata$turned.away.within.2[xdata$searched == "yes" &
  xdata$time.from.last.gaze.to.back.of.head.to.camera < 2] <- "no"

xdata <- xdata[order(xdata$name, xdata$day, xdata$trial), ]
xdata$trial.per.condition <-
  as.numeric(ave(xdata$name,
    list(xdata$name, xdata$condition),
    FUN = seq_along
  ))
xdata$trial.total <-
  as.numeric(ave(xdata$name,
    list(xdata$name),
    FUN = seq_along
  ))

ftable(searched ~ turned.away.within.2, xdata)
ftable(searched ~ name + condition, xdata)

# Interrater Reliability ------------------------------------------------
interrater <-
  xdata %>%
  filter(!is.na(searched.reli))
# percentage of reliability coding 50%
nrow(interrater) / nrow(xdata)
# number of ind for which the reliability coding was not the same
sum(interrater$searched.reli !=
  interrater$searched)
(sum(interrater$searched.reli ==
  interrater$searched)) / nrow(interrater)

# Cohen's Kappa
library(irr)
kappa_results <-
  kappa2(cbind(
    as.factor(interrater$searched),
    as.factor(interrater$searched.reli)
  ), "unweighted")
print(kappa_results)

# Chimp demographics------------------------------------------------
min(xdata$age)
max(xdata$age)
table(xdata$gender) / 54

# Descriptives
ftable(searched ~ condition + food.received, t.data)
round(tapply(
  (as.numeric(t.data$searched) - 1),
  list(t.data$condition, t.data$food.received), mean
), 3)

# Prepare data for model fitting------------------------------------
xx.fe.re <- fe.re.tab(
  fe.model = "searched ~
  condition*food.received + trial.total",
  re = "(1|name)",
  other.vars = c("age", "gender", "trial.per.condition"),
  data = xdata
) # maybe add age
xx.fe.re$summary
t.data <- xx.fe.re$data
str(t.data)

## Center dummy variables ----
## (necessary for the random effects in the model and potentially plotting)
t.data$food.received.low.code <-
  t.data$food.received.low - mean(t.data$food.received.low)
t.data$food.received.none.code <-
  t.data$food.received.none - mean(t.data$food.received.none)
t.data$condition.no.choice.code <-
  t.data$condition.no.choice - mean(t.data$condition.no.choice)
t.data$condition.choice.code <-
  t.data$condition.choice - mean(t.data$condition.choice)
t.data$z.age <- scale(t.data$age)

# Fitting models ----------------------------------------------------
## (as pre-registered)
contr <-
  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000000))

full <-
  glmer(
    searched ~ (food.received + condition)^2 +
      (1 + (condition.no.choice.code + condition.choice.code) +
        (food.received.low.code + food.received.none.code) | name),
    data = t.data, control = contr, family = binomial(link = "logit")
  )

main <-
  glmer(
    searched ~ condition + food.received +
      (1 + (condition.no.choice.code + condition.choice.code) +
        (food.received.low.code + food.received.none.code) | name),
    data = t.data, control = contr, family = binomial(link = "logit")
  )

null <-
  glmer(
    searched ~ 1 +
      (1 + (condition.no.choice.code + condition.choice.code) +
        (food.received.low.code + food.received.none.code) | name),
    data = t.data, control = contr, family = binomial(link = "logit")
  )

summary(full)$varcor
round(summary(full)$coefficients, 3)

## Test assumptions ----
overdisp.test(full) # no overdispersion
vif(main) # no colliniarity
ranef.diagn.plot(full)

## Checking model stability
m.stab.b <-
  glmm.model.stab(model.res = full, contr = contr, use = c("name"))
m.stab.b$detailed$warnings
as.data.frame(round(m.stab.b$summary[, -1], 3))
m.stab.plot(round(m.stab.b$summary[, -1], 3))

# Model Comparisons -------------------------------------------------------

## Full-null models comparison
round(anova(full, null, test = "Chisq"), 3)

## Reduced model comparisons
tests.full <- drop1p(
  model.res = full,
  para = F, data = NULL, contr = contr,
  n.cores = c("all-1", "all"), to.del = NULL
)
round(tests.full$drop1.res, 3)
tests.main <- drop1p(
  model.res = main,
  para = F, data = NULL, contr = contr,
  n.cores = c("all-1", "all"), to.del = NULL
)
round(tests.main$drop1.res, 3)

## First peek at effects
library("effects")
plot(effect("food.received:condition", full))
plot(effect(c("condition"), main))

## Pairwise comparisons
emm1 <- emmeans(full, ~ food.received * condition)
summary(emm1, type = "response")
emmeans(full, pairwise ~ food.received | condition, type = "response")

emm2 <- emmeans(main, ~condition)
summary(emm2, type = "response")
emmeans(main, pairwise ~ condition, type = "response")

# Bootstraps --------------------------------------------------------------
## Bootstraps of full model
boot.full <-
  boot.glmm.pred(
    model.res = full,
    excl.warnings = T,
    nboots = 1000, para = F,
    level = 0.95,
    use = c("condition", "food.received")
  )

as.data.frame(round(boot.full$ci.estimates, 3))
m.stab.plot(round(boot.full$ci.estimates, 3))
boot.full$ci.predicted

## Bootstraps of plot.model
plot.model <-
  glmer(
    searched ~ condition +
      (food.received.low.code + food.received.none.code) +
      (1 + (condition.no.choice.code + condition.choice.code) +
        (food.received.low.code + food.received.none.code) | name),
    data = t.data,
    control = contr,
    family = binomial(link = "logit")
  )

boot.plot <- boot.glmm.pred(
  model.res = plot.model, excl.warnings = T,
  nboots = 1000, para = F, level = 0.95, use = c("condition")
)

as.data.frame(round(boot.plot$ci.estimates, 3))
m.stab.plot(round(boot.plot$ci.estimates, 3))
boot.plot$ci.predicted

model.reward <-
  glmer(
    searched ~
      food.received +
      (condition.no.choice.code + condition.choice.code) +
      (1 + (condition.no.choice.code + condition.choice.code) +
        (food.received.low.code + food.received.none.code) | name),
    data = t.data, control = contr,
    family = binomial(link = "logit")
  )

boot.reward <-
  boot.glmm.pred(
    model.res = model.reward, excl.warnings = T,
    nboots = 1000, para = F, level = 0.95, use = c("food.received")
  )
round(boot.reward$ci.estimates, 3)
as.data.frame(round(boot.reward$ci.estimates, 3))
m.stab.plot(round(boot.reward$ci.estimates, 3))
boot.reward$ci.predicted

# Plotting----------------------------------------------------------------------
library(gghalves)
library(ggthemes)
library(cowplot)

which(xdata$searched == "")

xdata.agg <- xdata %>%
  mutate(searched.2 = as.numeric(searched) - 1) %>%
  group_by(name, food.received) %>%
  summarise(mean.resp = mean(searched.2, na.rm = T)) %>%
  ungroup()

xdata.agg$food.received2 <-
  jitter(as.numeric(as.factor(xdata.agg$food.received)), amount = 0.13)
xdata.agg$mean.resp2 <-
  jitter(xdata.agg$mean.resp, amount = 0.04)

exp1_plot_rewards <-
  ggplot(data = xdata.agg, aes(
    x = factor(food.received,
      levels = c("high", "low", "none")
    ),
    y = mean.resp2
  )) +
  geom_line(aes(x = food.received2, y = mean.resp2, group = name),
    color = "gray", lty = 1, alpha = .5
  ) +
  geom_point(
    data = xdata.agg %>% filter(food.received == "high"),
    aes(x = food.received2), color = "chartreuse4", size = 1.5,
    alpha = .6
  ) +
  geom_point(
    data = xdata.agg %>% filter(food.received == "low"),
    aes(x = food.received2), color = "gold3", size = 1.5,
    alpha = .6
  ) +
  geom_point(
    data = xdata.agg %>% filter(food.received == "none"),
    aes(x = food.received2), color = "red3", size = 1.5,
    alpha = .4
  ) +
  geom_violin(
    data = xdata.agg %>% filter(food.received == "high"),
    aes(x = food.received, y = mean.resp),
    position = position_nudge(x = 0.00),
    fill = "chartreuse4", alpha = .2
  ) +
  geom_violin(
    data = xdata.agg %>% filter(food.received == "low"),
    aes(x = food.received, y = mean.resp),
    position = position_nudge(x = 0.00),
    fill = "gold3", alpha = .2
  ) +
  geom_violin(
    data = xdata.agg %>% filter(food.received == "none"),
    aes(x = food.received, y = mean.resp),
    position = position_nudge(x = 0.00),
    fill = "red3", alpha = .2
  ) +
  geom_errorbar(
    data = boot.reward$ci.predicted %>%
      filter(food.received == "high"),
    aes(
      x = as.numeric(food.received) - 0.25,
      y = fitted,
      ymin = lower.cl, ymax = upper.cl
    ),
    color = "chartreuse4", width = 0.1, size = 1
  ) +
  geom_errorbar(
    data = boot.reward$ci.predicted %>%
      filter(food.received == "low"),
    aes(
      x = as.numeric(food.received) + 0.25,
      y = fitted, ymin = lower.cl, ymax = upper.cl
    ),
    color = "gold3", width = 0.1, size = 1
  ) +
  geom_errorbar(
    data = boot.reward$ci.predicted %>%
      filter(food.received == "none"),
    aes(
      x = as.numeric(food.received) + 0.25,
      y = fitted,
      ymin = lower.cl, ymax = upper.cl
    ),
    color = "red3", width = 0.1, size = 1
  ) +
  geom_point(
    data = boot.reward$ci.predicted %>%
      filter(food.received == "high"),
    aes(x = as.numeric(food.received) - 0.25, y = fitted),
    color = "darkgreen", size = 3.5
  ) +
  geom_point(
    data = boot.reward$ci.predicted %>%
      filter(food.received == "low"),
    aes(x = as.numeric(food.received) + 0.25, y = fitted),
    color = "gold4", size = 3.5
  ) +
  geom_point(
    data = boot.reward$ci.predicted %>%
      filter(food.received == "none"),
    aes(x = as.numeric(food.received) + 0.25, y = fitted),
    color = "darkred", size = 3.5
  ) +
  scale_x_discrete(
    limits = c("high", "low", "none"),
    name = "Reward received",
    labels = c(
      "High-value",
      "Low-value",
      "No reward"
    )
  ) +
  scale_y_continuous(
    name = "Post-decisional information search",
    labels = scales::percent
  ) +
  labs(title = "Effect of received reward value") +
  theme_classic() +
  theme(
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    strip.text.x = element_text(size = 13),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.title.y.left = element_text(vjust = 3),
    plot.title = element_text(color = "black", size = 15, face = "bold"),
    axis.title.x = element_text(
      margin =
        margin(t = 10, r = 0, b = 0, l = 0)
    )
  )

exp1_plot_rewards

# Condition Plot----------------------------------------------------------------
ftable(searched ~ name + condition, xdata)

xdata.agg <- xdata %>%
  mutate(searched2 = as.numeric(searched) - 1) %>%
  group_by(name, condition) %>%
  summarise(mean.resp = mean(searched2, na.rm = T)) %>%
  ungroup()

xdata.agg$condition2 <-
  jitter(as.numeric(as.factor(xdata.agg$condition)), amount = 0.07)

exp1_plot_condition <-
  ggplot() +
  geom_violin(
    data = xdata.agg,
    aes(
      x = factor(condition,
        levels = c("control", "choice", "no-choice")
      ),
      y = mean.resp, fill = condition
    ),
    position = position_nudge(x = 0.00), alpha = .2
  ) +
  scale_fill_manual(values = c("dodgerblue", "darkorange", "darkorange")) +
  geom_point(
    data = xdata.agg %>% filter(condition == "control"),
    aes(x = condition2, y = mean.resp), color = "dodgerblue", size = 1.5,
    alpha = .4
  ) +
  geom_point(
    data = xdata.agg %>% filter(condition == "choice"),
    aes(x = condition2, y = mean.resp), color = "darkorange", size = 1.5,
    alpha = .4
  ) +
  geom_point(
    data = xdata.agg %>% filter(condition == "no-choice"),
    aes(x = condition2, y = mean.resp), color = "darkorange", size = 1.5,
    alpha = .6
  ) +
  geom_line(
    data = xdata.agg, aes(x = condition2, y = mean.resp, group = name),
    color = "gray", lty = 1, alpha = .7
  ) +
  geom_errorbar(
    data = boot.plot$ci.predicted %>%
      filter(condition == "control"),
    aes(
      x = as.numeric(condition) - 0.25, y = fitted,
      ymin = lower.cl, ymax = upper.cl
    ), color = "dodgerblue",
    width = 0.1, size = 1
  ) +
  geom_errorbar(
    data = boot.plot$ci.predicted %>%
      filter(condition == "choice"),
    aes(
      x = as.numeric(condition) + 0.25, y = fitted,
      ymin = lower.cl, ymax = upper.cl
    ), color = "darkorange",
    width = 0.1, size = 1
  ) +
  geom_errorbar(
    data = boot.plot$ci.predicted %>%
      filter(condition == "no-choice"),
    aes(
      x = as.numeric(condition) + 0.25, y = fitted,
      ymin = lower.cl, ymax = upper.cl
    ), color = "darkorange",
    width = 0.1, size = 1
  ) +
  geom_point(
    data = boot.plot$ci.predicted %>% filter(condition == "control"),
    aes(x = as.numeric(condition) - 0.25, y = fitted),
    color = "dodgerblue4", size = 3.5
  ) +
  geom_point(
    data = boot.plot$ci.predicted %>% filter(condition == "choice"),
    aes(x = as.numeric(condition) + 0.25, y = fitted),
    color = "darkorange4", size = 3.5
  ) +
  geom_point(
    data = boot.plot$ci.predicted %>% filter(condition == "no-choice"),
    aes(x = as.numeric(condition) + 0.25, y = fitted),
    color = "darkorange4", size = 3.5
  ) +
  scale_x_discrete(
    limits = c("control", "choice", "no-choice"),
    name = "Condition",
    labels = c("Control", "Choice", "No-choice")
  ) +
  scale_y_continuous(
    name = "Post-decisional information search",
    labels = scales::percent
  ) +
  labs(title = "Effect of condition") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    strip.text.x = element_text(size = 13),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.title.y.left = element_text(vjust = 3),
    plot.title = element_text(color = "black", size = 15, face = "bold"),
    axis.title.x = element_text(
      margin =
        margin(t = 10, r = 0, b = 0, l = 0)
    )
  )
exp1_plot_condition

# Combine plots
library(ggpubr)
# theme_set(theme_pubr())
figure <- ggarrange(exp1_plot_condition, exp1_plot_rewards,
  labels = c("(a)", "(b)"),
  hjust = -1, # hjust = 0 for left alignment
  vjust = 2.5, # hjust = 0 for left alignment
  ncol = 2, nrow = 1
)
figure



# Outcome effect without control condition:  ------------------------------

without_contr <- subset(t.data, t.data$condition == "choice" | t.data$condition == "no-choice")

full <-
  glmer(
    searched ~ condition * food.received +
      (1 + (condition.choice.code) +
        (food.received.low.code + food.received.none.code) | name),
    data = without_contr, control = contr, family = binomial(link = "logit")
  )

main <-
  glmer(
    searched ~ condition + food.received +
      (1 + (condition.choice.code) +
        (food.received.low.code + food.received.none.code) | name),
    data = without_contr, control = contr, family = binomial(link = "logit")
  )

null <-
  glmer(
    searched ~ condition + food.received +
      (1 + (condition.choice.code) +
        (food.received.low.code + food.received.none.code) | name),
    data = without_contr, control = contr, family = binomial(link = "logit")
  )

anova(full, null, test = "Chisq")
drop1(full, test = "Chisq")
library(effects)

plot(effect("condition:food.received", full))

save.image("./study1/R_objects/analysis_1.RData")
