# Analysis of Countefactual Curiosity Study 1
# Authors: Hanna Schleihauf, Nika Ghavami, Lou Haux
# Julia Fischer, Esther Herrmann, Jan Engelmann
# Date: June 8th, 2023

# Loading packages and functions
library("lme4")
library("readxl")
library("tidyverse")
library("tidyselect")
library("parallel")
library("optimx")
library("emmeans")
library("emmeans")
library("car")

source("./study1/functions/diagnostic_fcns.r")
source("./study1/functions/glmm_stability.r")
source("./study1/functions/drop1_para.r")
source("./study1/functions/boot_glmm.r")


# load("./study1/images/Analysis_1.RData")

# Load data
xdata <-
  read_xlsx("./study1/data/Coding_Sheets_Study_1_Kids_06082023.xlsx",
    col_names = TRUE, na = c("NA")
  )

str(xdata)

levels(as.factor(xdata$condition))
xdata <- subset(xdata, xdata$condition.trial == "choice" |
  xdata$condition.trial == "no.choice" |
  xdata$condition.trial == "control")
xdata$condition <- droplevels(xdata$condition)

## Check data
xdata$searched.after <-
  as.factor(xdata$searched.after)
ftable(searched.after ~
  condition + reward.received, xdata)
# xdata$condition =
#  factor(xdata$condition,
#                levels = c("control", "choice", "no.choice"),
#                labels = c("control", "choice", "no.choice"))
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

## Prepare data for model fitting------------------------------------
xx.fe.re <- fe.re.tab(
  fe.model = "searched.after ~ age*condition*reward.received +
                   gender + trial.per.child",
  re = "(1|id)", other.vars = c("age.group"),
  data = xdata
)
xx.fe.re$summary
t.data <- xx.fe.re$data
str(t.data)

## Center dummy variables
## (necessary for the random effects in the model and potentially plotting)
t.data$reward.received.low.code <-
  t.data$reward.received.low - mean(t.data$reward.received.low)
t.data$reward.received.none.code <-
  t.data$reward.received.none - mean(t.data$reward.received.none)
t.data$condition.no.choice.code <-
  t.data$condition.no.choice - mean(t.data$condition.no.choice)
t.data$condition.choice.code <-
  t.data$condition.choice - mean(t.data$condition.choice)
t.data$z.age <- scale(t.data$age)
t.data$z.trial <- scale(t.data$trial.per.child)
t.data$age.group <- as.factor(t.data$age.group)

## Fitting the model (age had been forgotten in the pre-registration)
contr <-
  glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 10000000)
  )
full <-
  glmer(searched.after ~
    reward.received * condition * z.age + gender +
    (1 + (reward.received.low.code +
      reward.received.none.code) || id),
  data = t.data, control = contr,
  family = binomial(link = "logit")
  )

red <-
  glmer(searched.after ~
    (reward.received + condition + z.age)^2 + gender +
    (1 + (reward.received.low.code +
      reward.received.none.code) || id),
  data = t.data, control = contr,
  family = binomial(link = "logit")
  )

red2 <-
  glmer(searched.after ~
    reward.received + condition * z.age + gender +
    (1 + (reward.received.low.code +
      reward.received.none.code) || id),
  data = t.data, control = contr,
  family = binomial(link = "logit")
  )

main <-
  glmer(searched.after ~
    reward.received + condition + z.age + gender +
    (1 + (reward.received.low.code +
      reward.received.none.code) || id),
  data = t.data, control = contr,
  family = binomial(link = "logit")
  )

null <-
  glmer(searched.after ~ 1 +
    (1 + (reward.received.low.code +
      reward.received.none.code) || id),
  data = t.data, control = contr,
  family = binomial(link = "logit")
  )

summary(full)$varcor
round(summary(full)$coefficients, 3)

ftable(searched.after ~
  condition + reward.received, t.data)
ftable(searched.after ~
  condition + age.group, t.data)
ftable(searched.after ~
  reward.received, t.data)

## Test assumptions
overdisp.test(full) # no overdispersion
vif(main) # no colliniarity
diagnostics.plot(full)
ranef.diagn.plot(full)

## Checking model stability
m.stab.b <-
  glmm.model.stab(model.res = full, contr = contr, use = c("id"))
m.stab.b$detailed$warnings
xx <- as.data.frame(round(m.stab.b$summary[, -1], 3))
dev.off()
m.stab.plot(round(m.stab.b$summary[, -1], 3))
write.table(xx, "full.model.stab.txt", quote = FALSE, sep = "\t")
xx

## Full-null models comparison
round(anova(full, null, test = "Chisq"), 3)

## Reduced model comparisons
tests.full <- drop1p(
  model.res = full, para = F, data = NULL,
  contr = contr, n.cores = c("all-1", "all"), to.del = NULL
)
round(tests.full$drop1.res, 3)
tests.red <- drop1p(
  model.res = red, para = F, data = NULL,
  contr = contr, n.cores = c("all-1", "all"), to.del = NULL
)
round(tests.red$drop1.res, 3)
tests.red2 <- drop1p(
  model.res = red2, para = F, data = NULL,
  contr = contr, n.cores = c("all-1", "all"), to.del = NULL
)
round(tests.red2$drop1.res, 3)
tests.main <- drop1p(
  model.res = main, para = F, data = NULL,
  contr = contr, n.cores = c("all-1", "all"), to.del = NULL
)
round(tests.main$drop1.res, 3)

## First peek at effects
library("effects")
plot(effect("condition:z.age", red2))
plot(effect(c("reward.received"), red2))
plot(effect(c("gender"), red2))

## Pairwise comparisons
emm1 <- emmeans(red2, ~condition)
summary(emm1, type = "response")
emmeans(red2, pairwise ~ condition)
summary(pairs(emm1), type = "response")
summary(pairs(regrid(emm1)), type = "response")

emm2 <- emmeans(red2, ~reward.received)
summary(emm2, type = "response")
emmeans(red2, pairwise ~ reward.received)
summary(pairs(emm2), type = "response")
summary(pairs(regrid(emm2)), type = "response")


## Bootstraps of full model
# The bootstrap has already been run and is saved in the image
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
  glmer(searched.after ~ condition * z.age +
    (reward.received.low.code + reward.received.none.code) +
    (1 | id),
  data = t.data, control = contr,
  family = binomial(link = "logit")
  )
# The bootstrap has already been run and is saved in the image
boot.condition.age <- boot.glmm.pred(
  model.res = model.condition.age,
  excl.warnings = T,
  nboots = 1000, para = F, level = 0.95,
  use = c("condition", "z.age")
)
round(boot.condition.age$ci.estimates, 3)
as.data.frame(round(boot.condition.age$ci.estimates, 3))
m.stab.plot(round(boot.condition.age$ci.estimates, 3))
boot.condition.age$ci.predicted

model.reward <-
  glmer(searched.after ~ reward.received +
    z.age * (condition.no.choice.code + condition.choice.code) +
    (1 | id),
  data = t.data, control = contr,
  family = binomial(link = "logit")
  )
# The bootstrap has already been run and is saved in the image
boot.reward <-
  boot.glmm.pred(
    model.res = model.reward, excl.warnings = T,
    nboots = 1000, para = F,
    level = 0.95, use = c("reward.received")
  )
round(boot.reward$ci.estimates, 3)
as.data.frame(round(boot.reward$ci.estimates, 3))
m.stab.plot(round(boot.reward$ci.estimates, 3))
boot.reward$ci.predicted



# Plotting Reward model----------------------------------------------------------------
library(gghalves)
library(ggthemes)
library(cowplot)

xdata.agg <- xdata %>%
  mutate(searched.after.numeric = as.numeric(searched.after) - 1) %>%
  group_by(id, reward.received) %>%
  summarise(mean.resp = mean(searched.after.numeric, na.rm = T)) %>%
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
    color = "chartreuse4", size = 2.5
  ) +
  geom_point(
    data = boot.reward$ci.predicted %>%
      filter(reward.received == "low"),
    aes(x = as.numeric(reward.received) + 0.25, y = fitted),
    color = "gold3", size = 2.5
  ) +
  geom_point(
    data = boot.reward$ci.predicted %>%
      filter(reward.received == "none"),
    aes(x = as.numeric(reward.received) + 0.25, y = fitted),
    color = "red3", size = 2.5
  ) +
  scale_x_discrete(
    limits = c("high", "low", "none"),
    name = "Reward received",
    labels = c("High-value reward", "Low-value reward", "No reward")
  ) +
  scale_y_continuous(
    name = "Proportion of trials with information-search",
    labels = scales::percent
  ) +
  labs(title = "Effect of received reward value") +
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


# Plotting Condition model----------------------------------------------------------------
xdata$z.age <- scale(xdata$age)

xdata.agg.condition <- xdata %>%
  mutate(searched.after.numeric = as.numeric(searched.after) - 1) %>%
  group_by(id, z.age, condition) %>%
  summarise(mean.resp = mean(searched.after.numeric, na.rm = T)) %>%
  ungroup()

xdata.agg.condition$mean.resp2 <-
  jitter(xdata.agg.condition$mean.resp, amount = 0.02)

exp1_plot_condition <-
  ggplot() +
  geom_point(
    data = xdata.agg.condition,
    aes(x = z.age, y = mean.resp2, color = condition),
    size = 2.5, alpha = .4
  ) +
  scale_color_manual(values = c("dodgerblue", "darkorange", "darkorange")) +
  geom_ribbon(
    data = boot.condition.age$ci.predicted,
    aes(
      ymin = lower.cl, ymax = upper.cl, x = z.age,
      group = condition,
      fill = condition
    ),
    alpha = 0.3, color = "black"
  ) +
  scale_fill_manual(values = c("dodgerblue", "darkorange", "darkorange")) +
  geom_line(data = boot.condition.age$ci.predicted, aes(x = z.age, y = fitted)) +
  facet_wrap(~condition, switch = "x") +
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
    name = "Proportion of trials with information-search",
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

# save.image("./study1/images/Analysis_1.RData")
