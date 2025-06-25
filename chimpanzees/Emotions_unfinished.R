# Analysis Study 1 - Chimpanzees

# Load packages -------------------------------------------------------------
library("lme4")
library("readxl")
library("tidyverse")
library("tidyselect")
library("parallel")
library("optimx")
library("emmeans")
source("./study1/functions/diagnostic_fcns.r")
source("./study1/functions/glmm_stability.r")
source("./study1/functions/drop1_para.r")
source("./study1/functions/boot_glmm.r")
library("emmeans")
library("car")

# Load data -------------------------------------------------------------
xdata <-
  read.csv("./study1/data/Counterfactual_Curiosity_Study_1_Chimps.csv",
            header = TRUE, na = c("NA")
  )

xdata$emo.plus.turn <- NA
xdata$emo.plus.turn[xdata$searched == "yes" &
                      xdata$time.from.last.gaze.to.back.of.head.to.camera >= 2] <- "yes"
xdata$emo.plus.turn[xdata$searched == "yes" &
                      xdata$turned.away == "no"] <- "no"
xdata$emo.plus.turn[xdata$searched == "yes" &
                      xdata$time.from.last.gaze.to.back.of.head.to.camera < 2] <- "no"
xdata$emo.plus.turn[xdata$emotional.displays == "yes"] <- "yes"

# Analyse emotional displays and turn away within 2 sec ----------------------------------------------------------------
xx.fe.re <- fe.re.tab(
  fe.model = "emo.plus.turn ~ value.remaining.seen*condition",
  re = "(1|name)",
  other.vars = c(
    "age", "gender",
    "trial.per.condition", "trial.total"
  ),
  data = xdata
) # maybe add age
xx.fe.re$summary
t.data <- xx.fe.re$data
str(t.data)

## Center dummy variables
t.data$value.remaining.seen.better.code <-
  t.data$value.remaining.seen.better -
  mean(t.data$value.remaining.seen.better)
t.data$value.remaining.seen.same.code <-
  t.data$value.remaining.seen.same -
  mean(t.data$value.remaining.seen.same)

t.data$condition.no.choice.code <-
  t.data$condition.no.choice -
  mean(t.data$condition.no.choice)
t.data$condition.choice.code <-
  t.data$condition.choice -
  mean(t.data$condition.choice)

t.data$z.age <-
  scale(t.data$age)
t.data$z.trial.per.condition <-
  scale(t.data$trial.per.condition)
t.data$z.trial.total <-
  scale(t.data$trial.total)

## Fitting the model as pre-registered
contr <-
  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000000))
full.emo <-
  glmer(emo.plus.turn ~ value.remaining.seen * condition +
          (1 + (value.remaining.seen.better.code +
                  value.remaining.seen.same.code) +
             (condition.no.choice.code +
                condition.choice.code) | name),
        data = t.data, control = contr,
        family = binomial(link = "logit")
  )

summary(full)$varcor
round(summary(full)$coefficients, 3)

## Fit main effects model
main.emo <-
  glmer(emo.plus.turn ~
          value.remaining.seen +
          condition +
          (1 + (value.remaining.seen.better.code +
                  value.remaining.seen.same.code) +
             (condition.no.choice.code +
                condition.choice.code) | name),
        data = t.data, control = contr, family = binomial(link = "logit")
  )
## Reduced model comparisons
tests.full <- drop1p(
  model.res = full.emo, para = F, data = NULL,
  contr = contr, n.cores = c("all-1", "all"), to.del = NULL
)
round(tests.full$drop1.res, 3)
tests.main <- drop1p(
  model.res = main.emo, para = F, data = NULL,
  contr = contr, n.cores = c("all-1", "all"), to.del = NULL
)
round(tests.main$drop1.res, 3)




