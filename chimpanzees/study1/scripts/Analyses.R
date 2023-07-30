library("lme4")
library("readxl")
library("tidyverse")
library("tidyselect")
library("parallel")
library("optimx")
library("emmeans")
source("./functions/diagnostic_fcns.r")
source("./functions/glmm_stability.r")
source("./functions/drop1_para.r")
source("./functions/boot_glmm.r")
library("emmeans")
library("car")

getwd()

xdata <- 
  read_xlsx("/Users/hanna.schleihauf/Dropbox/Research/my_projects/Counterfactual_Curiousity_Chimps/Counterfactual_Curiosity_Study_1.xlsx", 
            col_names = TRUE, na= c("NA"))

## Check data 
str(xdata)
xdata$searched = as.factor(xdata$searched)
ftable(searched ~ condition + food.received, xdata)
xdata$condition = factor(xdata$condition, levels = c("control", "choice", "no-choice"))
xdata$value.remaining.seen = as.factor(xdata$value.remaining.seen)
xdata$value.remaining.seen = relevel(xdata$value.remaining.seen, ref = "worse")

levels(xdata$value.remaining.seen)


xdata$turned.away.within.2 = NA
xdata$turned.away.within.2[xdata$searched == 
                             "yes" & 
                             xdata$time.from.last.gaze.to.back.of.head.to.camera >= 2] = "yes"
xdata$turned.away.within.2[xdata$searched == 
                             "yes" & 
                             xdata$turned.away == "no"] = "no"
xdata$turned.away.within.2[xdata$searched == 
                             "yes" & 
                             xdata$time.from.last.gaze.to.back.of.head.to.camera < 2] = "no"

xdata$emo.plus.turn = NA 
xdata$emo.plus.turn[xdata$searched == 
                      "yes" & 
                      xdata$time.from.last.gaze.to.back.of.head.to.camera >= 2] = "yes"
xdata$emo.plus.turn[xdata$searched == 
                      "yes" & 
                      xdata$turned.away == "no"] = "no"
xdata$emo.plus.turn[xdata$searched == 
                      "yes" & 
                      xdata$time.from.last.gaze.to.back.of.head.to.camera < 2] = "no"
xdata$emo.plus.turn[xdata$emotional.displays == "yes"] = "yes"


xdata <- xdata[
  order( xdata$name, xdata$day, xdata$trial ),
]
xdata$trial.per.condition = 
  as.numeric(ave(xdata$name, 
                 list(xdata$name, xdata$condition), 
                 FUN=seq_along))
xdata$trial.total = 
  as.numeric(ave(xdata$name, 
                 list(xdata$name), 
                 FUN=seq_along))

ftable(searched ~ turned.away.within.2 , xdata)
ftable(searched ~ name + condition , xdata)

## Prepare data for model fitting------------------------------------
xx.fe.re=fe.re.tab(fe.model="searched ~ 
                   condition*food.received + trial.total",
                   re="(1|name)", 
                   other.vars= c("age", "gender", "trial.per.condition"), 
                   data=xdata)  #maybe add age
xx.fe.re$summary
t.data=xx.fe.re$data 
str(t.data) 

## Center dummy variables (necessary for the random effects in the model and potentially plotting)
t.data$food.received.low.code = 
  t.data$food.received.low - mean(t.data$food.received.low)
t.data$food.received.none.code = 
  t.data$food.received.none - mean(t.data$food.received.none)
t.data$condition.no.choice.code = 
  t.data$condition.no.choice - mean(t.data$condition.no.choice)
t.data$condition.choice.code = 
  t.data$condition.choice - mean(t.data$condition.choice)

t.data$z.age = scale(t.data$age)
t.data$z.trial.per.condition = scale(t.data$trial.per.condition)
t.data$z.trial.total = scale(t.data$trial.total)

## Fitting the model as pre-registered
contr <- 
  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000000))
full <- 
  glmer(searched ~ (food.received+condition)^2 +  
                (1 + (condition.no.choice.code + condition.choice.code) +
                   (food.received.low.code + food.received.none.code) | name), 
              data = t.data, control = contr, family = binomial(link = "logit"))

summary(full)$varcor
round(summary(full)$coefficients, 3)

ftable(searched ~ condition + food.received, t.data)
xx = tapply((as.numeric(t.data$searched)-1), 
            list(t.data$condition, t.data$food.received), mean)

## Test assumptions
overdisp.test(full) #no overdispersion
xx=glmer(searched ~ condition + food.received +  
           (1 + (condition.no.choice.code + condition.choice.code) +
              (food.received.low.code + food.received.none.code) | name), 
         data = t.data, control = contr, family = binomial(link = "logit"))
vif(xx) #no colliniarity

## Checking model stability
m.stab.b=glmm.model.stab(model.res=full, contr=contr, use=c("name")) 
m.stab.b$detailed$warnings
xx=as.data.frame(round(m.stab.b$summary[, -1], 3))
dev.off()
m.stab.plot(round(m.stab.b$summary[, -1], 3))
write.table(xx, "full.model.stab.txt",quote=FALSE, sep = "\t")
xx

## Fit null model
null <- 
  glmer(searched ~ 1 +  
                (1 + (condition.no.choice.code + condition.choice.code) + 
                   (food.received.low.code + food.received.none.code)  | name), 
              data = t.data, control = contr, family = binomial(link = "logit"))
round(anova(full, null, test="Chisq"), 3)
drop1
## Fit main effects model
main <- 
  glmer(searched ~ condition + food.received +
                (1 + (condition.no.choice.code + condition.choice.code) + 
                   (food.received.low.code + food.received.none.code) | name), 
              data = t.data, control = contr, family = binomial(link = "logit"))

## Reduced model comparisons
tests.full=drop1p(model.res=full, 
                  para=F, data=NULL, contr=contr, 
                  n.cores=c("all-1", "all"), to.del=NULL)
round(tests.full$drop1.res, 3) 
tests.main=drop1p(model.res=main, 
                  para=F, data=NULL, contr=contr, 
                  n.cores=c("all-1", "all"), to.del=NULL)
round(tests.main$drop1.res, 3)

## First peek at effects
library("effects")
plot(effect("food.received:condition", full))
plot(effect( c("condition") , main))


## Pairwise comparisons
emm1 = emmeans(full,  ~ food.received*condition)
summary(emm1, type = "response")
emmeans(full, pairwise ~ condition*food.received)
xx = summary(pairs(emm1), type = "response")
summary(pairs(regrid(emm1)), type = "response")

emm2 = emmeans(main,  ~ condition)
summary(emm2, type = "response")
emmeans(main, pairwise ~ condition)
xx = summary(pairs(emm2), type = "response")
summary(pairs(regrid(emm2)), type = "response")


## Bootstraps of full model
boot.full = 
  boot.glmm.pred(model.res=full, 
                 excl.warnings=T,
                 nboots=1000, para=F, 
                 level=0.95, 
                 use=c("condition", "food.received"))
round(boot.full$ci.estimates, 3)
as.data.frame(round(boot.full$ci.estimates, 3))
m.stab.plot(round(boot.full$ci.estimates, 3))
boot.full$ci.predicted

## Bootstraps of plot.model
plot.model <- 
  glmer(searched ~ condition + 
        (food.received.low.code + food.received.none.code) +  
         (1 + (condition.no.choice.code + condition.choice.code) +
                         (food.received.low.code + food.received.none.code) | name), 
        data = t.data, 
        control = contr, 
        family = binomial(link = "logit"))

boot.plot=boot.glmm.pred(model.res=plot.model, excl.warnings=T,
                         nboots=1000, para=F, level=0.95, use=c("condition"))
round(boot.plot$ci.estimates, 3)
as.data.frame(round(boot.plot$ci.estimates, 3))
m.stab.plot(round(boot.plot$ci.estimates, 3))
boot.plot$ci.predicted

model.reward <- 
  glmer(searched ~ 
          food.received + 
          (condition.no.choice.code + condition.choice.code) +
          (1 + (condition.no.choice.code + condition.choice.code) +
          (food.received.low.code + food.received.none.code) | name), 
          data = t.data, control = contr, 
        family = binomial(link = "logit"))
boot.reward = 
  boot.glmm.pred(model.res=model.reward, excl.warnings=T,
                           nboots=1000, para=F, level=0.95, use=c("food.received"))
round(boot.reward$ci.estimates, 3)
as.data.frame(round(boot.reward$ci.estimates, 3))
m.stab.plot(round(boot.reward$ci.estimates, 3))
boot.reward$ci.predicted

ftable(searched ~ turned.away.within.2, xdata)

# Analyse emotional displays and turn away within 2 sec ----------------------------------------------------------------
xx.fe.re=fe.re.tab(fe.model="emo.plus.turn ~ value.remaining.seen*condition",
                   re="(1|name)", 
                   other.vars= c("age", "gender", 
                                 "trial.per.condition", "trial.total"), 
                   data=xdata)  #maybe add age
xx.fe.re$summary
t.data=xx.fe.re$data 
str(t.data) 

## Center dummy variables
t.data$value.remaining.seen.better.code = 
  t.data$value.remaining.seen.better - 
  mean(t.data$value.remaining.seen.better)
t.data$value.remaining.seen.same.code = 
  t.data$value.remaining.seen.same - 
  mean(t.data$value.remaining.seen.same)

t.data$condition.no.choice.code = 
  t.data$condition.no.choice - 
  mean(t.data$condition.no.choice)
t.data$condition.choice.code = 
  t.data$condition.choice - 
  mean(t.data$condition.choice)

t.data$z.age = 
  scale(t.data$age)
t.data$z.trial.per.condition = 
  scale(t.data$trial.per.condition)
t.data$z.trial.total = 
  scale(t.data$trial.total)

## Fitting the model as pre-registered
contr <- 
  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000000))
full.emo <- 
  glmer(emo.plus.turn ~ value.remaining.seen *condition +  
                    (1 + (value.remaining.seen.better.code + 
                            value.remaining.seen.same.code) +
                       (condition.no.choice.code + 
                          condition.choice.code)  | name), 
                  data = t.data, control = contr, 
        family = binomial(link = "logit"))

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
          condition.choice.code)  | name), 
       data = t.data, control = contr, family = binomial(link = "logit"))
## Reduced model comparisons
tests.full=drop1p(model.res=full.emo, para=F, data=NULL, 
                  contr=contr, n.cores=c("all-1", "all"), to.del=NULL)
round(tests.full$drop1.res, 3) 
tests.main=drop1p(model.res=main.emo, para=F, data=NULL, 
                  contr=contr, n.cores=c("all-1", "all"), to.del=NULL)
round(tests.main$drop1.res, 3)


# Plotting----------------------------------------------------------------------
library(gghalves)
library(ggthemes)
library(cowplot)

xdata.agg <- xdata %>%
  mutate(searched.2=as.numeric(searched)-1)%>%
  group_by(name, food.received) %>%
  summarise(mean.resp = mean(searched.2, na.rm = T)) %>%
  ungroup()

xdata.agg$food.received2 <- 
  jitter(as.numeric(as.factor(xdata.agg$food.received)), amount = 0.13)
xdata.agg$mean.resp2 <- 
  jitter(xdata.agg$mean.resp, amount = 0.04)

exp1_plot_rewards <- 
  ggplot(data = xdata.agg, aes(x = factor(food.received, 
                                          levels=c("high", "low", "none")), 
                               y = mean.resp2)) +
  
  geom_line(aes(x = food.received2, y = mean.resp2, group = name), 
            color = "gray", lty = 1, alpha = .5) +
  
  geom_point(data = xdata.agg %>% filter(food.received == "high"), 
             aes(x = food.received2), color = "chartreuse4", size = 1.5, 
             alpha = .6) +
  geom_point(data = xdata.agg %>% filter(food.received == "low"), 
             aes(x = food.received2), color = "gold3", size = 1.5, 
             alpha = .6) +
  geom_point(data = xdata.agg %>% filter(food.received == "none"), 
             aes(x = food.received2), color = "red3", size = 1.5, 
             alpha = .4) +
  
  geom_violin(data = xdata.agg %>% filter(food.received == "high"), 
                    aes(x = food.received, y = mean.resp), 
              position = position_nudge(x = 0.00),  
                    fill = "chartreuse4", alpha = .2) +
  geom_violin(data = xdata.agg %>% filter(food.received == "low"), 
                    aes(x = food.received, y = mean.resp), 
              position = position_nudge(x = 0.00), 
                    fill = "gold3", alpha = .2) +
  geom_violin(data = xdata.agg %>% filter(food.received == "none"), 
                    aes(x = food.received, y = mean.resp), 
              position = position_nudge(x = 0.00), 
                    fill = "red3", alpha = .2) + 
  
  geom_errorbar(data = boot.reward$ci.predicted %>% 
                  filter(food.received == "high"), 
                aes(x = as.numeric(food.received) - 0.25, 
                    y = fitted, 
                    ymin = lower.cl, ymax = upper.cl), 
                    color = "chartreuse4", width = 0.1, size = 1) +
  geom_errorbar(data = boot.reward$ci.predicted %>% 
                  filter(food.received == "low"), 
                aes(x = as.numeric(food.received) + 0.25, 
                    y = fitted, ymin = lower.cl, ymax = upper.cl), 
                    color = "gold3", width = 0.1, size = 1) +
  geom_errorbar(data = boot.reward$ci.predicted %>% 
                  filter(food.received == "none"), 
                aes(x = as.numeric(food.received) + 0.25, 
                    y = fitted, 
                    ymin = lower.cl, ymax = upper.cl), 
                    color = "red3", width = 0.1, size = 1) +
  
  geom_point(data = boot.reward$ci.predicted %>% 
               filter(food.received == "high"), 
             aes(x = as.numeric(food.received) - 0.25, y = fitted), 
             color = "chartreuse4", size = 2.5) +
  geom_point(data = boot.reward$ci.predicted %>% 
               filter(food.received == "low"), 
             aes(x = as.numeric(food.received) + 0.25, y = fitted), 
             color = "gold3", size = 2.5) +
  geom_point(data = boot.reward$ci.predicted %>% 
               filter(food.received == "none"), 
             aes(x = as.numeric(food.received) + 0.25, y = fitted), 
             color = "red3", size = 2.5) +
  
  scale_x_discrete(limits = c("high", "low", "none"), 
                   name = "Reward received", 
                   labels = c("High-value reward", 
                              "Low-value reward", 
                              "No reward")) +
  scale_y_continuous(name = "Proportion of trials with information-search", 
                     labels = scales::percent) +
  
  labs(title="Effect of received reward value") +
  theme_classic() + 
  theme(axis.ticks.x = element_blank(),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),
        strip.text.x = element_text(size = 13),
        plot.margin=unit(c(1,1,1,1), "cm"),
        axis.title.y.left = element_text(vjust=3),
        plot.title = element_text(color="black", size=15, face="bold"), 
        axis.title.x = element_text(margin = 
                                      margin(t = 10, r = 0, b = 0, l = 0)))

exp1_plot_rewards



# Condition Plot----------------------------------------------------------------
ftable(searched ~ name + condition, xdata)

xdata.agg <- xdata %>%
  mutate(searched2=as.numeric(searched)-1)%>%
  group_by(name, condition) %>%
  summarise(mean.resp = mean(searched2, na.rm = T)) %>%
  ungroup()

xdata.agg$condition2 <- 
  jitter(as.numeric(as.factor(xdata.agg$condition)), amount = 0.07)

exp1_plot_condition <- 
  ggplot(data = xdata.agg, aes(x = factor(condition, 
                               levels=c("control", "choice", "no-choice")), 
                               y = mean.resp)) +
  
  geom_point(data = xdata.agg %>% filter(condition == "control"), 
             aes(x = condition2), color = "dodgerblue", size = 1.5, 
             alpha = .4) +
  geom_point(data = xdata.agg %>% filter(condition == "choice"), 
             aes(x = condition2), color = "darkorange", size = 1.5, 
             alpha = .4) +
  geom_point(data = xdata.agg %>% filter(condition == "no-choice"), 
             aes(x = condition2), color = "darkorange", size = 1.5, 
             alpha = .6) +
  
  geom_line(aes(x = condition2, y = mean.resp, group = name), 
            color = "gray", lty = 1, alpha = .7) +
  
  geom_violin(data = xdata.agg %>% filter(condition == "control"), 
                    aes(x = condition, y = mean.resp), 
                    position = position_nudge(x = 0.00),  
                    fill = "dodgerblue", alpha = .2) +
  geom_violin(data = xdata.agg %>% filter(condition == "choice"), 
                    aes(x = condition, y = mean.resp), 
                    position = position_nudge(x = 0.00), 
                    fill = "darkorange", alpha = .2) +
  geom_violin(data = xdata.agg %>% filter(condition == "no-choice"), 
                    aes(x = condition, y = mean.resp), 
                    position = position_nudge(x = 0.00), 
                    fill = "darkorange", alpha = .2) + 
  
  geom_errorbar(data = boot.plot$ci.predicted %>% 
                  filter(condition == "control"), 
                aes(x = as.numeric(condition) - 0.25, y = fitted, 
                    ymin = lower.cl, ymax = upper.cl), color = "dodgerblue", 
                width = 0.1, size = 1) +
  geom_errorbar(data = boot.plot$ci.predicted %>% 
                  filter(condition == "choice"), 
                aes(x = as.numeric(condition) + 0.25, y = fitted, 
                    ymin = lower.cl, ymax = upper.cl), color = "darkorange", 
                width = 0.1, size = 1) +
  geom_errorbar(data = boot.plot$ci.predicted %>% 
                  filter(condition == "no-choice"), 
                aes(x = as.numeric(condition) + 0.25, y = fitted, 
                    ymin = lower.cl, ymax = upper.cl), color = "darkorange", 
                width = 0.1, size = 1) +
  
  geom_point(data = boot.plot$ci.predicted %>% filter(condition == "control"), 
             aes(x = as.numeric(condition) - 0.25, y = fitted), 
             color = "dodgerblue", size = 2.5) +
  geom_point(data = boot.plot$ci.predicted %>% filter(condition == "choice"), 
             aes(x = as.numeric(condition) + 0.25, y = fitted), 
             color = "darkorange", size = 2.5) +
  geom_point(data = boot.plot$ci.predicted %>% filter(condition == "no-choice"), 
             aes(x = as.numeric(condition) + 0.25, y = fitted), 
             color = "darkorange", size = 2.5) +
  
  scale_x_discrete(limits = c("control", "choice", "no-choice"), 
                   name = "Condition", 
                   labels = c("Control", "Choice", "No-choice")) +
  scale_y_continuous(name = "Proportion of trials with information-search", 
                     labels = scales::percent) +
  
  labs(title="Effect of condition") +
  theme_classic() + 
  theme(axis.ticks.x = element_blank(),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),
        strip.text.x = element_text(size = 13),
        plot.margin=unit(c(1,1,1,1), "cm"),
        axis.title.y.left = element_text(vjust=3),
        plot.title = element_text(color="black", size=15, face="bold"), 
        axis.title.x = element_text(margin = 
                                      margin(t = 10, r = 0, b = 0, l = 0)))

exp1_plot_condition



save.image("analysis.RData")
