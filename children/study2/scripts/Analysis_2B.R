# Analysis Study 2 - Children

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
source("./study1/functions/diagnostic_fcns.r")
source("./study1/functions/glmm_stability.r")
source("./study1/functions/drop1_para.r")
source("./study1/functions/boot_glmm.r")
options(scipen = 9999)

# Load Data -----------------------------------------------------------------

# To load all the objects, including model results and bootstraps, you can run the follwing line:
# load("./study2/R_objects/Analysis_2.RData")

xdata <-
  read.csv("./study2/data/data_study2B_children.csv",
           header = TRUE, na = c("NA", "")
  )

xdata <-
  xdata %>%
  filter(trial == 4 | trial == 5 | trial == 6 | trial == 7)
xdata$searched.in.targeted.box <- as.factor(xdata$searched.in.targeted.box)

## Check data
xdata$searched.in.targeted.box

## Participant information ------------------------------------------
table(xdata$age.group)/4
table(xdata$gender)/4
tapply(xdata$age, xdata$age.group, mean)
tapply(xdata$age, xdata$age.group, min)
tapply(xdata$age, xdata$age.group, max)

mean(xdata$age)
sd(xdata$age)

round((table(xdata$ethnic.info, useNA = "always")/4)/48, 2)
xx <- rbind(xdata$education.parent1, xdata$education.parent2)
round((table(xx, useNA = "always")/8)/48, 2)
round((table(xdata$family.income, useNA = "always")/4)/48, 2)

# Interrater Reliability ------------------------------------------------
interrater <-
  xdata %>%
  filter(!is.na(searched.in.targeted.box.reliability))
# percentage of reliability coding 25%
nrow(interrater) / nrow(xdata)
interrater$searched.in.targeted.box.reliability <-
  as.factor(interrater$searched.in.targeted.box.reliability)

# number of dyads for which the reliability coding was not the same
sum(as.numeric(interrater$searched.in.targeted.box.reliability) !=
   (as.numeric(interrater$searched.in.targeted.box)))
(sum(as.numeric(interrater$searched.in.targeted.box.reliability) ==
     as.numeric(interrater$searched.in.targeted.box))) / nrow(interrater)

# Cohen's Kappa
library(irr)
kappa_results <-
  kappa2(cbind(
    as.factor(interrater$searched.in.targeted.box.reliability),
    interrater$searched.in.targeted.box
  ), "unweighted")
print(kappa_results)

# FITTING THE MODEL -----------------------------------------------------------------
xx.fe.re=fe.re.tab(fe.model="searched.in.targeted.box ~ trial",
                   re="(1|nr)", data=xdata)  #maybe add age
xx.fe.re$summary
t.data=xx.fe.re$data
str(t.data)
t.data$z.trial = scale(t.data$trial)

contr <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000000))
chance.test_exp2 <- glmer(searched.in.targeted.box ~ 1 +
                            (1 + z.trial | nr),
                          data = t.data, family = binomial, control = contr
)


# CHECKING MODEL RESULTS -----------------------------------------------------------------
summary(chance.test_exp2)
plogis(0.5022)   # --> 62.8 % chance they look into the one that was available

table(t.data$searched.in.targeted.box)
ftable(searched.in.targeted.box ~ nr, t.data)

# SETTING PLOT THEME -----------------------------------------------------------------
boot.reward <-
  boot.glmm.pred(
    model.res = chance.test_exp2, excl.warnings = T,
    nboots = 1000, para = F,
    level = 0.95)

conf.int <- as.data.frame(round(boot.reward$ci.estimates, 3))
m.stab.plot(round(boot.reward$ci.estimates, 3))

plogis(0.502)
plogis(0.231)
plogis(0.884)

my_theme <-
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification='left',
    legend.background = element_blank(),
    #legend.box.background = element_rect(colour = "black"),
    legend.text = element_text(size=15),

    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),

    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13), #, family = "Arial Narrow"),
    strip.text.x = element_text(size = 15),

    #text = element_text(family = "Arial Narrow"),

    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.title.y.left = element_text(vjust = 3),
    plot.title = element_text(color = "black", size = 15, face = "bold"),
    axis.title.x = element_text(
      margin =
        margin(t = 10, r = 0, b = 0, l = 0)
    ))

# PLOTTING -----------------------------------------------------------------
xdata$searched.in.targeted.box.nr <- as.numeric(xdata$searched.in.targeted.box)-1
xdata.agg <- xdata %>%
  group_by(nr) %>%
  summarise(mean.resp =
              mean(searched.in.targeted.box.nr, na.rm = T))

emm1 <- emmeans(chance.test_exp2, ~ 1)
conf.int <- as.data.frame(summary(emm1, type = "response"))

plot_study2_children <-
  ggplot(
    data = xdata.agg,
    aes(x = 1, y = mean.resp)
  ) +
  geom_point(
    data = xdata.agg, color = "gray30", size = 1.5,
    alpha = .4, position = position_jitter(w = 0.2, h = 0.025)
  ) +
  geom_violin(
    data = xdata.agg,
    aes(x = 1, y = mean.resp),
    position = position_nudge(x = 0.00),
    fill = "gray70", alpha = .2, color = NA
  ) +
  geom_errorbar(
    data = conf.int,
    aes(
      x = 1, y = plogis(boot.reward$ci.estimates$orig),
      ymin = plogis(boot.reward$ci.estimates$X2.5.),
      ymax = plogis(boot.reward$ci.estimates$X97.5.)
    ), color = "black",
    width = 0.1, linewidth = 1
  ) +
  geom_point(
    data = conf.int,
    aes(x = 1, y = plogis(boot.reward$ci.estimates$orig),
    color = "brown1", size = 4.5)
  )  +
  geom_hline(yintercept=0.5, linetype='dotted', col = 'blue') +
  scale_x_continuous(
    name = "Children"
  ) +
  scale_y_continuous(
    name = "Searches in the available box",
    labels = scales::percent
  ) +
  theme_light() +
  my_theme +
  theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none")

plot_study2_children

save.image("./study2/R_objects/Analysis_2.RData")
