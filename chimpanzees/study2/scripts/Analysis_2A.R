# Analysis Study 2 - Chimpanzees

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
# load("./study2/R_objects/analysis_2.RData")

xdata <-
  read.csv("./study2/data/data_study2A_chimpanzees.csv",
           header = TRUE, na = c("NA", "")
  )

# Interrater Reliability ------------------------------------------------
# number of dyads for which the reliability coding was not the same
sum(xdata$searched.in.available.reli !=
  xdata$searched.in.available, na.rm = TRUE)

# Cohen's Kappa
library(irr)
kappa_results <-
  kappa2(cbind(
    xdata$searched.in.available.reli,
    xdata$searched.in.available
  ), "unweighted")
print(kappa_results)

# Descriptives
prop.table(table(xdata$searched.in.available))

# Fitting the Model ------------------------------------------------
xx.fe.re <- fe.re.tab(
  fe.model = "searched.in.available ~ trial",
  re = "(1|name)", data = xdata
) # maybe add age
xx.fe.re$summary
t.data <- xx.fe.re$data
str(t.data)

t.data$z.trial <- scale(t.data$trial)

contr <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000000))
chance.test_exp2 <- glmer(
  searched.in.available ~ 1 +
    (1 + z.trial | name),
  data = t.data, family = binomial, control = contr
)

## Test assumptions ----
overdisp.test(chance.test_exp2) # no overdispersion
ranef.diagn.plot(full)

## Checking model stability
m.stab.b <-
  glmm.model.stab(model.res = chance.test_exp2, contr = contr, use = c("name"))
m.stab.b$detailed$warnings
as.data.frame(round(m.stab.b$summary[, -1], 3))
m.stab.plot(round(m.stab.b$summary[, -1], 3))

# CHECKING MODEL RESULTS ------------------------------------------------
round(summary(chance.test_exp2)$coefficient, 3)
plogis(0.4390) # --> 60.8 % chance they look into the one that was available
ftable(searched.in.available ~ name, t.data)

boot <-
  boot.glmm.pred(
    model.res = chance.test_exp2, excl.warnings = T,
    nboots = 1000, para = F,
    level = 0.95
  )

plogis(0.502)
plogis(0.231)
plogis(0.884)

as.data.frame(round(boot$ci.estimates, 3))
m.stab.plot(round(boot$ci.estimates, 3))

# SETTING PLOT THEME ------------------------------------------------
my_theme <-
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "left",
    legend.background = element_blank(),
    # legend.box.background = element_rect(colour = "black"),
    legend.text = element_text(size = 15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13), # , family = "Arial Narrow"),
    strip.text.x = element_text(size = 15),

    # text = element_text(family = "Arial Narrow"),

    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.title.y.left = element_text(vjust = 3),
    plot.title = element_text(color = "black", size = 15, face = "bold"),
    axis.title.x = element_text(
      margin =
        margin(t = 10, r = 0, b = 0, l = 0)
    )
  )

# PLOTTING ESTIMATES ------------------------------------------------

conf.int <- as.data.frame(round(boot$ci.estimates, 3))
xdata$searched.in.available.nr <- NA
xdata <-
  xdata %>%
  mutate(searched.in.available.nr = ifelse(searched.in.available == "yes", 1,
    ifelse(searched.in.available == "no", 0, NA)
  ))

xdata.agg <- xdata %>%
  group_by(name) %>%
  summarise(
    mean.resp =
      mean(searched.in.available.nr, na.rm = T)
  ) %>%
  na.omit()

plot_study2_chimpanzees <-
  ggplot(
    data = xdata.agg,
    aes(x = 1, y = mean.resp)
  ) +
  theme_light() +
  geom_point(
    data = xdata.agg, color = "gray30", size = 1.5,
    alpha = .4, position = position_jitter(w = 0.25, h = 0)
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
      x = 1, y = plogis(boot$ci.estimates$orig),
      ymin = plogis(boot$ci.estimates$X2.5.),
      ymax = plogis(boot$ci.estimates$X97.5.)
    ), color = "black",
    width = 0.1, linewidth = 1
  ) +
  geom_point(
    data = conf.int,
    aes(x = 1, y = plogis(boot$ci.estimates$orig), ),
    color = "brown1", size = 4.5
  ) +
  geom_hline(yintercept = 0.5, linetype = "dotted", col = "blue") +
  scale_x_continuous(
    name = "Chimpanzees"
  ) +
  scale_y_continuous(
    name = "Searches in the available box",
    limits = c(0, 1),
    labels = scales::percent,
  ) +
  my_theme +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

plot_study2_chimpanzees

# Plot for individual subjects ------------------------------------------------
xdataind <- xdata %>%
  group_by(name, searched.in.available) %>%
  count(searched.in.available)
xdataind <- pivot_wider(xdataind, names_from = searched.in.available, values_from = n)
xdataind$prob <- xdataind$yes / (xdataind$yes + xdataind$no)

xdata_plot <- ggplot(data = xdataind, aes(y = prob, x = name)) +
  geom_boxplot() +
  ylab("Proportion of information search in available box") +
  xlab("") +
  geom_hline(yintercept = 0.5, linetype = "dotted", col = "blue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text = element_text(size = 15))

xdata_plot

save.image("./study2/R_objects/analysis_2.RData")
