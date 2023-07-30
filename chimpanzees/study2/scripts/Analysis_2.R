library("lme4")
library("readxl")
library("tidyverse")
library("tidyselect")
library("parallel")
library("dfoptim")
library("optimx")
library("emmeans")
source("./functions/diagnostic_fcns.r")
source("./functions/glmm_stability.r")
source("./functions/drop1_para.r")
source("./functions/boot_glmm.r")
library("emmeans")
library("car")


xdata <- 
  read_xlsx("/Users/hanna.schleihauf/Dropbox/Research/my_projects/Counterfactual_Curiousity_Chimps/Counterfactual_Curiosity_Study_2.xlsx", col_names = TRUE, na= c("NA"))

## Check data 
xdata$searched.in.available

## Prepare data for model fitting------------------------------------
xx.fe.re=fe.re.tab(fe.model="searched.in.available ~ trial",
                   re="(1|name)", data=xdata)  #maybe add age
xx.fe.re$summary
t.data=xx.fe.re$data 
str(t.data) 

t.data$z.trial = scale(t.data$trial)

contr <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000000))
chance.test_exp2 <- glmer(searched.in.available ~ 1 +
                            (1 + z.trial | name),
                          data = t.data, family = binomial, control = contr
)

summary(chance.test_exp2)

plogis(0.4390)   # --> 60.8 % chance they look into the one that was available

ftable(searched.in.available ~ name, t.data)

#plot
# For individual subjects
xdataind <- xdata %>% group_by(name, searched.in.available)  %>% count(searched.in.available)
xdataind <- pivot_wider(xdataind, names_from= searched.in.available, values_from= n) 
xdataind$prob <- xdataind$yes / (xdataind$yes + xdataind$no)

xdata_plot <- ggplot(data = xdataind, aes(y = prob, x = name)) +  
  geom_boxplot() + ylab("Proportion of information search in target box") + 
  xlab("") + geom_hline(yintercept=0.5, linetype='dotted', col = 'blue')  + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust=1))  + 
  theme(axis.title = element_text(size = 15)) + theme(axis.text = element_text(size = 15))

xdata_plot


