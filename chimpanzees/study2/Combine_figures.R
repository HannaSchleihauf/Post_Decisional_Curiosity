load("/Users/hanna.schleihauf/Dropbox/Research/my_projects/Counterfactual_Curiousity/Counterfactual_Curiosity_shared/children/study2/images/Analysis_2.RData")

# Combine plots
library(ggpubr)
# theme_set(theme_pubr())
figure <- ggarrange(plot_study2, plot_study2_chimpanzees,
                    labels = c(
                      "(a)",
                      "(b)"
                    ),
                    # label.x = 0,    # X position 0 for left
                    # label.y = 1,    # Y position 1 for top
                    hjust = -1, # hjust = 0 for left alignment
                    vjust = 2.5, # hjust = 0 for left alignment
                    ncol = 2, nrow = 1
)
figure
