# calculate survival and plot graphs
library(haven)
install.packages(c("haven"))
library(survival)
library(survminer)
library(stringr)

setwd('survival_mdr/')
outcome.dat <- data.frame(read_sav('All CLL outcome.sav'))
# call gene function
MDR_gene('EIF3IP1',a)

all.algo.sample.list$Del_ARID1B <- '1'
#change name of ox sample ids
all.algo.sample.list$Sample[1] <- "ARC00245"
all.algo.sample.list$Sample[2] <- "ADM00110"
all.algo.sample.list$Sample[3] <- "ADM00398"

outcome.merged <- left_join(x=outcome.dat,y=all.algo.sample.list, by='Sample')
outcome.merged$Del_ARID1B[is.na(outcome.merged$Del_ARID1B)]  <- 0

# overall survival
fitos <- survfit(Surv(outcome.merged$OSdays/365.25, outcome.merged$OSstatus) ~ outcome.merged$Del_ARID1B, data = outcome.merged)
summary(fitos)

ggsurvplot(fitos,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

# TTFT survival
fit_ttft <- survfit(Surv(outcome.merged$TTFTdays/365.25, outcome.merged$TTFTstatus) ~ outcome.merged$Del_ARID1B, 
                    data = outcome.merged)
summary(fit_ttft)

ggsurvplot(fit_ttft,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

# OSPT survival
fit_pt <- survfit(Surv(outcome.merged$OS_PTdays/365.25, outcome.merged$OS_PTstatus) ~ outcome.merged$Del_ARID1B, 
                    data = outcome.merged) 

ggsurvplot(fit_pt,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))
