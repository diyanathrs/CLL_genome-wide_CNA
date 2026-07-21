# plot centre-wise survival

library(survival)
library(survminer)
library(dplyr)
library(gtsummary)

# load data
cll.outcome <- read.csv("../survival_mdr/All CLL outcome_spss.csv")
cll.outcome$OSdays <- cll.outcome$OSdays/365.25
cll.outcome$TTFTdays <- cll.outcome$TTFTdays/365.25
cll.outcome$OS_PTdays <- cll.outcome$OS_PTdays/365.25
table(cll.outcome$Study)

# add study2
cll.outcome$Study <- gsub('^Hull.*','Hull',cll.outcome$Study)
cll.outcome$Study <- gsub('^Newcastle.*','Newcastle',cll.outcome$Study)
cll.outcome$Study <- gsub('^Oxford.*','Oxford',cll.outcome$Study)
table(cll.outcome$Study)
str(cll.outcome$Study)
levels(cll.outcome$Study)

names(cll.outcome)

# median follow up by centre
cll.outcome %>%
  group_by(Study) %>%
  summarise(
    n = dplyr::n(),
    events = sum(OSstatus, na.rm = T),
    median_followup = median(OSdays, na.rm = T)
  ) %>% mutate(across(where(is.numeric), ~ round(.x, 2))) %>% gt()

##################
## summary stats
##################
names(cll.outcome)
dim(cll.outcome)
table(cll.outcome$Zap70)

cll.outcome <- cll.outcome %>%
  mutate(Sex = factor(Sex, levels = c(1, 2), labels = c("Female", "Male")),
         Binet = factor(Binet,
                        levels = c(1, 2, 3),
                        labels = c("A", "B", "C")),
         VH = factor(VH,
                     levels = c(0, 1),
                     labels = c("Unmutated", "Mutated")),
         CD38 = factor(CD38,
                       levels = c(0, 1),
                       labels = c("Negative", "Positive")),
         Zap70 = factor(Zap70,
                        levels = c(0, 1),
                        labels = c("Negative", "Positive")))

tbl1 <- cll.outcome %>%
  select(Age,
    Sex,    Binet,
    VH, CD38, Zap70) %>%
  tbl_summary(statistic = list(
      all_continuous() ~ "{median} ({p25}–{p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      Age ~ "Age (years)",
      Sex ~ "Sex",
      Binet ~ "Binet stage",
      VH ~ "IGHV mutated",
      CD38 ~ "CD38 positive",
      Zap70 ~ "Zap-70 positive"
    ),
    missing = "ifany",
    missing_text = "Missing"
  ) %>%
  bold_labels()

tbl1


#########################
## survival for vars
#########################
library(survminer)
#sex
fit <- survfit(Surv(OSdays, OSstatus) ~ Sex, data = cll.outcome)
sex <- ggsurvplot(fit, pval = TRUE, risk.table = F, palette = 'Set1', xlab="Time(years)")

#binet
fit <- survfit(Surv(OSdays, OSstatus)~ Binet, data = cll.outcome)
binet <- ggsurvplot(fit, pval = TRUE, risk.table = F, palette = 'Set1', xlab="Time(years)")

#IGHV
fit <- survfit(Surv(OSdays, OSstatus) ~ VH, data = cll.outcome)
vh <- ggsurvplot(fit, pval = TRUE, risk.table = F, palette = 'Set1', xlab="Time(years)")

#CD38
fit <- survfit(Surv(OSdays, OSstatus) ~ CD38, data = cll.outcome)
cd38 <- ggsurvplot(fit, pval = TRUE, risk.table = F, palette = 'Set1', xlab="Time(years)") 


#arrange
surv.plts <- sex$plot + binet$plot + vh$plot + cd38$plot +
  plot_layout(ncol = 2) + plot_annotation(tag_levels = 'A')
surv.plts

ggsave(filename = "plots/CLLsurv.plts_OS.png",
       plot = surv.plts,
       width = 11,
       height = 9,
       dpi = 300)


# multivariable model for clinical factors
cll.outcome$Study <- factor(cll.outcome$Study, levels = c("Oxford" ,"Bournemouth" ,"Cardiff", 
                                                         "Hull" ,"Leicester" ,"Newcastle", "Southampton"))
table(cll.outcome$Study)
names(cll.outcome)[names(cll.outcome) == "VH"] <- "IGHV"

fit.clinical <- coxph(Surv(OSdays, OSstatus) ~
    Age +
    Sex +
    Binet +
    IGHV +
    Study,
  data = cll.outcome)

summary(fit.clinical)

# prop hazard
ph.test <- cox.zph(fit.clinical)

ph.test
plot(ph.test)

library(forestmodel)

forest_model(fit.clinical)

p <- ggforest(fit.clinical, data = cll.outcome, fontsize = 0.8)
p

ggsave("plots/Clinical_multivariable_cox_forest.pdf",
  p,dpi = 300,
  width = 8,
  height = 8)

# diagnosis
cll.outcome |>
  dplyr::filter(
    complete.cases(OSdays, OSstatus, Age, Sex, Binet, IGHV, Study)
  ) |>
  dplyr::count(Study, OSstatus)

table(model.frame(fit.clinical)$Study,
      model.frame(fit.clinical)$OSstatus)

#table 
library(broom)
library(dplyr)

cox.table <- tidy(
  fit.clinical,
  exponentiate = TRUE,
  conf.int = TRUE
) %>%
  select(
    Variable = term,
    HR = estimate,
    Lower95 = conf.low,
    Upper95 = conf.high,
    P = p.value
  )

cox.table

cox.table$Variable <- recode(
  cox.table$Variable,
  "Age" = "Age (per year)",
  "SexMale" = "Male",
  "BinetB" = "Binet stage B",
  "BinetC" = "Binet stage C",
  "VHMutated" = "IGHV mutated",
  "StudyCardiff" = "Cardiff",
  "StudyHull" = "Hull",
  "StudyLeicester" = "Leicester",
  "StudyNewcastle" = "Newcastle",
  "StudySouthampton" = "Southampton"
)

###########################
# Survival
#########################
time_var <- 'OS_PTdays' #OS_PTdays
event_var <- 'OS_PTstatus' #OS_PTstatus
var <- "Study"
covars <- c(var, "Sex", "Binet", "Zap70")

#cox
fit <- coxph(reformulate(var,
              response = paste0("Surv(", time_var, ", ", event_var, ")")),
  data = cll.outcome)

summary(fit)

# cox adjusted
fit <- coxph(reformulate(covars,
                         response = paste0("Surv(", time_var, ", ", event_var, ")")),
             data = cll.outcome)

summary(fit)

# Kaplan- meier for centre
fit <- survfit(reformulate("Study", response = paste0("Surv(", time_var, ", ", event_var, ")")),
  data = cll.outcome)
centre <- ggsurvplot(fit, pval = TRUE, risk.table = F, palette = 'Set3', xlab='Time(years)',
                     legend = "right",, font.legend=c(10))

centre

# save 
ggsave(filename = "plots/CLLsurv.centre_OS.png",
       plot = centre$plot,
       width = 7.5,
       height = 5,
       dpi = 300)

# check more on centre 
coxph(Surv(OS_time, OS_event) ~ study +
        age + Binet + IGHV + CD38)

# ox vs others
time_var <- 'OSdays' #OS_PTdays
event_var <- 'OSstatus' #OS_PTstatus

cll.outcome <- cll.outcome %>% mutate(centre2 = if_else(grepl("^Oxford", Study), "Oxford", "Others"))
fit <- survfit(reformulate("centre2", response = paste0("Surv(", time_var, ", ", event_var, ")")),
               data = cll.outcome)
ggsurvplot(fit, pval = TRUE, risk.table = TRUE)

####################
## adjusted cox
########################
time_var <- 'OSdays' #OS_PTdays
event_var <- 'OSstatus' #OS_PTstatus
var <- "Study2"
covars <- c(var, "Age", "Sex", "Binet", "VH")

coxph(reformulate(covars, response = paste0("Surv(", time_var, ", ", event_var, ")")), data = cll.outcome)

# proportional hazard
fit <- coxph(reformulate(covars, response = paste0("Surv(", time_var, ", ", event_var, ")")), 
             data = cll.outcome)
plot(cox.zph(fit))

#Stratified Cox model
coxph(Surv(OSdays, OSstatus) ~ Age + Sex + strata(Study2), data = cll.outcome)
