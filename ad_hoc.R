# Ad hoc coding

fitos <- survfit(Surv(OSdays/365.25, OSstatus) ~Study, data = outcomes)
fitttfs <- survfit(Surv(TTFTdays/365.25, TTFTstatus) ~Study, data = outcomes)
fitpt <- survfit(Surv(OS_PTdays/365.25, OS_PTstatus) ~Study, data = outcomes)

unique(outcomes$study2)
outcomes$study2 <- gsub('Newcastle2','Newcastle',outcomes$study2)
stu.order <- c("Leicester", "Bournemouth","Hull","Newcastle","Southampton",
               "Cardiff", "Oxford-ARC","Oxford-ADM" )
outcomes$study2 <- factor(outcomes$study2,levels = stu.order )

ggsurvplot(fitos, pval = TRUE, conf.int = F,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", ncensor.plot = TRUE,# Specify median survival
           ggtheme = theme_bw())

ggsurvplot(fitttfs, pval = TRUE, conf.int = F,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata",
           ncensor.plot = TRUE,# Change line type by groups # Specify median survival
           ggtheme = theme_bw())

ggsurvplot(fitpt, pval = TRUE, conf.int = F,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata",# Change line type by groups
           ncensor.plot = TRUE, # Specify median survival
           ggtheme = theme_bw())


## check interval length of each algo###
#######################################
# use length instead of probes - > 500kb
gr.cll.cnv <- read.csv('All_cnvs_with_ox.csv') %>% 
  filter(!study %in% c('Oxford','Oxford-ARC','Oxford-ADM')) %>%  
  filter(!chr %in% c('X','Y','23')) %>% filter(CNV_type==type) %>% toGRanges() 

type <- 'Loss'
algo <- 'QuantiSNP'

if (!algo =='all') {
  gr.cll.fil <- gr.cll.cnv %>% filter(method==algo)
}
if (length(gr.cll.cnv) < 1) {stop("No CNAs left after algorithm filtering, check algo argument..")}

seqlevelsStyle(gr.cll.fil) <- 'UCSC'
head(gr.cll.fil)
range(gr.cll.fil$length)
range(gr.cll.fil$numSNP)
table(gr.cll.fil$method)
#plot
data.frame(gr.cll.fil) %>% ggplot(aes(length, fill = method)) + geom_histogram(bins=10)
# check min gap
red <- reduce(gr.cll.fil, min.gapwidth=100000)
min(data.frame(red)$width)
max(data.frame(red)$width)
data.frame(red) %>% ggplot(aes(width)) + geom_histogram(bins = 10)
