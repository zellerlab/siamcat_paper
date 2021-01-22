# ##############################################################################
#
## Combine results
#
# ##############################################################################

.libPaths('/g/scb2/zeller/SHARED/software/R/3.5')

library("tidyverse")
library("here")

n.tax.jobs <- read_tsv(here('parameter_space', 'job_info', 'job_info.tsv'),
                       col_types = cols(
                         filt.method = col_character(),
                         filt.cutoff = col_character(),
                         norm.method = col_character(),
                         log.n0 = col_double(),
                         sd.min.q = col_double(),
                         ml.method = col_character(),
                         fs = col_logical(),
                         fs.method = col_character(),
                         fs.cutoff = col_double(),
                         job.id = col_double()
                       )) %>% 
  nrow()
n.func.jobs <- read_tsv(here('parameter_space', 'job_info', 
                             'job_info_func.tsv'),
                        col_types = cols(
                          filt.method = col_character(),
                          filt.cutoff = col_character(),
                          norm.method = col_character(),
                          log.n0 = col_double(),
                          sd.min.q = col_double(),
                          ml.method = col_character(),
                          fs = col_logical(),
                          fs.method = col_character(),
                          fs.cutoff = col_double(),
                          job.id = col_double()
                        )) %>% 
  nrow()

all.tasks <- read_tsv(here('parameter_space', 'data_info', 'all_tasks.tsv'))

expected.number <- tibble(type=c('RDP', 'mOTUs2', 'metaphlan', 
                                 'eggNOG', 'humann2'),
                          n.exp=c(n.tax.jobs, n.tax.jobs,n.tax.jobs,
                                  n.func.jobs, n.func.jobs)) %>% 
  full_join(all.tasks)


# ##############################################################################
# collect results
res.files <- list.files(
  here("parameter_space", "results"),
  pattern = '(RDP)|(mOTUs2)|(metaphlan)|(eggNOG)|(humann2)',
  full.names = TRUE)
res.files <- res.files[which(str_detect(res.files, '.tsv'))]
res <- map(res.files, read_tsv, col_types=cols(
  filt.method = col_character(),
  filt.cutoff = col_character(),
  norm.method = col_character(),
  log.n0 = col_double(),
  sd.min.q = col_double(),
  ml.method = col_character(),
  fs = col_logical(),
  fs.method = col_character(),
  fs.cutoff = col_double(),
  job.id = col_double(),
  dataset.id = col_character(),
  type = col_character(),
  case = col_character(),
  job.id = col_double(),
  processed = col_logical(),
  auroc = col_double(),
  time = col_double(),
  ci.low = col_double(),
  ci.high = col_double()
))
all.results <- bind_rows(res)


n.completed.jobs <- all.results %>% 
  group_by(dataset.id, type, case) %>% 
  summarise(n=n())

missing <- full_join(n.completed.jobs, expected.number, 
                     by=c('dataset.id', 'case', 'type')) %>% 
  ungroup() %>% mutate(n=as.double(n)) %>% 
  mutate(n=case_when(is.na(n)~0,TRUE~n))

if (nrow(missing %>% filter(n!=n.exp)) > 0){
  stop("Some jobs are still missing for: \n\t",
          missing %>% 
            filter(n!=n.exp) %>% 
            select(dataset.id, type, case) %>% 
            mutate(x=paste0(dataset.id, '-', type, '-', case)) %>% 
            pull(x) %>% paste(collapse = '\n\t'),
          "\nPlease double-check!")
    
} else {
  message("All jobs have finished! Yay!!!")
}

# ##############################################################################
# adjust some stuff for nicer plotting later
disease.match <- c(
  'NASH'='LIV', 'NAFLD'='LIV', 'LIV'='LIV', 'cirrhosis'='LIV',
  'ASD'='Neuro', 'PAR'='Neuro', 'PD'='Neuro',
  'CRC'='CRC', 'adenoma'='CRC', 'ADA'='CRC',
  'HIV'='HIV',
  'pre-hypertension'='CVD', 'hypertension'='CVD', 'ACVD'='CVD', 'HTN'='CVD', 'pHTN'='CVD',
  'MS'='OB', 'OB'='OB', 'T2D'='OB', 'IGT'='OB',
  'RA'='ART', 'ART'='ART', 'AS'='ART', 'T1D'='ART',
  'CDI'='DIR', 'NONCDI'='DIR', 'EDD'='DIR',
  'IBD'='IBD', 'CD'='IBD', 'UC'='IBD')
all.results$case.split <- disease.match[all.results$case]

disease.match.nice <- c('adenoma'='ADA', 'hypertension'='HTN',
                        'pre-hypertension'='pHTN', 'cirrhosis'='LIV',
                        'PAR'='PD', 'NONCDI'='nCDI', 'RA'='ART')
all.results <- all.results %>% 
  mutate(case.nice=case) %>% 
  mutate(case.nice=case_when(
    case.nice%in%names(disease.match.nice)~disease.match.nice[case.nice],
    TRUE~case.nice))

# adjust dataset IDs
dataset.match <- c(
  "scher"='Scher 2013', "kang"='Kang 2013', "son"='Son 2015',
  "schubert"='Schubert 2014', "vincent"='Vincent 2013', "baxter"='Baxter 2016',
  "chen"='Chen 2012', "wang"='Wang 2012', "zeller"='Zeller 2014',
  "singh"='Singh 2015', "dinh"='Dinh 2015', "lozupone"='Lozupone 2013',
  "noguerajulian"='Noguera-Julian 2016', "gevers"='Gevers 2014',
  "morgan"='Morgan 2012', "papa"='Papa 2012', "willing"='Willing 2010',
  "zhang"='Zhang 2013', "wong"='Wong 2013', "zhu"='Zhu 2013',
  "goodrich"='Goodrich 2014', "ross"='Ross 2015', "turnbaugh"='Turnbaugh 2009',
  "zupancic"='Zupancic 2012', "scheperjans"='Scheperjans 2015',
  "alkanani"='Alkanani 2015', 
  "metaHIT"='Nielsen 2014', "HMP2"="Lloyd-Price 2019",
  "FengQ_2015"='Feng 2015', "HanniganGD_2017"='Hannigan 2017',
  "KarlssonFH_2013"='Karlsson 2013', "KosticAD_2015"='Kostic 2015',
  "LiJ_2017"='Li 2017', "NielsenHB_2014"='Nielsen 2014',
  "QinJ_2012"='Qin 2012', "QinN_2014"='Qin 2014',
  "ThomasAM_2018a"='Thomas 2019a',
  "ThomasAM_2018b"='Thomas 2019b', "VincentC_2016"='Vincent 2016',
  "VogtmannE_2016"='Vogtmann 2016', "YuJ_2015"='Yu 2015',
  "ZellerG_2014"='Zeller 2014')

all.results <- all.results %>% 
  mutate(dataset.nice=dataset.id) %>% 
  mutate(dataset.nice=case_when(
    dataset.nice%in%names(dataset.match)~dataset.match[dataset.nice],
    TRUE~dataset.nice)) %>% 
  mutate(dataset.nice=str_replace(dataset.nice, '_', ' '))

all.results <- all.results %>% 
  mutate(label.nice=paste0(dataset.nice, ' ', case.nice))


write_tsv(all.results, path = here("parameter_space", "job_info",
                                   "full_results.tsv"))
