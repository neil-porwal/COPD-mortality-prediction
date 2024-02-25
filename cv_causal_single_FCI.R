# Training different causal discoverly algorithms to figure out which one to actually run on the real data to find the best biomarkers for mortality rate prediction

library(rCausalMGM)
library(survival)
library(dplyr)
library(optparse)
library(readxl)

option_list = list(
  make_option(c("-k", "--fold"), type="numeric", default="0",
              help="CV Fold (1 to 10)"),
  make_option(c("-g", "--gold"), type="character", default=NA,
              help="All or II_IV"),
  make_option(c("-m", "--mgm"), action='store_true', default='FALSE',
              help="Learn an initial skeleton with MGM"),
  make_option(c("-f", "--fci"), action='store_true', default='FALSE',
              help="Learn causal graph with FCI"),
  make_option(c("-o", "--orientrule"), type="character", default="majority",
              help="Orientation rule: majority, maxp, or conservative"),
  make_option(c("-a", "--alpha"), type="numeric", default="0.1",
              help="Alpha for FDR control")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(opt)

k <- opt$fold
goldstage <- opt$gold
mgmFlag <- opt$mgm
fciFlag <- opt$fci
rankFlag <- TRUE
orientRule <- opt$orientrule
alpha <- opt$alpha

mgm_or_fci <- ""

if (mgmFlag) {
  mgm_or_fci <- "MGM"
} else if (fciFlag) {
  mgm_or_fci <- "FCI"
}

alpha_label <- alpha
if (alpha == 0.2) {
  alpha_label <- 0.2
} else if (alpha == 0.05) {
  alpha_label <- 0.05
}

gold_label <- ""

if (goldstage == 'All') {
  gold_label <- "All"
  train <- read.csv(paste('/net/talisker/home/benos/porwaln/__MACOSX/COPDGene/Cross-Validation-Matricies/All-GOLD-Stages/train/all_copd_model_data_train_', k, '.csv', sep = ""), header=T, row.names=1)  # Add absolute path
  train <- train[, colnames(train) != 'basophl_pct']
} else if (goldstage == 'II_IV') {
  gold_label <- "II_IV"
  train <- read.csv(paste('/net/talisker/home/benos/porwaln/__MACOSX/COPDGene/Cross-Validation-Matricies/2-4-GOLD-Stages/train/2_4_copd_model_data_train_', k, '.csv', sep = ""), header=T, row.names=1)  # Add absolute path
  train <- train[, colnames(train) != 'basophl_pct']
}

if (any(c('overall.time', 'overall.status') %in% colnames(train))) {
  train$overall <- Surv(train$overall.time, train$overall.status)
}

train <- train[,!colnames(train) %in% c('overall.time', 'overall.status')]

# Factor Vars
clin.dict <- read_excel('/net/talisker/home/benos/porwaln/__MACOSX/COPDGene/COPDGene_Data-P1P2P3.2021Aug/COPDGene_P1P2P3_Visitlevel_DataDict_Mar20_rev_16Aug21.xlsx', 1) %>% as.data.frame

rownames(clin.dict) <- clin.dict$VariableName

factor.vars <- c(setdiff(rownames(clin.dict)[clin.dict$CodedVariable=='Y'],
                         c('HealthStatus', 'SchoolCompleted', 'ATS_ERS',
                           'DrugCostCovered', 'LungDiseaseInformed')),
                 c('Internet', 'Insurance', 'CVD'))

train <- train %>% mutate_if(colnames(train) %in% factor.vars, factor)

ig <- NULL

knowledge24 <- readRDS('copd_2_4_prior_knowledge.rds')
knowledge24$tiers[[3]] <- knowledge24$tiers[[3]][-which(knowledge24$tiers[[3]] == "basophl_pct")]
knowledge <- knowledge24

if (mgmFlag) {
  
  ig.steps <- steps(train, rank=rankFlag, verbose=T)
  
  pdf(paste0('Causal-Discovery-Model-Testing/', gold_label, '/MGM/steps.instabs.metabric.rna.cv', k, '.pdf'), width=5, height=5)
  
  log10params <- log10(ig.steps$lambdas)
  plot(x=log10params, y=ig.steps$instability[,6], col='black', pch=19,
       xlab=expression(log10(lambda)),
       ylab="Edge instability across subsamples",
       ylim=c(0,min(0.5, 2*max(ig.steps$instab, na.rm=T))),
       cex=1)
  
  points(x=log10params, y=ig.steps$instability[,1], col='red', pch=19, cex=1)
  points(x=log10params, y=ig.steps$instability[,2], col='dodgerblue', pch=19, cex=1)
  points(x=log10params, y=ig.steps$instability[,3], col='purple', pch=19, cex=1)
  points(x=log10params, y=ig.steps$instability[,4], col='orange', pch=19, cex=1)
  points(x=log10params, y=ig.steps$instability[,5], col='green', pch=19, cex=1)
  
  abline(h=ig.steps$gamma, lty=5, col='gray', lwd=3, cex=1)
  
  allIdx <- c(1, which(ig.steps$instability[1:which.max(ig.steps$instability[,6]),6]<ig.steps$gamma))
  allIdx <- allIdx[length(allIdx)]
  
  abline(v=log10params[allIdx], col='black',  lty=2, lwd=3)
  abline(v=log10(ig.steps$graph$lambda[1]), col='red',  lty=2, lwd=3)
  abline(v=log10(ig.steps$graph$lambda[2]), col='dodgerblue',  lty=2, lwd=3)
  abline(v=log10(ig.steps$graph$lambda[3]), col='purple',  lty=2, lwd=3)
  abline(v=log10(ig.steps$graph$lambda[4]), col='orange',  lty=2, lwd=3)
  abline(v=log10(ig.steps$graph$lambda[5]), col='green',  lty=2, lwd=3)
  
  legend(x = "topleft", title="Edge Type",
         legend = c("All", "CC", "CD", "DD", "SC", "SD"),
         col = c("black","red", "dodgerblue", "purple", "orange", "green"),
         pch = 19, cex=1)
  
  dev.off()
  
  saveGraph(ig.steps$graph,
            paste0('Causal-Discovery-Model-Testing/', gold_label, '/MGM/out/mgmStEPS.metabric.rna',
                   ifelse(rankFlag, '.rank', '.linear'),
                   '.cv', k, '.txt'))
  saveGraph(ig.steps$graph,
            paste0('Causal-Discovery-Model-Testing/', gold_label, '/MGM/graph/mgmStEPS.metabric.rna',
                   ifelse(rankFlag, '.rank', '.linear'),
                   '.cv', k, '.sif'))
  
  ig.stars <- coxmgm(train, lambda=ig.steps$lambdas[allIdx], rank=rankFlag, verbose=T)
  
  saveGraph(ig.stars,
            paste0('Causal-Discovery-Model-Testing/', gold_label, '/MGM/out/mgmStARS.metabric.rna',
                   ifelse(rankFlag, '.rank', '.linear'),
                   '.cv', k, '.txt'))
  saveGraph(ig.stars,
            paste0('Causal-Discovery-Model-Testing/', gold_label, '/MGM/graph/mgmStARS.metabric.rna',
                   ifelse(rankFlag, '.rank', '.linear'),
                   '.cv', k, '.sif'))
  
  ig.path <- coxmgmPath(train, rank=rankFlag, verbose=T)
  
  p <- ncol(train)
  n <- nrow(train)
  kappa <- log(p) / log(n)
  gamma <- max(0, 1-1/(4*kappa)) / 2
  gamma <- 0.5
  penalty <- 1 + 4 * gamma * kappa
  ebic <- -2 *ig.path$loglik + penalty * log(nrow(train)) * ig.path$nParams
  
  pdf(paste0('Causal-Discovery-Model-Testing/', gold_label, '/MGM/solution.path.metabric.rna.cv', k, '.pdf'), width=5, height=5)
  plot(ig.path)
  points(log10(ig.path$lambdas), ebic / (2 * nrow(train)), col='purple', pch=19)
  abline(v=log10(ig.path$lambdas)[which.min(ebic)], lty=2, col='purple', lwd=2)
  dev.off()
  
  ig.path$graph.ebic <- ig.path$graphs[[which.min(ebic)]]
  
  saveGraph(ig.path$graph.bic,
            paste0('Causal-Discovery-Model-Testing/', gold_label, '/MGM/out/mgmBIC.metabric.rna',
                   ifelse(rankFlag, '.rank', '.linear'),
                   '.cv', k, '.txt'))
  saveGraph(ig.path$graph.bic,
            paste0('Causal-Discovery-Model-Testing/', gold_label, '/MGM/graph/mgmBIC.metabric.rna',
                   ifelse(rankFlag, '.rank', '.linear'),
                   '.cv', k, '.sif'))
  
  saveGraph(ig.path$graph.ebic,
            paste0('Causal-Discovery-Model-Testing/', gold_label, '/MGM/out/mgmEBIC.metabric.rna',
                   ifelse(rankFlag, '.rank', '.linear'),
                   '.cv', k, '.txt'))
  saveGraph(ig.path$graph.ebic,
            paste0('Causal-Discovery-Model-Testing/', gold_label, '/MGM/graph/mgmEBIC.metabric.rna',
                   ifelse(rankFlag, '.rank', '.linear'),
                   '.cv', k, '.sif'))
}

if (fciFlag) {
  ig <- loadGraph(paste0('Causal-Discovery-Model-Testing/', gold_label, '/MGM/out/mgmStARS.metabric.rna',
                         ifelse(rankFlag, '.rank', '.linear'),
                         '.cv', k, '.txt'))
  
  g <- fciStable(train, initialGraph=ig, alpha=alpha, orientRule=orientRule,
                 knowledge=knowledge, fdr=T, rank=rankFlag, verbose=T)
  
  saveGraph(g,
            paste0('Causal-Discovery-Model-Testing/', gold_label, '/FCI/', orientRule, '/', alpha, '/out/mgmStARS.metabric.rna',
                   ifelse(rankFlag, '.rank', '.linear'),'.FDR',
                   gsub('0[.]', '', as.character(alpha)),
                   '.cv', k, '.txt'))
  
  saveGraph(g,
            paste0('Causal-Discovery-Model-Testing/', gold_label, '/FCI/', orientRule, '/', alpha, '/graph/mgmStARS.metabric.rna',
                   ifelse(rankFlag, '.rank', '.linear'),'.FDR',
                   gsub('0[.]', '', as.character(alpha)),
                   '.cv', k, '.sif'))
  
  ig <- loadGraph(paste0('Causal-Discovery-Model-Testing/', gold_label, '/MGM/out/mgmStEPS.metabric.rna',
                         ifelse(rankFlag, '.rank', '.linear'),
                         '.cv', k, '.txt'))
  
  g <- fciStable(train, initialGraph=ig, alpha=alpha, orientRule=orientRule,
                 knowledge=knowledge, fdr=T, rank=rankFlag, verbose=T)
  
  saveGraph(g,
            paste0('Causal-Discovery-Model-Testing/', gold_label, '/FCI/', orientRule, '/', alpha, '/out/mgmStEPS.metabric.rna',
                   ifelse(rankFlag, '.rank', '.linear'),'.FDR',
                   gsub('0[.]', '', as.character(alpha)),
                   '.cv', k, '.txt'))
  
  saveGraph(g,
            paste0('Causal-Discovery-Model-Testing/', gold_label, '/FCI/', orientRule, '/', alpha, '/graph/', 'mgmStEPS.metabric.rna',
                   ifelse(rankFlag, '.rank', '.linear'),'.FDR',
                   gsub('0[.]', '', as.character(alpha)),
                   '.cv', k, '.sif'))
  
  ig <- loadGraph(paste0('Causal-Discovery-Model-Testing/', gold_label, '/MGM/out/mgmEBIC.metabric.rna',
                         ifelse(rankFlag, '.rank', '.linear'),
                         '.cv', k, '.txt'))
  
  g <- fciStable(train, initialGraph=ig, alpha=alpha, orientRule=orientRule,
                 knowledge=knowledge, fdr=T, rank=rankFlag, verbose=T)
  
  saveGraph(g,
            paste0('Causal-Discovery-Model-Testing/', gold_label, '/FCI/', orientRule, '/', alpha, '/out/mgmEBIC.metabric.rna',
                   ifelse(rankFlag, '.rank', '.linear'),'.FDR',
                   gsub('0[.]', '', as.character(alpha)),
                   '.cv', k, '.txt'))
  
  saveGraph(g,
            paste0('Causal-Discovery-Model-Testing/', gold_label, '/FCI/', orientRule, '/',alpha, '/graph/mgmEBIC.metabric.rna',
                   ifelse(rankFlag, '.rank', '.linear'),'.FDR',
                   gsub('0[.]', '', as.character(alpha)),
                   '.cv', k, '.sif'))
  
  ig <- loadGraph(paste0('Causal-Discovery-Model-Testing/', gold_label, '/MGM/out/mgmBIC.metabric.rna',
                         ifelse(rankFlag, '.rank', '.linear'),
                         '.cv', k, '.txt'))
  
  g <- fciStable(train, initialGraph=ig, alpha=alpha, orientRule=orientRule,
                 knowledge=knowledge, fdr=T, rank=rankFlag, verbose=T)
  
  saveGraph(g,
            paste0('Causal-Discovery-Model-Testing/', gold_label, '/FCI/', orientRule, '/', alpha, '/out/mgmBIC.metabric.rna',
                   ifelse(rankFlag, '.rank', '.linear'),'.FDR',
                   gsub('0[.]', '', as.character(alpha)),
                   '.cv', k, '.txt'))
  
  saveGraph(g,
            paste0('Causal-Discovery-Model-Testing/', gold_label, '/FCI/', orientRule, '/', alpha, '/graph/mgmBIC.metabric.rna',
                   ifelse(rankFlag, '.rank', '.linear'),'.FDR',
                   gsub('0[.]', '', as.character(alpha)),
                   '.cv', k, '.sif'))
  
  
}