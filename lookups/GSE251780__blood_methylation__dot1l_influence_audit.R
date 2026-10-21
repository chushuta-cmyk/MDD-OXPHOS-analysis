.libPaths(c('/Users/tan/Developer/envs/R/lib/R/library', .libPaths()))
library(sandwich)
out <- '/Users/tan/Documents/Codex/2026-09-19/ple/outputs/step31_step4_audit_2026-10-04'
root <- '/Users/tan/Documents/MDD-MultiOmics/analysis_framework/GSE251780'
x <- read.delim(file.path(out,'target15_processed_matrix.tsv'),check.names=FALSE)
m <- read.csv(file.path(root,'data/GSE251780_metadata.csv'))
c <- read.delim(file.path(root,'results/GSE251780_EpiDISH_cell_fractions.tsv'))
d <- merge(m,c,by='SampleID',sort=FALSE)
stopifnot(nrow(d)==147,!anyDuplicated(d$SampleID),all(as.character(d$SampleID)%in%names(x)))
d$Group <- factor(d$Group,levels=c('Control','MDD'))
d$Sex <- factor(d$Sex)
old <- read.delim(file.path(out,'Shen2025_crosscheck_corrected.tsv'))
checks <- list(); fits <- list()
for(i in seq_len(nrow(x))) {
  b <- as.numeric(x[i,as.character(d$SampleID)])
  d$M <- log2(pmin(pmax(b,1e-6),1-1e-6)/(1-pmin(pmax(b,1e-6),1-1e-6)))
  for(model in c('sex_adjusted','cell_adjusted')) {
    f <- if(model=='sex_adjusted') M~Group+Sex else M~Group+Sex+B+NK+CD4T+CD8T+Mono+Eosino
    fit <- lm(f,d)
    stopifnot(fit$rank==length(coef(fit)))
    target <- old[old$CpG==x$ID_REF[i]&old$Model==model,]
    delta <- unname(coef(fit)['GroupMDD'])-target$Local_M_difference
    stopifnot(abs(delta)<1e-10)
    checks[[length(checks)+1]] <- data.frame(CpG=x$ID_REF[i],Model=model,Coefficient_difference=delta)
    if(x$ID_REF[i]=='cg04173586') {
      se <- sqrt(vcovHC(fit,type='HC3')['GroupMDD','GroupMDD'])
      beta <- coef(fit)['GroupMDD']; p <- 2*pt(-abs(beta/se),df.residual(fit))
      loo <- sapply(seq_len(nrow(d)),function(j)coef(lm(f,d[-j,]))['GroupMDD'])
      fits[[model]] <- data.frame(Model=model,N=nrow(d),Beta_M=beta,SE_HC3=se,P_HC3_t=p,
        CI_low=beta-qt(.975,df.residual(fit))*se,CI_high=beta+qt(.975,df.residual(fit))*se,
        P_OLS=coef(summary(fit))['GroupMDD','Pr(>|t|)'],LOO_beta_min=min(loo),LOO_beta_max=max(loo),Max_Cooks_distance=max(cooks.distance(fit)))
      write.table(data.frame(SampleID=d$SampleID,LOO_beta=loo,Cooks_distance=cooks.distance(fit),Leverage=hatvalues(fit)),file.path(out,paste0('DOT1L_influence_',model,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
    }
  }
}
write.table(do.call(rbind,checks),file.path(out,'target14_coefficient_reproduction.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
res <- do.call(rbind,fits); print(res)
write.table(res,file.path(out,'DOT1L_diagnostics.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
# Diagnostic figure only: raw beta distributions are not adjusted disease effects.
b <- as.numeric(x[x$ID_REF=='cg04173586',as.character(d$SampleID)])
png(file.path(out,'DOT1L_diagnostic.png'),width=1500,height=700,res=140)
par(mfrow=c(1,2),mar=c(5,5,3,1))
boxplot(b~d$Group,ylab='Raw methylation beta',xlab='',main='cg04173586 (DOT1L): raw distribution',col=c('#A6CEE3','#FDBF6F'))
set.seed(20261004);stripchart(b~d$Group,vertical=TRUE,method='jitter',pch=16,add=TRUE,col='#33333366')
plot(res$Beta_M,1:2,xlim=range(res$CI_low,res$CI_high),ylim=c(.5,2.5),yaxt='n',ylab='',xlab='M-value difference (MDD - control)',main='OLS coefficient with HC3 95% CI',pch=19)
axis(2,at=1:2,labels=c('Sex','Sex + cells'),las=1);segments(res$CI_low,1:2,res$CI_high,1:2);abline(v=0,lty=2,col='gray50')
dev.off()
cat('Detection P range: ',range(as.numeric(x[x$ID_REF=='cg04173586',paste0(d$SampleID,'.Detection.Pval')])),'\n')
