.libPaths(c(.libPaths(), '/Users/tan/Developer/envs/R/lib/R/library'))
suppressPackageStartupMessages({library(data.table);library(limma);library(sandwich);library(AnnotationDbi);library(org.Hs.eg.db);library(ggplot2)})
root <- '/Users/tan/Documents/MDD-MultiOmics'
out <- '/Users/tan/Documents/Codex/2026-09-19/ple/outputs/Claude_review_2026-09-30'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
set.seed(20260930)
sc <- fread(file.path(root,'analysis_framework/results_strict/blood_OXPHOS_sample_scores.tsv'))
genes <- sort(unique(toupper(na.omit(select(org.Hs.eg.db,keys='GO:0006119',keytype='GOALL',columns='SYMBOL')$SYMBOL))))
fwrite(data.table(Gene=genes),file.path(out,'GO0006119_membership.tsv'),sep='\t')
paths <- c(GSE98793='transcriptomics/blood/GSE98793/results_fixed/GSE98793_expression_cleaned.rds',GSE39653='transcriptomics/blood/GSE39653/results_fixed/GSE39653_gene_expression.rds',GSE19738='analysis_framework/GSE19738/results/GSE19738_gene_expression.rds')
forms <- c(GSE98793='~Group+Age+Sex+Batch',GSE39653='~Group',GSE19738='~Group+Age+Sex')
cam <- list(); models <- list()
for(id in names(paths)) {
 d <- as.data.frame(sc[Cohort==id]); d$Group <- factor(d$Group,levels=c('Control','MDD'));d$Sex<-factor(d$Sex);d$Batch<-factor(d$Batch)
 x <- readRDS(file.path(root,paths[[id]])); x<-x[,match(d$Sample,colnames(x)),drop=FALSE]
 stopifnot(identical(colnames(x),d$Sample));rownames(x)<-toupper(trimws(rownames(x)))
 x<-x[!is.na(rownames(x)) & nzchar(rownames(x)) & rowSums(is.finite(x))==ncol(x),,drop=FALSE]
 if(anyDuplicated(rownames(x))) {o<-order(rownames(x),-matrixStats::rowVars(x));x<-x[o,,drop=FALSE];x<-x[!duplicated(rownames(x)),,drop=FALSE]}
 design<-model.matrix(as.formula(forms[[id]]),d); idx<-list(OXPHOS=which(rownames(x)%in%genes))
 for(cor in c(.01,NA_real_)) {
  y<-camera(x,idx,design,contrast='GroupMDD',inter.gene.cor=cor)
  cam[[paste(id,cor)]]<-data.table(Cohort=id,Method=ifelse(is.na(cor),'CAMERA estimated correlation','CAMERA fixed 0.01'),NGenes=y$NGenes,Correlation=if('Correlation'%in%names(y)) y$Correlation else cor,Direction=y$Direction,P=y$PValue)
 }
 y<-fry(x,idx,design,contrast='GroupMDD');cam[[paste0(id,'fry')]]<-data.table(Cohort=id,Method='FRY directional self contained',NGenes=y$NGenes,Correlation=NA_real_,Direction=y$Direction,P=y$PValue)
 for(f in unique(c('~Group',forms[[id]]))) {
  fit<-lm(as.formula(paste('Score_all',f)),d);j<-'GroupMDD';s<-sqrt(vcovHC(fit,type='HC3')[j,j])
  models[[paste(id,f)]]<-data.table(Cohort=id,Formula=f,N=nobs(fit),Controls=sum(d$Group=='Control'),Cases=sum(d$Group=='MDD'),Score_SD=sd(d$Score_all),Residual_SD=sigma(fit),R2=summary(fit)$r.squared,Beta=coef(fit)[j],SE_HC3=s,SE_OLS=coef(summary(fit))[j,2],P_HC3_normal=2*pnorm(-abs(coef(fit)[j]/s)),Max_leverage=max(hatvalues(fit)))
 }
}
cc<-rbindlist(cam);cc[,BH_within_method:=p.adjust(P,'BH'),by=Method];fwrite(cc,file.path(out,'CAMERA_FRY_recheck.tsv'),sep='\t');print(cc)
mm<-rbindlist(models);fwrite(mm,file.path(out,'score_SE_audit.tsv'),sep='\t');print(mm)
# Existing meta outputs are retained with their actual modified-HK implementation.
mt<-fread(file.path(root,'analysis_framework/results_strict/blood_OXPHOS_random_effects_meta.tsv'));fwrite(mt,file.path(out,'meta_verified_existing.tsv'),sep='\t')
# Residualized rank sensitivity and robust full-model regression for paired scores.
pbase<-file.path(root,'analysis_framework/GSE251780/results')
pd<-merge(fread(file.path(pbase,'GSE251780_paired_sample_scores.tsv')),fread(file.path(pbase,'GSE251780_EpiDISH_cell_fractions.tsv')),by='SampleID')
cv<-setdiff(names(fread(file.path(pbase,'GSE251780_EpiDISH_cell_fractions.tsv'))),c('SampleID','Neutro'))
rhs<-paste(c('Group','Sex',cv),collapse='+');form<-as.formula(paste('RNA_OXPHOS_all~DNAm_OXPHOS_promoter_z+',rhs))
fit<-lm(form,pd);rob<-MASS::rlm(form,pd,maxit=100)
pd$RNA_residual<-resid(lm(as.formula(paste('RNA_OXPHOS_all~',rhs)),pd));pd$DNAm_residual<-resid(lm(as.formula(paste('DNAm_OXPHOS_promoter_z~',rhs)),pd))
cs<-cor.test(pd$RNA_residual,pd$DNAm_residual,method='spearman',exact=FALSE)
j<-'DNAm_OXPHOS_promoter_z';hc<-sqrt(vcovHC(fit,type='HC3')[j,j])
res<-data.table(N=nrow(pd),OLS_beta=coef(fit)[j],HC3_SE=hc,HC3_t_P=2*pt(-abs(coef(fit)[j]/hc),df=df.residual(fit)),Residual_Spearman_rho=unname(cs$estimate),Residual_Spearman_P=cs$p.value,Huber_fullmodel_beta=coef(rob)[j])
fwrite(res,file.path(out,'paired_robust_sensitivity.tsv'),sep='\t');print(res)
p<-ggplot(pd,aes(DNAm_residual,RNA_residual))+geom_point(aes(color=Group),alpha=.7)+geom_smooth(method='lm',formula=y~x,color='black')+theme_bw(base_size=12)+labs(title='Paired blood promoter methylation and RNA pathway scores',subtitle=sprintf('n=%s; adjusted beta=%.3f; HC3 P=%.3f; residual Spearman rho=%.3f',nrow(pd),res$OLS_beta,res$HC3_t_P,res$Residual_Spearman_rho),x='Residual OXPHOS promoter methylation score',y='Residual OXPHOS RNA score')
ggsave(file.path(out,'Figure3_paired_revised.png'),p,width=8,height=5.4,dpi=220);ggsave(file.path(out,'Figure3_paired_revised.pdf'),p,width=8,height=5.4)
br<-fread(file.path(root,'analysis_framework/results/predefined_OXPHOS_across_transcriptomes.tsv'));br<-br[grepl('GSE53987|GSE54563|GSE80655',Dataset)];br[,BH_7_regions:=p.adjust(P,'BH')];fwrite(br,file.path(out,'brain_rank_BH7.tsv'),sep='\t')
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
labs<-c('GSE98793 whole blood','GSE39653 PBMC','GSE19738 whole blood','GSE251778 whole blood','Three cohort mixed specimen synthesis','Four cohort sensitivity')
fp<-data.table(Label=factor(labs,levels=rev(labs)),Beta=c(-.1259399,-.2776345,-.5380509,.106,-.243609552,-.147026804),Low=c(-.280,-.886,-1.002,-.173,-.824484207,-.548941974),High=c(.028,.331,-.074,.385,.337265102,.254888365))
p1<-ggplot(fp,aes(Beta,Label))+geom_vline(xintercept=0,linetype=2,color='grey')+geom_errorbar(aes(xmin=Low,xmax=High),width=.12,orientation='y')+geom_point(color='#215b75',size=2.5)+theme_bw(base_size=11)+labs(title='Blood pathway scores and exploratory synthesis',x='Standardized pathway score difference with 95% CI\nMDD minus control',y=NULL)
p2<-ggplot(br,aes(NES,reorder(Dataset,NES),color=BH_7_regions<.05))+geom_vline(xintercept=0,linetype=2,color='grey')+geom_point(size=3)+scale_color_manual(values=c('grey','#be4b38'))+theme_bw(base_size=11)+theme(legend.position='bottom')+labs(title='Exploratory brain OXPHOS ranks',color='BH q < 0.05 across seven regions',x='Normalized enrichment score',y=NULL)
gz<-'/Users/tan/Documents/Codex/2026-09-19/ple/outputs/MDD_Anxiety_online_followup_2026-09-23/GWAS_results/01_结构层_FUMA_MAGMA/gene_level_cross_trait_all.tsv.gz'
gt<-fread(gz);cut<--log10(.05/nrow(gt))
p4<-ggplot(gt,aes(-log10(P_MDD),-log10(P_ANX)))+geom_point(alpha=.3,size=.5,color='#78909c')+geom_vline(xintercept=cut,linetype=2)+geom_hline(yintercept=cut,linetype=2)+geom_text(data=gt[SYMBOL%in%c('DCC','SORCS3','DRD2')],aes(label=SYMBOL),nudge_y=.3,size=3)+theme_bw(base_size=11)+labs(title='Shared gene associations are contextual evidence',x='MDD MAGMA gene -log10(P)',y='Anxiety MAGMA gene -log10(P)')
for(nm in c('Figure1_scores','Figure2_brain','Figure4_MAGMA')){zz<-switch(nm,Figure1_scores=p1,Figure2_brain=p2,Figure4_MAGMA=p4);ggsave(file.path(out,paste0(nm,'.png')),zz,width=8,height=5.3,dpi=220);ggsave(file.path(out,paste0(nm,'.pdf')),zz,width=8,height=5.3)}
