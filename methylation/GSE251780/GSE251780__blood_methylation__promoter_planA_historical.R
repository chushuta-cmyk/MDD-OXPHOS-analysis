#!/usr/bin/env Rscript
.libPaths(c('<DEV_R_LIBRARY_PATH>',.libPaths()))
suppressPackageStartupMessages({library(readxl);library(minfi);library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19);library(org.Hs.eg.db);library(data.table);library(sandwich)})
# --- Path configuration: EDIT these placeholder roots for your machine ----------
DATA_ROOT    <- "<DATA_ROOT>"      # raw GEO / supplementary files (see README "Inputs")
WORK_ROOT    <- "<WORK_ROOT>"       # author analysis workspace (intermediate result tables)
ARCHIVE_ROOT <- "<ARCHIVE_ROOT>"    # author manuscript output archive
OUTPUT_ROOT  <- "./output"          # results are written here (relative to repository root)

out<-file.path(OUTPUT_ROOT,'methylation/secondary_planA_2026-10-05');dir.create(out,recursive=TRUE,showWarnings=FALSE)
root<-file.path(WORK_ROOT,'analysis_framework/GSE251780')
write_tsv<-function(x,f)fwrite(x,file.path(out,f),sep='\t',na='NA')
mode<-commandArgs(TRUE)[1]
if(mode=='prepare') {
 mc<-as.data.frame(read_excel(file.path(out,'sources/Human.MitoCarta3.0.xls'),sheet=2,col_types='text'))
 # Fixed subunit definition: exclude assembly factors, cytochrome C and non-complex regulators.
 iscore<-grepl('OXPHOS > Complex (I|II|III|IV|V) > C(I|II|III|IV|V) subunits',mc$MitoCarta3.0_MitoPathways)
 core<-mc[iscore,c('HumanGeneID','Symbol','hg19_Chromosome','MitoCarta3.0_MitoPathways')]
 names(core)<-c('EntrezID','MitoCartaSymbol','Chromosome','Pathways')
 core$Complex<-vapply(core$Pathways,function(p)regmatches(p,regexpr('Complex (III|II|IV|I|V)',p)),character(1))
 core$Nuclear<-!core$Chromosome%in%c('MT','chrM','M')
 write_tsv(core,'MitoCarta_core_frozen.tsv')
 ann<-as.data.frame(minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19));ann$Probe<-rownames(ann)
 retained<-fread(file.path(root,'results/GSE251780_MDD_DMP_all.tsv'),select='Probe')$Probe
 ann<-ann[ann$Probe%in%retained,]
 syms<-unique(unlist(strsplit(ann$UCSC_RefGene_Name,';',fixed=TRUE)));syms<-syms[!is.na(syms)&nzchar(syms)]
 canonical<-AnnotationDbi::keys(org.Hs.eg.db,keytype='SYMBOL');known<-intersect(syms,canonical)
 direct<-AnnotationDbi::select(org.Hs.eg.db,keys=known,keytype='SYMBOL',columns='ENTREZID')
 names(direct)[1]<-'OldSymbol';direct$MappingType<-'canonical'
 unknown<-setdiff(syms,known);aliases<-intersect(unknown,AnnotationDbi::keys(org.Hs.eg.db,keytype='ALIAS'))
 alias<-AnnotationDbi::select(org.Hs.eg.db,keys=aliases,keytype='ALIAS',columns='ENTREZID');names(alias)[1]<-'OldSymbol';alias$MappingType<-'alias'
 mp<-rbind(direct,alias);mp<-unique(mp[!is.na(mp$ENTREZID),]);nk<-table(mp$OldSymbol);mp$Ambiguous<-nk[mp$OldSymbol]>1
 write_tsv(mp,'EPIC_symbol_Entrez_audit.tsv');safe<-mp[!mp$Ambiguous,];idmap<-setNames(safe$ENTREZID,safe$OldSymbol)
 go<-AnnotationDbi::select(org.Hs.eg.db,keys='GO:0006119',keytype='GOALL',columns='ENTREZID')$ENTREZID;go<-unique(na.omit(go))
 targets<-unique(c(core$EntrezID[core$Nuclear],go))
 rows<-lapply(seq_len(nrow(ann)),function(i){
   sy<-strsplit(ann$UCSC_RefGene_Name[i],';',fixed=TRUE)[[1]]; gr<-strsplit(ann$UCSC_RefGene_Group[i],';',fixed=TRUE)[[1]]
   if(length(sy)!=length(gr))return(NULL)
   ids<-unname(idmap[sy]);ok<-gr%in%c('TSS1500','TSS200','1stExon')&!is.na(ids)&ids%in%targets
   if(!any(ok))return(NULL)
   unique(data.frame(Probe=ann$Probe[i],EntrezID=ids[ok],OldSymbol=sy[ok],Region=gr[ok]))
 })
 pm<-unique(do.call(rbind,rows));write_tsv(pm,'promoter_mapping_frozen.tsv')
 sets<-rbind(data.frame(Set='GO0006119_original_ID_repaired',EntrezID=go),data.frame(Set='MitoCarta_complex_subunits_nuclear',EntrezID=core$EntrezID[core$Nuclear]))
 for(k in unique(core$Complex))sets<-rbind(sets,data.frame(Set=gsub(' ','_',k),EntrezID=core$EntrezID[core$Nuclear&core$Complex==k]))
 write_tsv(unique(sets),'analysis_sets_frozen.tsv')
 saveRDS(list(ann=ann,core=core,pm=pm,sets=sets),file.path(out,'prepared_mapping.rds'))
 cat('core total',nrow(core),'nuclear',sum(core$Nuclear),'mapped promoter probes',length(unique(pm$Probe)),'\n')
 print(table(core$Complex,core$Nuclear))
 quit(save='no')
}
stopifnot(mode=='analyse')
x<-fread(file.path(out,'target_promoter_beta.tsv'),check.names=FALSE);meta<-fread(file.path(root,'data/GSE251780_metadata.csv'));cf<-fread(file.path(root,'results/GSE251780_EpiDISH_cell_fractions.tsv'))
d<-merge(meta,cf,by='SampleID',sort=FALSE);d$SampleID<-as.character(d$SampleID);d$Group<-factor(d$Group,levels=c('Control','MDD'));d$Sex<-factor(d$Sex)
stopifnot(nrow(d)==147,!anyDuplicated(d$SampleID),all(d$SampleID%in%names(x)))
beta<-as.matrix(x[,d$SampleID,with=FALSE]);rownames(beta)<-x[[1]];storage.mode(beta)<-'double';stopifnot(all(is.finite(beta)),all(beta>=0&beta<=1))
M<-log2(pmin(pmax(beta,1e-6),1-1e-6)/(1-pmin(pmax(beta,1e-6),1-1e-6)))
obj<-readRDS(file.path(out,'prepared_mapping.rds'));pm<-obj$pm;sets<-obj$sets
ids<-unique(pm$EntrezID);geneM<-t(vapply(ids,function(id)colMeans(M[unique(pm$Probe[pm$EntrezID==id]),,drop=FALSE]),numeric(nrow(d))));rownames(geneM)<-ids
sdg<-apply(geneM,1,sd);ok<-is.finite(sdg)&sdg>0;geneM<-geneM[ok,,drop=FALSE]
Z<-t(scale(t(geneM)))
res<-list();coverage<-list();sample_scores<-list();leverage<-list();loo<-list()
for(s in unique(sets$Set)) {
 allids<-unique(sets$EntrezID[sets$Set==s]);present<-intersect(allids,rownames(Z))
 coverage[[s]]<-data.frame(Set=s,Defined_genes=length(allids),Covered_genes=length(present),Probes=length(unique(pm$Probe[pm$EntrezID%in%present])))
 if(length(present)<3)next
 score<-colMeans(Z[present,,drop=FALSE]);d$Score<-as.numeric(scale(score))
 sample_scores[[s]]<-data.frame(SampleID=d$SampleID,Set=s,Score_SD=d$Score)
 for(model in c('sex_adjusted','cell_adjusted')) {
  f<-if(model=='sex_adjusted')Score~Group+Sex else Score~Group+Sex+B+NK+CD4T+CD8T+Mono+Eosino
  fit<-lm(f,d);stopifnot(fit$rank==length(coef(fit)),nobs(fit)==147)
  b<-coef(fit)['GroupMDD'];se<-sqrt(vcovHC(fit,type='HC3')['GroupMDD','GroupMDD']);ci<-b+c(-1,1)*qt(.975,df.residual(fit))*se
  res[[length(res)+1]]<-data.frame(Set=s,Model=model,N=nobs(fit),Beta_score_SD=b,SE_HC3=se,CI_low=ci[1],CI_high=ci[2],P_HC3=2*pt(-abs(b/se),df.residual(fit)),P_OLS=coef(summary(fit))['GroupMDD','Pr(>|t|)'],Max_leverage=max(hatvalues(fit)))
  if(model=='cell_adjusted') {
   if(length(leverage)==0)leverage[[1]]<-data.frame(SampleID=d$SampleID,Group=d$Group,Sex=d$Sex,Leverage=hatvalues(fit),Eosino=d$Eosino)
   high<-which.max(hatvalues(fit));f2<-lm(f,d[-high,]);s2<-sqrt(vcovHC(f2,type='HC3')['GroupMDD','GroupMDD']);b2<-coef(f2)['GroupMDD'];ci2<-b2+c(-1,1)*qt(.975,df.residual(f2))*s2
   loo[[s]]<-data.frame(Set=s,Model='cell_adjusted_six_fractions',Excluded_SampleID=d$SampleID[high],Reason='maximum design leverage; diagnostic only, not primary exclusion',N=nobs(f2),Beta_score_SD=b2,SE_HC3=s2,CI_low=ci2[1],CI_high=ci2[2],P_HC3=2*pt(-abs(b2/s2),df.residual(f2)),Max_leverage=max(hatvalues(f2)),InferenceStatus=if(is.finite(s2))'estimable'else'HC3 non-estimable: remaining sparse-cell sample has leverage 1')
   # A stability diagnostic, not a replacement primary model: remove sparse Eosino.
   f3<-lm(Score~Group+Sex+B+NK+CD4T+CD8T+Mono,d[-high,]);s3<-sqrt(vcovHC(f3,type='HC3')['GroupMDD','GroupMDD']);b3<-coef(f3)['GroupMDD'];ci3<-b3+c(-1,1)*qt(.975,df.residual(f3))*s3
   loo[[paste0(s,'_without_Eosino')]]<-data.frame(Set=s,Model='cell_adjusted_without_sparse_Eosino',Excluded_SampleID=d$SampleID[high],Reason='exploratory stability diagnostic; differs from full cell model',N=nobs(f3),Beta_score_SD=b3,SE_HC3=s3,CI_low=ci3[1],CI_high=ci3[2],P_HC3=2*pt(-abs(b3/s3),df.residual(f3)),Max_leverage=max(hatvalues(f3)),InferenceStatus='estimable')
  }
 }
}
res<-do.call(rbind,res);res$BH_all_reported_tests<-p.adjust(res$P_HC3,'BH')
write_tsv(res,'PlanA_local_score_results.tsv');write_tsv(do.call(rbind,coverage),'PlanA_gene_probe_coverage.tsv');write_tsv(do.call(rbind,sample_scores),'PlanA_sample_scores.tsv')
write_tsv(do.call(rbind,leverage),'PlanA_cell_model_leverage.tsv');loo<-do.call(rbind,loo);loo$BH_all_14_diagnostics<-p.adjust(loo$P_HC3,'BH',n=nrow(loo));write_tsv(loo,'PlanA_high_leverage_sensitivity.tsv')
pb<-data.frame(Probe=rownames(beta),Mean_beta=rowMeans(beta),SD_beta=apply(beta,1,sd),Fraction_beta_below_005=rowMeans(beta<.05));write_tsv(pb,'promoter_beta_distribution.tsv')
png(file.path(out,'PlanA_local_score_forest.png'),width=1400,height=1000,res=150);par(mar=c(5,17,3,1));yy<-rev(seq_len(nrow(res)))
plot(res$Beta_score_SD,yy,xlim=range(res$CI_low,res$CI_high),yaxt='n',ylab='',xlab='Adjusted score difference in within-cohort SD (MDD - control)',pch=19,main='GSE251780 exploratory measurement sensitivity')
axis(2,at=yy,labels=paste(res$Set,res$Model),las=1,cex.axis=.65);segments(res$CI_low,yy,res$CI_high,yy);abline(v=0,lty=2,col='grey');dev.off()
capture.output(sessionInfo(),file=file.path(out,'sessionInfo_PlanA.txt'));print(res)
