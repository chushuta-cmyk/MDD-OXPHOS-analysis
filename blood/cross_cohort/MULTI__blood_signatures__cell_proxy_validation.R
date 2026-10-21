options(warn=1)
.libPaths(c('/Users/tan/Developer/envs/R/lib/R/library',.libPaths()))
suppressPackageStartupMessages({library(limma);library(edgeR);library(org.Hs.eg.db);library(matrixStats);library(sandwich)})
ROOT <- '/Users/tan/Documents/MDD-MultiOmics'
OLD <- file.path(ROOT,'analysis_framework/sensitivity_2026-10-04/erady_validation')
OUT <- '/Users/tan/Documents/Codex/2026-09-19/ple/outputs/secondary_completion_2026-10-05/signatures'
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
set.seed(20261005)
gs <- read.delim(file.path(OLD,'gene_sets.tsv'))
canonical <- AnnotationDbi::keys(org.Hs.eg.db,keytype='SYMBOL')
all_alias <- AnnotationDbi::keys(org.Hs.eg.db,keytype='ALIAS')
mapping <- list()
map_symbols <- function(g,label) {
  targets <- setdiff(unique(g),canonical)
  av <- intersect(targets,all_alias)
  mp <- if(length(av)) suppressMessages(AnnotationDbi::select(org.Hs.eg.db,av,'SYMBOL','ALIAS')) else data.frame(ALIAS=character(),SYMBOL=character())
  choices <- split(mp$SYMBOL,mp$ALIAS)
  translated <- vapply(g,function(a){if(a%in%canonical)return(a); z<-unique(na.omit(choices[[a]]));if(length(z)==1)z else NA_character_},character(1))
  status <- ifelse(g%in%canonical,'canonical',ifelse(is.na(translated),'unresolved_or_ambiguous','unique_alias'))
  mapping[[length(mapping)+1]] <<- data.frame(Source=label,Original=g,Symbol=translated,Status=status)
  translated
}
gs$gene <- map_symbols(gs$gene,'signature_input')
gs <- unique(gs[!is.na(gs$gene),])
sets <- split(gs$gene,gs$set)
# Conservative somatic canonical reconstruction. Keep previous set as sensitivity.
sets$Histone_RD_legacy <- sets$Histone_RD_canonical
rd <- sets$Histone_RD_canonical
sets$Histone_RD_canonical <- rd[grepl('^H1-[1-5]$|^H2AC[0-9]+$|^H2BC[0-9]+L?$|^H3C[0-9]+$|^H4C[0-9]+$',rd)]
sets$Histone_RI_variant <- union(sets$Histone_RI_variant,'H1-10')
write.table(data.frame(Gene=setdiff(rd,sets$Histone_RD_canonical),Action='excluded_from_conservative_RD_not_automatically_RI'),file.path(OUT,'RD_classification_changes.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
frozen <- do.call(rbind,lapply(names(sets),function(s)data.frame(Signature=s,Gene=sets[[s]])))
write.table(frozen,file.path(OUT,'frozen_local_signatures.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(c('Frozen 2026-10-05 after prior data inspection; this is a correction analysis, not preregistration.',
 'Primary signatures: GO0006119 OXPHOS and conservative reconstructed canonical RD histones.',
 'Primary strata: overall in four cohorts; primary test: weighted FRY, BH over all eligible cohort x primary-signature tests.',
 'CAMERA estimated correlation, mean-z OLS/HC3, and matched-background resampling are secondary; separate BH families.',
 'Full signature panel, additional covariate models, legacy RD and sex strata are exploratory; BH within run and test/model.',
 'Original author gene lists not obtained; no exact Erady replication claim.',
 'Matched-background resampling: two-sided, 1000 draws, 0.1 expression bins, unique substitute genes, exclude all tested signatures; expand pools up to +/-4 log expression when needed.',
 'Marker adjustment uses explicitly named median log-expression columns; these are proxy scores, not proportions.',
 'No meta-analysis treats sex strata or models as independent cohorts.'),file.path(OUT,'analysis_plan_frozen.txt'))
collapse_expression <- function(E,label) {
  sym <- map_symbols(rownames(E),label)
  good <- !is.na(sym)&apply(E,1,function(x)all(is.finite(x)))
  E<-as.matrix(E[good,,drop=FALSE]);sym<-sym[good]
  # Selection independent of disease test statistics.
  o<-order(sym,-rowMeans(E));E<-E[o,,drop=FALSE];sym<-sym[o]
  E<-E[!duplicated(sym),,drop=FALSE];rownames(E)<-sym[!duplicated(sym)]
  E
}
hc3 <- function(fit) {
  stopifnot(fit$rank==length(coef(fit)))
  b<-coef(fit)['GroupMDD'];se<-sqrt(vcovHC(fit,type='HC3')['GroupMDD','GroupMDD']);df<-df.residual(fit)
  data.frame(Beta=unname(b),SE_HC3=se,P_HC3=2*pt(-abs(b/se),df),CI_low=b-qt(.975,df)*se,CI_high=b+qt(.975,df)*se,P_OLS=coef(summary(fit))['GroupMDD',4],N=nobs(fit),Max_leverage=max(hatvalues(fit)),Design_condition=kappa(model.matrix(fit)))
}
results<-scores<-coverage<-samples<-errors<-list()
run_one <- function(E,mm,rhs,cohort,stratum,technology,primary=FALSE) {
  # GEO metadata may have repeated descriptive characteristic columns.
  # Required model variables retain their names; only duplicate names are suffixed.
  names(mm) <- make.unique(names(mm))
  rownames(mm)<-as.character(mm$SampleID)
  stopifnot(!anyDuplicated(mm$SampleID),all(rownames(mm)%in%colnames(E)))
  E<-E[,rownames(mm),drop=FALSE]
  X<-model.matrix(as.formula(paste('~',rhs)),mm)
  saveRDS(list(metadata=mm,design=X),file.path(OUT,paste0(cohort,'_',stratum,'_model_inputs.rds')))
  stopifnot(nrow(X)==ncol(E),qr(X)$rank==ncol(X),sum(mm$Group=='MDD')>1,sum(mm$Group=='Control')>1)
  cid<-match('GroupMDD',colnames(X)); stopifnot(!is.na(cid))
  if(technology=='RNAseq') {
    y<-DGEList(E);keep<-filterByExpr(y,design=X);y<-calcNormFactors(y[keep,,keep.lib.sizes=FALSE]);v<-voom(y,X,plot=FALSE);expr<-v$E; fit<-eBayes(lmFit(v,X));inp<-v
  } else {expr<-E;fit<-eBayes(lmFit(expr,X),robust=TRUE,trend=TRUE);inp<-expr}
  id<-paste(cohort,stratum,sep='__');cat('\nRUN',id,'N',nrow(mm),'genes',nrow(expr),'\n')
  samples[[id]]<<-data.frame(Cohort=cohort,Stratum=stratum,N=nrow(mm),Cases=sum(mm$Group=='MDD'),Controls=sum(mm$Group=='Control'),Technology=technology,Model=rhs,Primary=primary)
  idx<-lapply(sets,function(g)which(rownames(expr)%in%g))
  coverage[[id]]<<-data.frame(Cohort=cohort,Stratum=stratum,Signature=names(sets),Defined=sapply(sets,length),Covered=sapply(idx,length))
  idx<-idx[sapply(idx,length)>=5]
  idx<-idx[intersect(names(idx),c('OXPHOS_GO0006119','Histone_RD_canonical'))]
  ca<-camera(inp,idx,X,contrast=cid,inter.gene.cor=NA)
  cf<-camera(inp,idx,X,contrast=cid,inter.gene.cor=.01)
  fr<-fry(inp,idx,X,contrast=cid)
  lfc<-setNames(fit$coefficients[,'GroupMDD'],rownames(expr))
  baseline<-rowMeans(expr[,mm$Group=='Control',drop=FALSE])
  excluded<-unique(unlist(sets));background<-setdiff(rownames(expr),excluded)
  stopifnot(length(background)>100)
  bins<-floor(baseline*10);names(bins)<-rownames(expr)
  bg_by_bin<-split(background,bins[background])
  Z<-t(scale(t(expr)))
  marker<-data.frame(row.names=rownames(mm))
  for(s in grep('^MCP_',names(sets),value=TRUE)) {
    ii<-which(rownames(expr)%in%sets[[s]])
    if(length(ii)>=3)marker[[paste0('CELL_',s)]]<-matrixStats::colMedians(expr[ii,,drop=FALSE])
  }
  stopifnot(!anyDuplicated(names(marker)))
  wanted<-paste0('CELL_MCP_',c('T_cells','Cytotoxic_lymphocytes','NK_cells','B_lineage','Monocytic_lineage'))
  if(cohort!='GSE39653')wanted<-c(wanted,'CELL_MCP_Neutrophils','CELL_MCP_Myeloid_dendritic_cells')
  ct<-intersect(wanted,names(marker));ct<-ct[vapply(marker[ct],function(x)sd(x)>1e-8,logical(1))]
  d<-cbind(mm,marker);stopifnot(!anyDuplicated(names(d)))
  for(s in names(idx)) {
    g<-rownames(expr)[idx[[s]]];obs<-median(lfc[g]);B<-1000
    groups<-split(g,bins[g]); pools<-list()
    for(k in names(groups)) {
      pool<-bg_by_bin[[k]];if(is.null(pool))pool<-character()
      if(length(pool)<length(groups[[k]])) {
        pool<-background[abs(baseline[background]-median(baseline[groups[[k]]]))<=4]
      }
      if(length(pool)<length(groups[[k]]))stop('Insufficient matched background for ',s)
      pools[[k]]<-pool
    }
    null<-rep(NA_real_,B)
    for(b in seq_len(B)) {
      chosen<-character()
      for(k in names(groups)) {
        avail<-setdiff(pools[[k]],chosen);n<-length(groups[[k]])
        if(length(avail)<n){chosen<-NULL;break}
        chosen<-c(chosen,sample(avail,n,replace=FALSE))
      }
      if(length(chosen)==length(g))null[b]<-median(lfc[chosen])
    }
    valid<-sum(is.finite(null));stopifnot(valid>=950)
    # Equal-tail two-sided probability; no tail selection based on observed sign.
    pleft<-(1+sum(null<=obs,na.rm=TRUE))/(valid+1);pright<-(1+sum(null>=obs,na.rm=TRUE))/(valid+1)
    results[[length(results)+1]]<<-data.frame(Cohort=cohort,Stratum=stratum,Primary_run=primary,Signature=s,NGenes=length(g),Median_logFC=obs,FRY_P=fr[s,'PValue'],FRY_direction=fr[s,'Direction'],CAMERA_P=ca[s,'PValue'],Correlation=ca[s,'Correlation'],CAMERA_fixed_P=cf[s,'PValue'],Background_P=min(1,2*min(pleft,pright)),Background_draws=valid,Background_low=quantile(null,.025,na.rm=TRUE),Background_high=quantile(null,.975,na.rm=TRUE))
    d$Score<-colMeans(Z[idx[[s]],,drop=FALSE]);stopifnot(all(is.finite(d$Score)))
    for(adj in c(FALSE,TRUE)) {
      if(adj && (startsWith(s,'MCP_')||length(ct)<3))next
      form<-as.formula(paste('Score ~',rhs,if(adj)paste('+',paste(ct,collapse='+'))else''))
      fit_sc<-lm(form,d)
      if(fit_sc$rank!=length(coef(fit_sc))) {
        errors[[length(errors)+1]]<<-data.frame(Run=id,Signature=s,Model='marker_adjusted',Issue='rank_deficient');next
      }
      z<-hc3(fit_sc);z$Cohort<-cohort;z$Stratum<-stratum;z$Primary_run<-primary;z$Signature<-s;z$Model<-if(adj)'marker_proxy_adjusted'else'base_covariates';z$Score_SD<-sd(d$Score)
      scores[[length(scores)+1]]<<-z
    }
  }
  write.table(do.call(rbind,results),file.path(OUT,'set_tests_checkpoint.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
}
# Microarrays
for(cohort in c('GSE98793','GSE39653','GSE19738')) {
  if(cohort=='GSE98793') {
    p<-file.path(ROOT,'transcriptomics/blood/GSE98793/results_fixed')
    E<-readRDS(file.path(p,'GSE98793_expression_cleaned.rds'));m<-readRDS(file.path(p,'GSE98793_metadata.rds'));m<-m[colnames(E),];m$SampleID<-colnames(E)
    m$Group<-factor(ifelse(grepl('case',m$source_name_ch1),'MDD','Control'),levels=c('Control','MDD'))
    m$Sex<-factor(ifelse(grepl('F',m$gender),'Female','Male'));m$Age<-as.numeric(m$age);m$Batch<-factor(m$batch);rhs<-'Group+Age+Sex+Batch'
  } else if(cohort=='GSE39653') {
    p<-file.path(ROOT,'transcriptomics/blood/GSE39653/results_fixed');E<-readRDS(file.path(p,'GSE39653_gene_expression.rds'));m<-readRDS(file.path(p,'GSE39653_metadata.rds'));m<-m[colnames(E),];m$SampleID<-colnames(E);m$Group<-factor(m$Diagnosis,levels=c('Control','MDD'));rhs<-'Group'
  } else {
    p<-file.path(ROOT,'analysis_framework/GSE19738/results');E<-readRDS(file.path(p,'GSE19738_gene_expression.rds'));m<-readRDS(file.path(p,'GSE19738_metadata.rds'));m<-m[colnames(E),];m$SampleID<-colnames(E);m$Group<-factor(ifelse(m$MDD==1,'MDD','Control'),levels=c('Control','MDD'));m$Sex<-factor(tolower(m$Sex));m$Age<-as.numeric(m$Age);rhs<-'Group+Age+Sex'
  }
  E<-collapse_expression(E,cohort)
  run_one(E,m,rhs,cohort,'overall','microarray',TRUE)
  if(FALSE)for(sx in levels(m$Sex)){
    sub<-droplevels(m[m$Sex==sx,,drop=FALSE]);r<-gsub('+Sex','',rhs,fixed=TRUE);run_one(E,sub,r,cohort,tolower(sx),'microarray',FALSE)
  }
}
# RNAseq
p<-file.path(ROOT,'transcriptomics/blood/GSE251778')
cnt<-read.csv(gzfile(file.path(p,'data/GSE251778_Mokhtari_rna_counts.csv.gz')),row.names=1,check.names=FALSE)
ens<-cnt$gene_ID;cnt<-cnt[,-1,drop=FALSE];colnames(cnt)<-sub('_rna$','',colnames(cnt))
emap<-read.delim(file.path(ROOT,'analysis_framework/sensitivity_2026-10-04/ENSG_to_SYMBOL.tsv'))
rownames(cnt)<-make.unique(ifelse(is.na(emap$symbol[match(ens,emap$ensembl_id)]),paste0('UNMAPPED_',ens),emap$symbol[match(ens,emap$ensembl_id)]))
# Avoid make.unique aliases: collapse true symbols explicitly.
sym<-map_symbols(emap$symbol[match(ens,emap$ensembl_id)],'GSE251778');keep<-!is.na(sym);cnt<-as.matrix(cnt[keep,,drop=FALSE]);sym<-sym[keep];o<-order(sym,-rowMeans(cnt));cnt<-cnt[o,,drop=FALSE];sym<-sym[o];cnt<-cnt[!duplicated(sym),,drop=FALSE];rownames(cnt)<-sym[!duplicated(sym)]
m<-read.csv(file.path(p,'results_fixed/GSE251778_metadata_aligned.csv'));m<-m[match(colnames(cnt),as.character(m$SampleID)),];m$SampleID<-as.character(m$SampleID);m$Group<-factor(m$Group,levels=c('Control','MDD'));m$Sex<-factor(m$Sex)
run_one(cnt,m,'Group+Sex','GSE251778','overall','RNAseq',TRUE)
if(FALSE)for(sx in levels(m$Sex))run_one(cnt,droplevels(m[m$Sex==sx,,drop=FALSE]),'Group','GSE251778',tolower(sx),'RNAseq',FALSE)
res<-do.call(rbind,results);sc<-do.call(rbind,scores)
primary<-res$Primary_run & res$Signature%in%c('OXPHOS_GO0006119','Histone_RD_canonical')
for(col in c('FRY_P','CAMERA_P','Background_P')) {
  res[[paste0(col,'_q_primary')]]<-NA_real_;res[[paste0(col,'_q_primary')]][primary]<-p.adjust(res[[col]][primary],'BH')
  res[[paste0(col,'_q_panel')]]<-ave(res[[col]],interaction(res$Cohort,res$Stratum),FUN=function(p)p.adjust(p,'BH'))
}
sc$P_HC3_q_panel<-ave(sc$P_HC3,interaction(sc$Cohort,sc$Stratum,sc$Model),FUN=function(p)p.adjust(p,'BH'))
ix<-sc$Primary_run&sc$Signature%in%c('OXPHOS_GO0006119','Histone_RD_canonical')&sc$Model=='base_covariates'
sc$P_HC3_q_primary<-NA_real_;sc$P_HC3_q_primary[ix]<-p.adjust(sc$P_HC3[ix],'BH')
sc$P_HC3_q_across_cohorts_by_model<-ave(sc$P_HC3,sc$Model,FUN=function(p)p.adjust(p,'BH'))
write.table(res,file.path(OUT,'corrected_set_tests.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(sc,file.path(OUT,'corrected_sample_scores.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,coverage),file.path(OUT,'signature_coverage.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,samples),file.path(OUT,'sample_manifest.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(unique(do.call(rbind,mapping)),file.path(OUT,'gene_symbol_mapping_audit.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
if(length(errors))write.table(do.call(rbind,errors),file.path(OUT,'model_failures.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(OUT,'sessionInfo.txt'))
cat('\nPRIMARY RESULTS\n');print(res[primary,c('Cohort','Signature','NGenes','Median_logFC','FRY_P','FRY_P_q_primary','CAMERA_P','Background_P')]);cat('\nDONE\n')
