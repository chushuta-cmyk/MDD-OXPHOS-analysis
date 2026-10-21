options(warn=1)
suppressPackageStartupMessages({library(limma);library(edgeR);library(org.Hs.eg.db)})

# --- Path configuration: EDIT these placeholder roots for your machine ----------
DATA_ROOT    <- "<DATA_ROOT>"      # raw GEO / supplementary files (see README "Inputs")
WORK_ROOT    <- "<WORK_ROOT>"       # author analysis workspace (intermediate result tables)
ARCHIVE_ROOT <- "<ARCHIVE_ROOT>"    # author manuscript output archive
OUTPUT_ROOT  <- "./output"          # results are written here (relative to repository root)

ROOT <- WORK_ROOT
OUT <- file.path(OUTPUT_ROOT, "brain/brain_boundary")
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
genes <- unique(read.delim(file.path(ARCHIVE_ROOT, "outputs/Claude_review_2026-09-30/GO0006119_membership.tsv"))$Gene)
writeLines(genes,file.path(OUT,'frozen_GO0006119.txt'))
results <- manifests <- samples <- failures <- sources <- list()
read_src <- function(path,rds=TRUE) {
  stopifnot(file.exists(path))
  sources[[length(sources)+1]] <<- data.frame(Path=path,Bytes=file.info(path)$size,MD5=unname(tools::md5sum(path)))
  if(rds)readRDS(path) else read.csv(path,check.names=FALSE)
}
clean <- function(E) {
  E <- as.matrix(E);storage.mode(E)<-'double'
  rownames(E)<-toupper(trimws(rownames(E)))
  keep<-nzchar(rownames(E)) & apply(E,1,function(z)all(is.finite(z)))
  E<-E[keep,,drop=FALSE]; vv<-apply(E,1,var)
  E<-E[order(rownames(E),-vv),,drop=FALSE]
  E<-E[!duplicated(rownames(E)),,drop=FALSE]
  E[apply(E,1,var)>0,,drop=FALSE]
}
test <- function(E,m,formula,cohort,region,family,technology='microarray',contrasts=NULL,note='') {
  vars<-all.vars(formula);valid<-complete.cases(m[,vars,drop=FALSE])
  m<-droplevels(m[valid,,drop=FALSE]);stopifnot(!anyDuplicated(m$SampleID),all(m$SampleID%in%colnames(E)))
  E<-E[,m$SampleID,drop=FALSE];X<-model.matrix(formula,m)
  stopifnot(qr(X)$rank==ncol(X),nrow(X)>ncol(X)+2)
  saveRDS(list(metadata=m,design=X),file.path(OUT,paste0(cohort,'_',gsub('[^A-Za-z0-9_]','_',region),'_model_inputs.rds')))
  if(is.null(contrasts)) {
    cid<-match('DiagnosisMDD',colnames(X));stopifnot(!is.na(cid))
    C<-matrix(0,nrow=ncol(X),ncol=1,dimnames=list(colnames(X),'MDD_vs_Control'));C[cid,1]<-1
  } else C<-makeContrasts(contrasts=contrasts,levels=X)
  if(technology=='RNAseq') {
    y<-DGEList(E);y<-y[filterByExpr(y,design=X),,keep.lib.sizes=FALSE];y<-calcNormFactors(y)
    inp<-voom(y,X,plot=FALSE);expr<-inp$E
  } else {inp<-E;expr<-E}
  ix<-which(rownames(expr)%in%genes);stopifnot(length(ix)>=5)
  fit<-eBayes(contrasts.fit(lmFit(inp,X),C),robust=TRUE,trend=technology!='RNAseq')
  for(j in seq_len(ncol(C))) {
    id<-paste(cohort,region,colnames(C)[j],sep='__');cat('RUN',id,'N',nrow(m),'GO genes',length(ix),'\n')
    ca<-camera(inp,list(OXPHOS=ix),X,contrast=C[,j],inter.gene.cor=NA)
    cf<-camera(inp,list(OXPHOS=ix),X,contrast=C[,j],inter.gene.cor=.01)
    fr<-fry(inp,list(OXPHOS=ix),X,contrast=C[,j])
    results[[length(results)+1]]<<-data.frame(ID=id,Cohort=cohort,Region=region,Contrast=colnames(C)[j],Family=family,N=nrow(m),NGenes=length(ix),Direction=fr$Direction,Median_coefficient=median(fit$coefficients[ix,j]),CAMERA_est_P=ca$PValue,Correlation=ca$Correlation,FRY_P=fr$PValue,FRY_mixed_P=fr$PValue.Mixed,CAMERA_fixed_P=cf$PValue,Note=note)
  }
  id<-paste(cohort,region,sep='__')
  manifests[[length(manifests)+1]]<<-data.frame(Cohort=cohort,Region=region,Family=family,N=nrow(m),Excluded_missing_covariates=sum(!valid),Design_columns=ncol(X),Rank=qr(X)$rank,Model=paste(deparse(formula),collapse=''),Technology=technology,Note=note)
  samples[[length(samples)+1]]<<-data.frame(Cohort=cohort,Region=region,SampleID=m$SampleID,Group=if('Diagnosis'%in%names(m))as.character(m$Diagnosis)else if('ClinicalGroup'%in%names(m))as.character(m$ClinicalGroup)else if('Group'%in%names(m))as.character(m$Group)else if('AnxietyGroup'%in%names(m))as.character(m$AnxietyGroup)else 'continuous')
  write.table(do.call(rbind,results),file.path(OUT,'tests_checkpoint.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
}
safe <- function(label,code)tryCatch(force(code),error=function(e){cat('FAIL',label,conditionMessage(e),'\n');failures[[length(failures)+1]]<<-data.frame(Analysis=label,Error=conditionMessage(e))})
safe('GSE53987',{
 p<-file.path(ROOT,'transcriptomics/brain/GSE53987/results_fixed/MOFA_inputs')
 E<-clean(read_src(file.path(p,'GSE53987_expression_cleaned.rds')));m<-read_src(file.path(p,'GSE53987_metadata.rds'))
 m$SampleID<-as.character(m$GSM_ID);m$Diagnosis<-factor(m$Diagnosis,levels=c('Control','MDD'))
 for(r in unique(m$Brain_Area))safe(paste('53987',r),{
   d<-m[m$Brain_Area==r,,drop=FALSE];d$Tetrad_f<-factor(d$Tetrad)
   tab<-table(d$Tetrad_f,d$Diagnosis);paired<-rownames(tab)[rowSums(tab>0)==2];d<-droplevels(d[d$Tetrad_f%in%paired,,drop=FALSE])
   for(v in c('Age','PMI','pH'))d[[v]]<-as.numeric(as.character(d[[v]]))
   test(E,d,~Diagnosis+Age+PMI+pH+Tetrad_f,'GSE53987',r,'Brain_regions',note='Matched tetrads; available age/PMI/pH adjusted. RIN not included in historical model.')
 })
})
safe('GSE54563',{
 p<-file.path(ROOT,'transcriptomics/brain/GSE54563/results_fixed/MOFA_inputs')
 E<-clean(read_src(file.path(p,'GSE54563_expression_cleaned.rds')));m<-read_src(file.path(p,'GSE54563_metadata.rds'))
 m$SampleID<-as.character(m$sample_id);m$Diagnosis<-factor(m$disease,levels=c('Control','MDD'));m$pair<-factor(m$pair)
 test(E,m,~Diagnosis+pair,'GSE54563','ACC','Brain_regions',note='Paired model; other postmortem covariates absent from saved metadata.')
})
safe('GSE80655',{
 path<-file.path(DATA_ROOT, "Psychology_study/datasets/GSE80655_raw_counts_GRCh38.p13_NCBI.tsv")
 sources[[length(sources)+1]]<-data.frame(Path=path,Bytes=file.info(path)$size,MD5=unname(tools::md5sum(path)))
 E<-as.matrix(read.delim(path,row.names=1,check.names=FALSE));stopifnot(all(E>=0),all(E==round(E)))
 mp<-AnnotationDbi::mapIds(org.Hs.eg.db,keys=rownames(E),keytype='ENTREZID',column='SYMBOL',multiVals='first')
 k<-!is.na(mp);E<-E[k,,drop=FALSE];sym<-toupper(unname(mp[k]));o<-order(sym,-rowSums(E));E<-E[o,,drop=FALSE];sym<-sym[o];E<-E[!duplicated(sym),,drop=FALSE];rownames(E)<-sym[!duplicated(sym)]
 m<-read_src(file.path(ROOT,'transcriptomics/brain/GSE80655/results_fixed/MOFA_inputs/GSE80655_metadata.rds'))
 m$SampleID<-as.character(m$GSM_ID);m$Diagnosis<-factor(m$Diagnosis,levels=c('Control','MDD'))
 for(r in unique(m$Brain_Region))safe(paste('80655',r),test(E,m[m$Brain_Region==r,,drop=FALSE],~Diagnosis,'GSE80655',r,'Brain_regions','RNAseq',note='NEW TMM/voom weighted reanalysis of raw counts; not numerically equivalent to historical DESeq2/prerank. No postmortem covariates in saved metadata.'))
})
safe('GSE98793_GAD',{
 p<-file.path(ROOT,'transcriptomics/blood/GSE98793/results_fixed')
 E<-clean(read_src(file.path(p,'GSE98793_expression_cleaned.rds')));m<-read_src(file.path(p,'GSE98793_metadata.rds'));m<-m[colnames(E),];m$SampleID<-colnames(E)
 iscase<-grepl('case',m$source_name_ch1);gad<-as.character(m$GAD)%in%c('1','TRUE')
 m$ClinicalGroup<-factor(ifelse(!iscase,'Control',ifelse(gad,'MDD_GAD','MDD_only')),levels=c('Control','MDD_only','MDD_GAD'))
 stopifnot(all(table(m$ClinicalGroup)==64));m$Age<-as.numeric(m$age);m$Sex<-factor(m$gender);m$Batch<-factor(m$batch)
 test(E,m,~0+ClinicalGroup+Age+Sex+Batch,'GSE98793','GAD_subgroups','GAD_contrasts',contrasts=c('ClinicalGroupMDD_GAD-ClinicalGroupControl','ClinicalGroupMDD_only-ClinicalGroupControl','ClinicalGroupMDD_GAD-ClinicalGroupMDD_only'),note='No pure GAD group; within-MDD comparison is not an independent GAD effect.')
})
safe('GSE61672',{
 p<-file.path(ROOT,'analysis_framework/GSE61672/data');E<-clean(read_src(file.path(p,'GSE61672_gene_expression.rds')));m<-read_src(file.path(p,'GSE61672_metadata.rds'));m$SampleID<-as.character(m$Sample)
 for(v in c('GAD7','Age','BMI','RIN'))m[[paste0(v,'_z')]]<-as.numeric(scale(m[[v]]))
 m$AnxietyGroup<-factor(m$AnxietyGroup,levels=c('Control','Case'));m$Ethnicity3<-factor(m$Ethnicity3);m$Batch<-factor(m$Batch)
 test(E,m,~GAD7_z+Age_z+BMI_z+Ethnicity3+Batch+RIN_z,'GSE61672','GAD7','Anxiety_boundary',contrasts='GAD7_z',note='Original covariate specification, continuous symptom coefficient.')
 test(E,m,~AnxietyGroup+Age_z+BMI_z+Ethnicity3+Batch+RIN_z,'GSE61672','CaseControl','Anxiety_boundary',contrasts='AnxietyGroupCase',note='Original case/control covariate specification.')
})
safe('GSE228702',{
 p<-file.path(ROOT,'analysis_framework/GSE228702/data');path<-file.path(p,'GSE228702_adjusted_expression_Cellcounts_granulatorAbis0_nnls.csv.gz')
 sources[[length(sources)+1]]<-data.frame(Path=path,Bytes=file.info(path)$size,MD5=unname(tools::md5sum(path)))
 E<-read.csv(gzfile(path),row.names=1,check.names=FALSE);colnames(E)<-sub('^X(?=[0-9])','',colnames(E),perl=TRUE);E<-clean(E)
 E<-t(apply(E,1,function(z)qnorm((rank(z,ties.method='average')-.5)/length(z))))
 m<-read_src(file.path(p,'GSE228702_metadata_parsed.csv'),FALSE);m$SampleID<-m$Sample;m$Group<-factor(m$Group,levels=c('noSAD_noELA','noSAD_ELA','SAD_noELA','SAD_ELA'));m$Gender<-factor(m$Gender)
 test(E,m,~0+Group+Age_z+Gender,'GSE228702','SAD_ELA','SAD_ELA_boundary',contrasts=c('(GroupSAD_noELA+GroupSAD_ELA-GroupnoSAD_noELA-GroupnoSAD_ELA)/2','(GroupnoSAD_ELA+GroupSAD_ELA-GroupnoSAD_noELA-GroupSAD_noELA)/2','GroupSAD_ELA-GroupSAD_noELA-GroupnoSAD_ELA+GroupnoSAD_noELA'),note='Deposited cell-adjusted processed values; phenotype-blind within-gene rank INT, NOT raw-count voom.')
})
res<-do.call(rbind,results)
for(v in c('CAMERA_est_P','FRY_P','FRY_mixed_P','CAMERA_fixed_P')) {
 res[[paste0(v,'_BH_family')]]<-ave(res[[v]],res$Family,FUN=function(p)p.adjust(p,'BH'))
 res[[paste0(v,'_BH_all')]]<-p.adjust(res[[v]],'BH')
}
write.table(res,file.path(OUT,'correlation_aware_tests.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
for(n in c('manifests','samples','sources','failures'))if(length(get(n)))write.table(do.call(rbind,get(n)),file.path(OUT,paste0(n,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(OUT,'sessionInfo.txt'))
print(res[,c('Cohort','Region','N','NGenes','Direction','CAMERA_est_P','FRY_P','FRY_P_BH_family')]);cat('DONE\n')
