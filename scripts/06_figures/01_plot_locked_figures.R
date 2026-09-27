# Publication-ready figure set v1. Visual refinement only; all biological values read from locked outputs.
# Usage: Rscript Publication_ready_v1_plot_source.R <existing output dir> [all|1|2|3|4|S1|S3]
args <- commandArgs(trailingOnly=TRUE)
ROOT <- if (length(args)>=3) args[[3]] else getwd()
OUT <- if(length(args)>=1) args[[1]] else file.path(ROOT,"figures","regenerated")
ONLY <- if(length(args)>=2) args[[2]] else "all"
if(!dir.exists(OUT)) stop("Output directory does not exist: ",OUT)
options(stringsAsFactors=FALSE)
neg <- "#345F8A"; pos <- "#B97937"; ink <- "#222B35"; muted <- "#596570"
gray <- "#D9DFE4"; bg <- "#FFFFFF"; accent <- "#476B76"; neutral <- "#F2F4F5"
times6 <- c("0.5h","4h","24h","72h","7d","28d")
times4 <- c("7d","14d","1m","2m")
programs <- c("ComplexI_assembly","Aerobic_respiration")
program_label <- c("Complex I assembly","Aerobic respiration")
cells <- c("Neurons","Astrocytes","Mature/myelin-forming oligodendroglial cells","Microglia")
short_cells <- c("Neurons","Astrocytes","Mature OLs","Microglia")
read_locked <- function(path) read.csv(file.path(ROOT,path),check.names=FALSE,fileEncoding="UTF-8-BOM")
save_figure <- function(stem,width,height,draw) {
  subdir <- if(stem %in% c("FigureS1","FigureS3")) "Supplementary_Figures" else "Main_Figures"
  for(fmt in c("svg","pdf","png","tiff")) {
    target <- file.path(OUT,subdir,paste0(stem,"_final.",fmt))
    if(fmt=="svg") svg(target,width=width,height=height,pointsize=10,bg=bg,family="sans")
    if(fmt=="pdf") pdf(target,width=width,height=height,pointsize=10,bg=bg,useDingbats=FALSE,family="Helvetica")
    if(fmt=="png") png(target,width=width*600,height=height*600,res=600,bg=bg,type="cairo")
    if(fmt=="tiff") tiff(target,width=width*600,height=height*600,res=600,bg=bg,type="cairo",compression="lzw")
    draw(); dev.off()
  }
}
panel <- function(letter,title_text,line=1,cex=1.05) title(main=paste0(letter,"  ",title_text),adj=0,font.main=2,col.main=ink,cex.main=cex,line=line)
text_box <- function(x1,y1,x2,y2,heading,body,fill="#F5F8FB",border=gray,cex=.78) {
  rect(x1,y1,x2,y2,col=fill,border=border,lwd=1.5)
  text(x1+.012,y2-.025,heading,adj=c(0,1),font=2,cex=1.0,col=ink)
  text(x1+.012,y2-.085,body,adj=c(0,1),cex=cex,col=muted,leading=1.15)
}

figure1 <- function() save_figure("Figure1",7.1,4.8,function(){
  par(mar=c(0,0,0,0),family="sans")
  plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab="",ylab="",xaxs="i",yaxs="i")
  text(.035,.955,"A  Parallel cross-resolution analyses",adj=0,font=2,cex=1.13,col=ink)
  rect(.03,.47,.305,.83,col="#F8FAFB",border=gray,lwd=1)
  text(.05,.79,"GSE5296  Whole tissue",adj=0,font=2,cex=.91,col=ink)
  text(.05,.70,"96 microarray arrays\nSCI arrays pool four mice\nTemporal / regional scores",adj=c(0,.5),cex=.75,col=muted)
  rect(.33,.45,.97,.855,col="#FBFCFC",border="#AEBAC2",lwd=1)
  text(.35,.82,"GSE234774  Shared atlas",adj=0,font=2,cex=.95,col=accent)
  rect(.35,.51,.64,.77,col="#F8FAFB",border=gray,lwd=.8)
  rect(.66,.51,.95,.77,col="#F8FAFB",border=gray,lwd=.8)
  text(.37,.735,"snRNA-seq",adj=0,font=2,cex=.86,col=ink)
  text(.37,.67,"20 libraries\nLibrary-by-cell-type pseudobulk\nCell-resolved transcription",adj=c(0,.5),cex=.75,col=muted)
  text(.68,.735,"Spatial transcriptomics",adj=0,font=2,cex=.86,col=ink)
  text(.68,.67,"9 biological samples, 36 sections\nUI, SCI 7d and SCI 2m\nSpot-level spatial context",adj=c(0,.5),cex=.75,col=muted)
  text(.65,.475,"One atlas; snRNA and spatial are not independent replication",cex=.70,col=muted)
  segments(c(.17,.49,.81),c(.47,.45,.45),c(.17,.49,.81),c(.38,.38,.38),col=accent,lwd=1)
  rect(.23,.275,.77,.38,col="#F1F4F5",border="#B9C5CB",lwd=.9)
  text(.5,.35,"Shared locked respiratory programs",font=2,cex=.92,col=ink)
  text(.5,.306,"Complex I assembly (GO:0032981)  /  Aerobic respiration (GO:0009060)",cex=.68,col=ink)
  segments(.29,.275,.20,.23,col="#9CAAB3",lwd=1,lty=2)
  rect(.03,.085,.41,.23,col="#FAFBFC",border=gray,lwd=.8)
  text(.05,.198,"GSE319931  Secondary regional analysis",adj=0,font=2,cex=.83,col=ink)
  text(.05,.145,"Four SCI animals; Above, Epicenter, Below\nrepeated within each animal",adj=c(0,.5),cex=.73,col=muted)
  text(.46,.218,"Central question",adj=0,font=2,cex=.88,col=ink)
  text(.46,.15,"Does tissue-level respiratory suppression imply\nuniform cell-type-resolved suppression?",adj=c(0,.5),cex=.81,col=ink)
  text(.5,.032,"Shared programs connect parallel datasets; no matched samples or causal sequence is implied.",cex=.69,col=muted)
})

figure2 <- function(){
  s3 <- read_locked("supplementary_source_tables/Table_S3_GSE5296_program_results_Program_results.csv")
  s3 <- s3[s3$Region=="impact (I)",]
  save_figure("Figure2",7.1,5.25,function(){
    layout(matrix(1:4,nrow=2,byrow=TRUE))
    par(oma=c(1.5,.4,.3,.2),mar=c(3.2,3.8,2.1,.7),family="sans",cex.axis=.82,cex.lab=.82)
    for(pi in 1:2) for(mi in 1:2){
      p <- programs[[pi]]; method <- c("GSE5296 ssGSEA","GSE5296 Eigengene")[[mi]]
      d <- s3[s3$Program==p & s3[["Score method"]]==method,]
      if(nrow(d)!=6) stop("Expected six locked impact timepoints")
      d <- d[match(times6,d$Time),]
      y <- as.numeric(d[["Existing SCI_minus_sham delta (mean SCI - mean sham)"]])
      if(any(!is.finite(y)) || any(y>=0)) stop("Locked impact direction mismatch")
      ylim <- range(c(0,y)); ylim <- ylim+c(-.08,.06)*diff(ylim)
      plot(1:6,y,pch=19,cex=1.1,col=neg,type="p",xaxt="n",xlab="",ylab="SCI - sham delta",ylim=ylim,bty="l")
      axis(1,at=1:6,labels=times6)
      abline(h=0,col="#B7C0C8",lty=2,lwd=.8)
      panel(LETTERS[(pi-1)*2+mi],paste(program_label[[pi]],"-",if(mi==1)"ssGSEA" else "Eigengene"),cex=.88,line=.7)
    }
    mtext("72 h primary points retained; flagged-array sensitivity not evaluable",side=1,outer=TRUE,line=.5,cex=.73,col=muted)
  })
}

direction_heat <- function(s4,program,letter,title_text,numeric_size=.72){
  d <- s4[s4$Program==program & s4$Time %in% times4 & s4[["Analysis phase"]]=="Primary",]
  plot(c(0,4),c(0,8),type="n",axes=FALSE,xlab="",ylab="",xaxs="i",yaxs="i")
  for(ci in 1:4) for(mi in 1:2) for(ti in 1:4){
    score <- c("snRNA MeanZ","snRNA Eigengene")[[mi]]
    one <- d[d[["Cell type"]]==cells[[ci]] & d$Score==score & d$Time==times4[[ti]],]
    if(nrow(one)!=1) stop("Missing/duplicate corrected S4 tile")
    value <- as.numeric(one[["Score difference (delta)"]] )
    yy <- 8-((ci-1)*2+mi)
    rect(ti-1,yy,ti,yy+1,col=if(value<0)neg else pos,border=bg,lwd=2)
    text(ti-.5,yy+.5,sprintf(if(mi==1)"%.3f" else "%.2f",value),col="white",font=2,cex=numeric_size)
  }
  axis(1,at=(1:4)-.5,labels=times4,tick=FALSE)
  axis(2,at=7.5:0.5,labels=rep(c("MeanZ","Eigengene"),4),las=1,tick=FALSE,cex.axis=.75)
  mtext(short_cells,side=2,at=c(7,5,3,1),line=5.3,las=1,cex=.79,col=ink)
  box(col=gray);panel(letter,title_text,cex=1.05)
}

figure3 <- function(){
  s4 <- read_locked("supplementary_source_tables/Table_S5_snRNA_program_results__Program_results.csv")
  scores <- read_locked("derived_data/GSE234774_snRNA/corrected_program_scores.csv")
  genes <- read_locked("supplementary_source_tables/Table_S6_Astrocyte_Aerobic_gene_level__Gene_level.csv")
  breadth <- read_locked("supplementary_source_tables/Table_S6_Astrocyte_Aerobic_gene_level__Summary_by_time.csv")
  if(nrow(genes)!=49 || !identical(as.numeric(breadth[["Positive genes"]]),c(40,39,38,33))) stop("Locked 49-gene breadth mismatch")
  save_figure("Figure3",7.1,11.4,function(){
    layout(matrix(c(1,1,1,1,2,2,3,4,5,5,5,5),nrow=3,byrow=TRUE),heights=c(2.8,2.5,6.1))
    par(oma=c(1.6,.15,1.0,.15),family="sans",cex.axis=.78,cex.lab=.78)
    par(mar=c(2.2,11.5,2.2,.5));direction_heat(s4,programs[[2]],"A","Aerobic respiration - injury - UI",numeric_size=.73)
    par(mar=c(2.2,10.5,2.2,.5));direction_heat(s4,programs[[1]],"B","Complex I assembly - injury - UI",numeric_size=.60)
    plot_library <- function(method,tag){
      d <- scores[scores$celltype=="Astrocytes" & scores$program==programs[[2]] & scores$method==method & scores$time %in% c("uninjured",times4),]
      x <- match(d$time,c("uninjured",times4))
      if(nrow(d)!=14 || anyNA(x)) stop("Astrocyte locked library score selection mismatch")
      plot(NA,xlim=c(.6,5.4),ylim=range(d$score),xaxt="n",xlab="",ylab="Raw library score",bty="l")
      axis(1,at=1:5,labels=c("UI","7d","14d","1m","2m"))
      for(i in 1:5){idx <- which(x==i);points(i+seq(-.08,.08,length.out=length(idx)),as.numeric(d$score[idx]),pch=19,col=if(i==1)ink else accent,cex=1.25)}
      panel(tag,if(method=="A")"Astrocytes - MeanZ" else "Astrocytes - Eigengene",cex=.83,line=.6)
    }
    par(mar=c(2.6,3.0,2.1,.3));plot_library("A","C1")
    par(mar=c(2.6,3.0,2.1,.3));plot_library("B","C2")
    par(mar=c(2.0,4.3,2.9,.4))
    first <- c("Ndufs4","Cox10","Ndufs1","Fxn")
    if(!all(first %in% genes$Gene)) stop("Persistent negative gene missing")
    genes <- genes[match(c(first,sort(setdiff(genes$Gene,first))),genes$Gene),]
    values <- as.matrix(genes[paste0(times4," diff_logCPM")]);storage.mode(values) <- "numeric"
    if(any(!is.finite(values)) || any(values[1:4,]>=0)) stop("Locked negative-gene values mismatch")
    lim <- max(abs(values)) # plotting-only symmetric limit; original values unchanged
    pal <- colorRampPalette(c(neg,"#F5F6F6",pos))(101)
    plot(c(0,4.85),c(0,49),type="n",axes=FALSE,xlab="",ylab="",xaxs="i",yaxs="i")
    for(i in 1:49) for(j in 1:4){
      z <- values[i,j]; ci <- round((z+lim)/(2*lim)*100)+1
      rect(j-1,49-i,j,50-i,col=pal[[ci]],border=NA)
    }
    axis(1,at=(1:4)-.5,labels=times4,tick=FALSE)
    axis(2,at=49-(1:49)+.5,labels=genes$Gene,las=1,tick=FALSE,cex.axis=.65)
    rect(0,0,4,49,border=gray,lwd=.5)
    title(main="D  Astrocyte aerobic genes",adj=0,line=1.7,cex.main=.97,col.main=ink,font.main=2)
    mtext(paste0(as.numeric(breadth[["Positive genes"]]),"/49"),side=3,at=(1:4)-.5,line=.3,cex=.72,col=ink)
    for(k in 1:100) rect(4.21,11+(k-1)*.27,4.46,11+k*.27,col=pal[[k]],border=NA)
    text(4.52,c(11,24.5,38),sprintf("%+.3f",c(-lim,0,lim)),adj=0,cex=.67,col=ink)
    text(4.20,41.5,"diff_logCPM",adj=0,cex=.65,font=2,col=ink)
    mtext("Direction colors are not cross-method effect-size scales.",side=3,outer=TRUE,line=.1,cex=.72,col=muted)
    mtext("All post-injury comparisons use the same three UI libraries.",side=1,outer=TRUE,line=.15,cex=.72,col=muted)
  })
}

figure4 <- function(){
  s6 <- read_locked("supplementary_source_tables/Table_S7_Spatial_results__Replicate_level_summary.csv")
  primary <- read_locked("supplementary_source_tables/Table_S7_Spatial_results__Primary_comparisons_sensitivity.csv")
  sensitivity <- read_locked("supplementary_source_tables/Table_S7_Spatial_results__Sensitivity.csv")
  gz <- gzfile(file.path(ROOT,"derived_data","GSE234774_spatial","spatial_program_scores.csv.gz"),open="rt")
  spots <- read.csv(gz);close(gz)
  chosen <- c("105_A_uninjured_1","105_A_7days_1","105_A_2months_1")
  if(!all(chosen %in% spots$section_id)) stop("Locked representative section absent")
  maps <- spots[spots$section_id %in% chosen & spots$method=="A",]
  pal <- colorRampPalette(c("#E8EFF1","#88A8AD","#2D5978"))(101)
  save_figure("Figure4",7.1,6.8,function(){
    layout(matrix(c(1,1,1,2,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,9,9,9,9),nrow=4,byrow=TRUE),heights=c(1.5,.76,.76,.68))
    par(oma=c(.7,.2,.25,.2),family="sans",cex.axis=.78,cex.lab=.76)
    for(pi in 1:2){
      p <- programs[[pi]]
      d <- s6[s6$Program==p & s6$Score=="Spatial Mean",]
      if(nrow(d)!=9) stop("Spatial sample summary must have 9 rows per program")
      x <- match(d$Condition,c("uninjured","7d","2m"))
      par(mar=c(2.6,3.6,2.0,.4))
      plot(NA,xlim=c(.6,3.4),ylim=range(d[["Median score"]]),xaxt="n",xlab="Condition",ylab="Biological-sample median score",bty="l")
      axis(1,at=1:3,labels=c("UI n=3","7d n=3","2m n=3"))
      for(j in 1:3){z <- d[x==j,]; yy <- as.numeric(z[["Median score"]]);points(j+seq(-.07,.07,length.out=nrow(z)),yy,pch=19,cex=1.05,col=if(pi==1)neg else pos);segments(j-.17,median(yy),j+.17,median(yy),col=ink,lwd=1.1)}
      panel(if(pi==1)"A1" else "A2",paste(program_label[[pi]],"- Spatial Mean"),line=.6,cex=.89)
    }
    map_panel <- function(p,j,letter){
      d <- maps[maps$program==p & maps$section_id==chosen[[j]],]
      if(nrow(d)==0) stop("Selected section has no spots")
      values <- maps$score[maps$program==p];lim <- range(values,finite=TRUE)
      idx <- pmax(1,pmin(101,round((d$score-lim[1])/diff(lim)*100)+1))
      par(mar=c(1.25,.3,1.35,.3))
      xr <- range(d$x);yr <- range(d$y); dx <- diff(xr)
      plot(d$x,d$y,pch=15,cex=.52,col=pal[idx],asp=1,axes=FALSE,xlab="",ylab="",xlim=c(xr[1]-dx*.03,xr[2]+dx*.25),ylim=yr)
      panel(letter,paste(program_label[[match(p,programs)]],c("UI","7d","2m")[[j]],sep=" - "),line=.1,cex=.74)
      mtext(chosen[[j]],side=1,line=.1,cex=.63,col=ink)
      if(j==3){
        bx1 <- xr[2]+dx*.08;bx2 <- xr[2]+dx*.12
        for(k in 1:100)rect(bx1,yr[1]+(k-1)/100*diff(yr),bx2,yr[1]+k/100*diff(yr),col=pal[[k]],border=NA)
        text(bx2+dx*.02,c(yr[1],mean(yr),yr[2]),sprintf("%.2f",c(lim[1],mean(lim),lim[2])),adj=0,cex=.56,col=ink)
        text(bx1,yr[2]+diff(yr)*.07,"score",adj=0,cex=.57,col=ink,xpd=TRUE)
      }
    }
    for(j in 1:3) map_panel(programs[[1]],j,if(j==1)"B" else "")
    for(j in 1:3) map_panel(programs[[2]],j,if(j==1)"C" else "")
    par(mar=c(.3,.3,1.4,.3))
    plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab="",ylab="")
    text(.01,.92,"D  Locked robustness",adj=0,font=2,cex=.88,col=ink)
    facts <- c("8/8 primary comparisons down","3/3 injured samples concordant","Low-UMI: no reversal","Leave-one-section-out: 0 reversals")
    for(i in 1:4){xx <- c(.01,.50,.01,.50)[i]; yy <- c(.60,.60,.20,.20)[i];rect(xx,yy-.035,xx+.018,yy+.035,col=neg,border=NA);text(xx+.032,yy,facts[i],adj=0,cex=.73,col=ink)}
  })
}

figureS1 <- function(){
  s3 <- read_locked("supplementary_source_tables/Table_S3_GSE5296_program_results_Program_results.csv")
  s8 <- read_locked("supplementary_source_tables/Table_S4_Quality_control_and_sensitivity_analysis_exclusions__GSE5296_flagged_GSM.csv")
  impact <- s3[s3$Region=="impact (I)",]
  save_figure("FigureS1",7.1,6.2,function(){
    layout(matrix(1:4,nrow=2,byrow=TRUE),heights=c(1,1.1))
    par(oma=c(.8,.2,.2,.2),family="sans",mar=c(3,3.6,2.0,.5),cex.axis=.78,cex.lab=.77)
    for(pi in 1:2){
      p <- programs[[pi]]
      d <- s3[s3$Program==p & s3[["Score method"]]=="GSE5296 ssGSEA",]
      plot(NA,xlim=c(1,6),ylim=range(as.numeric(d[["Existing SCI_minus_sham delta (mean SCI - mean sham)"]])),xaxt="n",xlab="Categorical time point",ylab="SCI - sham ssGSEA delta",bty="l")
      axis(1,at=1:6,labels=times6);abline(h=0,col=gray,lty=2)
      for(i in 1:3){
        region <- c("rostral (A)","impact (I)","caudal (B)")[[i]]
        z <- d[d$Region==region,];z <- z[match(times6,z$Time),]
        points(1:6+c(-.08,0,.08)[[i]],as.numeric(z[["Existing SCI_minus_sham delta (mean SCI - mean sham)"]]),pch=c(15,19,17)[[i]],col=c("#8A969E","#35424A","#B1B9BE")[[i]],cex=.95)
      }
      legend("bottomleft",legend=c("Rostral","Impact","Caudal"),pch=c(15,19,17),col=c("#8A969E","#35424A","#B1B9BE"),bty="n",cex=.69)
      panel(if(pi==1)"A" else "B",paste(program_label[[pi]],"- regional ssGSEA"),line=.5,cex=.85)
    }
    draw_grid <- function(qc=FALSE){
      rows <- expand.grid(p=programs,m=c("GSE5296 ssGSEA","GSE5296 Eigengene"),stringsAsFactors=FALSE)
      par(mar=c(2.5,9.7,2.1,.4))
      plot(c(0,6),c(0,4),type="n",axes=FALSE,xlab="",ylab="",xaxs="i",yaxs="i")
      for(i in 1:4) for(j in 1:6){
        one <- impact[impact$Program==rows$p[[i]] & impact[["Score method"]]==rows$m[[i]] & impact$Time==times6[[j]],]
        if(nrow(one)!=1) stop("S3 impact method/QC grid mismatch")
        status <- if(qc){
          if(one[["Sensitivity evaluable"]]=="No")"NE" else if(one[["Sensitivity direction_changed"]]=="TRUE")"Changed" else "Unchanged"
        } else as.character(one$Direction)
        fill <- if(status=="NE")gray else if(status=="Changed" || status=="up")pos else if(qc)accent else neg
        yy <- 4-i
        rect(j-1,yy,j,yy+1,col=fill,border=bg,lwd=2)
        text(j-.5,yy+.5,if(qc && status=="Unchanged")"U" else status,col=if(status=="NE")ink else "white",font=2,cex=.66)
      }
      axis(1,at=(1:6)-.5,labels=times6,tick=FALSE)
      axis(2,at=3.5:0.5,labels=c("Complex I ssGSEA","Aerobic ssGSEA","Complex I Eigengene","Aerobic Eigengene"),las=1,tick=FALSE,cex.axis=.78)
      box(col=gray)
      panel(if(qc)"D" else "C",if(qc)"QC - sensitivity evaluability" else "Impact-site direction",line=.5,cex=.86)
    }
    draw_grid(FALSE);draw_grid(TRUE)
    invisible(nrow(s8)) # Five flagged arrays are explained in the figure legend.
  })
}

figureS3 <- function(){
  scores <- read_locked("derived_data/GSE319931/GSE319931_program_scores.csv")
  omni <- read_locked("supplementary_source_tables/Table_S8_GSE319931_mixed_models__Omnibus.csv")
  pair <- read_locked("supplementary_source_tables/Table_S8_GSE319931_mixed_models__Pairwise.csv")
  sci <- scores[scores$condition=="SCI",]
  animals <- sort(unique(as.character(sci$animal)))
  if(!identical(animals,c("1","2","3","4"))) stop("SCI animal identities mismatch")
  colors <- c("#253A48","#52707C","#82939A","#ADB8BE")
  save_figure("FigureS3",7.1,5.4,function(){
    layout(matrix(c(1,2,5,3,4,5),nrow=2,byrow=TRUE),widths=c(1,1,1.05))
    par(oma=c(1.3,.2,.3,.2),mar=c(3,3.1,2,.4),family="sans",cex.axis=.74,cex.lab=.72)
    for(pi in 1:2) for(mi in 1:2){
      p <- programs[[pi]];m <- c("A","B")[[mi]]
      d <- sci[sci$program==p & sci$method==m,]
      if(nrow(d)!=12) stop("Expected 4 animals x 3 regions")
      plot(NA,xlim=c(1,3),ylim=range(d$score),xaxt="n",xlab="Repeated region within animal",ylab="Locked score",bty="l")
      axis(1,at=1:3,labels=c("Above","Epicenter","Below"))
      for(a in 1:4){z <- d[as.character(d$animal)==animals[[a]],];z <- z[match(c("Above","Epicenter","Below"),z$region),];lines(1:3,z$score,type="l",col=colors[[a]],lwd=1.05);points(1:3,z$score,pch=c(15,16,17,18)[a],col=colors[[a]],cex=.9)}
      panel(LETTERS[(pi-1)*2+mi],paste(program_label[[pi]],if(mi==1)"MeanZ" else "Eigengene"),line=.5,cex=.77)
    }
    par(mar=c(.4,.2,1.8,.2))
    plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab="",ylab="")
    text(.02,.98,"E  Mixed-model statistics",adj=0,font=2,cex=.91,col=ink)
    text(.02,.90,"score ~ region + (1 | animal)",adj=0,cex=.72,col=ink)
    text(.02,.84,"n = 4 SCI animals",adj=0,cex=.72,col=ink)
    legend("topleft",inset=c(.02,.23),legend=paste("Animal",animals),col=colors,lty=1,pch=c(15,16,17,18),bty="n",cex=.70,ncol=2)
    text(.03,.58,"Program / method",adj=0,font=2,cex=.67,col=ink)
    text(.56,.58,"Omnibus",adj=0,font=2,cex=.67,col=ink)
    text(.78,.58,"A-B Holm",adj=0,font=2,cex=.67,col=ink)
    for(i in 1:nrow(omni)){
      o <- omni[i,]; q <- pair[pair$Program==o$Program & pair$Score==o$Score & pair$Contrast=="A-B",]
      if(nrow(q)!=1) stop("Locked A-B Holm row missing")
      lab <- paste(if(o$Program==programs[[1]])"Complex I" else "Aerobic",if(grepl("MeanZ",o$Score))"MeanZ" else "Eigengene")
      yy <- .51-(i-1)*.11
      text(.03,yy,lab,adj=0,cex=.67,col=ink)
      text(.58,yy,sprintf("%.5f",as.numeric(o[["Omnibus region P"]])),adj=0,cex=.67,col=ink)
      text(.80,yy,sprintf("%.5f",as.numeric(q[["Holm-adjusted P"]])),adj=0,cex=.67,col=ink)
    }
    text(.03,.048,"Within injured cords only; region-matched\nuninjured controls unavailable.",adj=0,cex=.72,col=ink)
    mtext("Four paired Above-Epicenter-Below trajectories per score",side=1,outer=TRUE,line=.3,cex=.69,col=muted)
  })
}

if(ONLY %in% c("all","1")) figure1()
if(ONLY %in% c("all","2")) figure2()
if(ONLY %in% c("all","3")) figure3()
if(ONLY %in% c("all","4")) figure4()
if(ONLY %in% c("all","S1")) figureS1()
if(ONLY %in% c("all","S3")) figureS3()
