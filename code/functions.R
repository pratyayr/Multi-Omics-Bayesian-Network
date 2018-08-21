
##### Code for user-defined functions for the analysis #####




###### The required Fsher (and related other) 
####  transformation functions ######

fisher=function(r,n){
  z=myforward(r)
  se=1/sqrt(n-3)
  return(z/se)
}

fisherinv=function(z,n){
  se=1/sqrt(n-3)
  z=z*se
  return(myinv(z))
}

myforward=function(s){
  w=.5*log((1+s)/(1-s))
  return(w)}

myinv=function(w){
  s=(exp(2*w)-1)/(exp(2*w)+1)
  return(s)}

logit=function(x){return(log(x/(1-x)))}
expit=function(x){return(exp(x)/(1+exp(x)))}





#### Function for gene-eqtl analysis using Z-test on 
### Fisher transformed correlations ####

geneeqtl.fisher=function(exp,geno){
  n=dim(exp)[2]
  cormat=cor(t(exp),t(geno))
  z=fisher(cormat,n)
  pval=pchisq(z^2,1,lower.tail=F)
  return(pval) #Returns #genes x #SNPs p-value matrix
}



#### Function to generate blacklist matrix for Bayesian Network Analysis (BNA) ###
blackmatgen=function(n1,n2,n3){
  blackmat=matrix(c(rep(c(paste0("gene",1:n2),paste0("mir",1:n3),"phe"),n1),
                    rep("phe",n2+n3),
                    rep(paste0("snp",1:n1),each=n2+n3+1),
                    c(paste0("gene",1:n2),paste0("mir",1:n3))),nrow=(n2+n3+1)*n1+n2+n3)
  temp=expand.grid(paste0("snp",1:n1),paste0("snp",1:n1))
  colnames(blackmat)=c("Var1","Var2")
  blackmat=rbind(blackmat,temp)
  return(blackmat)
}



###### Function to run exhaustive BNA with 4 variables #####

myhillclimb.exhaust=function(candidrow,candid,counts,genos,genes,pheno,all.models,plot=T){
  print(candid[candidrow,])
  smalldata=data.frame(snp=as.numeric(genos[candid[candidrow,2],]),
                       gene=genes[candid[candidrow,3],],
                       mir=counts[candid[candidrow,1],],pheno=pheno)
  
  scorevec=numeric(length(all.models))
  for (i in 1:length(all.models)){
    scorevec[i]=score(model2network(all.models[i]),smalldata)
  }
  bn.hc=model2network(all.models[which.max(scorevec)])
  
  interesting=c(2,5,7,11,12,17,19,21,25,28,31,39,
                41,42,46,47,52,55,61,63,64,70,72,
                76,78,79,81,83,84,90,93,94)#32
  interesting.g=c(2,5,17,19,21,25,31,39,41,42,46,52,55,61,63,64,76,78,79,81,83,90,93,94)#24
  bicgap1=sort(scorevec,decreasing = T)[1]-sort(scorevec,decreasing = T)[2]
  bicgap2=max(scorevec[interesting])-max(scorevec[-interesting])
  bicgap3=max(scorevec[interesting.g])-max(scorevec[setdiff(interesting,interesting.g)])
  
  if (plot==T){
    graphviz.plot(bn.hc)
  }
  
  return(list(bn.hc=bn.hc,scorevec=scorevec,bicgap=c(bicgap1,bicgap2,bicgap3)))
}





#### Function to generate QQ plots ####
qqunif.plot<-function(pvalues, 
                      should.thin=T, thin.obs.places=2, thin.exp.places=2, 
                      xlab=expression(paste("Expected (",-log[10], " p-value)")),
                      ylab=expression(paste("Observed (",-log[10], " p-value)")), 
                      title="QQ-plot", yrange=NULL,
                      draw.conf=TRUE, conf.points=1000, conf.col="lightgray", conf.alpha=.05,
                      already.transformed=FALSE, pch=20, aspect="iso", prepanel=prepanel.qqunif,
                      par.settings=list(superpose.symbol=list(pch=pch)), ...) {
  
  
  #error checking
  if (length(pvalues)==0) stop("pvalue vector is empty, can't draw plot")
  if(!(class(pvalues)=="numeric" || 
       (class(pvalues)=="list" && all(sapply(pvalues, class)=="numeric"))))
    stop("pvalue vector is not numeric, can't draw plot")
  if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
  if (already.transformed==FALSE) {
    if (any(unlist(pvalues)==0)) noprob=1
    #stop("pvalue vector contains zeros, can't draw plot")
  } else {
    if (any(unlist(pvalues)<0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
  }
  
  
  grp<-NULL
  n<-1
  exp.x<-c()
  if(is.list(pvalues)) {
    nn<-sapply(pvalues, length)
    rs<-cumsum(nn)
    re<-rs-nn+1
    n<-min(nn)
    if (!is.null(names(pvalues))) {
      grp=factor(rep(names(pvalues), nn), levels=names(pvalues))
      names(pvalues)<-NULL
    } else {
      grp=factor(rep(1:length(pvalues), nn))
    }
    pvo<-pvalues
    pvalues<-numeric(sum(nn))
    exp.x<-numeric(sum(nn))
    for(i in 1:length(pvo)) {
      if (!already.transformed) {
        pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
        exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method="first")-.5)/nn[i])
      } else {
        pvalues[rs[i]:re[i]] <- pvo[[i]]
        exp.x[rs[i]:re[i]] <- -log10((nn[i]+1-rank(pvo[[i]], ties.method="first")-.5)/(nn[i]+1))
      }
    }
  } else {
    n <- length(pvalues)+1
    if (!already.transformed) {
      exp.x <- -log10((rank(pvalues, ties.method="first")-.5)/n)
      pvalues <- -log10(pvalues)
    } else {
      exp.x <- -log10((n-rank(pvalues, ties.method="first")-.5)/n)
    }
  }
  
  
  #this is a helper function to draw the confidence interval
  panel.qqconf<-function(n, conf.points=1000, conf.col="gray", conf.alpha=.05, ...) {
    require(grid)
    conf.points = min(conf.points, n-1);
    mpts<-matrix(nrow=conf.points*2, ncol=2)
    for(i in seq(from=1, to=conf.points)) {
      mpts[i,1]<- -log10((i-.5)/n)
      mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
      mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
      mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
    }
    grid.polygon(x=mpts[,1],y=mpts[,2], gp=gpar(fill=conf.col, lty=0), default.units="native")
  }
  
  #reduce number of points to plot
  if (should.thin==T) {
    if (!is.null(grp)) {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places),
                                grp=grp))
      grp = thin$grp
    } else {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places)))
    }
    pvalues <- thin$pvalues
    exp.x <- thin$exp.x
  }
  gc()
  
  prepanel.qqunif= function(x,y,...) {
    A = list()
    #A$xlim = range(x, y)*1.02
    if (is.null(yrange)==T) {yrange=range(x, y)*1.02}
    A$xlim=yrange
    A$xlim[1]=0
    A$ylim = A$xlim
    return(A)
  }
  
  #draw the plot
  noninf=!is.infinite(pvalues) & !is.infinite(pvalues)
  pvalues=pvalues[noninf]
  exp.x=exp.x[noninf]
  xyplot(pvalues~exp.x, groups=grp, xlab=xlab, ylab=ylab, main=title, aspect=aspect,
         prepanel=prepanel, 
         scales=list(axs="i"), pch=pch,
         panel = function(x, y, ...) {
           if (draw.conf) {
             panel.qqconf(n, conf.points=conf.points, 
                          conf.col=conf.col, conf.alpha=conf.alpha)
           };
           panel.xyplot(x,y, ...);
           panel.abline(0,1);
         }, par.settings=par.settings, ...
  )
}





























