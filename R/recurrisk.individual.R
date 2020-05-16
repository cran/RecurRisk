recurrisk.individual <-
function(data,stratum,covar,timevar,eventvar,stagevar,stage.dist.value,link="Log-logistic",adj.r=1){
  data<-data.frame(data)
  RR <- adj.r
  if(link=="Weibull"){
    distribution<-"weibull"
  }
  if(link=="Log-logistic"){
    distribution<-"llogis"
  }
  if(is.numeric(RR)== F){
    warning.str<-"Warning: Adjustment Factor r should be numeric."
    print(warning.str)    
  }
  stage.dist.name <- stagevar
  nstratum<-length(stratum)
  ncovar<-length(covar)
  
  if(nstratum>0 & ncovar>0){
    stratum.nm.string<-paste("!is.na(data$",stratum,")", collapse=" & ",sep="")
    covar.nm.string<-paste("!is.na(data$",covar,")", collapse=" & ",sep="")  
    allvar.nm.string<-paste(stratum.nm.string," & ",covar.nm.string,sep="")
  }
  
  if(ncovar==0 & nstratum>0){
    stratum.nm.string<-paste("!is.na(data$",stratum,")", collapse=" & ",sep="")
    allvar.nm.string<-stratum.nm.string
  }
  if(nstratum==0 & ncovar>0){
    covar.nm.string<-paste("!is.na(data$",covar,")", collapse=" & ",sep="")  
    allvar.nm.string<-covar.nm.string
  }
  if(nstratum==0 & ncovar==0){
    stratum<-"nostratum"
    data[,"nostratum"]<-0
  }
  
  if(nstratum>0 | ncovar>0){
    data.nm <- eval(parse(text=paste("data[",allvar.nm.string,",]")))
    data<-data.nm
  }
  
  
  int.max<-max(data[,timevar],na.rm=T) 
  
  allvar<-c(covar,stratum)
  nallvar<-length(allvar)
  allvar.group.nvalue <- rep(NA,nallvar)
  allvar.group <- vector("list",nallvar)
  
  for (iallvar in 1:nallvar){
    if(nstratum==0 & ncovar==0){
      allvar.group[[iallvar]] <- 0
    }
    if(nstratum>0 | ncovar>0){
      iallvar.group<-unique(data[,allvar[iallvar]])
      iallvar.group<-sort(iallvar.group, na.last = NA)
      iallvar.group<-iallvar.group[which(!is.na(iallvar.group))]
      allvar.group[[iallvar]] <- iallvar.group
    }
    allvar.group.nvalue[iallvar] <- length(allvar.group[[iallvar]])
  }
  intervals<-list(1:int.max)
  out.combo.list <- append(intervals,allvar.group) 
  out.combo <- expand.grid(out.combo.list)
  colnames(out.combo)<-c("fup",allvar)
  
  if(ncovar>0){
    covar.char<-paste(covar,".char",sep="")
    out.combo[,covar.char]<-NA
    for(icovar in 1:length(covar)){
      out.combo[,covar.char[icovar]]<-paste(covar[icovar],out.combo[,covar[icovar]],sep="")
    }
  }
  if(stage.dist.name %in% allvar){
    out.combo<-out.combo[-which(out.combo[,stage.dist.name]==stage.dist.value),]
  }
  addcols<-c("cure","lambda","k","r","theta","theta.se","obs_surv","obs_dist_surv","surv_curemodel",
             "surv_notcured","median_surv_notcured","s1_analytical","se_CI_analytical",
             "d_c_int","d_u_int","d_s_int","d_theta","d_cure_int","d_lambda_int","d_k_int")
  out.combo[,addcols]<-NA
  x.combo<-function(x){
    paste(x,collapse="_",sep="")
  }
  
  if(nstratum>0){
    stratum.group<-allvar.group[which(allvar %in% stratum)] 
    stratum.group.nvalue<-allvar.group.nvalue[which(allvar %in% stratum)] 
    stratum.all <- expand.grid(stratum.group)
    colnames(stratum.all)<-stratum
    if(!stagevar %in% stratum){
      stratum.nodist<-stratum.all
    }
    if(stagevar %in% stratum){
      stratum.nodist<-stratum.all[which(stratum.all[,stage.dist.name]!=stage.dist.value),] 
    }
    if(nstratum==1 & (stagevar %in% stratum)){
      dist.pos<-which(stratum.all[,stage.dist.name]==stage.dist.value)
      stratum.nodist<-expand.grid(stratum.group[[1]][c(-dist.pos)])
      colnames(stratum.nodist)<-stratum
    }
    if(nstratum==1 & (!stagevar %in% stratum)){
      stratum.nodist<-expand.grid(stratum.group[[1]])
      colnames(stratum.nodist)<-stratum
    }
    out.combo[,"stratum.combo"]<-apply(data.frame(out.combo[,stratum]), 1, x.combo)
  }
  if(ncovar>0){
    covar.group<-allvar.group[which(allvar %in% covar)] 
    covar.group.nvalue<-allvar.group.nvalue[which(allvar %in% covar)] 
    covar.all <- expand.grid(covar.group)
    colnames(covar.all)<-covar
    
    c.cols<-paste("c_",covar,sep="")
    u.cols<-paste("u_",covar,sep="")
    d.c.cols<-gsub("c_","d_c_",c.cols)
    d.u.cols<-gsub("u_","d_u_",u.cols)
    covar.all[,c("cure","lambda",c.cols,u.cols)]<-NA
    out.combo[,c(d.c.cols,d.u.cols)]<-NA
    
    all.c.cols<-vector("list", ncovar)
    all.u.cols<-vector("list", ncovar)
    cov.c.cols<-vector("list", ncovar)
    cov.u.cols<-vector("list", ncovar)
    
    data[,covar.char]<-NA
    for(icovar in 1:ncovar){
      data[,covar.char[icovar]]<-paste(covar[icovar],data[,covar[icovar]],sep="")
      all.c.cols[[icovar]]<-paste(c.cols[icovar],"_",covar.group[[icovar]],sep="")
      all.u.cols[[icovar]]<-paste(u.cols[icovar],"_",covar.group[[icovar]],sep="")
      cov.c.cols[[icovar]]<-paste(c.cols[icovar],"_",covar.group[[icovar]][-1],sep="")
      cov.u.cols[[icovar]]<-paste(u.cols[icovar],"_",covar.group[[icovar]][-1],sep="")
    }
    out.combo[,"covar.combo"]<-apply(data.frame(out.combo[,covar]), 1, x.combo)
    if(stage.dist.name %in% covar){
      covar.all<-covar.all[-which(covar.all[,stage.dist.name]==stage.dist.value),]
    }
  }
  out.combo$r<-RR
  str.km<-1
  str<-1
  if(ncovar>0){
    str.km<-paste(covar,collapse="+",sep="")
    str.scale<-paste("scale(",covar.char,")",collapse=" + ", sep="")
    str.cure<-paste(covar.char,collapse="+",sep="")
    str<-paste(str.cure,"+",str.scale,sep="")
  }
  str.km<-paste("Surv(",timevar,",", eventvar,")~",str.km,sep="")
  form.km<-as.formula(str.km)
  
  str<-paste("Surv(",timevar,",", eventvar,")~",str,sep="")
  form<-as.formula(str)
  
  time.count0<-proc.time()
  if(nstratum==0){
    nstratumgroup<-1
  }
  if(nstratum>0){
    nstratumgroup<-dim(stratum.nodist)[1]
  }
  
  for(istratumgroup in 1:nstratumgroup){
    if(nstratum>0){
      stratum.value<-stratum.nodist[istratumgroup,]
      stratum.combo<-x.combo(stratum.value)
      stratum.dist.value<-stage.dist.value
      if(nstratum>1 & (stagevar %in% stratum)){
        stratum.dist.value<-stratum.value
        stratum.dist.value[which(stratum==stage.dist.name)]<-stage.dist.value
      }
      condition.string<-paste("data$",stratum,"==",stratum.value, collapse=" & ",sep="")
      data.km <- eval(parse(text=paste("data[",condition.string,",]")))
      dist.string<-paste("data$",stratum,"==",stratum.dist.value, collapse=" & ",sep="")   
      distdata.km <- eval(parse(text=paste("data[",dist.string,",]")))
      if(!stagevar %in% stratum){
        dist.string<-paste(condition.string, " & data$", stage.dist.name, "==",stage.dist.value,sep="")   
        distdata.km <- eval(parse(text=paste("data[",dist.string,",]")))
      }
    }
    
    if(nstratum==0){
      data.km<-data
      #data.km<-data[which(data[,stage.dist.name]!=stage.dist.value),]
      distdata.km<-data[which(data[,stage.dist.name]==stage.dist.value),]
    }
    
    km.model<-survfit(form.km, data=data.km,type="kaplan-meier")
    fit.km<-summary(km.model)
    nobs<-km.model$strata
    
    km.model.dist<-survfit(form.km, data=distdata.km,type="kaplan-meier")
    fit.km.dist<-summary(km.model.dist)
    nobs.dist<-km.model.dist$strata
    
    cure.model<-flexsurvcure(form, data=data.km,dist=distribution, mixture=T)
    res<-cure.model$res
    cov<-cure.model$cov
    if(ncovar==0){
      row.names(res)<-c("cure_int","lambda_int","k_int")
    }
    if(ncovar>0){
      row.names(res)<-c("cure_int","lambda_int","k_int",unlist(cov.c.cols),unlist(cov.u.cols))
      scale.pos<-which(substr(row.names(res),1,2)=="u_")
      cure.pos<-c(4:(scale.pos[1]-1))
      c.covar.vec<-vector("list", ncovar)
      u.covar.vec<-vector("list", ncovar)
      for(icovar in 1:ncovar){
        covar.i<-covar[icovar]
        if(icovar==1){
          covar.i.curepos<-cure.pos[1:(covar.group.nvalue[1]-1)]
          covar.i.scalepos<-scale.pos[1:(covar.group.nvalue[1]-1)]
        }
        if(icovar>1){
          covar.i.curepos<-cure.pos[(sum(covar.group.nvalue[1:(icovar-1)]-1)+1):sum(covar.group.nvalue[1:(icovar)]-1)]
          covar.i.scalepos<-scale.pos[(sum(covar.group.nvalue[1:(icovar-1)]-1)+1):sum(covar.group.nvalue[1:(icovar)]-1)]
        }
        c.covar.vec[[icovar]]<-c(0,res[covar.i.curepos,1])
        u.covar.vec[[icovar]]<-c(0,res[covar.i.scalepos,1])
      }
    }
    colnames(cov)<-row.names(res)
    rownames(cov)<-colnames(cov)
    
    cure.int<-res[1,1]
    k.int<-res[2,1]
    lambda.int<-res[3,1]
    
    c.int<-log(cure.int/(1-cure.int))
    l.int<-log(lambda.int)
    
    if(nstratum==0){
      rows.istratum<-1:dim(out.combo)[1]
    }
    if(nstratum>0){
      rows.istratum<-which(out.combo$stratum.combo==stratum.combo)
    }
    
    out.combo[rows.istratum,"k"]<-k.int
    s.int<--log(k.int)
    out.combo[rows.istratum,"c_int"]<-c.int
    out.combo[rows.istratum,"u_int"]<-l.int
    out.combo[rows.istratum,"s_int"]<-s.int
    
    theta0 <- 0.1
    
    if(ncovar==0){
      ncovargroup<-1
      out.combo[rows.istratum,"cure"]<-cure.int
      out.combo[rows.istratum,"lambda"]<-lambda.int
    }  
    if(ncovar>0){
      ncovargroup<-dim(covar.all)[1]
    }
    
    for(icovargroup in 1:ncovargroup){
      
      d.vec.keep<-c("d_cure_int","d_lambda_int","d_k_int","d_theta")  
      cols.keep<-c("cure_int","lambda_int","k_int")
      ### vector dataframe value change
      if(ncovar>0){
        covar.value<-covar.all[icovargroup,covar]
        covar.char.value<-paste(covar,covar.value,sep="")
        covar.string<-paste(covar,"=",covar.value, collapse=", ",sep="")
        covar.char.string<-paste(covar.char,"=",covar.char.value, collapse=",",sep="")
        d.vec.keep<-c("d_cure_int","d_lambda_int","d_k_int",d.c.cols,d.u.cols,"d_theta")  
        cols.keep.c<-vector("list", ncovar)
        cols.keep.u<-vector("list", ncovar)
        
        for(icovar in 1:ncovar){
          covar.all[icovargroup,c.cols[icovar]]<-c.covar.vec[[icovar]][which(covar.group[[icovar]]==covar.all[icovargroup,covar[[icovar]]])]
          covar.all[icovargroup,u.cols[icovar]]<-u.covar.vec[[icovar]][which(covar.group[[icovar]]==covar.all[icovargroup,covar[[icovar]]])]
          cols.keep.c[[icovar]]<-all.c.cols[[icovar]][which(covar.group[[icovar]]==covar.all[icovargroup,covar[[icovar]]])]
          cols.keep.u[[icovar]]<-all.u.cols[[icovar]][which(covar.group[[icovar]]==covar.all[icovargroup,covar[[icovar]]])]
          if(covar.all[icovargroup,covar[[icovar]]]==covar.group[[icovar]][1]){
            cols.keep.c[[icovar]]<-"cure_int"
            cols.keep.u[[icovar]]<-"lambda_int"
            d.vec.keep[which(d.vec.keep==d.c.cols[icovar])]<-NA
            d.vec.keep[which(d.vec.keep==d.u.cols[icovar])]<-NA
          }
        }
        d.vec.keep<-d.vec.keep[which(!is.na(d.vec.keep))]
        cols.keep<-c("cure_int","lambda_int","k_int",unlist(cols.keep.c),unlist(cols.keep.u))
        covar.all[icovargroup,"cure"]<-1/(1+exp(-(c.int+sum(covar.all[icovargroup,c.cols]))))
        covar.all[icovargroup,"lambda"]<-exp(l.int+sum(covar.all[icovargroup,u.cols]))
      }
      
      cols.keep.pos<-which(colnames(cov) %in% cols.keep)
      cov.icovargroup<-cov[cols.keep.pos,cols.keep.pos]
      
      if(ncovar==0){
        surv<-fit.km$surv
        surv.se<-fit.km$std.err
        distsurv<-fit.km.dist$surv
        distsurv.se<-fit.km.dist$std.err
        if(length(fit.km$surv)-length(fit.km.dist$surv)>0){
          surv<-c(surv,rep(NA,int.max-length(surv)))
          surv.se<-c(surv.se,rep(NA,int.max-length(surv.se)))
          distsurv<-c(distsurv,rep(distsurv[length(distsurv)],length(fit.km$surv)-length(distsurv)),rep(NA,int.max-length(fit.km$surv)))
          distsurv.se<-c(distsurv.se,rep(distsurv.se[length(distsurv.se)],length(fit.km$std.err)-length(distsurv.se)),rep(NA,int.max-length(fit.km$std.err)))
        }
        if(length(fit.km$surv)-length(fit.km.dist$surv)<=0){
          distsurv<-c(distsurv,rep(NA,int.max-length(distsurv)))
          distsurv.se<-c(distsurv.se,rep(NA,int.max-length(distsurv.se)))
          surv<-c(surv,rep(surv[length(surv)],length(fit.km.dist$surv)-length(surv)),rep(NA,int.max-length(fit.km.dist$surv)))
          surv.se<-c(surv.se,rep(surv.se[length(surv.se)],length(fit.km.dist$std.err)-length(surv.se)),rep(NA,int.max-length(fit.km.dist$std.err)))
          
        }
        surv.cure<-summary(cure.model, t=seq(from=1,to=int.max,by=1), type="survival")[[1]][,"est"]
      }
      
      if(ncovar>0){
        surv<-fit.km$surv[which(fit.km$strata==covar.string)]
        surv.se<-fit.km$std.err[which(fit.km$strata==covar.string)]
        surv<-c(surv,rep(surv[length(surv)],int.max-length(surv)))
        surv.se<-c(surv.se,rep(surv.se[length(surv.se)],int.max-length(surv.se)))
        strata.pos<-which(names(km.model$strata)==covar.string)
        if(nobs[strata.pos]<int.max){
          surv[c((nobs[strata.pos]+1):length(surv))]<-NA
          surv.se[c((nobs[strata.pos]+1):length(surv.se))]<-NA
        }
        
        if((stage.dist.name %in% covar) & ncovar==1){
          distsurv<-fit.km.dist$surv
          distsurv.se<-fit.km.dist$std.err
          distsurv<-c(distsurv,rep(distsurv[length(distsurv)],nobs[strata.pos]-length(distsurv)),rep(NA,int.max-nobs[strata.pos]))
          distsurv.se<-c(distsurv.se,rep(distsurv.se[length(distsurv.se)],nobs[strata.pos]-length(distsurv.se)),rep(NA,int.max-nobs[strata.pos]))
        }
        if((stage.dist.name %in% covar) & ncovar>1){
          dist.covar.value<-covar.value
          dist.covar.value[which(covar==stage.dist.name)]<-stage.dist.value
          dist.covar.string<-paste(covar,"=",dist.covar.value, collapse=", ",sep="")
          distsurv<-fit.km.dist$surv[which(fit.km.dist$strata==dist.covar.string)]
          distsurv.se<-fit.km.dist$std.err[which(fit.km.dist$strata==dist.covar.string)]
          distsurv<-c(distsurv,rep(distsurv[length(distsurv)],int.max-length(distsurv)))
          distsurv.se<-c(distsurv.se,rep(distsurv.se[length(distsurv.se)],int.max-length(distsurv.se)))
          strata.dist.pos<-which(names(km.model.dist$strata)==dist.covar.string)
          if(nobs.dist[strata.dist.pos]<int.max){
            distsurv[c((nobs.dist[strata.dist.pos]+1):length(distsurv))]<-NA
            distsurv.se[c((nobs.dist[strata.dist.pos]+1):length(distsurv.se))]<-NA
          }
        }
        if(!stage.dist.name %in% covar){
          distsurv<-fit.km.dist$surv[which(fit.km.dist$strata==covar.string)]
          distsurv.se<-fit.km.dist$std.err[which(fit.km.dist$strata==covar.string)]
          distsurv<-c(distsurv,rep(distsurv[length(distsurv)],int.max-length(distsurv)))
          distsurv.se<-c(distsurv.se,rep(distsurv.se[length(distsurv.se)],int.max-length(distsurv.se)))
          strata.dist.pos<-which(names(km.model.dist$strata)==covar.string)
          if(nobs.dist[strata.dist.pos]<int.max){
            distsurv[c((nobs.dist[strata.dist.pos]+1):length(distsurv))]<-NA
            distsurv.se[c((nobs.dist[strata.dist.pos]+1):length(distsurv.se))]<-NA
          }
        }
        
        surv.cure<-summary(cure.model, t=seq(from=1,to=int.max,by=1), type="survival")[[strata.pos]][,"est"]
      }
      
      if(istratumgroup==1){
        rows.add<-(int.max*(icovargroup-1)+1):(icovargroup*int.max)
      }
      if(istratumgroup>1){
        rows.add<-(int.max*(icovargroup-1)+1):(icovargroup*int.max)+sum(table(out.combo$stratum.combo)[1:(istratumgroup-1)])
      }
      out.combo[rows.add,"obs_surv"]<-surv
      out.combo[rows.add,"obs_dist_surv"]<-distsurv
      out.combo[rows.add,"surv_curemodel"]<-surv.cure
      
      if(ncovar>0){
        out.combo[rows.add,"cure"]<-covar.all[icovargroup,"cure"]
        out.combo[rows.add,"lambda"]<-covar.all[icovargroup,"lambda"]
      }
      
      fup <- 1:int.max
      fit<-NULL
      y<-distsurv[which(distsurv!=0 & distsurv.se!=0)]
      se<-distsurv.se[which(distsurv!=0 & distsurv.se!=0)]
      x<-fup[which(distsurv!=0 & distsurv.se!=0)]
      
      if(dim(distdata.km)[1]>2 & length(which(!is.na(y)))>2 & length(table(y))>2){
        try(fit <- nls(y~exp(-theta*x),start=list(theta=theta0),control=nls.control(maxiter=2000),weights=1/se^2))
        if(is.null(fit)){
          error.str<-"Theta cannot be estimated due to the nls function. Please check the data for the subgroup."
          print(error.str)
          out.combo[rows.add,"theta"] <- NA
          out.combo[rows.add,"theta.se"] <- NA
        }
        if(!is.null(fit)){
          theta <- summary(fit)$coefficients[1,1]
          theta.se <- summary(fit)$coefficients[1,2]
          out.combo[rows.add,"theta"] <- RR*theta
          out.combo[rows.add,"theta.se"] <- theta.se
        }
      }
      
      if (distribution=="weibull"){
        out.combo[rows.add,"surv_notcured"] <- with(out.combo[rows.add,],exp(-(fup/lambda)**k))
        out.combo[rows.add,"s1_analytical"] <- with(out.combo[rows.add,],surv_notcured*(1-k/theta/lambda*(fup/lambda)**(k-1)))
        out.combo[rows.add,"d_u_int"] <- with(out.combo[rows.add,],(1-cure)*(k*s1_analytical*(fup/lambda)**k+exp(-(fup/lambda)**k)*k*k/theta/lambda*(fup/lambda)**(k-1)))
        out.combo[rows.add,"d_s_int"]  <- with(out.combo[rows.add,],(1-cure)*(s1_analytical*k*(fup/lambda)**k*log(fup/lambda)+exp(-(fup/lambda)**k)*k/theta*lambda**(-k)*fup**(k-1)*(1+k*log(fup/lambda))))
        out.combo[rows.add,"d_theta"]  <- with(out.combo[rows.add,],(1-cure)*exp(-(fup/lambda)**k)* k/theta**2/lambda*(fup/lambda)**(k-1))
        out.combo[rows.add,"median_surv_notcured"]<-with(out.combo[rows.add,],lambda * log(2)**(1/k))
        
      }
      if (distribution=="llogis"){
        out.combo[rows.add,"surv_notcured"] <- with(out.combo[rows.add,],1/(1+(fup/lambda)**k))
        out.combo[rows.add,"s1_analytical"] <- with(out.combo[rows.add,],(1+(fup/lambda)**k-k/theta/lambda*(fup/lambda)**(k-1))/(1+(fup/lambda)**k)**2)
        out.combo[rows.add,"d_u_int"] <- with(out.combo[rows.add,],(1-cure)*k*(s1_analytical*2*(fup/lambda)**k/(1+(fup/lambda)**k)-s1_analytical+1/(1+(fup/lambda)**k)**2))
        out.combo[rows.add,"d_s_int"]  <- with(out.combo[rows.add,],(1-cure)*k*(fup/lambda)**k/(1+(fup/lambda)**k)*((1+k*log(fup/lambda)-theta*fup*log(fup/lambda))/theta/fup/(1+(fup/lambda)**k)+2*s1_analytical*log(fup/lambda)))
        out.combo[rows.add,"d_theta"]  <- with(out.combo[rows.add,],(1-cure)*k/lambda/theta**2*(fup/lambda)**(k-1)/(1+(fup/lambda)**k)**2)
        out.combo[rows.add,"median_surv_notcured"]<-with(out.combo[rows.add,],lambda)
      }
      
      out.combo[rows.add,"d_c_int"]  <- with(out.combo[rows.add,],cure*(1-cure)*(1-s1_analytical))  ## different format from sas formula
      out.combo[rows.add,"d_cure_int"]<- out.combo[rows.add,"d_c_int"]*(1/cure.int-1/(1-cure.int))
      out.combo[rows.add,"d_lambda_int"]<- out.combo[rows.add,"d_u_int"]*1/lambda.int
      out.combo[rows.add,"d_k_int"]<- out.combo[rows.add,"d_s_int"]*(-1/k.int)
      
      if(ncovar>0){
        out.combo[rows.add,d.c.cols] <- out.combo[rows.add,"d_c_int"]
        out.combo[rows.add,d.u.cols] <- out.combo[rows.add,"d_u_int"]
      }
      
      cov.sub<-cov.icovargroup
      row.0 <- rep(0,dim(cov.sub)[1])
      cov.sub <- rbind(cov.sub,row.0)
      col.0 <- rep(0,dim(cov.sub)[2]+1)
      cov.sub <- cbind(cov.sub,col.0)
      colnames(cov.sub)[dim(cov.icovargroup)[1]+1]<-"theta"
      rownames(cov.sub)[dim(cov.icovargroup)[1]+1]<-"theta"
      cov.sub[dim(cov.sub)[1],dim(cov.sub)[1]]<-theta.se^2
      
      ## currently work for model with 1 covariate, need to test for multiple covars, could add columns 0's for group 0 
      for(k in 1:length(rows.add)){
        d.vec.k<-out.combo[rows.add[k],d.vec.keep]
        d.vec.k<-as.numeric(d.vec.k)
        varcov.k<-t(d.vec.k)%*%as.matrix(cov.sub)%*%d.vec.k
        out.combo[rows.add[k],"se_CI_analytical"]<-sqrt(varcov.k[1,1])
      }
    }
  }
  del.cols<-c("c_int","u_int","s_int","d_c_int","d_u_int","d_s_int","d_theta",
              "d_cure_int","d_lambda_int","d_k_int","theta.se")
  
  if(nstratum>0 & ncovar>0){
    del.cols<-c(del.cols,d.c.cols,d.u.cols, covar.char,"covar.combo","stratum.combo")
  }
  if(nstratum==0 & ncovar>0){
    del.cols<-c(del.cols,d.c.cols,d.u.cols, covar.char,"covar.combo")
  }
  if(nstratum>0 & ncovar==0){
    del.cols<-c(del.cols,"stratum.combo")
  }
  
  out.combo[,del.cols]<-NULL
  out<-out.combo
  
  t_surv  <- out$surv_notcured
  t2_surv <- out$obs_dist_surv
  
  rows.fup1  <- which(out$fup==1)
  vec.fup1   <- c(rows.fup1,dim(out)[1]+1)
  diff.fup1  <- diff(vec.fup1)
  old_t_surv <- c(1,t_surv[-length(t_surv)])
  old_t_surv[rows.fup1] <- 1
  old_t2_surv <- c(1,t2_surv[-length(t2_surv)])
  old_t2_surv[rows.fup1] <- 1
  
  
  t_pdf  <- old_t_surv-t_surv
  t2_pdf <- old_t2_surv^RR-t2_surv^RR
  t2_pdf[rows.fup1] <- 1-t2_surv[rows.fup1]^RR
  temp_pdf <- cbind(t_pdf,t2_pdf)
  
  
  volterra_continuous <- function(fstar,f2,delta=1){
    f1      <- rep(NA,length(fstar))
    f1_pdf  <- rep(NA,length(fstar))
    f1_surv <- rep(NA,length(fstar))
    f1[1]   <- fstar[1]/(delta*f2[1])
    f1_pdf[1]  <- f1[1]
    f1_surv[1] <- 1-f1[1]*delta
    for (i in 2:length(fstar)){
      f1[i] <- fstar[i]
      for (j in 1:(i-1)){
        f1[i] <- f1[i]-f1[j]*f2[i-j+1]*delta
      }
      f1[i] <- f1[i]/(delta*f2[1])
      f1_pdf[i] <- f1[i]
      f1_surv[i] <- f1_surv[i-1]-f1[i]*delta
    }
    return(f1_surv)
  }
  
  out[,"s1_numerical"]<-NA
  
  for(ifup1 in 1:length(rows.fup1)){
    
    istart <- rows.fup1[ifup1]
    iend   <- (rows.fup1+diff.fup1-1)[ifup1]
    fstar  <- t_pdf[istart:iend]
    f2     <- t2_pdf[istart:iend]
    delta  <- 1
    f1_surv_ifup1 <- volterra_continuous(fstar,f2)
    out$s1_numerical[istart:iend] <- f1_surv_ifup1
  }
  
  out$s1_numerical[which(out$s1_numerical<0)] <- 0
  out$s1_analytical[which(out$s1_analytical<0)] <- 0
  out[,"G_numerical"]   <- with(out,cure+(1-cure)*s1_numerical)
  out[,"CI_numerical"]  <- 1-out$G_numerical
  out[,"G_analytical"]  <- with(out,cure+(1-cure)*s1_analytical)
  out[,"CI_analytical"] <- 1-out$G_analytical
  out[,"link"]<-link
  
  #out <- out[which(out$fup<=int.max.out),]
  if(nstratum==0 & ncovar==0){
    stratum<-NULL
  }
  out <- out[,c(stratum,covar,"fup","link","r","cure","lambda","k","theta",
                "surv_curemodel","surv_notcured","median_surv_notcured",
                "s1_numerical","G_numerical","CI_numerical",
                "s1_analytical","G_analytical","CI_analytical","se_CI_analytical",
                "obs_surv","obs_dist_surv")]
  colnames(out)[which(colnames(out)=="fup")]<-"followup"
  return(out)
}
