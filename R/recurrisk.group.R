recurrisk.group <-
function(data,data.cansurv,stagevar,stage.dist.value,adj.r=1){
    seerdata<-data
    ## clean headers
    header<-colnames(seerdata)
    header.clean = gsub("[,()<>={}!@#$%^&*+-]", "", header);
    header.seer = gsub(" ", "_", header.clean);
    colnames(seerdata)<-header.seer
    csdata <- data.cansurv
    cnames.seer <- colnames(seerdata)
    cnames.cs <- colnames(csdata)
    RR <- adj.r
    if(max(nchar(cnames.seer))>32){
      warning.str<-"Warning: SEER*Stat data contains a variable with name exceeding 32 characters which may cause problems in CanSurv output."
      print(warning.str)    
    }
    if(is.numeric(RR)== F){
      warning.str<-"Warning: Adjustment Factor r should be numeric."
      stop(warning.str)    
    }
    
    stage.dist.name <- stagevar
    
    if(length(which(grepl("Cause", cnames.seer)==T))>1){
      surv.name <- cnames.seer[which(grepl("Survival", cnames.seer)==T & grepl("Cum", cnames.seer)==T & grepl("Cause", cnames.seer))]
      survse.name <- cnames.seer[which(grepl("SE", cnames.seer)==T & grepl("Cum", cnames.seer)==T & grepl("Cause", cnames.seer))]
    }
    if(length(which(grepl("Relative", cnames.seer)==T))>1){
      surv.name <- cnames.seer[which(grepl("Survival", cnames.seer)==T & grepl("Cum", cnames.seer)==T & grepl("Relative", cnames.seer))]
      survse.name <- cnames.seer[which(grepl("SE", cnames.seer)==T & grepl("Cum", cnames.seer)==T & grepl("Relative", cnames.seer))]
    }
    ## convert to proportions if the survival rates are percentages
    if(length(which(seerdata[,surv.name]>1))>0){
      seerdata[,surv.name]<-seerdata[,surv.name]*0.01
      seerdata[,survse.name]<-seerdata[,survse.name]*0.01
    }
    
    int.max.seer <- max(seerdata$Interval,na.rm=T)
    allvar.seer<-cnames.seer[(which(cnames.seer!="Page_type")[1]):(which(cnames.seer=="Interval")-1)]
    cols.keep.seer<-c(allvar.seer,"Interval",surv.name)
    seerdata.keep <- seerdata[which(seerdata[,stage.dist.name]!=stage.dist.value),]   
    seerdata.out <- seerdata.keep[,cols.keep.seer]
    seerdata.dist <- seerdata[which(seerdata[,stage.dist.name]==stage.dist.value),]   
    seerdata.dist.out <- seerdata.dist[,cols.keep.seer]
    seerdata.dist.out[,stage.dist.name]<-NULL
    colnames(seerdata.out)[which(colnames(seerdata.out)==surv.name)]<-"obs_surv"
    colnames(seerdata.dist.out)[which(colnames(seerdata.dist.out)==surv.name)]<-"obs_dist_surv"
    
    no.stratum <- 0
    no.covar   <- 0
    no.covar.c <- 0
    no.covar.u <- 0
    no.cure    <- 0
    
    if(!"c_int" %in% cnames.cs){
      no.cure<-1 
    }  
    
    if(no.cure==0){
      cint.pos <- which(substr(cnames.cs,1,5)=="c_int")
    }
    uint.pos <- which(substr(cnames.cs,1,5)=="u_int")
    sint.pos <- which(substr(cnames.cs,1,5)=="s_int")
    
    if(no.cure==1){
      csdata[,"c_int"]<-NA
      csdata<-csdata[,c(1:(uint.pos-1),(sint.pos+1),uint.pos:sint.pos)]
      cnames.cs<-colnames(csdata)
      cint.pos <- which(substr(cnames.cs,1,5)=="c_int")
      uint.pos <- which(substr(cnames.cs,1,5)=="u_int")
      sint.pos <- which(substr(cnames.cs,1,5)=="s_int")
    }
    
    
    if(uint.pos-cint.pos==1){
      no.covar.c<-1
      covar.c<-NA 
    }  
    if(sint.pos-uint.pos==1){
      no.covar.u<-1
      covar.u<-NA 
    }  
    
    if(no.covar.c==1 & no.covar.u==1){
      no.covar<-1
    }
    
    if (cint.pos==3){
      no.stratum<-1
      csdata[,"Stratum_1__NoStratum"]<-0
      csdata<-csdata[,c(1,(sint.pos+1),2:sint.pos)]
      cnames.cs<-colnames(csdata)
      cint.pos <- which(substr(cnames.cs,1,5)=="c_int")
      uint.pos <- which(substr(cnames.cs,1,5)=="u_int")
      sint.pos <- which(substr(cnames.cs,1,5)=="s_int")
    }
    
    
    stratum.pos <- which(substr(cnames.cs,1,7)=="Stratum")
    nstratum <- length(stratum.pos)
    stratum.header <- cnames.cs[stratum.pos]
    stratum.str <- strsplit(stratum.header, "__")
    stratum <- sapply(stratum.str, "[", 2)
    ## stratum.vec  Stratum_1 Stratum_2
    stratum.vec <- sapply(stratum.str, "[", 1)
    
    if(length(which(!stratum %in% cnames.seer))>0){
      stratum<-gsub("[^[:alnum:][:space:]_]", "", stratum)
      stratum<-gsub("__+","_",stratum)
      cnames.cs[stratum.pos]<-paste(stratum.vec,"__",stratum,sep="")
    }
    
    if(no.covar.c==0){
      covar.header.c <- cnames.cs[(cint.pos+1):(uint.pos-1)]
      covar.str.c <- strsplit(covar.header.c, "__")
      covar.c <- unique(sapply(covar.str.c, "[", 2))
      ## covar.vec  covar_1 covar_2
      covar.vec.c <- unique(sapply(covar.str.c, "[", 1))
      covar.vec.c <- gsub("c","c_covar",covar.vec.c)
      ncovar.c <- length(covar.c)
      if(length(which(!covar.c %in% cnames.seer))>0){
        covar.c<-gsub("[^[:alnum:][:space:]_]", "", covar.c) 
        covar.c<-gsub("__+","_", covar.c)
        
        covar.c.prefix <- sapply(covar.str.c, "[", 1)
        covar.c.middle <- sapply(covar.str.c, "[", 2)
        covar.c.suffix <- sapply(covar.str.c, "[", 3)
        covar.c.middle <- gsub("[^[:alnum:][:space:]_]", "", covar.c.middle)
        covar.c.middle <- gsub("__+","_", covar.c.middle)
        cnames.cs[(cint.pos+1):(uint.pos-1)]<-paste(covar.c.prefix,"__",covar.c.middle,"__",covar.c.suffix,sep="")
      }
    }
    
    
    if(no.covar.u==0){
      covar.header.u <- cnames.cs[(uint.pos+1):(sint.pos-1)]
      covar.str.u <- strsplit(covar.header.u, "__")
      covar.u <- unique(sapply(covar.str.u, "[", 2))
      ## covar.vec  covar_1 covar_2
      covar.vec.u <- unique(sapply(covar.str.u, "[", 1))
      covar.vec.u <- gsub("u","u_covar",covar.vec.u)
      ncovar.u <- length(covar.u)
      if(length(which(!covar.u %in% cnames.seer))>0){
        covar.u<-gsub("[^[:alnum:][:space:]_]", "", covar.u) 
        covar.u<-gsub("__+","_", covar.u)
        
        covar.u.prefix <- sapply(covar.str.u, "[", 1)
        covar.u.middle <- sapply(covar.str.u, "[", 2)
        covar.u.suffix <- sapply(covar.str.u, "[", 3)
        covar.u.middle <- gsub("[^[:alnum:][:space:]_]", "", covar.u.middle)
        covar.u.middle <- gsub("__+","_", covar.u.middle)
        cnames.cs[(uint.pos+1):(sint.pos-1)]<-paste(covar.u.prefix,"__",covar.u.middle,"__",covar.u.suffix,sep="")
      }
    }
    
    if(no.covar.c==0 | no.covar.u==0){
      covar  <- union(covar.c,covar.u)
      if(no.covar.c==1){
        covar  <- covar.u
      }
      if(no.covar.u==1){
        covar  <- covar.c
      }
      ncovar <- length(covar)
      covar.vec<-paste("covar_", 1:ncovar, sep="")
    }
    
    cnames.cs.new <- cnames.cs
    for (istratum in 1:nstratum){
      stratum.i <- stratum[istratum]
      stratum.i.str <- paste("__",stratum.i,sep="")
      cnames.cs.new <- gsub(stratum.i.str,"",cnames.cs.new)
      if(no.covar==0){
        for (icovar in 1:ncovar){
          covar.i <- covar[icovar]
          covar.i.str <- paste("__",covar.i,"_",sep="")
          cnames.cs.new <- gsub(covar.i.str,"",cnames.cs.new)
        }
      }
    }
    colnames(csdata) <- cnames.cs.new
    
    stratum.level.nvalue <- rep(NA,nstratum)
    stratum.level <- vector("list", nstratum)
    for (istratum in 1:nstratum){
      stratum.col<-paste("Stratum_",istratum,sep="")
      stratum.level[[istratum]] <- unique(csdata[,stratum.col])
      stratum.level.nvalue[istratum] <- length(stratum.level[[istratum]])
    }
    
    level.combo <- expand.grid(stratum.level)
    nlevel <- dim(level.combo)[1]
    
    ## create stratum.combo
    ## for example, stratum_1.value=1, stratum_2.value=0, then stratum.combo=1_0
    
    x.combo<-function(x){
      paste(x,collapse="_",sep="")
    }
    
    stratum.combo <- rep(NA,dim(csdata)[1])
    stratum.combo <- apply(data.frame(csdata[,stratum.vec]), 1, x.combo)
    
    if(no.cure==0){
      cus.str.percombo<-cnames.cs.new[cint.pos:sint.pos]
    }
    if(no.cure==1){
      cus.str.percombo<-cnames.cs.new[uint.pos:sint.pos]
    }
    cus.str <- rep(cus.str.percombo,nlevel)
    rnames.cs<-paste(stratum.combo,"__",cus.str,sep="")
    
    
    ########################################
    ### if long covar name not in seer var
    if(no.covar==0){
      for(icovar in 1:length(covar)){
        covar.i0 <- covar[icovar]
        covar.i<-covar[icovar]
        if(!covar.i %in% allvar.seer){
          seer.check<-gsub("_","",allvar.seer)
          covar.i<-gsub("_","",covar.i)
          covar.i<-allvar.seer[which(grepl(covar.i,seer.check))]
          covar[icovar]<-covar.i
          if(no.covar.c==0){
            covar.c[which(covar.c==covar.i0)]<-covar.i
          }
          if(no.covar.u==0){
            covar.u[which(covar.u==covar.i0)]<-covar.i
          }
        }
      }
    }
    
    if(no.covar==0 & no.stratum==0){
      covar.group.nvalue <- rep(NA,ncovar)
      covar.group <- vector("list", ncovar)
      for (icovar in 1:ncovar){
        #covar.group[[icovar]] <- unique(seerdata.dist[,covar[icovar]])
        covar.group[[icovar]] <- unique(seerdata[,covar[icovar]])
        covar.group.nvalue[icovar] <- length(covar.group[[icovar]])
      }
      allvar<-c(stratum,covar)
      allvar.vec<-c(stratum.vec,covar.vec)
    }
    
    
    if(no.stratum==1 & no.covar==0){
      allvar<-covar
      allvar.vec<- covar.vec
    }
    if(no.covar==1){
      allvar<-stratum
      allvar.vec<- stratum.vec
    }
    
    nallvar<-length(allvar)
    allvar.group.nvalue <- rep(NA,nallvar)
    allvar.group <- vector("list",nallvar)
    
    for (iallvar in 1:nallvar){
      if(no.stratum==1 & no.covar==1){
        allvar.group[[iallvar]] <- 0
      }
      if(no.stratum==0 | no.covar==0){
        allvar.group[[iallvar]] <- unique(seerdata[,allvar[iallvar]])
      }
      allvar.group.nvalue[iallvar] <- length(allvar.group[[iallvar]])
    }
    
    
    nallvar.seer<-length(allvar.seer)
    allvar.seer.group.nvalue <- rep(NA,nallvar.seer)
    allvar.seer.group <- vector("list",nallvar.seer)
    
    for (iallvar.seer in 1:nallvar.seer){
      allvar.seer.group[[iallvar.seer]] <- unique(seerdata[,allvar.seer[iallvar.seer]])
      allvar.seer.group.nvalue[iallvar.seer] <- length(allvar.seer.group[[iallvar.seer]])
    }
    
    group.combo <- expand.grid(allvar.group)
    ngroup <- dim(group.combo)[1]
    
    seer.group.combo <- expand.grid(allvar.seer.group)
    nseergroup <- dim(seer.group.combo)[1]
    
    ##### estimate exponential distribution parameter ####
    
    theta0 <- 0.1
    theta.sum <- array(NA,dim=c(nseergroup,nallvar.seer+3))
    colnames(theta.sum)<-c(allvar.seer,"r","theta","theta.se")
    
    for (iseergroup in 1:nseergroup){
      seer.group.combo.igroup<-seer.group.combo[iseergroup,]
      seer.group.combo.igroup[which(allvar.seer==stage.dist.name)]<-stage.dist.value
      condition.string<-paste("seerdata.dist$",allvar.seer,"==",seer.group.combo.igroup, collapse=" & ",sep="")    
      data.sub <- eval(parse(text=paste("seerdata.dist[",condition.string,",]")))
      y0 <- data.sub[,surv.name] 
      x0 <- data.sub$Interval 
      se0 <- data.sub[,survse.name]
      y<-y0[which(y0!=0 & se0!=0)]
      se<-se0[which(y0!=0 & se0!=0)]
      x<-x0[which(y0!=0 & se0!=0)]
      fit<-NULL
      if(dim(data.sub)[1]>2 & length(which(!is.na(y)))>2 & length(table(y))>2){
        try(fit <- nls(y~exp(-theta*x),start=list(theta=theta0),control=nls.control(maxiter=2000),weights=1/se^2))
        if(is.null(fit)){
          errorgroup.str<-paste(allvar.seer,"==",seer.group.combo.igroup, collapse=" & ",sep="")    
          error.str<-"Theta cannot be estimated due to the nls function error. Please check the data for the below group:"
          error.print<-paste(error.str,errorgroup.str,sep="")
          print(error.print)
          nls.error.ind[iseergroup]<-1
          theta.sum[iseergroup,] <- c(as.numeric(seer.group.combo[iseergroup,]),RR,NA,NA)
        }
        if(!is.null(fit)){
          theta <- summary(fit)$coefficients[1,1]
          theta.se <- summary(fit)$coefficients[1,2]
          theta.sum[iseergroup,] <- c(as.numeric(seer.group.combo[iseergroup,]),RR,RR*theta,theta.se)
        }
      }
    }
    
    theta.sum<-data.frame(theta.sum)
    theta.sum[,"allvar.seer.combo"] <- NA
    theta.sum$allvar.seer.combo <- apply(data.frame(theta.sum[,allvar.seer]), 1, x.combo)
    
    
    ##### create data frame for output
    intervals<-list(1:int.max.seer)
    out.combo.list <- append(intervals,allvar.group) 
    out.combo <- expand.grid(out.combo.list)
    colnames(out.combo)<-c("fup",allvar.vec)
    
    seer.out.combo.list <- append(intervals,allvar.seer.group) 
    seer.out.combo <- expand.grid(seer.out.combo.list)
    colnames(seer.out.combo)<-c("Interval",allvar.seer)
    seer.out.combo<-seer.out.combo[which(seer.out.combo[,stage.dist.name]!=stage.dist.value),]
    
    if(no.stratum==1){
      out.combo[,stratum.vec] <- 0
    }
    
    if(no.covar==1){
      out.combo[,c("Link","allvar.combo","stratum.combo","c_int","u_int","s_int","cure","lambda","k")] <- NA
    }
    
    if(no.covar==0){
      if(no.covar.c==0 & no.covar.u==1){
        out.combo[,c("Link","allvar.combo","stratum.combo","covar.combo","c_int",covar.vec.c,"u_int","s_int","cure","lambda","k")] <- NA
      }
      if(no.covar.c==1 & no.covar.u==0){
        out.combo[,c("Link","allvar.combo","stratum.combo","covar.combo","c_int","u_int",covar.vec.u,"s_int","cure","lambda","k")] <- NA
      }
      if(no.covar.c==0 & no.covar.u==0){  
        out.combo[,c("Link","allvar.combo","stratum.combo","covar.combo","c_int",covar.vec.c,"u_int",covar.vec.u,"s_int","cure","lambda","k")] <- NA
      }
    }
    
    out.combo$Link <- csdata$Link[1]
    ## create matrix for estimates 
    est.matrix <- csdata[,c(2:(cint.pos-1))]
    est.matrix[,"stratum.combo"] <- stratum.combo
    est.matrix[,"rnames.cs"] <- rnames.cs
    
    
    out.combo$stratum.combo <- apply(data.frame(out.combo[,stratum.vec]), 1, x.combo)
    out.combo$allvar.combo <- apply(data.frame(out.combo[,allvar.vec]), 1, x.combo)
    
    for (krow in 1:dim(out.combo)[1]){
      
      start.row <- match(out.combo$stratum.combo[krow],est.matrix$stratum.combo)
      if(no.cure==0){
        cus.vec <- est.matrix$estimate[start.row:(start.row+sint.pos-cint.pos)]
      }
      if(no.cure==1){
        cus.vec <- est.matrix$estimate[start.row:(start.row+sint.pos-uint.pos)]
        cus.vec <- c(0,cus.vec)
      }
      cus.order <- cnames.cs.new[cint.pos:sint.pos]
      out.combo$c_int[krow] <- cus.vec[which(cus.order=="c_int")]
      out.combo$u_int[krow] <- cus.vec[which(cus.order=="u_int")]
      out.combo$s_int[krow] <- cus.vec[which(cus.order=="s_int")]
      
      c.cols <- c("c_int")
      u.cols <- c("u_int")  
      if(no.covar==0){
        out.combo$covar.combo[krow] <- paste(out.combo[krow,covar.vec],collapse="_", sep="")
        if (no.covar.c==0){
          for(icovar.c in 1:ncovar.c){
            covar.i.c.match <- match(covar.c[icovar.c],covar)
            covar.i.c.match.str <- paste("covar_",covar.i.c.match,sep="")  
            covar.i.c <- covar.vec.c[icovar.c]
            covar.i.c.str <- paste(covar.i.c,"_",out.combo[krow,covar.i.c.match.str],sep="")
            covar.i.c.str <- gsub("covar_","",covar.i.c.str)
            out.combo[krow,covar.i.c]<-cus.vec[which(cus.order==covar.i.c.str)]
          }
          c.cols <- c("c_int",covar.vec.c)
        }
        if (no.covar.u==0){
          for(icovar.u in 1:ncovar.u){
            covar.i.u.match <- match(covar.u[icovar.u],covar)
            covar.i.u.match.str <- paste("covar_",covar.i.u.match,sep="")  
            covar.i.u <- covar.vec.u[icovar.u]
            covar.i.u.str <- paste(covar.i.u,"_",out.combo[krow,covar.i.u.match.str],sep="")
            covar.i.u.str <- gsub("covar_","",covar.i.u.str)
            out.combo[krow,covar.i.u]<-cus.vec[which(cus.order==covar.i.u.str)]
          }
          u.cols <- c("u_int",covar.vec.u)
        }
        
      }
      
      if(no.cure==0){
        out.combo[krow,"cure"]<-1/(1+exp(-sum(out.combo[krow,c.cols])))
      }
      if(no.cure==1){
        out.combo[krow,"cure"]<-0
      }
      out.combo[krow,"lambda"]<-exp(sum(out.combo[krow,u.cols]))
      out.combo[krow,"k"]<-exp(-out.combo[krow,"s_int"])
    }
    
    cols.merge.seer<- cols.keep.seer[which(!cols.keep.seer %in% c(stage.dist.name,surv.name))]
    cols.merge.int <- c(allvar,"Interval")
    cols.merge.seerint <- c(allvar.seer,"Interval")
    
    seer.out.combo[,"allvar.seer.combo"] <- apply(data.frame(seer.out.combo[,allvar.seer]), 1, x.combo)
    seer.out.combo[,"seerint.combo"] <- apply(data.frame(seer.out.combo[,cols.merge.seerint]), 1, x.combo)
    if(no.stratum==0 | no.covar==0){
      seer.out.combo[,"int.combo"] <- apply(data.frame(seer.out.combo[,cols.merge.int]), 1, x.combo)
    }
    if(no.stratum==1 & no.covar==1){
      seer.out.combo[,"int.combo"]<-paste(0,"_",seer.out.combo[,"Interval"],sep="")
    } 
    
    seerdata.out[,"merge.combo"] <- apply(data.frame(seerdata.out[,cols.merge.seer]), 1, x.combo)
    seerdata.out[,"seerint.combo"] <- apply(data.frame(seerdata.out[,cols.merge.seerint]), 1, x.combo)
    seerdata.dist.out[,"merge.combo"] <- apply(data.frame(seerdata.dist.out[,cols.merge.seer]), 1, x.combo)
    
    seerdata.out[,"obs_dist_surv"]<-NA
    rows.match <- match(seerdata.out$merge.combo, seerdata.dist.out$merge.combo)
    seerdata.out$obs_dist_surv[which(!is.na(rows.match))]<-seerdata.dist.out$obs_dist_surv[rows.match[which(!is.na(rows.match))]]
    
    seer.out.combo[,"obs_surv"]<-NA
    seer.out.combo[,"obs_dist_surv"]<-NA
    rows.match <- match(seer.out.combo$seerint.combo, seerdata.out$seerint.combo)
    seer.out.combo$obs_surv[which(!is.na(rows.match))]<-seerdata.out$obs_surv[rows.match[which(!is.na(rows.match))]]
    seer.out.combo$obs_dist_surv[which(!is.na(rows.match))]<-seerdata.out$obs_dist_surv[rows.match[which(!is.na(rows.match))]]
    
    out.combo[,"int.combo"] <- paste(out.combo$allvar.combo,"_",out.combo$fup,sep="")
    out.combo[,"ID.out"] <- 1:dim(out.combo)[1]
    seer.out.combo[,"ID"] <- 1:dim(seer.out.combo)[1]
    
    out.merge <- merge(out.combo,seer.out.combo,by="int.combo")
    out.merge.sort <- out.merge[order(as.numeric(out.merge$ID)),]
    
    theta.sum[,allvar.seer]<-NULL
    out.theta <- merge(out.merge.sort,theta.sum,by="allvar.seer.combo")
    out <- out.theta[order(as.numeric(out.theta$ID)),]
    
    ###################################################################
    
    ## calculate survival 
    if (out$Link[1]=="Weibull"){
      out[,"surv_notcured"] <- with(out,exp(-(fup/lambda)**k))
      out[,"s1_analytical"] <- with(out,surv_notcured*(1-k/theta/lambda*(fup/lambda)**(k-1)))
      out[,"d_mu_int"] <- with(out,(1-cure)*(k*s1_analytical*(fup/lambda)**k+exp(-(fup/lambda)**k)*k*k/theta/lambda*(fup/lambda)**(k-1)))
      out[,"d_sigma"]  <- with(out,(1-cure)*(s1_analytical*k*(fup/lambda)**k*log(fup/lambda)+exp(-(fup/lambda)**k)*k/theta*lambda**(-k)*fup**(k-1)*(1+k*log(fup/lambda))))
      out[,"d_theta"]  <- with(out,(1-cure)*exp(-(fup/lambda)**k)* k/theta**2/lambda*(fup/lambda)**(k-1))
      out[,"median_surv_notcured"]<-with(out,lambda * log(2)**(1/k))
    }
    if (out$Link[1]=="Loglogistic"){
      out[,"surv_notcured"] <- with(out,1/(1+(fup/lambda)**k))
      out[,"s1_analytical"] <- with(out,(1+(fup/lambda)**k-k/theta/lambda*(fup/lambda)**(k-1))/(1+(fup/lambda)**k)**2)
      out[,"d_mu_int"] <- with(out,(1-cure)*k*(s1_analytical*2*(fup/lambda)**k/(1+(fup/lambda)**k)-s1_analytical+1/(1+(fup/lambda)**k)**2))
      out[,"d_sigma"]  <- with(out,(1-cure)*k*(fup/lambda)**k/(1+(fup/lambda)**k)*((1+k*log(fup/lambda)-theta*fup*log(fup/lambda))/theta/fup/(1+(fup/lambda)**k)+2*s1_analytical*log(fup/lambda)))
      out[,"d_theta"]  <- with(out,(1-cure)*k/lambda/theta**2*(fup/lambda)**(k-1)/(1+(fup/lambda)**k)**2)
      out[,"median_surv_notcured"]<-with(out,lambda)
    }
    
    out[,"surv_curemodel"] <- with(out,cure+(1-cure)*surv_notcured)
    d.c.cols  <- gsub("c_","d_c_",c.cols)
    d.mu.cols <- gsub("u_","d_mu_",u.cols)
    out[,d.c.cols]  <- with(out,cure*(1-cure)*(1-s1_analytical))  ## different format from sas formula
    out[,d.mu.cols] <- out$d_mu_int 
    
    cov.matrix <- csdata[,c(cint.pos:sint.pos)]
    cov.matrix[,"stratum.combo"] <- stratum.combo
    cov.matrix[,"estimate"] <- NULL
    
    for (k in 1:dim(out)[1]){
      stratum.combo.k <- out$stratum.combo[k]
      stratum.value.k <- out[k,stratum.vec]
      
      cols.keep <- c("c_int","u_int","s_int")
      
      if(no.covar==0){
        covar.value.k <- out[k,covar.vec]
        if(no.covar.c==0){
          covar.k.c.match <- match(covar.c,covar)
          covar.value.k.c <- covar.value.k[covar.k.c.match]
          covar.cols.c <- paste(covar.vec.c,"_",covar.value.k.c,sep="")
          covar.cols.c <- gsub("covar_","",covar.cols.c)
        }
        if(no.covar.u==0){
          covar.k.u.match <- match(covar.u,covar)
          covar.value.k.u <- covar.value.k[covar.k.u.match]
          covar.cols.u <- paste(covar.vec.u,"_",covar.value.k.u,sep="")
          covar.cols.u <- gsub("covar_","",covar.cols.u)
        }
        if(no.covar.c==0 & no.covar.u==0){
          cols.keep <- c("c_int",covar.cols.c,"u_int",covar.cols.u,"s_int")
        }
        if(no.covar.c==1 & no.covar.u==0){
          cols.keep <- c("c_int","u_int",covar.cols.u,"s_int")
        }
        if(no.covar.c==0 & no.covar.u==1){
          cols.keep <- c("c_int",covar.cols.c,"u_int","s_int")
        }
      }
      
      if(no.cure==1){
        cols.keep <- cols.keep[-which(cols.keep=="c_int")]
      }
      rows.k <- which(cov.matrix$stratum.combo==stratum.combo.k)
      if(no.cure==0){
        cov.matrix.k <- cov.matrix[rows.k,1:(dim(cov.matrix)[2]-1)]
      }
      if(no.cure==1){
        cov.matrix.k <- cov.matrix[rows.k,2:(dim(cov.matrix)[2]-1)]
      }
      cols.keep.pos <- which(colnames(cov.matrix.k) %in% cols.keep)
      covariance <- cov.matrix.k[cols.keep.pos,cols.keep.pos]
      row.0 <- rep(0,dim(covariance)[1])
      covariance.k <- rbind(covariance,row.0)
      col.0 <- rep(0,dim(covariance.k)[2]+1)
      covariance.k <- cbind(covariance.k,col.0)
      covariance.k[dim(covariance.k)[1],dim(covariance.k)[1]]<-out$theta.se[k]^2
      if(no.cure==0){
        d.vec.k <- as.numeric(out[k,c(d.c.cols,d.mu.cols,"d_sigma","d_theta")])
      }
      if(no.cure==1){
        d.vec.k <- as.numeric(out[k,c(d.mu.cols,"d_sigma","d_theta")])
      }
      varcov.k<-t(d.vec.k)%*%as.matrix(covariance.k)%*%d.vec.k
      out[k,"var_CI_analytical"] <- varcov.k
      out[k,"se_CI_analytical"]  <- sqrt(varcov.k)
    }
    
    out[,c(c.cols,u.cols,"s_int","allvar.combo","stratum.combo","covar.combo","int.combo","ID","ID.out")]<-NULL
    out[,c(d.c.cols,d.mu.cols,"d_sigma","d_theta")]<-NULL
    
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
    
    colnames(out)[which(colnames(out)=="Link")] <- "link"
    colnames(out)[which(colnames(out)=="fup")] <- "followup"
    othervar<-allvar.seer[which(!allvar.seer %in% allvar)]
    
    #out <- out[which(out$Interval<=int.max.out),]
    out.cols.keep<-c("followup","link","r","cure","lambda","k","theta",
                     "surv_curemodel","surv_notcured","median_surv_notcured",
                     "s1_numerical","G_numerical","CI_numerical",
                     "s1_analytical","G_analytical","CI_analytical","se_CI_analytical",
                     "obs_surv","obs_dist_surv")
    
    if(no.stratum==0 | no.covar==0){
      out <- out[,c(allvar,othervar,out.cols.keep)]
    }
    if(no.stratum==1 & no.covar==1){
      out <- out[,c(allvar.seer,out.cols.keep)]
    }  
    
    return(out)
  }
