########################################################################################
# Function to estimate pop'n @ size
########################################################################################

est_SzPop_strata<-function(stratum,lengths,strata_dat,cpue_dat,lfreq_dat){
  
  # Set up results matrices
  SZPOP_M_st<-matrix(nrow=length(stratum),ncol=length(lengths))
  rownames(SZPOP_M_st)<-stratum
  colnames(SZPOP_M_st)<-lengths
  SZPOP_F_st<-matrix(nrow=length(stratum),ncol=length(lengths))
  rownames(SZPOP_F_st)<-stratum
  colnames(SZPOP_F_st)<-lengths
  SZPOP_U_st<-matrix(nrow=length(stratum),ncol=length(lengths))
  rownames(SZPOP_U_st)<-stratum
  colnames(SZPOP_U_st)<-lengths
  
  strata_vec<-NULL
  summ_area_vec<-NULL
  
  #########################
  # Loop thru strata
  for(st in 1:length(stratum)){
    
    # Subset data to strata level
    strata_st<-subset(strata_dat,strata_dat$STRATUM==stratum[st])           
    cpue_st<-subset(cpue_dat,cpue_dat$STRATUM==stratum[st])
    hls_cpue<-unique(cpue_st$HAULJOIN)
    lfreq_st<-subset(lfreq_dat,lfreq_dat$HAULJOIN %in% hls_cpue)
    
    # Subset data to sex-specific (M=males, F=females, U=unsexed)
    lfreq_M_st<-subset(lfreq_st,lfreq_st$SEX==1)
    lfreq_F_st<-subset(lfreq_st,lfreq_st$SEX==2)
    lfreq_U_st<-subset(lfreq_st,lfreq_st$SEX==3)
    
    # Determine number of hauls of catch with lengths
    count<-length(unique(c(lfreq_M_st$HAULJOIN,lfreq_F_st$HAULJOIN,lfreq_U_st$HAULJOIN)))
    
    # Identify hauls with catch but no lengths
    hls_l<-unique(c(lfreq_M_st$HAULJOIN,lfreq_F_st$HAULJOIN,lfreq_U_st$HAULJOIN))
    hls_c<-cpue_st$HAULJOIN[which(is.na(cpue_st$CATCHJOIN)==FALSE)]
    hls_nol<-hls_c[which(is.na(match(hls_c,hls_l)==TRUE))]
    
    # Calc pop'n #'s in strata
    st_num<-mean(cpue_st$NUMCPUE)*strata_st$AREA
    
    # Calc CPUE ratio among hauls
    cprat<-tapply(cpue_st$NUMCPUE,cpue_st$HAULJOIN,mean)/sum(cpue_st$NUMCPUE)
    
    # Calc Total lengths sampled by haul
    n_st<-tapply(lfreq_st$FREQUENCY,lfreq_st$HAULJOIN,sum)
    
    # Calc sex-specific numbers at length
    n_h_M<-tapply(lfreq_M_st$FREQUENCY,list(lfreq_M_st$HAULJOIN,lfreq_M_st$LENGTH),sum)
    n_h_M[is.na(n_h_M)] <- 0
    n_h_F<-tapply(lfreq_F_st$FREQUENCY,list(lfreq_F_st$HAULJOIN,lfreq_F_st$LENGTH),sum)
    n_h_F[is.na(n_h_F)] <- 0
    n_h_U<-tapply(lfreq_U_st$FREQUENCY,list(lfreq_U_st$HAULJOIN,lfreq_U_st$LENGTH),sum)
    n_h_U[is.na(n_h_U)] <- 0
    
    # Sex-specific ratio of total
    ratio_h_M<-n_h_M/as.vector(n_st[match(as.numeric(rownames(n_h_M)),as.numeric(names(n_st)))])
    ratio_h_F<-n_h_F/as.vector(n_st[match(as.numeric(rownames(n_h_F)),as.numeric(names(n_st)))])
    ratio_h_U<-n_h_U/as.vector(n_st[match(as.numeric(rownames(n_h_U)),as.numeric(names(n_st)))])
    
    # Estimate size comp for hauls with catch that did not sample lengths
    if(length(hls_nol)>0){
      
      ratio_h_M_unk<-colSums(ratio_h_M)/count
      ratio_h_F_unk<-colSums(ratio_h_F)/count
      ratio_h_U_unk<-colSums(ratio_h_U)/count
      total<-sum(ratio_h_M_unk,ratio_h_F_unk,ratio_h_U_unk)
      ratio_h_M_unk<-ratio_h_M_unk/total
      ratio_h_F_unk<-ratio_h_F_unk/total
      ratio_h_U_unk<-ratio_h_U_unk/total
      
      # Add unkown size com hauls to sex-specific ratio of total
      ratio_h_M_unk_add<-matrix(ratio_h_M_unk,nrow=length(hls_nol),ncol=length(ratio_h_M_unk),byrow=TRUE)
      rownames(ratio_h_M_unk_add)<-hls_nol
      colnames(ratio_h_M_unk_add)<-colnames(ratio_h_M)
      ratio_h_M<-rbind(ratio_h_M,ratio_h_M_unk_add)
      ratio_h_F_unk_add<-matrix(ratio_h_F_unk,nrow=length(hls_nol),ncol=length(ratio_h_F_unk),byrow=TRUE)
      rownames(ratio_h_F_unk_add)<-hls_nol
      colnames(ratio_h_F_unk_add)<-colnames(ratio_h_F)
      ratio_h_F<-rbind(ratio_h_F,ratio_h_F_unk_add)
      ratio_h_U_unk_add<-matrix(ratio_h_U_unk,nrow=length(hls_nol),ncol=length(ratio_h_U_unk),byrow=TRUE)
      rownames(ratio_h_U_unk_add)<-hls_nol
      colnames(ratio_h_U_unk_add)<-colnames(ratio_h_U)
      ratio_h_U<-rbind(ratio_h_U,ratio_h_U_unk_add)
      
      # End if statements 
    }
    
    # Put it all together to get numbers-at-sex-at-length by strata, and put it in results matrix
    szpop_M<-round(colSums(ratio_h_M*as.vector(cprat[match(as.numeric(rownames(ratio_h_M)),as.numeric(names(cprat)))])*st_num),digits=0)
    SZPOP_M_st[st,match(as.numeric(names(szpop_M)),lengths)]<-szpop_M
    szpop_F<-round(colSums(ratio_h_F*as.vector(cprat[match(as.numeric(rownames(ratio_h_F)),as.numeric(names(cprat)))])*st_num),digits=0)
    SZPOP_F_st[st,match(as.numeric(names(szpop_F)),lengths)]<-szpop_F
    szpop_U<-round(colSums(ratio_h_U*as.vector(cprat[match(as.numeric(rownames(ratio_h_U)),as.numeric(names(cprat)))])*st_num),digits=0)
    SZPOP_U_st[st,match(as.numeric(names(szpop_U)),lengths)]<-szpop_U
    
    strata_vec<-c(strata_vec,rep(stratum[st],length.out=length(lengths)))
    summ_area_vec<-c(summ_area_vec,rep(strata_dat$SUMMARY_AREA[st],length.out=length(lengths)))
    
    # End stratum loop
  }
  
  
  # Now put it all together
  SZPOP_M_st[is.na(SZPOP_M_st)] <- 0    
  SZPOP_F_st[is.na(SZPOP_F_st)] <- 0  
  SZPOP_U_st[is.na(SZPOP_U_st)] <- 0  
  SzPop<-cbind(strata_vec,summ_area_vec,rep.int(lengths,times=length(stratum)),
               c(t(SZPOP_M_st)),
               c(t(SZPOP_F_st)),
               c(t(SZPOP_U_st)))
  colnames(SzPop)<-c("STRATUM","SUMMARY_AREA","LENGTH","MALES","FEMALES","UNSEXED")
  SzPop<-as.data.frame(SzPop)
  
  # End function
  SzPop}
