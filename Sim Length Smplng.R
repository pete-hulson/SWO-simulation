#########################
# Script to estimate population #'s at size for GOA/AI species

#########################
# Set directories and create sub-folders
path<-getwd()
pathD<-paste0(path,"/Data")
pathR<-paste0(path,"/Results")
pathF<-paste0(path,"/Functions")
source(paste0(pathF,"/est_SzPop.R"))

#########################
# Read in data
lfreq_all<-read.csv(paste0(pathD,"/lfreq.csv"))
sizepop_RACE_all<-read.csv(paste0(pathD,"/sizepop.csv"))
cpue_all<-read.csv(paste0(pathD,"/CPUE.csv"))
strata_all<-read.csv(paste0(pathD,"/strata.csv"))

# Define region and subsequent data to use
region<-"GOA"
#region<-"AI"

lfreq_reg<-subset(lfreq_all,lfreq_all$SURVEY==region)
sizepop_RACE_reg<-subset(sizepop_RACE_all,sizepop_RACE_all$SURVEY==region)
cpue_reg<-subset(cpue_all,cpue_all$SURVEY==region)
strata_reg<-subset(strata_all,strata_all$SURVEY==region)

species<-sort(unique(lfreq_reg$SPECIES_CODE))
yrs<-sort(unique(lfreq_reg$YEAR))
stratum<-sort(unique(strata_reg$STRATUM))

# Create a folder for species-specific results
if(file.exists(paste0(pathR,"/",region))==FALSE){
  dir.create(paste0(pathR,"/",region))}

# Set desired max sex sample sizes
max_sx<-c(50,75,100,125,150)

# Set number of iterations
I<-10

#########################
# Loop thru species
for(sp in 1:length(species)){
# Loop for just pollock
#for(sp in 8:8){
# Loop for just POP
#for(sp in length(species):length(species)){
  
  # Start time
  st_time<-Sys.time()
  
  # Create a folder for species-specific results
  if(file.exists(paste0(pathR,"/",region,"/",species[sp]))==FALSE){
    dir.create(paste0(pathR,"/",region,"/",species[sp]))}
  
  # Subset to species-specific data and parameters
  lfreq_sp<-subset(lfreq_reg,lfreq_reg$SPECIES_CODE==species[sp])
  lengths<-sort(unique(lfreq_sp$LENGTH)) 
  sizepop_RACE_sp<-subset(sizepop_RACE_reg,sizepop_RACE_reg$SPECIES_CODE==species[sp])
  cpue_sp<-subset(cpue_reg,cpue_reg$SPECIES_CODE==species[sp])
  
  # Set up results matrices
  SZPOP_M<-matrix(nrow=length(yrs),ncol=length(lengths))
  colnames(SZPOP_M)<-as.character(lengths)
  rownames(SZPOP_M)<-as.character(yrs)
  SZPOP_F<-matrix(nrow=length(yrs),ncol=length(lengths))
  colnames(SZPOP_F)<-as.character(lengths)
  rownames(SZPOP_F)<-as.character(yrs)
  SZPOP_U<-matrix(nrow=length(yrs),ncol=length(lengths))
  colnames(SZPOP_U)<-as.character(lengths)
  rownames(SZPOP_U)<-as.character(yrs)
  SZPOP_TOT<-matrix(nrow=length(yrs),ncol=length(lengths))
  colnames(SZPOP_TOT)<-as.character(lengths)
  rownames(SZPOP_TOT)<-as.character(yrs)
  
  # Set up results matrices for simulation statistics
  Stats_F_iter<-array(NA,dim=c(I,length(yrs),length(max_sx)))
  colnames(Stats_F_iter)<-yrs
  rownames(Stats_F_iter)<-seq(1,I)
  Stats_M_iter<-array(NA,dim=c(I,length(yrs),length(max_sx)))
  colnames(Stats_M_iter)<-yrs
  rownames(Stats_M_iter)<-seq(1,I)
  Stats_TOT_iter<-array(NA,dim=c(I,length(yrs),length(max_sx)))
  colnames(Stats_TOT_iter)<-yrs
  rownames(Stats_TOT_iter)<-seq(1,I) 
  
  # Set up results matrices for other interesting output
  hl_prop<-matrix(nrow=length(yrs),ncol=length(max_sx))
  colnames(hl_prop)<-max_sx
  rownames(hl_prop)<-yrs
  sampl_redux<-matrix(nrow=length(yrs),ncol=length(max_sx))
  colnames(sampl_redux)<-max_sx
  rownames(sampl_redux)<-yrs
  
  #########################
  # Loop thru max sex sample size scenarios
  for(ss in 1:length(max_sx)){  
    
    #########################
    # Loop thru all years
    #for(y in 1:length(yrs)){
    # Loop thru last 3 years
    for(y in (length(yrs)-2):length(yrs)){
      
      # Subset to yearly data and parameters
      lfreq_sp_y<-subset(lfreq_sp,lfreq_sp$YEAR==yrs[y])
      cpue_sp_y<-subset(cpue_sp,cpue_sp$YEAR==yrs[y])
      sizepop_RACE_sp_y<-subset(sizepop_RACE_sp,sizepop_RACE_sp$YEAR==yrs[y])
      
      # Estimate expanded length comps with full dataset
      sp_SzPop_yr<-est_SzPop(stratum,lengths,strata_reg,cpue_sp_y,lfreq_sp_y)
      SZPOP_M[y,]<-sp_SzPop_yr$MALES
      SZPOP_F[y,]<-sp_SzPop_yr$FEMALES
      SZPOP_U[y,]<-sp_SzPop_yr$UNSEXED
      SZPOP_TOT[y,]<-sp_SzPop_yr$MALES+sp_SzPop_yr$FEMALES+sp_SzPop_yr$UNSEXED
      
      # Set up results matrices for resampled estimated expanded length comps
      SZPOP_M_iter<-matrix(nrow=I,ncol=length(lengths))
      colnames(SZPOP_M_iter)<-as.character(lengths)
      rownames(SZPOP_M_iter)<-seq(1,I)
      SZPOP_F_iter<-matrix(nrow=I,ncol=length(lengths))
      colnames(SZPOP_F_iter)<-as.character(lengths)
      rownames(SZPOP_F_iter)<-seq(1,I)
      SZPOP_U_iter<-matrix(nrow=I,ncol=length(lengths))
      colnames(SZPOP_U_iter)<-as.character(lengths)
      rownames(SZPOP_U_iter)<-seq(1,I)
      SZPOP_TOT_iter<-matrix(nrow=I,ncol=length(lengths))
      colnames(SZPOP_TOT_iter)<-as.character(lengths)
      rownames(SZPOP_TOT_iter)<-seq(1,I)
      
      # Determine which hauls exceeded desired max sex sample size and keep track of proportion
      hls<-as.numeric(names(which(tapply(subset(lfreq_sp_y$FREQUENCY,lfreq_sp_y$SEX!=3),subset(lfreq_sp_y$HAULJOIN,lfreq_sp_y$SEX!=3),sum)>max_sx[ss])))
      hls_tot<-unique(lfreq_sp_y$HAULJOIN)
      hl_prop[y,ss]<-length(hls)/length(hls_tot)
      
      # Determine sexed sample reduction
      sampl_redux[y,ss]<-sum(tapply(subset(lfreq_sp_y,lfreq_sp_y$HAULJOIN %in% hls)$FREQUENCY,subset(lfreq_sp_y,lfreq_sp_y$HAULJOIN %in% hls)$HAULJOIN,sum)-max_sx[ss])
      
      #########################
      # Loop thru iterations
      for(i in 1:I){
        
        # keep track of where simulation is at
        print(paste0("yr = ",y," ss = ",ss," iter = ",i))
        
        # Setup dummy data frame
        lfreq_sp_y_hl_smpld<-NULL
        
        #########################
        # Loop thru hauls w/ >100 samples
        for (h in 1:length(hls)){
          
          # Subset length freq to haul level
          lfreq_sp_y_hl<-subset(lfreq_sp_y,lfreq_sp_y$HAULJOIN==hls[h])
          
          # Start of code to sample #'s based on binomial and multinomial (keeping this in here for now in case I want to take a looksie later)
          #p_F<-sum(lfreq_sp_y_hl$FREQUENCY[which(lfreq_sp_y_hl$SEX==2)])/sum(lfreq_sp_y_hl$FREQUENCY[which(lfreq_sp_y_hl$SEX!=3)])
          #SampSz_F<-rbinom(I,max_sx,p_F)
          #SampSz_M<-max_sx-SampSz_F
          #M_n<-tapply(subset(lfreq_sp_y_hl$FREQUENCY,lfreq_sp_y_hl$SEX==1),subset(lfreq_sp_y_hl$LENGTH,lfreq_sp_y_hl$SEX==1),sum)
          #M_N<-sum(subset(lfreq_sp_y_hl$FREQUENCY,lfreq_sp_y_hl$SEX==1))
          #M_comp<-M_n/M_N
          #F_n<-tapply(subset(lfreq_sp_y_hl$FREQUENCY,lfreq_sp_y_hl$SEX==2),subset(lfreq_sp_y_hl$LENGTH,lfreq_sp_y_hl$SEX==2),sum)
          #F_N<-sum(subset(lfreq_sp_y_hl$FREQUENCY,lfreq_sp_y_hl$SEX==2))
          #F_comp<-F_n/F_N
          #LengthSamp_F<-rmultinom(I,SampSz_F,F_comp)
          #LengthSamp_M<-rmultinom(I,SampSz_M,M_comp)     
          
          #########################
          # Resample haul at max # for sex determination
          
          # Subset out data that has sexes determined and that don't
          lfreq_sp_y_hl_MF<-subset(lfreq_sp_y_hl,lfreq_sp_y_hl$SEX!=3)
          lfreq_sp_y_hl_U<-subset(lfreq_sp_y_hl,lfreq_sp_y_hl$SEX==3)
          
          # Set up vector of sampled lengths to then resample
          Ln_vec<-rep(lfreq_sp_y_hl_MF$LENGTH,times=lfreq_sp_y_hl_MF$FREQUENCY)
          Sx_vec<-rep(lfreq_sp_y_hl_MF$SEX,times=lfreq_sp_y_hl_MF$FREQUENCY)
          LnSx<-as.data.frame(cbind(Ln_vec,Sx_vec))
          
          # Resample sexed lengths at max sample size desired
          rws<-sample.int(length(LnSx$Sx_vec),max_sx[ss],replace=FALSE)
          
          # With those that are sampled and sexed, get organized into table
          LnSx_smpl<-LnSx[rws,]
          LnSx_smpl_M<-subset(LnSx_smpl,LnSx_smpl$Sx_vec==1)
          LnSx_smpl_F<-subset(LnSx_smpl,LnSx_smpl$Sx_vec==2)
          lfreq_smpl_M<-tapply(LnSx_smpl_M$Sx_vec,LnSx_smpl_M$Ln_vec,length)
          lfreq_smpl_F<-tapply(LnSx_smpl_F$Sx_vec,LnSx_smpl_F$Ln_vec,length)
          
          # For those left over that weren't sampled to be sexed, put into 'unsexed' table
          LnSx_notsmpl<-LnSx[-rws,]
          lfreq_smpl_U<-tapply(LnSx_notsmpl$Ln_vec,LnSx_notsmpl$Ln_vec,length)
          
          # Put the F/M/U data together (with if statements if there's a haul that only samples one sex)
          lfreq_smpl_raw<-as.matrix(rbind(
            if(length(lfreq_smpl_F)){cbind(lfreq_smpl_F,rep(2,len.out=length(lfreq_smpl_F)))},
            if(length(lfreq_smpl_M)){cbind(lfreq_smpl_M,rep(1,len.out=length(lfreq_smpl_M)))},
            cbind(lfreq_smpl_U,rep(3,len.out=length(lfreq_smpl_U)))))
          
          # Put it all together into a data frame that mimics RACE length frequency dataset
          for(j in 1:length(lfreq_smpl_raw[,1])){
            dummy_row<-lfreq_sp_y_hl[1,]
            dummy_row$LENGTH<-as.numeric(rownames(lfreq_smpl_raw)[j])
            dummy_row$FREQUENCY<-lfreq_smpl_raw[j,1]
            dummy_row$SEX<-lfreq_smpl_raw[j,2]
            lfreq_sp_y_hl_smpld<-rbind(lfreq_sp_y_hl_smpld,dummy_row)
          }
          
          # Don't forget to add back in the samples that were originally unsexed
          lfreq_sp_y_hl_smpld<-rbind(lfreq_sp_y_hl_smpld,lfreq_sp_y_hl_U)
          
          #########################
          # End resampling haul loop
        }
        
        # Replace actual hauls with resampled hauls
        lfreq_sp_y_smpld<-lfreq_sp_y[-which(lfreq_sp_y$HAULJOIN %in% hls),]
        lfreq_sp_y_smpld<-rbind(lfreq_sp_y_smpld,lfreq_sp_y_hl_smpld)
        
        # Estimate iteration's pop'n #'s at size
        sp_SzPop_iter<-est_SzPop(stratum,lengths,strata_reg,cpue_sp_y,lfreq_sp_y_smpld)
        SZPOP_M_iter[i,]<-sp_SzPop_iter$MALES
        SZPOP_F_iter[i,]<-sp_SzPop_iter$FEMALES
        SZPOP_U_iter[i,]<-sp_SzPop_iter$UNSEXED
        SZPOP_TOT_iter[i,]<-sp_SzPop_iter$MALES+sp_SzPop_iter$FEMALES+sp_SzPop_iter$UNSEXED
        
        # Compute statistics (ESS/SSQ/etc)
        full_comp_F<-SZPOP_F[y,]/sum(SZPOP_F[y,])
        iter_comp_F<-SZPOP_F_iter[i,]/sum(SZPOP_F_iter[i,])
        Stats_F_iter[i,y,ss]<-sum(iter_comp_F*(1-iter_comp_F))/sum((iter_comp_F-full_comp_F)^2)
        
        full_comp_M<-SZPOP_M[y,]/sum(SZPOP_M[y,])
        iter_comp_M<-SZPOP_M_iter[i,]/sum(SZPOP_M_iter[i,])
        Stats_M_iter[i,y,ss]<-sum(iter_comp_M*(1-iter_comp_M))/sum((iter_comp_M-full_comp_M)^2)
        
        full_comp_TOT<-SZPOP_TOT[y,]/sum(SZPOP_TOT[y,])
        iter_comp_TOT<-SZPOP_TOT_iter[i,]/sum(SZPOP_TOT_iter[i,])
        Stats_TOT_iter[i,y,ss]<-sum(iter_comp_TOT*(1-iter_comp_TOT))/sum((iter_comp_TOT-full_comp_TOT)^2)
        
        #########################
        # End iteration loop
      }
      
      #########################
      # End year loop
    }
    # Write stats results matrices
    write.csv(Stats_F_iter[,,ss],paste0(pathR,"/",region,"/",species[sp],"/Stats_F_iter_mx",max_sx[ss],"_",species[sp],".csv"))
    write.csv(Stats_M_iter[,,ss],paste0(pathR,"/",region,"/",species[sp],"/Stats_M_iter_mx",max_sx[ss],"_",species[sp],".csv"))
    write.csv(Stats_TOT_iter[,,ss],paste0(pathR,"/",region,"/",species[sp],"/Stats_TOT_iter_mx",max_sx[ss],"_",species[sp],".csv"))
    
    #########################
    # End max sex sample size loop
  }
  
  # Write results matrices
  write.csv(hl_prop,paste0(pathR,"/",region,"/",species[sp],"/hl_prop_",species[sp],".csv"))
  write.csv(sampl_redux,paste0(pathR,"/",region,"/",species[sp],"/sampl_redux_",species[sp],".csv"))
  
  # End time
  end_time<-Sys.time()
  
  # Compute species-spefic timet o run sim
  run_time<-end_time-st_time
  run_time
  
  #########################
  # End species loop 
}









