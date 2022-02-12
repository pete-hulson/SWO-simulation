#########################
# Script to estimate population #'s at size for GOA/AI species

#########################
# Set directories and create sub-folders
path<-getwd()
pathD<-paste0(path,"/Data")
pathR<-paste0(path,"/Results")
pathF<-paste0(path,"/Functions")
source(paste0(pathF,"/est_SzPop.R"))
source(paste0(pathF,"/est_SzPop_strata.R"))

#########################
# Read in data
lfreq_all<-read.csv(paste0(pathD,"/lfreq.csv"))
cpue_all<-read.csv(paste0(pathD,"/CPUE.csv"))
strata_all<-read.csv(paste0(pathD,"/strata.csv"))

# Define region and subsequent data to use
region<-"GOA"
#region<-"AI"

lfreq_reg<-subset(lfreq_all,lfreq_all$SURVEY==region)
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

# Get results matrices set up
# Raw length frequencies for example hauls to compare simulation to actual
Examp_hls_nms<-c("SURVEY","SPECIES_CODE","MAX_SEX_SS","YEAR","ITERATION","HAULJOIN","SEX","LENGTH","FREQUENCY")
Examp_hls<-matrix(nrow=1,ncol=length(Examp_hls_nms))
colnames(Examp_hls)<-Examp_hls_nms
Examp_hls<-Examp_hls[-1,]
# Simulation estimates of pop'n at size at management area scale
Sizepop_ma_nms<-c("SURVEY","SPECIES_CODE","MAX_SEX_SS","YEAR","ITERATION","LENGTH","MALES","FEMALES","UNSEXED","TOTAL")
Sizepop_ma<-matrix(nrow=1,ncol=length(Sizepop_ma_nms))
colnames(Sizepop_ma)<-Sizepop_ma_nms
Sizepop_ma<-Sizepop_ma[-1,]
# Simulation estimates of pop'n at size at strata scale
Sizepop_strata_nms<-c("SURVEY","SPECIES_CODE","MAX_SEX_SS","YEAR","STRATUM","SUMMARY_AREA","ITERATION","LENGTH","MALES","FEMALES","UNSEXED","TOTAL")
Sizepop_strata<-matrix(nrow=1,ncol=length(Sizepop_strata_nms))
colnames(Sizepop_strata)<-Sizepop_strata_nms
Sizepop_strata<-Sizepop_strata[-1,]
# Statistics to store from simulation
Sim_stats_nms<-c("SURVEY","SPECIES_CODE","MAX_SEX_SS","YEAR","ITERATION","MALE_ESS","FEMALE_ESS","TOTAL_ESS")
Sim_stats<-matrix(nrow=1,ncol=length(Sim_stats_nms))
colnames(Sim_stats)<-Sim_stats_nms
Sim_stats<-Sim_stats[-1,]
# Statistics to store from simulation at strata level
Sim_stats_st_nms<-c("SURVEY","SPECIES_CODE","MAX_SEX_SS","YEAR","STRATUM","ITERATION","MALE_ESS","FEMALE_ESS","TOTAL_ESS")
Sim_stats_st<-matrix(nrow=1,ncol=length(Sim_stats_st_nms))
colnames(Sim_stats_st)<-Sim_stats_st_nms
Sim_stats_st<-Sim_stats_st[-1,]
# Other results to keep track of (% of hauls reduced, # of samples reduced)
Other_res_nms<-c("SURVEY","SPECIES_CODE","MAX_SEX_SS","YEAR","HAUL_PROP","SAMP_REDUX")
Other_res<-matrix(nrow=1,ncol=length(Other_res_nms))
colnames(Other_res)<-Other_res_nms
Other_res<-Other_res[-1,]

#########################
# Loop thru species
#for(sp in 1:length(species)){
# Loop for just pollock
#for(sp in 8:8){
# Loop for just POP
for(sp in length(species):length(species)){
  
  # Start time
  st_time<-Sys.time()
  
  # Subset to species-specific data and parameters
  lfreq_sp<-subset(lfreq_reg,lfreq_reg$SPECIES_CODE==species[sp])
  lengths<-sort(unique(lfreq_sp$LENGTH)) 
  sizepop_RACE_sp<-subset(sizepop_RACE_reg,sizepop_RACE_reg$SPECIES_CODE==species[sp])
  cpue_sp<-subset(cpue_reg,cpue_reg$SPECIES_CODE==species[sp])
  
  #########################
  # Loop thru all years
  #for(y in 1:length(yrs)){
  # Loop thru last 3 years
  #for(y in (length(yrs)-2):length(yrs)){
  # Loop thru last year
  for(y in length(yrs):length(yrs)){
    
    # Subset to yearly data and parameters
    lfreq_sp_y<-subset(lfreq_sp,lfreq_sp$YEAR==yrs[y])
    cpue_sp_y<-subset(cpue_sp,cpue_sp$YEAR==yrs[y])
    sizepop_RACE_sp_y<-subset(sizepop_RACE_sp,sizepop_RACE_sp$YEAR==yrs[y])
    
    # Estimate expanded length comps with full dataset at management area scale
    sp_SzPop_yr<-est_SzPop(stratum,lengths,strata_reg,cpue_sp_y,lfreq_sp_y)
    Sizepop_ma_yr<-cbind(
      rep(region,times=length(sp_SzPop_yr$LENGTH)),
      rep(species[sp],times=length(sp_SzPop_yr$LENGTH)),
      rep(0,times=length(sp_SzPop_yr$LENGTH)),
      rep(yrs[y],times=length(sp_SzPop_yr$LENGTH)),
      rep(0,times=length(sp_SzPop_yr$LENGTH)),
      sp_SzPop_yr,
      sp_SzPop_yr$MALES+sp_SzPop_yr$FEMALES+sp_SzPop_yr$UNSEXED)
    colnames(Sizepop_ma_yr)<-Sizepop_ma_nms
    Sizepop_ma<-rbind(Sizepop_ma,Sizepop_ma_yr)
    
    # Estimate expanded length comps with full dataset at strata scale
    sp_SzPop_yr_strata<-est_SzPop_strata(stratum,lengths,strata_reg,cpue_sp_y,lfreq_sp_y)
    Sizepop_strata_yr<-cbind(
      rep(region,times=length(sp_SzPop_yr_strata$STRATUM)),
      rep(species[sp],times=length(sp_SzPop_yr_strata$STRATUM)),
      rep(0,times=length(sp_SzPop_yr_strata$STRATUM)),
      rep(yrs[y],times=length(sp_SzPop_yr_strata$STRATUM)),
      sp_SzPop_yr_strata$STRATUM,
      sp_SzPop_yr_strata$SUMMARY_AREA,
      rep(0,times=length(sp_SzPop_yr_strata$STRATUM)),
      sp_SzPop_yr_strata$LENGTH,
      sp_SzPop_yr_strata$MALES,
      sp_SzPop_yr_strata$FEMALES,
      sp_SzPop_yr_strata$UNSEXED,
      sp_SzPop_yr_strata$MALES+sp_SzPop_yr_strata$FEMALES+sp_SzPop_yr_strata$UNSEXED)
    colnames(Sizepop_strata_yr)<-Sizepop_strata_nms
    Sizepop_strata<-rbind(Sizepop_strata,Sizepop_strata_yr)
    
    #########################
    # Loop thru max sex sample size scenarios
    for(ss in 1:length(max_sx)){  
      
      # Determine which hauls exceeded desired max sex sample size and keep track of proportion
      hls<-as.numeric(names(which(tapply(subset(lfreq_sp_y$FREQUENCY,lfreq_sp_y$SEX!=3),subset(lfreq_sp_y$HAULJOIN,lfreq_sp_y$SEX!=3),sum)>max_sx[ss])))
      
      # Stupid long way to get the hauls with the smallest and largest samples to look example of sampling sensitivity at haul level
      hls_minmax<-c(as.numeric(names(which(tapply(subset(subset(lfreq_sp_y,lfreq_sp_y$HAULJOIN %in% hls)$FREQUENCY,subset(lfreq_sp_y,lfreq_sp_y$HAULJOIN %in% hls)$SEX!=3),
                                                  subset(subset(lfreq_sp_y,lfreq_sp_y$HAULJOIN %in% hls)$HAULJOIN,subset(lfreq_sp_y,lfreq_sp_y$HAULJOIN %in% hls)$SEX!=3),sum)==
                                             min(tapply(subset(subset(lfreq_sp_y,lfreq_sp_y$HAULJOIN %in% hls)$FREQUENCY,subset(lfreq_sp_y,lfreq_sp_y$HAULJOIN %in% hls)$SEX!=3),
                                                        subset(subset(lfreq_sp_y,lfreq_sp_y$HAULJOIN %in% hls)$HAULJOIN,subset(lfreq_sp_y,lfreq_sp_y$HAULJOIN %in% hls)$SEX!=3),sum))))),
                    as.numeric(names(which(tapply(subset(subset(lfreq_sp_y,lfreq_sp_y$HAULJOIN %in% hls)$FREQUENCY,subset(lfreq_sp_y,lfreq_sp_y$HAULJOIN %in% hls)$SEX!=3),
                                                  subset(subset(lfreq_sp_y,lfreq_sp_y$HAULJOIN %in% hls)$HAULJOIN,subset(lfreq_sp_y,lfreq_sp_y$HAULJOIN %in% hls)$SEX!=3),sum)==
                                             max(tapply(subset(subset(lfreq_sp_y,lfreq_sp_y$HAULJOIN %in% hls)$FREQUENCY,subset(lfreq_sp_y,lfreq_sp_y$HAULJOIN %in% hls)$SEX!=3),
                                                        subset(subset(lfreq_sp_y,lfreq_sp_y$HAULJOIN %in% hls)$HAULJOIN,subset(lfreq_sp_y,lfreq_sp_y$HAULJOIN %in% hls)$SEX!=3),sum))))))
      
      # Store actual sample frequencies for example hauls
      examp_hls_raw<-subset(lfreq_sp_y,lfreq_sp_y$HAULJOIN %in% hls_minmax & lfreq_sp_y$SEX!=3)
      examp_res<-cbind(
        examp_hls_raw$REGION,
        examp_hls_raw$SPECIES_CODE,
        rep(0,times=length(examp_hls_raw$HAULJOIN)),
        examp_hls_raw$YEAR,
        rep(0,times=length(examp_hls_raw$HAULJOIN)),   
        examp_hls_raw$HAULJOIN,
        examp_hls_raw$SEX,
        examp_hls_raw$LENGTH,
        examp_hls_raw$FREQUENCY)
      colnames(examp_res)<-Examp_hls_nms
      Examp_hls<-rbind(Examp_hls,examp_res)
      
      # Keep track of proportion of hauls with samples > max sex SS
      hls_tot<-unique(lfreq_sp_y$HAULJOIN)
      hl_prop<-length(hls)/length(hls_tot)
      
      # Determine sexed sample reduction
      sampl_redux<-sum(tapply(subset(subset(lfreq_sp_y,lfreq_sp_y$HAULJOIN %in% hls)$FREQUENCY,subset(lfreq_sp_y,lfreq_sp_y$HAULJOIN %in% hls)$SEX!=3),
                              subset(subset(lfreq_sp_y,lfreq_sp_y$HAULJOIN %in% hls)$HAULJOIN,subset(lfreq_sp_y,lfreq_sp_y$HAULJOIN %in% hls)$SEX!=3),sum)-max_sx[ss])
      
      # Place haul % and SS redux in results data
      Other_res<-rbind(Other_res,c(region,species[sp],max_sx[ss],yrs[y],hl_prop,sampl_redux))
      
      #########################
      # Loop thru iterations
      for(i in 1:I){
        
        # keep track of where simulation is at
        print(paste0("species = ",species[sp],"yr = ",yrs[y]," ss = ",max_sx[ss]," iter = ",i))
        
        # Setup dummy data frame
        lfreq_sp_y_hl_smpld<-NULL
        
        #########################
        # Loop thru hauls w/ >100 samples
        for (h in 1:length(hls)){
          
          # Subset length freq to haul level
          lfreq_sp_y_hl<-subset(lfreq_sp_y,lfreq_sp_y$HAULJOIN==hls[h])
          
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
        
        # Keep track of resampled samples for example hauls
        examp_hls_iter<-subset(lfreq_sp_y_hl_smpld,lfreq_sp_y_hl_smpld$HAULJOIN %in% hls_minmax & lfreq_sp_y_hl_smpld$SEX!=3)
        examp_res_iter<-cbind(
          examp_hls_iter$REGION,
          examp_hls_iter$SPECIES_CODE,
          rep(max_sx[ss],times=length(examp_hls_iter$HAULJOIN)),
          examp_hls_iter$YEAR,
          rep(i,times=length(examp_hls_iter$HAULJOIN)),   
          examp_hls_iter$HAULJOIN,
          examp_hls_iter$SEX,
          examp_hls_iter$LENGTH,
          examp_hls_iter$FREQUENCY)
        colnames(examp_res_iter)<-Examp_hls_nms
        Examp_hls<-rbind(Examp_hls,examp_res_iter)   
        
        # Replace actual hauls with resampled hauls
        lfreq_sp_y_smpld<-lfreq_sp_y[-which(lfreq_sp_y$HAULJOIN %in% hls),]
        lfreq_sp_y_smpld<-rbind(lfreq_sp_y_smpld,lfreq_sp_y_hl_smpld)
        
        # Estimate iteration's pop'n #'s at size
        # Management area scale
        sp_SzPop_iter<-est_SzPop(stratum,lengths,strata_reg,cpue_sp_y,lfreq_sp_y_smpld)
        iter_res<-cbind(
          rep(region,times=length(sp_SzPop_iter$LENGTH)),
          rep(species[sp],times=length(sp_SzPop_iter$LENGTH)),
          rep(max_sx[ss],times=length(sp_SzPop_iter$LENGTH)),
          rep(yrs[y],times=length(sp_SzPop_iter$LENGTH)),
          rep(i,times=length(sp_SzPop_iter$LENGTH)),
          sp_SzPop_iter,
          sp_SzPop_iter$MALES+sp_SzPop_iter$FEMALES+sp_SzPop_iter$UNSEXED)
        colnames(iter_res)<-Sizepop_ma_nms
        Sizepop_ma<-rbind(Sizepop_ma,iter_res)
        
        # Strata scale
        sp_SzPop_strata_iter<-est_SzPop_strata(stratum,lengths,strata_reg,cpue_sp_y,lfreq_sp_y_smpld)
        iter_res_strata<-cbind(
          rep(region,times=length(sp_SzPop_yr_strata$STRATUM)),
          rep(species[sp],times=length(sp_SzPop_yr_strata$STRATUM)),
          rep(max_sx[ss],times=length(sp_SzPop_yr_strata$STRATUM)),
          rep(yrs[y],times=length(sp_SzPop_yr_strata$STRATUM)),
          sp_SzPop_yr_strata$STRATUM,
          sp_SzPop_yr_strata$SUMMARY_AREA,
          rep(i,times=length(sp_SzPop_yr_strata$STRATUM)),
          sp_SzPop_yr_strata$LENGTH,
          sp_SzPop_yr_strata$MALES,
          sp_SzPop_yr_strata$FEMALES,
          sp_SzPop_yr_strata$UNSEXED,
          sp_SzPop_yr_strata$MALES+sp_SzPop_yr_strata$FEMALES+sp_SzPop_yr_strata$UNSEXED)
        colnames(iter_res_strata)<-Sizepop_strata_nms
        Sizepop_strata<-rbind(Sizepop_strata,iter_res_strata)
        
        # Compute and store management area scale statistics (ESS)
        # Females
        full_comp_F<-sp_SzPop_yr$FEMALES/sum(sp_SzPop_yr$FEMALES)
        iter_comp_F<-sp_SzPop_iter$FEMALES/sum(sp_SzPop_iter$FEMALES)
        ESS_F_iter<-sum(iter_comp_F*(1-iter_comp_F))/sum((iter_comp_F-full_comp_F)^2)
        # Males
        full_comp_M<-sp_SzPop_yr$MALES/sum(sp_SzPop_yr$MALES)
        iter_comp_M<-sp_SzPop_iter$MALES/sum(sp_SzPop_iter$MALES)
        ESS_M_iter<-sum(iter_comp_M*(1-iter_comp_M))/sum((iter_comp_M-full_comp_M)^2)
        # Total
        full_comp_TOT<-(sp_SzPop_yr$FEMALES+sp_SzPop_yr$MALES+sp_SzPop_yr$UNSEXED)/sum(sp_SzPop_yr$FEMALES+sp_SzPop_yr$MALES+sp_SzPop_yr$UNSEXED)
        iter_comp_TOT<-(sp_SzPop_iter$FEMALES+sp_SzPop_iter$MALES+sp_SzPop_iter$UNSEXED)/sum(sp_SzPop_iter$FEMALES+sp_SzPop_iter$MALES+sp_SzPop_iter$UNSEXED)
        ESS_TOT_iter<-sum(iter_comp_TOT*(1-iter_comp_TOT))/sum((iter_comp_TOT-full_comp_TOT)^2)
        # Store stats
        stats_iter<-cbind(
          region,
          species[sp],
          max_sx[ss],
          yrs[y],
          i,
          ESS_M_iter,
          ESS_F_iter,
          ESS_TOT_iter)
        colnames(stats_iter)<-Sim_stats_nms
        Sim_stats<-rbind(Sim_stats,stats_iter)
        
        # Compute and store strata scale statistics (ESS)
        for(st in 1:length(stratum)){
          # Females
          full_comp_F_st<-subset(sp_SzPop_yr_strata,sp_SzPop_yr_strata$STRATUM==stratum[st])$FEMALES/sum(subset(sp_SzPop_yr_strata,sp_SzPop_yr_strata$STRATUM==stratum[st])$FEMALES)
          iter_comp_F_st<-subset(sp_SzPop_strata_iter,sp_SzPop_strata_iter$STRATUM==stratum[st])$FEMALES/sum(subset(sp_SzPop_strata_iter,sp_SzPop_strata_iter$STRATUM==stratum[st])$FEMALES)
          ESS_F_iter_st<-sum(iter_comp_F_st*(1-iter_comp_F_st))/sum((iter_comp_F_st-full_comp_F_st)^2)
          # Males
          full_comp_M_st<-subset(sp_SzPop_yr_strata,sp_SzPop_yr_strata$STRATUM==stratum[st])$MALES/sum(subset(sp_SzPop_yr_strata,sp_SzPop_yr_strata$STRATUM==stratum[st])$MALES)
          iter_comp_M_st<-subset(sp_SzPop_strata_iter,sp_SzPop_strata_iter$STRATUM==stratum[st])$MALES/sum(subset(sp_SzPop_strata_iter,sp_SzPop_strata_iter$STRATUM==stratum[st])$MALES)
          ESS_M_iter_st<-sum(iter_comp_M_st*(1-iter_comp_M_st))/sum((iter_comp_M_st-full_comp_M_st)^2)
          # Total
          full_comp_TOT_st<-(subset(sp_SzPop_yr_strata,sp_SzPop_yr_strata$STRATUM==stratum[st])$FEMALES+subset(sp_SzPop_yr_strata,sp_SzPop_yr_strata$STRATUM==stratum[st])$MALES+subset(sp_SzPop_yr_strata,sp_SzPop_yr_strata$STRATUM==stratum[st])$UNSEXED)/sum(subset(sp_SzPop_yr_strata,sp_SzPop_yr_strata$STRATUM==stratum[st])$FEMALES+subset(sp_SzPop_yr_strata,sp_SzPop_yr_strata$STRATUM==stratum[st])$MALES+subset(sp_SzPop_yr_strata,sp_SzPop_yr_strata$STRATUM==stratum[st])$UNSEXED)
          iter_comp_TOT_st<-(subset(sp_SzPop_strata_iter,sp_SzPop_strata_iter$STRATUM==stratum[st])$FEMALES+subset(sp_SzPop_strata_iter,sp_SzPop_strata_iter$STRATUM==stratum[st])$MALES+subset(sp_SzPop_strata_iter,sp_SzPop_strata_iter$STRATUM==stratum[st])$UNSEXED)/sum(subset(sp_SzPop_strata_iter,sp_SzPop_strata_iter$STRATUM==stratum[st])$FEMALES+subset(sp_SzPop_strata_iter,sp_SzPop_strata_iter$STRATUM==stratum[st])$MALES+subset(sp_SzPop_strata_iter,sp_SzPop_strata_iter$STRATUM==stratum[st])$UNSEXED)
          ESS_TOT_iter_st<-sum(iter_comp_TOT_st*(1-iter_comp_TOT_st))/sum((iter_comp_TOT_st-full_comp_TOT_st)^2)
          # Store stats
          stats_iter_st<-cbind(
            region,
            species[sp],
            max_sx[ss],
            yrs[y],
            stratum[st],
            i,
            ESS_M_iter_st,
            ESS_F_iter_st,
            ESS_TOT_iter_st)
          colnames(stats_iter_st)<-Sim_stats_st_nms
          Sim_stats_st<-rbind(Sim_stats_st,stats_iter_st)
        }
        
        #########################
        # End iteration loop
      }
      
      #########################
      # End max sex sample size loop
    }
    
    #########################
    # End year loop
  }
  
  # Write results matrices
  write.csv(Examp_hls,paste0(pathR,"/",region,"/Examp_hls.csv"))
  write.csv(Sizepop_ma,paste0(pathR,"/",region,"/Sizepop_ma.csv"))
  write.csv(Sizepop_strata,paste0(pathR,"/",region,"/Sizepop_strata.csv"))
  write.csv(Sim_stats,paste0(pathR,"/",region,"/Sim_stats.csv"))
  write.csv(Sim_stats_st,paste0(pathR,"/",region,"/Sim_stats_st.csv"))
  write.csv(Other_res,paste0(pathR,"/",region,"/Other_res.csv"))
  
  # End time
  end_time<-Sys.time()
  
  # Compute species-specific time to run sim
  run_time<-end_time-st_time
  run_time
  
  #########################
  # End species loop 
}









