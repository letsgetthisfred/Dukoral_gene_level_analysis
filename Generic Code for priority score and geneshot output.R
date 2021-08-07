#Generic Code

#Bring in data from geneshot/AMGMA output or existing csv files

#If starting from geneshot, first setup r environment per here: https://github.com/Golob-Minot/geneshot/wiki/Reading-HDF5-with-R
        #Can skip if bringing in data back in from csv's after already reading hdf5 during prior session (code for this further down)  
        #relevant code from geneshot wiki pulled out directly below as of 4/27/21
        
                # Install the reticulate package
                install.packages("devtools")
                devtools::install_github("rstudio/reticulate")
                library("reticulate")
                # Use the path to Python on your system
                use_python("/usr/path/to/python")
                use_python("/Users/fjh2/miniconda3/pkgs/rpy2-3.4.2-py38r40ha1b04c9_1/lib/python3.8")
                # Create and use a new virtual enviornment named "r-reticulate"
                use_virtualenv("r-reticulate") 
                # Install packages required to read hdf5 data
                py_install("pandas")
                py_install("scipy")
                py_install("pytables")
                py_install("h5py")
                py_install("rpy2")
                # Import those packages into the environment
                h5py <- import("h5py")
                rpy2 <- import("rpy2")
                rpy2$robjects
                rpy2_ro <- import("rpy2.robjects")
                rpy2_pandas2ri <- import("rpy2.robjects.pandas2ri")
                rpy2_pandas2ri$py2rpy
                pd <- import("pandas", convert = FALSE)
                np <- import("numpy", convert = FALSE)


        #Now that environment has been set up, bring data into R directly from geneshot/ AMGMA output (hdf5 files)
                #cag associations
                key_CAGstats <- pd$read_hdf("/Volumes/Fred_backup/Third_Year_Data/after_alignment_mishap/2021-03-12-Weil-IRSV.results.hdf5", key = "/stats/cag/corncob", mode ='r+') 
                CAGstats_df <- py_to_r(key_CAGstats)
                        CAGstats_df$cag<-as.character(CAGstats_df$cag) #want CAG to be recognized as not a number
                
                #pull in CAG basic info (eg size) (new file)
                key_CAG <- pd$read_hdf("/Volumes/Fred_backup/Third_Year_Data/after_alignment_mishap/2021-03-12-Weil-IRSV.results.hdf5", key = "/annot/cag/all", mode ='r+') 
                CAG_info_df <- py_to_r(key_CAG)
                
                #bring in CAG-genome overlap info (new file)
                key_containment2 <- pd$read_hdf("/Volumes/Fred_backup/Third_Year_Data/2021-03-17-Weil-IRSV-AMGMA.hdf5", key = "/genomes/cags/containment", mode ='r+')
                containment_df2 <- py_to_r(key_containment2)
                
                #bring in genome associations (calculated using corncob) (new file)
                key_associations2 <- pd$read_hdf("/Volumes/Fred_backup/Third_Year_Data/2021-03-17-Weil-IRSV-AMGMA.hdf5", key = "/stats/genome/corncob", mode ='r+')
                associations_df2 <- py_to_r(key_associations2)
                
                #create manifest relating genome ID and name (new AMGMA file)
                key_manifest2 <- pd$read_hdf("/Volumes/Fred_backup/Third_Year_Data/2021-03-17-Weil-IRSV-AMGMA.hdf5", key = "/genomes/manifest", mode ='r+') #may need ID from manifest to translate from containment to association tables
                manifest_df2 <- py_to_r(key_manifest2)
                
                #Interrelate values from genome associations table with IDs in manifest
                #genomes associations table (key = "/stats/genome/corncob") does not have genome ID as a column, so this adds it in
                        associations_df2['pandas_index']<-NA 
                        for (i in 1:(nrow(associations_df2))){
                                associations_df2[i,"pandas_index"]<-as.character(attributes(associations_df2)$pandas.index[i-1]) #subtract 1 i bc pandas index starts at "row 0" instead of row 1 
                        }
                        
                        manifest_df2['pandas_index']<-NA
                        for (i in 0:(nrow(manifest_df2)-1)){ #manifest table starts at row zero
                                manifest_df2[i,"pandas_index"]<-as.character(attributes(manifest_df2)$pandas.index[i-1]) #subtract 1 i bc pandas index starts at "row 0" instead of row 1 
                        }
                        #index has decimal in associations table but not in manifest; fix to match
                        for (i in 0:nrow(manifest_df2)){ 
                                manifest_df2[i,"pandas_index"]<-paste(manifest_df2[i,"pandas_index"],".0",sep="")
                        }
                        
                        genome_associations<-merge(associations_df2, manifest_df2,  by.x = "pandas_index", sort=FALSE)
        
                #Export from R to have for reference:
                #update path for local device and names as you see fit
                        write.csv(containment_df2.2,"PATH/genome_alignment_all_SECOND_HALF.csv", row.names = FALSE)
                        write.csv(genome_associations,"PATH/genome_associations_all_w_names.csv", row.names = FALSE) 
                        write.csv(manifest_df2,"PATH/genomes_manifest.csv", row.names = FALSE)
                        write.csv(CAGstats_df,"PATH/CAG_associations_all.csv", row.names = FALSE)
                        write.csv(CAG_info_df,"PATH/CAG_info_all.csv", row.names = FALSE)
                        #containment_df needs to be exported as 2 files bc it had too many rows for one csv file
                                containment_df2.1 <-containment_df2[1:766029,]
                                containment_df2.2 <-containment_df2[766030:1532058,]
                                write.csv(containment_df2.1,"PATH/genome_alignment_all_FIRST_HALF.csv", row.names = FALSE)
                                write.csv(containment_df2.2,"PATH/genome_alignment_all_SECOND_HALF.csv", row.names = FALSE)
                                
#If uploading raw output from existing csv's instead of reading in data from hdf5 files
        genome_associations <-read.csv(file="PATH/genome_associations_all_w_names.csv", header=TRUE) 
        manifest_df2 <- read.csv(file="PATHgenomes_manifest.csv", header=TRUE) 
        CAGstats_df <-read.csv(file="PATH/CAG_associations_all.csv", header=TRUE) 
        CAGstats_df$cag<-as.character(CAGstats_df$cag) #want CAGs to be read as characters
        CAG_info_df <- read.csv(file="PATH/CAG_info_all.csv", row.names=1, header=TRUE)
        #recall, had to export containment_df as 2 separate csv files due to size, so import here and combine
                containment_df2.11<- read.csv(file="PATH/genome_alignment_all_FIRST_HALF.csv", header=TRUE) 
                containment_df2.22<- read.csv(file="PATH/genome_alignment_all_SECOND_HALF.csv", header=TRUE) 
                containment_df2<- rbind(containment_df2.11,containment_df2.22)
                

#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
#CREATING PRIORITY SCORE LISTS 
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
                
#PARAMETERS TO CHANGE IN BELOW CODE AS YOU SEE FIT
        #CAG q value threshold- line 120 (eg 0.1)
        #CAG size threshold- line 122 (eg 2)
        #experimental parameter of interest, lines 121 and 174 (eg "OSP_IgA_Responder")
                
#create data frame to receive info on top CAGs and their associations 
        topCAGs<-data.frame(est_coeff=as.numeric(), #create empty receptacle for top CAGs
                            wald=as.numeric(),
                            weighting=as.numeric(), #create column for their weighting term
                            q_value=as.numeric(),
                            size=as.numeric()) 
                
#pull in top CAGs and their associations; through "if" statements, set desired significance threshold, exp parameter, and CAG min size
        for (i in 1:nrow(CAGstats_df)){ #initialize loop through rows of CAG association table
                if (CAGstats_df[i, "q_value"] <= NUMBER){ #set min significance threshold as a value outside of quotes (eg 0.31622776601 which approximates -log(q)>=0.5, or 0.03162277 which approximates -log(q)≥1.5)
                        if (CAGstats_df[i, "parameter"] == "PARAMETER") {      #set experimental parameter of interest within quotes
                                if (CAG_info_df[CAGstats_df[i,"cag"],"size"]>=NUMBER) { #set size threshold as a number outside of quotes; note this info is stored separately from association info (in CAG_info_df instead of CAGstats_df)
                                        rownames(topCAGs_CT_IgA[nrow(topCAGs) + 1,]) <- CAGstats_df[i, "cag"] #pull out CAG as row name
                                        topCAGs_CT_IgA[CAGstats_df[i, "cag"],"est_coeff"] <- CAGstats_df[i,"estimate"]
                                        topCAGs_CT_IgA[CAGstats_df[i, "cag"],"wald"] <- CAGstats_df[i,"wald"]
                                        topCAGs_CT_IgA[CAGstats_df[i, "cag"],"q_value"] <- CAGstats_df[i,"q_value"]
                                        topCAGs_CT_IgA[CAGstats_df[i, "cag"],"size"] <- CAG_info_df[CAGstats_df[i, "cag"],"size"]
                                }}}}
                
        #add weighting term to CAG associations data frame 
                for (i in 1:nrow(topCAGs)){
                        topCAGs[i,"weighting"] <- sqrt(topCAGs[i, "est_coeff"]^2+ topCAGs[i, "wald"]^2)
                }
        #if coefficient negative, make weighting term negative
                for (i in 1:nrow(topCAGs)){
                        if (topCAGs[i,"est_coeff"]<0){
                                topCAGs[i,"weighting"]<-(-1*topCAGs[i,"weighting"])
                        }
                }
                
#Create data frame with top CAGs as columns to prepare for creation of a CAG gene vs genome alignment table (this method was a bit roundabout)
        top_CAGs<-rownames(topCAGs) #create list of top CAGs 
        strain_map<- data.frame(rbind(top_CAGs)) #assign top CAGs as row 1
        colnames(strain_map)<-top_CAGs #convert top CAGs from row 1 to column headings
        rownames(strain_map)<-"No alignments found" #make row 1, which formally had top CAGs, a counter for CAG lack of alignment
        strain_map[,1:ncol(strain_map)] <- sapply(strain_map[,ncol(strain_map)], as.numeric) #make that first row numeric class instead of character
        for (i in 1:length(top_CAGs)) {strain_map["No alignments found",i]<-as.numeric(0)} #set values of that first row at 0 across the board 
                
        #loop through genome gene alignment table by CAG, pulling out genomes
        for (t in 1:length(top_CAGs))  {  #go through top CAG list, can also use 2:length(OSP_IgA_topCAGs)
                for (row in 1:nrow(containment_df2)) { #for each top CAG, look for it in containment file
                        if (containment_df2[row,"CAG"] == top_CAGs[t]) {    #if CAG present in containment file...maybe need to do containment_df[,2]
                                if (containment_df2[row,"genome"] %in% rownames(strain_map)) {  #check if genome already in col 1 of strains table
                                        strain_map[containment_df2[row, "genome"],t]<-containment_df2[row,"n_genes"]  #add number of genes aligned to appropriate cell
                                }
                                else {
                                        rownames(strain_map[nrow(strain_map) + 1,]) <- containment_df2[row, "genome"]  #add row with strain name in first column
                                        strain_map[containment_df2[row, "genome"],t] <- containment_df2[row,"n_genes"]  #add number of genes aligned to appropriate cell
                                }}
                        else {
                                (strain_map["No alignments found",t]<-strain_map["No alignments found",t]+1) #this should create a count of genomes to which none of the genes aligned?  
                        }}}
                
#for genome-level associations
        #create receptacle data frame for relevant genome association info (only those genomes with genes aligning from top CAGs)
                topGenomes<-as.data.frame(matrix(0, ncol=4, nrow=length(rownames(strain_map)))) 
                rownames(topGenomes)<-rownames(strain_map) #pull in all genomes with top CAG genes aligning
                colnames(topGenomes)<-c("est_coeff","wald","weighting","q_value")
        
        #loop through association values for all genomes, pulling into topGenomes dataframe association values for those with genes aligning from top CAGs
                for (i in 1:nrow(topGenomes)){ #initialize loop through rows of genomes with top CAG genes
                        for (t in 1:nrow(genome_associations)){ #loop also through genome association table rows
                                if (rownames(topGenomes[i,])==genome_associations[t,"id"]){#when you've found genome of interest
                                        if (genome_associations[t,"parameter"]=="PARAMETER"){ #narrow to only experimental parameter of interest, expressed in quotes
                                                topGenomes_CT_IgA[rownames(topGenomes[i,]),"est_coeff"]<-genome_associations[t,"estimate"]
                                                topGenomes_CT_IgA[rownames(topGenomes[i,]),"wald"]<-genome_associations[t,"wald"]
                                                topGenomes_CT_IgA[rownames(topGenomes[i,]),"q_value"]<-genome_associations[t,"q_value"]
                                        } 
                                }
                        }
                }
                

        #add weighting term to genome associations data frame 
                for (i in 1:nrow(topGenomes)){
                        topGenomes[i,"weighting"] <- sqrt(topGenomes[i, "est_coeff"]^2+ topGenomes[i, "wald"]^2)
                }
        #if coefficient negative, make weighting term negative... Do I want to do this? Unlike w/ CAG term, this might change whole score's sign.... Decided yes, can manipulate signs in excel afterwards
                for (i in 1:nrow(topGenomes_)){
                        if (topGenomes[i,"est_coeff"]<0){
                                topGenomes[i,"weighting"]<-(-1*topGenomes[i,"weighting"])
                        }
                }
                
#calculate priority scores 
        #create data frame for priority score list
                PS<-as.data.frame(matrix(0, ncol=5, nrow=length(rownames(strain_map)))) 
                rownames(PS)<-rownames(strain_map)
                colnames(PS)<-c("priority score","genomes","CAG term numerator", "CAG term denominator","strain term")
                blacklist<- c() #if i want to exclude the influence of any CAGs...
        
        #calculate priority score for top genomes, adding that info to the PS data frame
                for (i in 1:nrow(strain_map)) {  #iterate through relevant genomes (rows of strainmap)
                        CAG_num <- 0 #these variables will be reset to zero for each genome at the start of each loop
                        CAG_denom <- 0
                        strain<- topGenomes[rownames(strain_map[i,]),"weighting"] #pull in strain weighting term specific to genome
                        for (t in 1:ncol(strain_map)){ #iterate through top CAGs (columns of strainmap); where genome and CAG have gene overlapp, add number of genes to denominator of CAG term and weighted number to numerator of CAG term
                                if (is.na(strain_map[i,t])==FALSE && as.numeric(strain_map[i,t])>0 && colnames(strain_map[t]) %in% blacklist == FALSE) { #check for overlap in genes between CAG and genome; excluding NA values and blacklisted CAGs
                                        CAG_denom <- CAG_denom + as.numeric(strain_map[i,t]) #where overlap is found, add number of genes aligning from CAG to denominator of CAG term
                                        CAG_num <- CAG_num + as.numeric(strain_map[i,t])*topCAGs_CT_IgA[colnames(strain_map)[t],"weighting"] #where overlap is found, add **weighted** number of genes alignign from CAG to numerator
                                }
                        }
                        PS_CT_IgA[rownames(strain_map)[i],"priority score"]<-(CAG_num/CAG_denom)*strain #now we have looped through all CAGs for this genome, assign the genome a priority score
                        PS_CT_IgA[rownames(strain_map)[i],"CAG term numerator"]<-CAG_num #note genome's CAG term numerator
                        PS_CT_IgA[rownames(strain_map)[i],"CAG term denominator"]<-CAG_denom #note genome's CAG term denominator
                        PS_CT_IgA[rownames(strain_map)[i],"strain term"]<-strain #note genome's strain term 
                        CAG_num <- 0  #now reset value prior to the loop moving onto the next genome next genome
                        CAG_denom <- 0
                        strain <- 0
                }
                
        #assign genome its name in addition to its AMGMA id 
                for (i in 1:nrow(strain_map)){    #go through mapped strain IDs
                        for (t in 1:nrow(manifest_df2)){  #go through manifest containing names for genomes
                                if (rownames(PS)[i]==manifest_df2[t,"id"]){ #identify when mapped strain IDs are found in manifest
                                        PS[i,"genomes"]<-manifest_df2[t,"name"] #add genome name alongside it's strain ID in priority score output
                                }
                        }
                }
                
        #sort descending, so top PS scores are at top (can easily be done in output csv/xls)
                PS_sorted <- PS[order(-PS$`priority score`), ]
                
        #Export ranked priority scores as well as intermediary files for documentation purposes
                #update path to desired output location and name as appropriate (eg to indicate parameters used)
                write.csv(PS_CT_IgA_sorted,"PATH/PS.csv", row.names = TRUE)
                write.csv(strain_mapCTIgA,"PATH/genome_alignment.csv", row.names = TRUE)
                write.csv(topCAGs_CT_IgA,"PATH/top_CAGs.csv", row.names = TRUE)
                write.csv(topGenomes_CT_IgA,"PATH/top_Genomes.csv", row.names = TRUE)
                
        #IMPORTANT NOTE ON FINAL PRIORITY SCORES:
                #After export of Priority Score CSV, the sign of priority score was altered so that it matched that of CAG term numerator 
                        #The R code above multiplies strain term by CAG term, leading to situations where Priority Score was positive if both terms were negative; that informative directionality was therefore obscured
                #This was done by creating a "final priority score" column and using the below code
                        # =IF(OR(AND(D2<0,B2>0),AND(D2>0,B2<0)),B2*-1,B2)
                        # where column B had initial priority score output and column D had CAG numerator
                        #this formula was pulled down the column and values generated were used as final, sign-appropriate priority scores
                
       