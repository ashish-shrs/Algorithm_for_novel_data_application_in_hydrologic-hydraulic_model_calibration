# Authors: Ashish Shrestha*, Margaret Garcia
# *Contact Email: ashres15@asu.edu (Please feel free to contact me with any questions, if you want to apply these algorithms for your study/application)

# Part 1 of 2
# Note: This algorithm extract time series data of nodes's flooding for one dimensional SWMM model and overland flood water depth for one-& two-dimensional coupled version of SWMM model, and establish statistical relationship between two models.

# Installing packages/loading libraries (Tools Dependencies)
x = c('swmmr', 'ggplot2', 'xts', 'zoo', 'lubridate', 'dplyr', 
      'grid', 'gridExtra', 'reshape2', 'scales', 'RColorBrewer','tidyr')
#lapply(x, install.packages, character.only = TRUE) # Install all packages
#lapply(x, library, character.only = TRUE) #load all packages

# Assigning filepath for SWMM 1D result file. Its content will be extracted for each nodes in loop below
base_file = "~/filedirectory/swmm1D.out"

# Importing .csv file containing corresponding 1D node, 1D and 2D coupled nodes IDs, generated using spatial analysis tools of ArcGIS
allnodes = read.csv("~/filedirectory/all1D2Dnodes.csv", header = TRUE)

# Loop to extract surface flooding for each coupled nodes
mat_ts = as.data.frame(matrix(0 ,nrow = 714, ncol = 1)) # Creating empty matrix for surface flooding #nrow is variable to time series report length
for (i in 1:length(allnodes$Node1d)){ # Then using read_out function of swmmr R package this loop reads flooding for every node ID
  swmm_run <- read_out(file=base_file, iType=1, object_name = as.character(allnodes[i,3]), vIndex=5)  
  names(swmm_run) <- "J" # Since the col name is assigned node ID by default, here it is manually changed to "J" to extract easily
  swmm_run_xts <- as.xts(swmm_run$J$surface_flooding) # This is extracted time series
  swmm_run_ts <- fortify.zoo(swmm_run_xts)  # and time series is converted to data base
  mat_ts <- cbind(mat_ts, swmm_run_ts$swmm_run_xts) # That time series data base is combined to empty matrix above
  names(mat_ts)[i+1]<-as.character(allnodes$Node1d[i]) # Each new column is given corresponding node ID as col name
}

# Assigning filepath for SWMM 1D-2D result file. Its content will be extracted for each nodes in loop below
base_file_2d = "~/filedirectory/swmm2D.out"

# Loop to extract 2D depth time series for each coupled nodes
mat_ts_2d = as.data.frame(matrix(0 ,nrow = 714, ncol = 1)) # Creating empty matrix for 2D depth #nrow is variable to time series report length
for (i in 1:length(allnodes$Node2d)){ # Then using read_out function of swmmr package this loop reads depth for every node ID
  swmm_run_2d <- read_out(file = base_file_2d, iType=1, object_name = as.character(allnodes[i,4]), vIndex=0)
  names(swmm_run_2d) <- "JJ" # Since the col name is assigned node ID by default, here it is manually changed to "JJ" to extract easily
  swmm_run_xts_2d <- as.xts(swmm_run_2d$JJ$water_depth) # This is extracted time series
  swmm_run_ts_2d <- fortify.zoo(swmm_run_xts_2d)  # and time series is converted to data base
  mat_ts_2d <- cbind(mat_ts_2d, swmm_run_ts_2d$swmm_run_xts_2d) # That time series data base is combined to empty matrix mat_ts_2d
  names(mat_ts_2d)[i+1]<-as.character(allnodes$Node2d[i]) # Each new column is given corresponding node id as col name, similar as above
}

# Area under the curve for flooding in SWMM 1D
for (i in 1:length(mat_ts$V1)){ # Changing first column into time step 
  mat_ts$V1[i] = i
}
mat_integration_flood <- as.data.frame(matrix(0, nrow=713, ncol =1)) # Creates empty matrix 
for (j in 2:ncol(mat_ts)){ # This loop creates time integration of flooding for each 1D node, one at a time
  col_num = j # Consider column (or node id) from column 2, as column 1 is empty
  mat1 = matrix(0 ,nrow = nrow(mat_ts) -1, ncol = 1) # Empty matrix for width which is uniform for each small rectangles
  mat2 = matrix(0 ,nrow = nrow(mat_ts) -1, ncol = 1) # Empty matrix for height of each rectangle, as average or mid point of two data points
  mat3 = matrix(0 ,nrow = nrow(mat_ts) -1, ncol = 1) # Empty matrix for area, width x height x 60 sec
  mat4 = matrix(0 ,nrow = nrow(mat_ts) -1, ncol = 1) # Empty matrix for area under the curve, cumulative area to the left of the floodvstime curve
  for (i in 1:length(mat_ts[,col_num])-1){
    mat1[i] <- mat_ts[i+1,1] - mat_ts[i,1]
    mat2[i] <- 0.5*(mat_ts[i+1,col_num]+mat_ts[i,col_num])
    mat3[i] <- mat1[i] * mat2[i] * 60
    mat4[i] <- ifelse (i == 1, mat3[1], mat4[i-1]+mat3[i])
  }
  mat_integration_flood <- cbind(mat_integration_flood, mat4) # This step combines area under the curve to empty matrix mat_integration_flood
  names(mat_integration_flood)[j] <- as.character(allnodes$Node1d[j-1]) # And assign node ID as corresponding col names
}

# Matrix to find out index of starting point, end point and maximum point (SEM) of hydrographs in three rows
mat_sem <- as.data.frame(matrix(0, nrow = 3, ncol =1)) 
for (m in 2:ncol(mat_ts)){
  mat_st_end <- matrix(0, nrow = 2, ncol =1)
  if (sum(mat_ts[,m])>0){ # This condition only applies to flooding nodes. 
    tempvect <- which(mat_ts[,m]>0) 
    mat_st_end[1] <- min(tempvect)
    mat_st_end[2] <- max(tempvect)
    mat_st_end[3] <- ifelse(length(which(mat_ts[,m] == max(mat_ts[,m])))>1,which(mat_ts[,m] == max(mat_ts[,m]))[[(length(which(mat_ts[,m] == max(mat_ts[,m]))))]],which(mat_ts[,m] == max(mat_ts[,m]))[[1]])
  }else{ # For non flooding nodes which has zeros, SEM is assigned as zeros
    mat_st_end[1] <- 0
    mat_st_end[2] <- 0
    mat_st_end[3] <- 0
  }
  mat_sem <- cbind(mat_sem, mat_st_end)
  names(mat_sem)[m] <- as.character(colnames(mat_ts)[m])
}

# Linear regression
# x and y are considered from time index from starting point till peak
mat_lm_coeff <- as.data.frame(matrix(0, nrow=2, ncol =1)) # Empty matrix for lm coefficients 
mat_lm_rsquared <- as.data.frame(matrix(0, nrow=1, ncol = ncol(mat_ts)-1)) # Empty matrix for Rsquare value for each model
for (k in 2:ncol(mat_ts)){
  if (sum(mat_sem[,k]) > 0 & sd(mat_sem[,k])> 0){ 
    y_all <- mat_ts_2d[,k] # These are temporary vector, updated each time step
    y <- y_all[mat_sem[1,k]:mat_sem[3,k]]
    x_all <- mat_integration_flood[,k]
    x <- x_all[mat_sem[1,k]:mat_sem[3,k]]
    lm.fit = lm(y ~ x)
    lm.fit.coeff <- as.data.frame(coef(lm.fit)) # Extract coefficients for each y~x fit
    mat_lm_coeff <- cbind(mat_lm_coeff, lm.fit.coeff) # Combines new column into empty matrix 
    names(mat_lm_coeff)[k] <-as.character(colnames(mat_ts[k])) # Assign corresponding node id
    mat_lm_rsquared[1,k-1]<- summary(lm.fit)$r.squared # Extract R square value of each y~x fit
    names(mat_lm_rsquared)[k-1] <- as.character(colnames(mat_ts[k])) # Assign R square to each node column
  }else{
    lm.fit.coeff <- as.data.frame(c(0,0)) # For non-flooding nodes zeros are assigned as coefficient and Rsquared values
    mat_lm_coeff <- cbind(mat_lm_coeff, lm.fit.coeff) # Combines new column into empty matrix 
    names(mat_lm_coeff)[k] <-as.character(colnames(mat_ts[k])) # Assign corresponding node id
    mat_lm_rsquared[1,k-1]<- 0 # Extract R square value of each y~x fit
    names(mat_lm_rsquared)[k-1] <- as.character(colnames(mat_ts[k])) # Assign R square to each node column
  }
}

# Final outputs of Part 1
Linear_Regression_Intercepts <- mat_lm_coeff
Linear_Regression_Intercepts <- Linear_Regression_Intercepts[,-1]
Linear_Regression_Result <- rbind(Linear_Regression_Intercepts, mat_lm_rsquared)
rownames(Linear_Regression_Result) <- c("Coeffcient a - RL", "Coefficient b - RL", "Rsquared - RL")
## End of part 1 of 2


# Part 2 of 2
# Note: This algorithm parameterize SWMM 1D version using "Genetic Algorithm", with single objective optimization in parallel computing nodes

# Installing packages/loading libraries (Tools Dependencies)
x = c('swmmr', 'ggplot2', 'knitr', 'dplyr', 'purrr', 'readr', 
      'Rcpp', 'tibble', 'utils', 'xts', 'zoo','kimisc', 'qpcR', 'reshape2','GA', 'parallel','doParallel')
#lapply(x, install.packages, character.only = TRUE) # Install all packages
#lapply(x, library, character.only = TRUE) #load all packages

input  <- "~/filedirectory/swmmfile.inp"
report <- "~/filedirectory/swmmfile.rpt"
output <- "~/filedirectory/swmmfile.out"
dcores <- detectCores()


# Run swmm using "swmmr" tool
swmm_files <- run_swmm(
  inp = input,
  rpt = report,
  out = tmp_out_file,
  #exec = "/packages/7x/swmm5/51013/swmm5"
)


# Read info files for conduits and subcatchments for clustering
conduit_info1 <- read.csv("~/filedirectory/conduit_info1.csv", header=TRUE) # conduit_info1.csv file is database for missing conduit types information
conduit_info2 <- read.csv("~/filedirectory/conduit_info2.csv", header = TRUE) # conduit_info2.csv file is database of conduit attributes of un-calibrated model
subcatchment_info <- read.csv("~/filedirectory/Subcatchment_info.csv", header = TRUE) # Subcatchment_info.csv file is database of subcatchment attributes of un-calibrated model

# Observed event time series
obs <- as.data.frame(read.csv("~/filedirectory/observed.csv", header = TRUE))

# Initial set of parameters values
x <- c(0.009, 0.009, 0.011, 0.011, 0.019, 0.019, 0.011, 0.011, 
       0.01, 0.01,
       0.06, 0.05, 0.1, 0.4,
       0.27, 0.24, 0.25, 0.19,
       0.05, 0.1,
       0.05, 0.1,
       8.27, 8.27, 7, 3,
       1, 1) # List below introduce parameters 

objectv <- function(x) {
  # Conduits' Roughness 
  r_a <- x[1]
  r_b <- x[2]
  r_c <- x[3]
  r_d <- x[4]
  r_e <- x[5]
  r_f <- x[6]
  r_g <- x[7]
  r_h <- x[8]
  # Mannings n for surface
  nperv <- x[9]
  nimperv <- x[10]
  # Sub-catchments' Conductivity 
  c_abc <- x[11]
  c_def <- x[12]
  c_ghi <- x[13]
  c_jkl <- x[14]
  # Sub-catchments' Initial Deficit
  i_abc <- x[15]
  i_def <- x[16]
  i_ghi <- x[17]
  i_jkl <- x[18]
  # Sub-catchments' Depression Storage for Pervious Surface
  dp_abdeghjk <- x[19]
  dp_cfil <- x[20]
  # Sub-catchments' Depression Storage for Impervious Surface
  dip_abdeghjk <- x[21]
  dip_cfil <- x[22]
  # Sub-catchments' Suction Head 
  s_abc <- x[23]
  s_def <- x[24]
  s_ghi <- x[25]
  s_jkl <- x[26]
  # Conduit Geometry for Unknown Material
  geom_na1 <- x[27]
  geom_na2 <- x[28]
  
  # Read input file
  input1 = read_inp(input)
  
  # Parameterizing algorithm - Conduits' Roughness 
  for (i in 1: length(input1$conduits$Name)){
    input1$conduits$Roughness[i] <- ifelse(conduit_info22$Category[which(conduit_info22$NAME == as.character(input1$conduits$Name[i]))] == "Cat A", r_a,
                                           ifelse(conduit_info22$Category[which(conduit_info22$NAME == as.character(input1$conduits$Name[i]))] == "Cat B",r_b,
                                                  ifelse(conduit_info22$Category[which(conduit_info22$NAME == as.character(input1$conduits$Name[i]))] == "Cat C",r_c,
                                                         ifelse(conduit_info22$Category[which(conduit_info22$NAME == as.character(input1$conduits$Name[i]))] == "Cat D",r_d,
                                                                ifelse(conduit_info22$Category[which(conduit_info22$NAME == as.character(input1$conduits$Name[i]))] == "Cat E",r_e,
                                                                       ifelse(conduit_info22$Category[which(conduit_info22$NAME == as.character(input1$conduits$Name[i]))] == "Cat F",r_f,
                                                                              ifelse(conduit_info22$Category[which(conduit_info22$NAME == as.character(input1$conduits$Name[i]))] == "Cat G",r_g,r_h
                                                                              )))))))
  }
  # Parameterizing algorithm - Mannings n for surface
  for (i in 1:length(input1$subareas$`N-Perv`)){
    input1$subareas$`N-Perv`[i]  <- nperv
    input1$subareas$`N-Imperv`[i] <- nimperv
  }
  # Parameterizing algorithm - Sub-catchments' Conductivity 
  for (i in 1:length(input1$infiltration$Subcatchment)){
    input1$infiltration$Ksat[i]<- ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 1", c_abc,
                                         ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 2", c_abc,
                                                ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 3", c_abc,
                                                       ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 4", c_def,
                                                              ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 5", c_def,
                                                                     ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 6", c_def,
                                                                            ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 7", c_ghi,
                                                                                   ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 8", c_ghi,
                                                                                          ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 9", c_ghi,
                                                                                                 ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 10", c_jkl,
                                                                                                        ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 11", c_jkl,c_jkl
                                                                                                        )))))))))))
  }
  
  # Parameterizing algorithm - Sub-catchments' Initial Deficit
  for (i in 1:length(input1$infiltration$Subcatchment)){
    input1$infiltration$IMD[i]<- ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 1", i_abc,
                                        ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 2", i_abc,
                                               ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 3", i_abc,
                                                      ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 4", i_def,
                                                             ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 5", i_def,
                                                                    ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 6", i_def,
                                                                           ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 7", i_ghi,
                                                                                  ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 8", i_ghi,
                                                                                         ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 9", i_ghi,
                                                                                                ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 10", i_jkl,
                                                                                                       ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 11", i_jkl,i_jkl
                                                                                                       )))))))))))
  }
  
  # Parameterizing algorithm - Sub-catchments' Depression Storage for Pervious Surface
  for (i in 1:length(input1$subareas$Subcatchment)){
    input1$subareas$`S-Perv`[i] <- ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$subareas$Subcatchment[i]))] == "Cat 1", dp_abdeghjk,
                                          ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$subareas$Subcatchment[i]))] == "Cat 2", dp_abdeghjk,
                                                 ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$subareas$Subcatchment[i]))] == "Cat 3", dp_cfil,
                                                        ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$subareas$Subcatchment[i]))] == "Cat 4", dp_abdeghjk,
                                                               ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$subareas$Subcatchment[i]))] == "Cat 5", dp_abdeghjk,
                                                                      ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$subareas$Subcatchment[i]))] == "Cat 6", dp_cfil,
                                                                             ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$subareas$Subcatchment[i]))] == "Cat 7", dp_abdeghjk,
                                                                                    ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$subareas$Subcatchment[i]))] == "Cat 8", dp_abdeghjk,
                                                                                           ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$subareas$Subcatchment[i]))] == "Cat 9", dp_cfil,
                                                                                                  ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$subareas$Subcatchment[i]))] == "Cat 10", dp_abdeghjk,
                                                                                                         ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$subareas$Subcatchment[i]))] == "Cat 11", dp_abdeghjk,dp_cfil
                                                                                                         )))))))))))
  }
  # Parameterizing algorithm - Sub-catchments' Depression Storage for Impervious Surface
  for (i in 1:length(input1$subareas$Subcatchment)){
    input1$subareas$`S-Imperv`[i] <- ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$subareas$Subcatchment[i]))] == "Cat 1", dip_abdeghjk,
                                            ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$subareas$Subcatchment[i]))] == "Cat 2", dip_abdeghjk,
                                                   ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$subareas$Subcatchment[i]))] == "Cat 3", dip_cfil,
                                                          ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$subareas$Subcatchment[i]))] == "Cat 4", dip_abdeghjk,
                                                                 ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$subareas$Subcatchment[i]))] == "Cat 5", dip_abdeghjk,
                                                                        ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$subareas$Subcatchment[i]))] == "Cat 6", dip_cfil,
                                                                               ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$subareas$Subcatchment[i]))] == "Cat 7", dip_abdeghjk,
                                                                                      ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$subareas$Subcatchment[i]))] == "Cat 8", dip_abdeghjk,
                                                                                             ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$subareas$Subcatchment[i]))] == "Cat 9", dip_cfil,
                                                                                                    ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$subareas$Subcatchment[i]))] == "Cat 10", dip_abdeghjk,
                                                                                                           ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$subareas$Subcatchment[i]))] == "Cat 11", dip_abdeghjk,dip_cfil
                                                                                                           )))))))))))
  }
  # Parameterizing algorithm - Sub-catchments' Suction Head
  for (i in 1:length(input1$infiltration$Suction)){
    input1$infiltration$Suction[i]<- ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 1", s_abc,
                                            ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 2", s_abc,
                                                   ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 3", s_abc,
                                                          ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 4", s_def,
                                                                 ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 5", s_def,
                                                                        ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 6", s_def,
                                                                               ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 7", s_ghi,
                                                                                      ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 8", s_ghi,
                                                                                             ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 9", s_ghi,
                                                                                                    ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 10", s_jkl,
                                                                                                           ifelse(subcatchment_info$Category[which(subcatchment_info$NAME == as.character(input1$infiltration$Subcatchment[i]))] == "Cat 11", s_jkl,s_jkl
                                                                                                           )))))))))))
    }
  # Parameterizing algorithm - Conduit Geometry for Unknown Material
  for (i in 1: length(input1$xsections$Link)){
    input1$xsections$Geom1[i] <- ifelse(conduit_info1$TAG[which(conduit_info1$NAME == as.character(input1$xsections$Link[i]))] == "ND", geom_na1, #ND = No Data / Unknown Material
                                        ifelse(conduit_info1$TAG[which(conduit_info1$NAME == as.character(input1$xsections$Link[i]))] == "NA",geom_na2, #NA, NA2, NA3 = Not Available / Unknown Material
                                               ifelse(conduit_info1$TAG[which(conduit_info1$NAME == as.character(input1$xsections$Link[i]))] == "NA2",geom_na2,
                                                      ifelse(conduit_info1$TAG[which(conduit_info1$NAME == as.character(input1$xsections$Link[i]))] == "NA3",geom_na2, input1$xsections$Geom1[i]
                                                      ))))
    
  }
  input_new <- tempfile()
  output_new <- tempfile()
  report_new <- tempfile()
  write_inp(input1, input_fil)
  swmm_files <- run_swmm( 
    inp = input_new,
    rpt = report_new,
    out = output_new,
    #exec = "/packages/7x/swmm5/51013/swmm5" 
  )
  swmmoutt <- read_out(out_fil, iType = 2, object_name = "conduit_name_X", vIndex = 1)
  swmmoutt_ts0 <- fortify.zoo(swmmoutt)
  swmmoutt_ts <- as.data.frame(swmmoutt_ts0$swmmoutt)
  colnames(swmmoutt_ts)<-"sim"
  result <- (sqrt(sum((obs$Observed - swmmoutt_ts$sim)^2)) / sqrt(sum((obs$Observed - mean(swmmoutt_ts$sim))^2)))*(-1) # RSR, Equation 1
  return(result)
}

# Single objective optimization using genetic algorithm
Ga_run <- ga(type = "real-valued", fitness = objectv, lower = c(0.009, 0.009, 0.011, 0.011, 0.019, 0.019, 0.011, 0.011, 
                                                                 0.01, 0.01,
                                                                 0.06, 0.02, 0.1, 0.33,
                                                                 0.27, 0.24, 0.25, 0.19,
                                                                 0.05, 0.1,
                                                                 0.05, 0.1,
                                                                 8.27, 8.27, 3.5, 2.4,
                                                                 1, 1), 
              upper = c(0.014, 0.017, 0.015, 0.017, 0.021, 0.033, 0.08, 0.08, 
                        0.2, 0.24,
                        0.07, 0.3, 0.3, 0.46,
                        0.35, 0.34, 0.347, 0.368,
                        0.4, 0.5,
                        0.5, 0.5,
                        14.2, 12.94, 13.55, 5,
                        5, 5), 
              seed = 123, maxiter = 200, popSize = 280, run = 100, parallel = dcores)
## End of part 2 of 2


