##############################################################################################
## R code written by Connie Li Wai Suen.                                                    ##
## Date: 28 May 2018                                                                        ##
##                                                                                          ##
## Description:                                                                             ##
## This R script fits a 5-parameter logistic standard curve                                 ##
## to the dilutions of the positive controls for each protein and                           ##
## converts the median fluorescence intensity (MFI) values into relative antibody units.    ##
##                                                                                          ##
## Input:                                                                                   ##
## - Luminex-200 or Magpix output file (required)                                           ##
## - Plate layout with bleedcode information (optional)                                     ##
##                                                                                          ##
## Output:                                                                                  ##
## - csv file of relative antibody unit for each protein,                                   ##
##   its reciprocal, the minimum and maximum standard,                                      ## 
##   and the minimum and maximum dilution for that protein.                                 ##
##                                                                                          ##
##############################################################################################

rm(list=ls())
## R packages are installed only if required before being loaded them.
Rpackages <- c("readxl", "drc", "plyr")
lapply(Rpackages, function(x){ require(x, character.only=TRUE) || install.packages(x, dependencies=TRUE)})
library(readxl)
library(drc)
library(plyr)



#######################################
## CHANGE FILE NAMES IN THIS SECTION ##
#######################################

#setwd("~/Desktop/test")
#setwd("C:/Users/liwaisuen/Documents/Rhea_Longley/Luminex")
# Make sure these files are in the same working directory as the R script
# Or set working directory with full path of file location above.
setwd("~/Documents/Protein_Expression_Taka/Experimental_data/Negative_Control/raw_data/")
Luminex.ouput_filename <- "negative_aus_redcross_2016_24112020"
Plate.layout_filename <- "Albinama_Run_2_13102020_layout"


#######################################################################################################################


fname1 <- paste0(Luminex.ouput_filename, ".xlsx")
fname2 <- paste0(Plate.layout_filename, ".xlsx")

## The first 41 rows are not relevant and are not imported into R.
L <- as.data.frame(read_excel(fname1, skip = 41))
## Find all blank rows (i.e. rows that are all NA).
## Then keep rows preceding the first blank row.
blank.row.number <- which(rowSums(is.na(L)) == length(names(L)))[1]
L <- L[1:(blank.row.number-1),]
## Exclude column that corresponds to "Total events"
L <- L[, !(colnames(L) %in% c("Total Events"))]
dim(L)

## Plate layout
layout <- as.data.frame(read_excel(fname2)); layout

## Protein name list obtained after removing variable names: "Location" and "Sample"
proteins <- names(L[, -c(1:2)]); proteins
## Add new variable to data frame to indicate the first letter of sample type ("B", "C", "S", "U")
## "B" = blank, "C"=control, "S"=standard (dilution of the pool), "U"=sample
L$type.letter <- substr(L$Sample, start=1, stop=1)
dilution <- c(1/50, 1/100, 1/200, 1/400, 1/800, 1/1600, 1/3200, 1/6400, 1/12800, 1/25600)
dilution.scaled <- dilution*25600; dilution.scaled

results.df.wide <- NULL
for (i in proteins){
  results.df <- NULL
  ## Taking the mean of duplicates for each standard
  ## and storing in object std in the following order: S1, S2, S3, ..., S9, S10.
  std <- NULL
  b <- c <- d <- e <- NULL
  for (r in 1:nrow(L)){
    if (L$type.letter[r]=="S"){
      std <- c(std, as.numeric(L[r,i])) 
    }
  }
  
  ## Log-log model to obtain a more linear relationship
  ## and therefore make it easier to interpolate around the lower asymptote.
  log.std <- log(as.numeric(std)); log.std
  
  ## Five-parameter logistic function is given by the expression:
  ## f(x) = c + \frac{d-c}{(1+\exp(b(\log(x)-\log(e))))^f}
  model1 <- drm(log.std ~ dilution, fct=LL.5(names=c("b", "c", "d", "e", "f")))
  summary(model1)
  
  Sys.sleep(0.1)  ## Suspends execution for 0.1 second to prevent RStudio errors when plotting within the loop.
  plot(model1, main=i)
  
  ## F(x) = ((A-D)/(1+((x/C)^B))) + D    ## where A=minimum asymptote, B=Hill slope, C=ED50, D=Maximum asymptote
  ## x = C*(((A-D)/(F(x)-D))-1)^(1/B) = e*(((c-d)/(log(mfi.X)-d))-1)^(1/b)
  b <- coef(summary(model1))[1]; b  ## slope
  c <- coef(summary(model1))[2]; c  ## lower asymptote
  d <- coef(summary(model1))[3]; d  ## upper asymptote
  e <- coef(summary(model1))[4]; e  ## ED50
  f <- coef(summary(model1))[5]; f  ## asymmetry parameter (f=1 for 4PL curves)
  
  for (r in 1:nrow(L)){
    results <- NULL
    if (L$type.letter[r]=="U"){
      mfi.X <- as.numeric(L[r, i])
      y <- log(mfi.X)
      
      if (y > max(log.std)) {
        dil.X <- max(dilution)
      } else {
        dil.X <- e*(( ((d-c)/(y-c))^(1/f) - 1 )^(1/b) )
      }
      dil.X <- ifelse(dil.X > 0.02, 0.02, dil.X)
      dil.X <- ifelse((is.na(dil.X) & y>log.std[2]), 0.02, dil.X)       ## Setting observations with very high MFI to 1/50.
      dil.X <- ifelse(dil.X < 1/51200, 1/51200, dil.X)
      dil.X <- ifelse((is.na(dil.X) & y<max(log.std)), 1/51200, dil.X)  ## Setting observations with very low MFI to 1/51200.
      location.X  <- L[r, "Location"]
      sample.X  <- L[r, "Sample"]
      results <- cbind(Location=location.X, Sample=sample.X, MFI=mfi.X, Dilution=dil.X, DilutionReciprocal=1/dil.X, MinStd=min(std), MaxDilution=min(dilution), MaxStd=max(std), MinDilution=max(dilution))
      results.colnames <- c("Location", "Sample", paste0(i, "_", c("MFI", "Dilution", "DilutionReciprocal", "MinStd", "MaxDilution", "MaxStd", "MinDilution")))
      colnames(results) <- results.colnames
    }
    results.df <- rbind(results.df, results) 
  }
  if (is.null(results.df.wide)){
    results.df.wide <- results.df
  } else { results.df.wide <- merge(results.df.wide, results.df, by=c("Location", "Sample")) }
}
results.df.wide <- as.data.frame(results.df.wide)
results.location <- matrix(unlist(strsplit(as.character(results.df.wide$Location), ",")), ncol=2, byrow=T)[,2]
results.location <- substr(results.location, 1, nchar(results.location)-1)
results.df.wide <- cbind(Location.2=results.location, results.df.wide)

## Matching bleedcode from plate layout to corresponding sample.
location.1 <- matrix(unlist(strsplit(L$Location, ",")), ncol=2, byrow=T)[,2]
location.1 <- substr(location.1, 1, nchar(location.1)-1)
location.2 <- data.frame(Location.2=location.1, alpha=gsub("[[:digit:]]", "", location.1), numeric=gsub("[^[:digit:]]", "", location.1), Bleedcode=NA, stringsAsFactors = FALSE)
for (i in location.2[, "Location.2"]){
  location.2[location.2$Location.2==i, "Bleedcode"] <- layout[layout$Plate==location.2[location.2$Location.2==i, "alpha"], colnames(layout)==location.2[location.2$Location.2==i, "numeric"]]
}

## Using join() from plyr package to add bleedcode information to results.df.wide.
if (length(Plate.layout_filename)!=0){
  results.df.wide <- join(results.df.wide, location.2[,c("Location.2", "Bleedcode")], by="Location.2", type="left") 
  ## Move bleedcode to first column
  results.df.wide <- results.df.wide[, c("Bleedcode", colnames(results.df.wide)[!(colnames(results.df.wide) %in% "Bleedcode")])]
  head(results.df.wide)
  ## Write out to csv file.
  write.csv(results.df.wide, file=paste(Luminex.ouput_filename, "5PL_results_bleed.csv", sep="_"), row.names=F)
} else write.csv(results.df.wide, file=paste(Luminex.ouput_filename, "5PL_results.csv", sep="_"), row.names=F)


#######################################################################################################################

