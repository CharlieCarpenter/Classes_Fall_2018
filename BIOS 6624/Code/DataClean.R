#############################################################
## Program Name: DataClean.R							                 ##
## Purpose: To read in the dataset and label the variables ##
## 			   Create any new variables.					             ##
## Output: A Clean .csv file is exported for reading in to ##
##         all other analysis files.					             ##
## Created by: Nichole Carlson (R Code by: Kevin Josey)		 ##
#############################################################

# READ IN THE RAW DATA AND CREATE VARIABLE TRANSFORMS
# Note: This is the place where I would come back to
#       for any data cleaning, outlier removal more
#       transforms.  This keeps all data cleaning   
#		    and variable decisions in one place.


# Dependencies
library(readr)
library(dplyr)
library(magrittr)

# psa <- read_csv("C:/Repositories/Bios6623ClassExamples/PSA-Project/G1Analysis/DataRaw/prostate.csv")
psa <- read_csv("~/Downloads/prostate.csv")

colnames(psa) <- tolower(colnames(psa))

psaclean <- psa %>%
  mutate(lpsa = log(psa),
         grade6 = ifelse(gleason == 6, 1, 0),
         grade7 = ifelse(gleason == 7, 1, 0),
         grade8 = ifelse(gleason == 8, 1, 0)) %>%
  mutate(gd6ageint = age*grade6,
         gd7ageint = age*grade7,
         gd8ageint = age*grade8)
         
# export the clean dataset for use in other analysis programs
# write_csv(psaclean, "C:/Repositories/Bios6623ClassExamples/PSA-Project/G1Analysis/Data/psaclean.csv")
write_csv(psaclean, "~/Downloads/psaclean.csv")

# further data cleaning based on reviewing descriptive statistics
psaclean2 <- psaclean %>%
  filter(wt <= 400) %>%
  mutate(cavolsviint = cavol*svi,
         wtsviint = wt*svi,
         agesviint = age*svi,
         bphsviint = bph*svi,
         cappensviint = cappen*svi,
         gd6sviint = grade6*svi,
         gd7sviint = grade7*svi,
         gd8sviint = grade8*svi)

# export the clean dataset for use in other analysis programs
# write_csv(psaclean2, "C:/Repositories/Bios6623ClassExamples/PSA-Project/G1Analysis/Data/psaclean2.csv")
write_csv(psaclean2, "~/Downloads/psaclean2.csv")
