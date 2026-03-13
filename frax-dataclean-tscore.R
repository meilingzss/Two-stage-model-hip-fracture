library(readr)
fraxicd <- read_csv("~/Desktop/2026twostage/2026newy/fraxicd.csv")
library(dplyr)
library(stringr)
draxicd <- fraxicd
library(dplyr)
library(stringr)
#######Region
#ICD-10 Codes
#Hip: S72.0, S72.1, S72.2
#S72.0: Fracture of neck of femur (femoral neck)
#S72.1: Pertrochanteric fracture
#S72.2: Subtrochanteric fracture
#Proximal Humerus (Upper Arm): S42.2
#Other Humerus: S42.3–S42.4
#•    S42.2 → Fracture of upper end of humerus (proximal humerus)
#•    S42.3 → Fracture of shaft of humerus
#•    S42.4 → Fracture of lower end of humerus
#Wrist (Distal Radius/Ulna:S52.5, S52.6
#•    S52.0 → Upper end of ulna
#•    S52.1 → Upper end of radius
#•    S52.2 → Shaft of ulna
#•    S52.3 → Shaft of radius
#•    S52.4 → Shaft of both bones

########
draxicd <- draxicd %>%
  mutate(
    SLDFX = if_else(
      str_detect(p41270, "S42\\.2|S42\\.3|S42\\.4"),
      1, 0
    ),
    WRSTFX = if_else(
      str_detect(p41270, "S52\\.5|S52\\.6"),
      1, 0
    )
  )
table(draxicd$SLDFX)
table(draxicd$WRSTFX)

draxicd <- draxicd %>%
  mutate(p41270 = if_else(str_detect(p41270, "S72\\.0|S72\\.1|S72\\.2"),
                          1,
                          0))
library(dplyr)
library(stringr)

table(draxicd$p41270)

colnames(draxicd) <- c(
  "ID",
  "height",
  "sex",
  "smoke",
  "weight",
  "age",
  "leftneckbmd",
  "rightneckbmd",
  "hip"
)
#UK <- draxicd[, -1]
UK <- draxicd
# Remove the least frequent level for 'walk'


# Remove the least frequent level for 'smoke'
UK <- UK[UK$smoke != "Prefer not to answer", ]


table(UK$smoke)
table(UK$hip)
## mros use right hip: MrOS DXA Quality Assurance Manual for Hologic QDR-4500 Bone Densitometers
## sof use left hip:https://pmc.ncbi.nlm.nih.gov/articles/PMC4388249/
SOF_clinical <- read.table("/Users/mac/Desktop/dissertation/cooperation/SOF_data/sof.txt",
                           header = TRUE, sep = "\t", na.strings = c("", "NA", " "))

MrOS_clinical <- read.table("/Users/mac/Desktop/dissertation/cooperation/MrOS_data/mros1.txt",
                            header = TRUE, sep = "\t", na.strings = c("", "NA", " "))
summary(SOF_clinical$HA_WLKSPED)
summary(MrOS_clinical$HA_WLKSPED)
### Walking speed Approach 1: Based on quantiles

# categorize HA_WLKSPED into three groups:
#•    Slow pace: Below the 1st quartile (Q1)
#•    Steady average pace: Between the 1st quartile (Q1) and the 3rd quartile (Q3)
#•    Brisk pace: Above the 3rd quartile (Q3)

# Define thresholds based on quantiles of HA_WLKSPED
Q1 <- quantile(SOF_clinical$HA_WLKSPED, 0.25, na.rm = TRUE)
Q3 <- quantile(SOF_clinical$HA_WLKSPED, 0.75, na.rm = TRUE)

# Create categorical variable 'walk'
SOF_clinical$walk <- cut(SOF_clinical$HA_WLKSPED,
                         breaks = c(-Inf, Q1, Q3, Inf),
                         labels = c("Slow pace", "Steady average pace", "Brisk pace"))

# Check distribution
table(SOF_clinical$walk)

# Define thresholds based on quantiles of HA_WLKSPED
Q1 <- quantile(MrOS_clinical$HA_WLKSPED, 0.25, na.rm = TRUE)
Q3 <- quantile(MrOS_clinical$HA_WLKSPED, 0.75, na.rm = TRUE)

# Create categorical variable 'walk'
MrOS_clinical$walk <- cut(MrOS_clinical$HA_WLKSPED,
                          breaks = c(-Inf, Q1, Q3, Inf),
                          labels = c("Slow pace", "Steady average pace", "Brisk pace"))

# Check distribution
table(MrOS_clinical$walk)


## in MrOS_clinical, change level 1.25 to 1,  change level 2.5 to 2, change level 3.75 to 3, change levle 5 to 4.   in SOF_clinical, change level 5 to 4.
# Load dplyr if not already loaded
library(dplyr)

# Recode levels for MrOS_clinical
MrOS_clinical$HA_IADL51 <- recode(MrOS_clinical$HA_IADL51,
                                  `1.25` = 1,
                                  `2.5` = 2,
                                  `3.75` = 3,
                                  `5` = 4)

# Recode level 5 to 4 for SOF_clinical
library(dplyr)

SOF_clinical$HA_IADL51[SOF_clinical$HA_IADL51 == 5] <- 4



UK_male   <- subset(UK, sex == "Male")
UK_female <- subset(UK, sex == "Female")
UK <- UK_male

UK$newneck <- ifelse(UK$sex == "Male", UK$rightneckbmd, UK$leftneckbmd)

### T score for UK
# Function to compute T-score
compute_t_score <- function(value, mean_val, sd_val) {
  (value - mean_val) / sd_val
}

UK$newneck <- ifelse(UK$sex == "Male", UK$rightneckbmd, UK$leftneckbmd)
# Initialize column for T-scores
UK$neckbmd_Tscore <- NA

# Define age groups
age_groups <- list("65_71" = c(65, 71), "71_86" = c(71, 86))


# Loop through each age group
for (age_label in names(age_groups)) {
  age_range <- age_groups[[age_label]]
  
  # Subset data for the current age group and HA_HIPFX == 0
  subset_data <- UK[UK$age >= age_range[1] &
                      UK$age <= age_range[2] &
                      UK$hip == 0, ]
  
  # Compute mean and standard deviation for spinebmd
  mean_val <- mean(subset_data$newneck, na.rm = TRUE)
  sd_val <- sd(subset_data$newneck, na.rm = TRUE)
  
  # Apply T-score transformation for matching individuals in the UK dataset
  match_idx <- which(UK$age >= age_range[1] &
                       UK$age <= age_range[2])
  
  UK$neckbmd_Tscore[match_idx] <- compute_t_score(UK$newneck[match_idx], mean_val, sd_val)
}




colnames(UK)
colnames(MrOS_clinical)


UK$smoke[UK$smoke == "Never"] <- 0
UK$smoke[UK$smoke == "Previous"] <- 1
UK$smoke[UK$smoke == "Current"] <- 2

library(dplyr)

library(dplyr)

# Specify only columns that exist in UK
keep_cols <- c("ID",
"height",
"sex",
"smoke",
"weight",
"age",
"neckbmd_Tscore",
"hip"
)

# Keep only those columns
UK <- UK[, keep_cols]

write.csv(UK, "/Users/mac/Desktop/2026twostage/2026newy/idUKfrax_male_newtscore.8var.dataset.csv", row.names = FALSE)


