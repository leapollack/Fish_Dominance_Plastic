####Script for how daily Elo score was calculated (how Elo score was turned into daily ranks can be found in Elodailyrankmodels.R)

library(EloRating)
library(tidyverse)

###First read in all chasing data###
chases <- read.csv("Ch3_ChasingData.csv") %>%
  mutate(across(.cols = c(Date, Winner, Loser), .fns = as.factor),
         rowid = row_number())

chases_s <- chases %>% 
  filter(str_detect(Date, "/")) %>% 
  mutate(Date = lubridate::mdy(Date)) 

chases_d <- chases %>% 
  filter(str_detect(Date, "[a-z]")) %>% 
  mutate(Date = lubridate::dmy(Date))

chases <- bind_rows(chases_s, chases_d)

#example:
#Group 604
Group604 <- subset(chases, GroupID == 604, select = c("Date", "Winner", "Loser"))

#confirm you got the right data frame
Group604 

#Now use elo.seq to calculate daily elo rankings and assign it a name
Group604_elo <- elo.seq(winner=Group604$Winner, loser=Group604$Loser, Date=Group604$Date, runcheck=FALSE)

summary(Group604_elo)

#repeat for all groups

# removing groups 210 and 220 for this calculation, as each of these groups only had a single chase. This single chase will be used to determine their Elo scores manually (1050 for winner, 950 for loser)
elo_all <- chases %>% 
  filter(!GroupID %in% c(210, 220)) %>% 
  split(., .$GroupID) %>% 
  map(function(x) {
    if(nrow(x) < 2) return("Cannot calculate Elo for a single date")
    elo.seq(x$Winner, x$Loser, x$Date, runcheck = FALSE, startvalue = 1000)
  })


elo_all$`631`$lmat

make_tidy_elo <- function(x){
  
  scores <- x[["lmat"]]
  dates <- x[["truedates"]]
  
  scores <- as.data.frame(scores)
  scores[["Date"]] <- dates
  
  stability <- stab_elo(x)
  scores[["Elo_stability"]] <- stability
  
  pivot_longer(scores, cols = -c("Date", "Elo_stability"), 
               names_to = "FishID", values_to = "Elo_daily")
}



map(elo_all, make_tidy_elo) %>% 
  bind_rows(.id = "GroupID") %>% 
  write_csv("Ch3_Elo_daily.csv")
