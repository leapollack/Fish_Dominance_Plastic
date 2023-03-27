####Script for how daily Elo score was calculated (how Elo score was turned into daily ranks can be found in Elodailyrankmodels.R)

library(EloRating)

###First read in all chasing data###
chases <-read.csv("Ch3_ChasingData.csv") %>%
  mutate(across(.cols = c(Date, Winner, Loser), .fns = as.factor),)
 
#example:
#Group 604
Group604 <- subset(chases, GroupID == 604, select = c("Date", "Winner", "Loser"))

#confirm you got the right data frame
Group604 

#Now use elo.seq to calculate daily elo rankings and assign it a name
Group604_elo <- elo.seq(winner=Group604$Winner, loser=Group604$Loser, Date=Group604$Date, runcheck=FALSE)

summary(Group604_elo)

#repeat for all groups