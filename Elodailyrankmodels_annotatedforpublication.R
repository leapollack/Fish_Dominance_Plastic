###Elo rank and foraging behavior main manuscript models##

library("ggplot2")
library("tidyverse")
library(dplyr)
library(tidyr)
library(forcats)
library("brms")
library("performance")
library("rptR")
library(tidybayes)
library("patchwork")

theme_set(theme_minimal() + 
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))


elodata <-read.csv("Ch3_Data_Clean.csv")%>%
  mutate(across(.cols = c(Group, FishID, Observer, Tank), .fns = as.factor),)

####Calculate daily rank based on daily Elo Score ####
elodata <- elodata %>% 
  arrange(Elo_daily) %>%
  group_by(Group, Date_Str) %>%
  mutate(daily_maxrank = rank(desc(Elo_daily), ties.method = "max"),
         daily_minrank = rank(desc(Elo_daily), ties.method = "min"),
         n_Elos = n_distinct(Elo_daily),
         has_tie = n_Elos != Group_Size) %>% 
  ungroup()


#check for any errant groups 
badgroups <- elodata %>%
  group_by(Group) %>%
  summarise(group_size = mean(Group_Size),
            trial = sum(Trial)) %>% 
  arrange(group_size) %>% 
  print(n=100)

##cacaluate group rank stability
elodata %>%
  group_by(Group) %>%
  summarise(Elo_stability = mean(Elo_stability),
            group_size = mean(Group_Size))%>% 
                    arrange(group_size) %>% 
                    print(n=100)

####Format variables####
elodata$Weight <-as.numeric(elodata$Weight)
elodata$Length <-as.numeric(elodata$Length)
elodata$daily_maxrank <- as.factor(elodata$daily_maxrank)
elodata$Weight_z <- scale(elodata$Weight, center=TRUE, scale = TRUE)
elodata$Length_z <- scale(elodata$Length, center=TRUE, scale = TRUE)

####Differences in likelihood to sample food first - familiar food####

#make a "first" column for outcome
elodata <- elodata %>% 
  mutate(first = if_else(Pellet_Order==1, T, F))

##model for groups of 2
m2_first <- brm(data = elodata %>% 
                  filter(Group_Size == 2), family = bernoulli,
                bf(first ~ Length_z + daily_maxrank + Trial + (1 | Group/FishID)),
                prior = c(set_prior("normal(0, 1)", class = "Intercept"),
                          set_prior("normal(0, 1)", class = "b")),
                iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(),
                control = list(adapt_delta = 0.99, max_treedepth = 15))

###contrast between ranks
m2_first_emm <- m2_first %>% 
  emmeans::emmeans(~ daily_maxrank) %>% 
  emmeans::contrast(method = "pairwise") %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = exp(.value)) 
###values
m2_first_emm  %>% 
  median_hdci()
###plot
pf1 <- m2_first %>% 
  emmeans::emmeans(~ daily_maxrank) %>% 
  emmeans::contrast(method = "pairwise") %>% 
  gather_emmeans_draws()%>% 
  mutate(.value = exp(.value))%>% 
  mutate(contrast = recode(contrast, "daily_maxrank1 - daily_maxrank2" = "1 / 2")) %>% 
  mutate(contrast = factor(contrast, levels = c("1 / 2"))) %>% 
  group_by(contrast) %>% 
  mutate(`q2.5` = quantile(.value, probs = 0.025),
         `q97.5` = quantile(.value, probs = 0.975),
         touch_1 = (`q2.5` <= 1 & `q97.5` >= 1)) %>% 
  ggplot(aes(y = contrast, x = .value, color = touch_1)) +
  stat_halfeye(.width = 0.95) +
  geom_vline(xintercept = 1, lty = 2, alpha = 0.3) +
  xlim(c(0,15)) +
  xlab("") +
  ylab("") +
  scale_color_manual("95% CI\noverlaps 1", values = c("orange")) +
  theme(legend.position = "none")

##model for groups of 3
m3_first <- brm(data = elodata %>% 
                  filter(Group_Size == 3), family = bernoulli,
                bf(first ~ Length_z + daily_maxrank + Trial + (1 | Group/FishID)),
                prior = c(set_prior("normal(0, 1)", class = "Intercept"),
                          set_prior("normal(0, 1)", class = "b")),
                iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(),
                control = list(adapt_delta = 0.99, max_treedepth = 15))

###contrast between ranks
m3_first_emm <- m3_first %>% 
  emmeans::emmeans(~ daily_maxrank) %>% 
  emmeans::contrast(method = "pairwise") %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = exp(.value)) 
###values
m3_first_emm  %>% 
  median_hdci()
###plot
pf2 <- m3_first %>% 
  emmeans::emmeans(~ daily_maxrank) %>% 
  emmeans::contrast(method = "pairwise") %>% 
  gather_emmeans_draws()%>% 
  mutate(.value = exp(.value))%>% 
  mutate(contrast = recode(contrast, 
                           "daily_maxrank1 - daily_maxrank2" = "1 / 2", "daily_maxrank1 - daily_maxrank3" = "1 / 3", "daily_maxrank2 - daily_maxrank3" = "2 / 3")) %>% 
  mutate(contrast = factor(contrast, levels = c("2 / 3", "1 / 3", "1 / 2"))) %>% 
  group_by(contrast) %>% 
  mutate(`q2.5` = quantile(.value, probs = 0.025),
         `q97.5` = quantile(.value, probs = 0.975),
         touch_1 = (`q2.5` <= 1 & `q97.5` >= 1)) %>% 
  ggplot(aes(y = contrast, x = .value, color = touch_1)) +
  stat_halfeyeh(.width = 0.95) +
  geom_vline(xintercept = 1, lty = 2, alpha = 0.3) +
  xlim(c(0,15)) +
  xlab("Odds ratio of eating familiar food first") +
  ylab("") +
  scale_color_manual("95% CI\noverlaps 1", values = c("springgreen4", "orange")) +
  theme(legend.position = "none")

##model for groups of 4
m4_first <- brm(data = elodata %>% 
                  filter(Group_Size == 4), family = bernoulli,
                bf(first ~ Length_z + daily_maxrank + Trial + (1 | Group/FishID)),
                prior = c(set_prior("normal(0, 1)", class = "Intercept"),
                          set_prior("normal(0, 1)", class = "b")),
                iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(),
                control = list(adapt_delta = 0.99, max_treedepth = 15))

###contrast between ranks
m4_first_emm <- m4_first %>% 
  emmeans::emmeans(~ daily_maxrank) %>% 
  emmeans::contrast(method = "pairwise") %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = exp(.value)) 
###values
m4_first_emm  %>% 
  median_hdci()
###plot
pf3 <- m4_first %>% 
  emmeans::emmeans(~ daily_maxrank) %>% 
  emmeans::contrast(method = "pairwise") %>% 
  gather_emmeans_draws()%>% 
  mutate(.value = exp(.value))%>% 
  mutate(contrast = recode(contrast, 
                           "daily_maxrank1 - daily_maxrank2" = "1 / 2", "daily_maxrank1 - daily_maxrank3" = "1 / 3", "daily_maxrank2 - daily_maxrank3" = "2 / 3", "daily_maxrank1 - daily_maxrank4"= "1 / 4", "daily_maxrank2 - daily_maxrank4" = "2 / 4", "daily_maxrank3 - daily_maxrank4" = "3 / 4")) %>% 
  mutate(contrast = factor(contrast, levels = c("3 / 4","2 / 4", "2 / 3", "1 / 4", "1 / 3", "1 / 2"))) %>% 
  group_by(contrast) %>% 
  mutate(`q2.5` = quantile(.value, probs = 0.025),
         `q97.5` = quantile(.value, probs = 0.975),
         touch_1 = (`q2.5` <= 1 & `q97.5` >= 1)) %>% 
  ggplot(aes(y = contrast, x = .value, color = touch_1)) +
  stat_halfeye(.width = 0.95) +
  geom_vline(xintercept = 1, lty = 2, alpha = 0.3) +
  xlim(c(0,15)) +
  xlab("") +
  ylab("") +
  scale_color_manual("95% CI\noverlaps 1", values = c("springgreen4", "orange")) +
  theme(legend.position = "none")

####Differences in likelihood to sample food first - novel food####

#make a novel "first" column for outcome
elodata <- elodata %>% 
  mutate(first_nov = if_else(NovelOrder==1, T, F))

##model for groups of 2
m2_first_nov <- brm(data = elodata %>% 
                      filter(Group_Size == 2, Trial > 5), family = bernoulli,
                    bf(first_nov ~ Length_z + daily_maxrank + as_factor(Trial)+ (1 | Group/FishID)),
                    prior = c(set_prior("normal(0, 1)", class = "Intercept"),
                              set_prior("normal(0, 1)", class = "b")),
                    iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(),
                    control = list(adapt_delta = 0.99, max_treedepth = 15))

###contrast between ranks
m2_first_nov_emm <- m2_first_nov %>% 
  emmeans::emmeans(~ daily_maxrank) %>% 
  emmeans::contrast(method = "pairwise") %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = exp(.value)) 
###values
m2_first_nov_emm  %>% 
  median_hdci()
###plot
pf4 <- m2_first_nov %>% 
  emmeans::emmeans(~ daily_maxrank) %>% 
  emmeans::contrast(method = "pairwise") %>% 
  gather_emmeans_draws()%>% 
  mutate(.value = exp(.value))%>% 
  mutate(contrast = recode(contrast, "daily_maxrank1 - daily_maxrank2" = "1 / 2")) %>% 
  mutate(contrast = factor(contrast, levels = c("1 / 2"))) %>% 
  group_by(contrast) %>% 
  mutate(`q2.5` = quantile(.value, probs = 0.025),
         `q97.5` = quantile(.value, probs = 0.975),
         touch_1 = (`q2.5` <= 1 & `q97.5` >= 1)) %>% 
  ggplot(aes(y = contrast, x = .value, color = touch_1)) +
  stat_halfeye(.width = 0.95) +
  geom_vline(xintercept = 1, lty = 2, alpha = 0.3) +
  xlim(c(0,15)) +
  xlab("") +
  ylab("") +
  scale_color_manual("95% CI\noverlaps 1", values = c("orange")) +
  theme(legend.position = "none")

##model for groups of 3
m3_first_nov <- brm(data = elodata %>% 
                      filter(Group_Size == 3, Trial > 5), family = bernoulli,
                    bf(first_nov ~ Length_z + daily_maxrank + + as_factor(Trial)+(1 | Group/FishID)),
                    prior = c(set_prior("normal(0, 1)", class = "Intercept"),
                              set_prior("normal(0, 1)", class = "b")),
                    iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(),
                    control = list(adapt_delta = 0.99, max_treedepth = 15))

###contrast between ranks
m3_first_nov_emm <- m3_first_nov %>% 
  emmeans::emmeans(~ daily_maxrank) %>% 
  emmeans::contrast(method = "pairwise") %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = exp(.value)) 
###values
m3_first_nov_emm  %>% 
  median_hdci()
###plot
pf5 <- m3_first_nov %>% 
  emmeans::emmeans(~ daily_maxrank) %>% 
  emmeans::contrast(method = "pairwise") %>% 
  gather_emmeans_draws()%>% 
  mutate(.value = exp(.value))%>% 
  mutate(contrast = recode(contrast, 
                           "daily_maxrank1 - daily_maxrank2" = "1 / 2", "daily_maxrank1 - daily_maxrank3" = "1 / 3", "daily_maxrank2 - daily_maxrank3" = "2 / 3")) %>% 
  mutate(contrast = factor(contrast, levels = c("2 / 3", "1 / 3", "1 / 2"))) %>% 
  group_by(contrast) %>% 
  mutate(`q2.5` = quantile(.value, probs = 0.025),
         `q97.5` = quantile(.value, probs = 0.975),
         touch_1 = (`q2.5` <= 1 & `q97.5` >= 1)) %>% 
  ggplot(aes(y = contrast, x = .value, color = touch_1)) +
  stat_halfeye(.width = 0.95) +
  geom_vline(xintercept = 1, lty = 2, alpha = 0.3) +
  xlim(c(0,15)) +
  xlab("Odds ratio of eating novel food first") +
  ylab("") +
  scale_color_manual("95% CI\noverlaps 1", values = c("orange")) +
  theme(legend.position = "none")

##model for groups of 4
m4_first_nov <- brm(data = elodata %>% 
                    filter(Group_Size == 4, Trial > 5), family = bernoulli,
                    bf(first_nov ~ Length_z + daily_maxrank + as_factor(Trial)+ (1 | Group/FishID)),
                    prior = c(set_prior("normal(0, 1)", class = "Intercept"),
                              set_prior("normal(0, 1)", class = "b")),
                    iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(),
                    control = list(adapt_delta = 0.99, max_treedepth = 15))

###contrast between ranks
m4_first_nov_emm <- m4_first_nov %>% 
  emmeans::emmeans(~ daily_maxrank) %>% 
  emmeans::contrast(method = "pairwise") %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = exp(.value)) 
###values
m4_first_nov_emm  %>% 
  median_hdci()
###plot
pf6 <- m4_first_nov %>% 
  emmeans::emmeans(~ daily_maxrank) %>% 
  emmeans::contrast(method = "pairwise") %>% 
  gather_emmeans_draws()%>% 
  mutate(.value = exp(.value))%>% 
  mutate(contrast = recode(contrast, 
                           "daily_maxrank1 - daily_maxrank2" = "1 / 2", "daily_maxrank1 - daily_maxrank3" = "1 / 3", "daily_maxrank2 - daily_maxrank3" = "2 / 3", "daily_maxrank1 - daily_maxrank4"= "1 / 4", "daily_maxrank2 - daily_maxrank4" = "2 / 4", "daily_maxrank3 - daily_maxrank4" = "3 / 4")) %>% 
  mutate(contrast = factor(contrast, levels = c("3 / 4","2 / 4", "2 / 3", "1 / 4", "1 / 3", "1 / 2"))) %>% 
  group_by(contrast) %>% 
  mutate(`q2.5` = quantile(.value, probs = 0.025),
         `q97.5` = quantile(.value, probs = 0.975),
         touch_1 = (`q2.5` <= 1 & `q97.5` >= 1)) %>% 
  ggplot(aes(y = contrast, x = .value, color = touch_1)) +
  stat_halfeye(.width = 0.95) +
  geom_vline(xintercept = 1, lty = 2, alpha = 0.3) +
  xlim(c(0,15)) +
  xlab("") +
  ylab("") +
  scale_color_manual("95% CI\noverlaps 1", values = c("springgreen4", "orange")) +
  theme(legend.position = "none")

####Fig 1: combining 6 plots ####
ppp <- pylab + ((pf1 + pf2 + pf3) / (pf4 + pf5 + pf6)) 

ppp[[1]] <- ppp[[1]] + plot_layout(tag_level = "new")

ppp + plot_annotation(tag_levels = "a") + plot_layout(widths = c(1,25)) & theme(axis.title.x = element_text(size = 14))

ggsave("fig1_new.pdf", width = 6, height = 6)

####Differences in number of bites - familiar food####

#overall model form
mform_zinb <- bf(Pellet_Bites ~ daily_maxrank + Trial + Length_z + (1 | Group/FishID),
                 zi ~ daily_maxrank + Trial + Length_z + (1 | Group/FishID))

##model for groups of 2
m2_zinb <- brm(mform_zinb, 
                data = elodata %>% 
                  filter(Group_Size == 2), 
                family = zero_inflated_negbinomial(),
                prior = c(set_prior("normal(0,3)", class = "b"),
                          set_prior("cauchy(0,2)", class = "sd")),
                cores = 3, chains = 3, iter = 5000, warmup = 1000,
                control = list(adapt_delta = 0.99, max_treedepth = 20))


###looking at contrasts between ranks
m2_emm <- m2_zinb %>% 
  emmeans::emmeans(~ daily_maxrank) %>% 
  emmeans::contrast(method = "pairwise") %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = exp(.value)) # put this on odds ratio scale instead of log odds
###values
m2_emm %>% 
  median_hdci()
###plot
pb1 <- m2_emm %>% 
mutate(contrast = recode(contrast, "daily_maxrank1 - daily_maxrank2" = "1 / 2")) %>% 
  mutate(contrast = factor(contrast, levels = c("1 / 2"))) %>% 
  group_by(contrast) %>% 
  mutate(`q2.5` = quantile(.value, probs = 0.025),
         `q97.5` = quantile(.value, probs = 0.975),
         touch_1 = (`q2.5` <= 1 & `q97.5` >= 1)) %>% 
  ggplot(aes(y = contrast, x = .value, color = touch_1)) +
  stat_halfeye(.width = 0.95) +
  geom_vline(xintercept = 1, lty = 2, alpha = 0.3) +
  xlim(c(0.4,2.5)) +
  xlab("") +
  ylab("Contrasts between dominance ranks") +
  scale_color_manual("95% CI\noverlaps 1", values = c("springgreen4")) +
  theme(legend.position = "none")

##model for groups of 3
m3_zinb <- brm(mform_zinb, 
                data = elodata %>% 
                filter(Group_Size == 3), 
                family = zero_inflated_negbinomial(),
                prior = c(set_prior("normal(0,3)", class = "b"),
                          set_prior("cauchy(0,2)", class = "sd")),
                cores = 3, chains = 3, iter = 5000, warmup = 1000,
                control = list(adapt_delta = 0.99, max_treedepth = 20))

### looking at contrasts between ranks
m3_emm <- m3_zinb %>% 
  emmeans::emmeans(~ daily_maxrank) %>% 
  emmeans::contrast(method = "pairwise") %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = exp(.value)) # put this on odds ratio scale instead of log odds
###values
m3_emm %>% 
  mutate(contrast = "higher / lower") %>% 
  median_qi()
###plot
pb2 <- m3_emm %>% 
  mutate(contrast_simple = "higher / lower") %>%
  mutate(contrast = recode(contrast, 
                           "daily_maxrank1 - daily_maxrank2" = "1 / 2", "daily_maxrank1 - daily_maxrank3" = "1 / 3", "daily_maxrank2 - daily_maxrank3" = "2 / 3")) %>% 
  mutate(contrast = factor(contrast, levels = c("2 / 3", "1 / 3", "1 / 2"))) %>% 
  group_by(contrast) %>% 
  mutate(`q2.5` = quantile(.value, probs = 0.025),
         `q97.5` = quantile(.value, probs = 0.975),
         touch_1 = (`q2.5` <= 1 & `q97.5` >= 1)) %>% 
  ggplot(aes(y = contrast, x = .value, color = touch_1)) +
  stat_halfeye(.width = 0.95) +
  geom_vline(xintercept = 1, lty = 2, alpha = 0.3) +
  xlim(c(0.4,2.5)) +
  xlab("Odds ratio of familiar bites taken") +
  ylab("") +
  scale_color_manual("95% CI\noverlaps 1", values = c("springgreen4", "orange")) +
  theme(legend.position = "none")

##model for group of 4
m4_zinb <- brm(mform_zinb, 
                data = elodata %>% 
                filter(Group_Size == 4), 
                family = zero_inflated_negbinomial(),
                prior = c(set_prior("normal(0,3)", class = "b"),
                          set_prior("cauchy(0,2)", class = "sd")),
                cores = 3, chains = 3, iter = 6000, warmup = 1000,
                control = list(adapt_delta = 0.9999, max_treedepth = 20))

### looking at contrasts between ranks
m4_emm <- m4_zinb %>% 
  emmeans::emmeans(~ daily_maxrank) %>% 
  emmeans::contrast(method = "pairwise") %>% 
  gather_emmeans_draws()%>% 
 mutate(.value = exp(.value))
###values
m4_emm %>% 
  median_hdci()
###plot 
pb3 <- m4_emm %>% 
  mutate(contrast = recode(contrast, 
                           "daily_maxrank1 - daily_maxrank2" = "1 / 2", "daily_maxrank1 - daily_maxrank3" = "1 / 3", "daily_maxrank2 - daily_maxrank3" = "2 / 3", "daily_maxrank1 - daily_maxrank4"= "1 / 4", "daily_maxrank2 - daily_maxrank4" = "2 / 4", "daily_maxrank3 - daily_maxrank4" = "3 / 4")) %>% 
  mutate(contrast = factor(contrast, levels = c("3 / 4","2 / 4", "2 / 3", "1 / 4", "1 / 3", "1 / 2"))) %>% 
  group_by(contrast) %>% 
  mutate(`q2.5` = quantile(.value, probs = 0.025),
         `q97.5` = quantile(.value, probs = 0.975),
         touch_1 = (`q2.5` <= 1 & `q97.5` >= 1)) %>% 
  ggplot(aes(y = contrast, x = .value, color = touch_1)) +
  stat_halfeye(.width = 0.95) +
  geom_vline(xintercept = 1, lty = 2, alpha = 0.3) +
  xlim(c(0.4,2.5)) +
  xlab("") +
  ylab("") +
  scale_color_manual("95% CI\noverlaps 1", values = c("springgreen4", "orange")) +
  theme(legend.position = "none")

####Differences in number of bites - novel food####

#overall model form
mform_nov <- bf(NovelBites ~ daily_maxrank + as_factor(Trial)+ Length_z + (1 | Group/FishID),
                 zi ~ daily_maxrank + as_factor(Trial) + Length_z + (1 | Group/FishID))  

##model for groups of 2  
m2_nov <- brm(mform_nov, 
            data = elodata %>% 
            filter(Group_Size == 2, Trial > 5), 
            family = zero_inflated_poisson(),
            prior = c(set_prior("normal(0,3)", class = "b"),
                       set_prior("cauchy(0,2)", class = "sd")),
           cores = 3, chains = 3, iter = 5000, warmup = 1000,
           control = list(adapt_delta = 0.999, max_treedepth = 20))

### looking at contrasts between ranks
m2_novrank_emm <- m2_nov %>% 
  emmeans::emmeans(~ daily_maxrank) %>% 
  emmeans::contrast(method = "pairwise") %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = exp(.value))
###values
m2_novrank_emm  %>% 
  median_hdci()
###plot 
pb4 <- m2_novrank_emm %>% 
  mutate(contrast = recode(contrast, "daily_maxrank1 - daily_maxrank2" = "1 / 2")) %>% 
  mutate(contrast = factor(contrast, levels = c("1 / 2"))) %>% 
  group_by(contrast) %>% 
  mutate(`q2.5` = quantile(.value, probs = 0.025),
         `q97.5` = quantile(.value, probs = 0.975),
         touch_1 = (`q2.5` <= 1 & `q97.5` >= 1)) %>% 
  ggplot(aes(y = contrast, x = .value, color = touch_1)) +
  stat_halfeye(.width = 0.95) +
  geom_vline(xintercept = 1, lty = 2, alpha = 0.3) +
  xlim(c(0.4,2.5)) +
  xlab("") +
  ylab("") +
  scale_color_manual("95% CI\noverlaps 1", values = c("springgreen4", "orange")) +
  theme(legend.position = "none")


##model for group of 3 
m3_nov <- brm(mform_nov, data = elodata %>% 
               filter(Group_Size == 3, Trial > 5), 
               family = zero_inflated_poisson(),
               prior = c(set_prior("normal(0,3)", class = "b"),
                         set_prior("cauchy(0,2)", class = "sd")),
               cores = 3, chains = 3, iter = 5000, warmup = 1000,
               control = list(adapt_delta = 0.999, max_treedepth = 20))

### looking at contrasts between ranks
m3_novrank_emm <- m3_nov %>% 
  emmeans::emmeans(~ daily_maxrank) %>% 
  emmeans::contrast(method = "pairwise") %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = exp(.value)) 
###values
m3_novrank_emm  %>% 
  median_hdci()
###plot 
pb5 <- m3_novrank_emm %>% 
  mutate(contrast = recode(contrast, "daily_maxrank1 - daily_maxrank2" = "1 / 2", "daily_maxrank1 - daily_maxrank3" = "1 / 3", "daily_maxrank2 - daily_maxrank3" = "2 / 3")) %>% 
  mutate(contrast = factor(contrast, levels = c("2 / 3", "1 / 3", "1 / 2"))) %>% 
  group_by(contrast) %>% 
  mutate(`q2.5` = quantile(.value, probs = 0.025),
         `q97.5` = quantile(.value, probs = 0.975),
         touch_1 = (`q2.5` <= 1 & `q97.5` >= 1)) %>% 
  ggplot(aes(y = contrast, x = .value, color = touch_1)) +
  stat_halfeye(.width = 0.95) +
  geom_vline(xintercept = 1, lty = 2, alpha = 0.3) +
  xlim(c(0.4,2.5)) +
  xlab("Odds ratio of novel bites taken") +
  ylab("") +
  scale_color_manual("95% CI\noverlaps 1", values = c("orange")) +
  theme(legend.position = "none")

##model for groups of 4
m4_nov <- brm(mform_nov, data = elodata %>% 
               filter(Group_Size == 4, Trial > 5), 
               family = zero_inflated_poisson(),
               prior = c(set_prior("normal(0,3)", class = "b"),
                         set_prior("cauchy(0,2)", class = "sd")),
               cores = 3, chains = 3, iter = 5000, warmup = 1000,
               control = list(adapt_delta = 0.999, max_treedepth = 20))

### looking at contrasts between ranks
m4_novrank_emm <- m4_nov %>% 
  emmeans::emmeans(~ daily_maxrank) %>% 
  emmeans::contrast(method = "pairwise") %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = exp(.value)) 
###values
m4_novrank_emm  %>% 
  median_hdci()
###plot 
pb6 <- m4_novrank_emm %>% 
  mutate(contrast = recode(contrast, 
                           "daily_maxrank1 - daily_maxrank2" = "1 / 2", "daily_maxrank1 - daily_maxrank3" = "1 / 3", "daily_maxrank2 - daily_maxrank3" = "2 / 3", "daily_maxrank1 - daily_maxrank4"= "1 / 4", "daily_maxrank2 - daily_maxrank4" = "2 / 4", "daily_maxrank3 - daily_maxrank4" = "3 / 4")) %>% 
  mutate(contrast = factor(contrast, levels = c("3 / 4","2 / 4", "2 / 3", "1 / 4", "1 / 3", "1 / 2"))) %>% 
  group_by(contrast) %>% 
  mutate(`q2.5` = quantile(.value, probs = 0.025),
         `q97.5` = quantile(.value, probs = 0.975),
         touch_1 = (`q2.5` <= 1 & `q97.5` >= 1)) %>% 
  ggplot(aes(y = contrast, x = .value, color = touch_1)) +
  stat_halfeye(.width = 0.95) +
  geom_vline(xintercept = 1, lty = 2, alpha = 0.3) +
  xlim(c(0.4,2.5)) +
  xlab("") +
  ylab("") +
  scale_color_manual("95% CI\noverlaps 1", values = c("springgreen4", "orange")) +
  theme(legend.position = "none")

pylab <- ggplot(data.frame(l = "Contrasts between dominance ranks", x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), angle = 90, size = 5) +
  theme_void() +
  coord_cartesian(clip = "off")

####Fig 2: combining 6 plots ####
pp <- pylab + ((pb1 + pb2 + pb3) / (pb4 + pb5 + pb6)) 

pp[[1]] <- pp[[1]] + plot_layout(tag_level = "new")

pp + plot_annotation(tag_levels = "a") + plot_layout(widths = c(1,25)) & theme(axis.title.x = element_text(size = 14))

ggsave("fig2_new.pdf", width = 6, height = 6)

####Differences in foraging across novel food types ####
#using models from above, just looking at differences between trials instead of ranks

## model for groups of 2
### looking at contrasts between trials
m2_nov_emm <- m2_nov %>% 
  emmeans::emmeans(~ Trial) %>% 
  emmeans::contrast(method = "pairwise") %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = exp(.value)) 
###values
m2_nov_emm  %>% 
  median_hdci()
###plot
m2_nov_emm %>% 
  mutate(contrast = recode(contrast, "6 - 7" = "shrimp/bead", "7 - 8" = "bead/pine", "8 - 9"= "pine/virgin", "9 - 10"="virgin/biofouled", "6 - 8" = "shrimp/pine", "6 - 9" = "shrimp/virgin", "6 - 10" = "shrimp/biofouled", "7 - 9" = "bead/virgin", "7 - 10"= "bead/biofouled", "8 - 10" = "pine/biofouled" )) %>% 
  group_by(contrast) %>% 
  mutate(`q2.5` = quantile(.value, probs = 0.025),
         `q97.5` = quantile(.value, probs = 0.975),
         touch_1 = (`q2.5` <= 1 & `q97.5` >= 1)) %>% 
  ggplot(aes(y = contrast, x = .value, color = touch_1)) +
  stat_halfeye(.width = 0.95) +
  geom_vline(xintercept = 1, lty = 2, alpha = 0.3) +
  theme_minimal() +
  xlab("Ratio of novel bites taken") +
  ylab("Contrasts between novel foods") +
  scale_color_manual("95% CI\noverlaps 1", values = c("springgreen4", "orange")) +
  theme(legend.position = "none")

##model for groups of 3
### looking at contrasts between trials
m3_nov_emm <- m3_nov %>% 
  emmeans::emmeans(~ Trial) %>% 
  emmeans::contrast(method = "pairwise") %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = exp(.value)) 
###values
m3_nov_emm  %>% 
  median_hdci()
###plot 
m3_nov_emm %>% 
  mutate(contrast = recode(contrast, "6 - 7" = "shrimp/bead", "7 - 8" = "bead/pine", "8 - 9"= "pine/virgin", "9 - 10"="virgin/biofouled", "6 - 8" = "shrimp/pine", "6 - 9" = "shrimp/virgin", "6 - 10" = "shrimp/biofouled", "7 - 9" = "bead/virgin", "7 - 10"= "bead/biofouled", "8 - 10" = "pine/boifouled" )) %>% 
  group_by(contrast) %>% 
  mutate(`q2.5` = quantile(.value, probs = 0.025),
         `q97.5` = quantile(.value, probs = 0.975),
         touch_1 = (`q2.5` <= 1 & `q97.5` >= 1)) %>% 
  ggplot(aes(y = contrast, x = .value, color = touch_1)) +
  stat_halfeye(.width = 0.95) +
  geom_vline(xintercept = 1, lty = 2, alpha = 0.3) +
  theme_minimal() +
  xlab("Ratio of novel bites taken") +
  ylab("Contrasts between novel foods") +
  scale_color_manual("95% CI\noverlaps 1", values = c("springgreen4", "orange")) +
  theme(legend.position = "none")

##model for groups of 4
### looking at contrasts between trials
m4_nov_emm <- m4_nov %>% 
  emmeans::emmeans(~ Trial) %>% 
  emmeans::contrast(method = "pairwise") %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = exp(.value)) 
###values
m4_nov_emm  %>% 
  median_hdci()
###plot
m4_nov_emm %>% 
  mutate(contrast = recode(contrast, "6 - 7" = "shrimp/bead", "7 - 8" = "bead/pine", "8 - 9"= "pine/virgin", "9 - 10"="virgin/biofouled", "6 - 8" = "shrimp/pine", "6 - 9" = "shrimp/virgin", "6 - 10" = "shrimp/biofouled", "7 - 9" = "bead/virgin", "7 - 10"= "bead/biofouled", "8 - 10" = "pine/boifouled" )) %>% 
  group_by(contrast) %>% 
  mutate(`q2.5` = quantile(.value, probs = 0.025),
         `q97.5` = quantile(.value, probs = 0.975),
         touch_1 = (`q2.5` <= 1 & `q97.5` >= 1)) %>% 
  ggplot(aes(y = contrast, x = .value, color = touch_1)) +
  stat_halfeye(.width = 0.95) +
  geom_vline(xintercept = 1, lty = 2, alpha = 0.3) +
  theme_minimal() +
  xlab("Ratio of novel bites taken") +
  ylab("Contrasts between novel foods") +
  scale_color_manual("95% CI\noverlaps 1", values = c("springgreen4", "orange"))