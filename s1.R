#author: Long Bui, Nguyen Thanh Truong
#mail: longbui189@gmail.com, truong2nt@gmail.com
#date 25/03/2020
#description a draft

library(tidyverse)
library(incidence)
library(EpiEstim)
library(zoo)
library(lubridate)
library(epitrix)
library(R0)
library(lubridate)
library(earlyR)
library(projections)
library(ggplotify)
library(magrittr)

set.seed(12345)

# read data
vietnam <- read_excel("vietnam-2403.xlsx") %>% janitor::clean_names()
first16 <- read_csv("first16.csv")
si <- read_csv("si.csv") # serial intervel

#data preprocessing
vietnam$type = ""
vietnam$type <- ifelse(vietnam$been_to_wuhan_hubei == "Yes" | vietnam$been_to_covid_19_countries_not_include_china== "Yes", "imported" ,"local")
vietnam$place_of_detect <- na.locf(vietnam$place_of_detect) #filling NA
vietnam$place_of_treatment <- na.locf(vietnam$place_of_treatment) #filling NA
vietnam$date<- na.locf(vietnam$date) #filling NA
first16$onset <- as.Date(first16$onset, format = "%d/%m/%Y")
vietnam$date <- as.Date(vietnam$date, format = "%d/%m/%Y")
secondp <-vietnam[!vietnam$case %in% first16$id,] # create secon phase

#tablulate
prop.table(table(secondp$type))

#create incidence object

onset = c(first16$onset, vietnam[17:123,]$date)
total.inc <- incidence(onset)
plot(total.inc) #plot epicurve

epicurve <- plot(imported.inc,  color ="Red") + 
  theme_classic() +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))  +
  scale_x_date(date_breaks = "4 days", date_labels = "%d %b") 


#epicurve by gender
gender.inc <- incidence(onset, groups =vietnam$gender)
plot(gender.inc)

#epicurve by imported
imported.inc <-incidence(onset, groups = vietnam$type)
plot(imported.inc)


#epicurve by province
prov.inc <- incidence(vietnam$date, groups = vietnam$place_of_detect)
plot(prov.inc)

# fit generation time
vsi <- c(3,1,7,11,6,8)
si_fit <- fit_disc_gamma(vsi)
si_fit

params <- gamma_mucv2shapescale(si_fit$mu, si_fit$cv)

GT.ncov = generation.time("gamma", c(params$shape, params$scale))

#Overall Reproductive number
totalR = get_R(total.inc, si_mean = si_fit$mu, si_sd = si_fit$sd)
totalR
plot(totalR)

##run 1000 simulation
tR_val <- sample_R(totalR, 1000)
quantile(tR_val, c(0.025, 0.975))
quantile(tR_val)
hist(tR_val, border = "grey", col = "navy",
     xlab = "Values of R",
     main = "Sample of likely R values")


plot(totalR, "lambdas", scale = length(onset) + 1)
abline(v = onset, lwd = 3, col = "grey")
abline(v = today, col = "blue", lty = 2, lwd = 2)
points(onset, seq_along(onset), pch = 1, cex = 1)

#First phase reproductive number
fR = get_R(subset(total.inc, from = min(onset), to = max(first16$onset)),
           si_mean = si_fit$mu, si_sd = si_fit$sd)
fR
fR_val <- sample_R(fR, 1000)
quantile(fR_val, c(0.025, 0.975))
quantile(fR_val)

plot(fR, "lambdas", scale = length(first16$onset) + 1)
abline(v = first16$onset, lwd = 3, col = "grey")
#abline(v = today, col = "blue", lty = 2, lwd = 2)
points(first16$onset, seq_along(1:length(first16$onset)), pch = 20, cex = 3)

# Second phase reproductive number
sR = get_R(subset(total.inc, from = max(first16$onset), to = max(onset)),
           si_mean = si_fit$mu, si_sd = si_fit$sd)
sR_val <- sample_R(sR, 1000)
quantile(sR_val, c(0.025, 0.975))
quantile(sR_val)

res_parametric_si <- estimate_R(imported.inc, 
                                method="parametric_si",
                                config = make_config(list(
                                  mean_si = si_fit$mu, 
                                  std_si = si_fit$sd))
)


plot(sR, "lambdas", scale = length(vietnam[17:nrow(vietnam),]$date) + 1)
abline(v = vietnam[17:nrow(vietnam),]$date, lwd = 3, col = "grey")
abline(v = today, col = "blue", lty = 2, lwd = 2)
points(vietnam[17:123,]$date, seq_along(1:nrow(vietnam[17:nrow(vietnam),])), pch = 20, cex = 3)

plot(res_parametric_si, add_imported_cases=TRUE)
#projection
si <- distcrete::distcrete("gamma", shape = params$shape,
                scale = params$scale,
                interval = 1, w = 0)

pred <- project(total.inc, totalR$R_ml, si = si, n_days = 10, n_sim = 1000)

plot(pred)
df <- as.data.frame(pred, long = TRUE)
head(df)

p <- ggplot(df, aes(x = date, y = incidence)) +
  geom_jitter(alpha = .3) 
p


plot(total.inc) %>% add_projections(pred, boxplots = TRUE)
pred_cum



