devtools::document()
devtools::install()
library(BJM)

#####
##### With competing risks
#####

data(pbc3)

############################################################
# CMT plot
############################################################
data.plot.all = pbc3[!is.na(pbc3$status4),]
pbc.cmt.cr <- cmtPlot(data.plot.all, condi_time2event = 5, 
                      event_type_variable = 'status4', 
                      event_type = c("0", "1"),
                      bio_variable = "albumin", 
                      time_variable = "year", 
                      survival_variable = "years", 
                      interval_time = 1/4,
                      id_variable = "id")

############################################################
# Survival model
############################################################

## extract unique data
data.survival.fitting =  data.raw.sim.1[!duplicated(data.raw.sim.1$id), ]
## survival model structure
formMarginalSurv = Surv(years, status3) ~ age + sex
## status4 is the competing risk event type indicator can be 0 or 1
formConditionalCR = status4 ~ years + age + sex

## fitting survival submodel
survival_fit_all = survivalSub(data.survival.fitting, formMarginalSurv, formConditionalCR)

############################################################
# Longitudinal model
############################################################

## longitudinal model structure
LongSubFixed = list(
  "long1" = serBilir ~ year + age + sex +  (years) + (years) * year + 
    status4 + status4 * year  + status4 * year * years, 
  "long2" = prothrombin ~ year + age + sex + (years) + (years) * year + 
    status4 + status4 * year  + status4 * year * years, 
  "long3" = albumin ~ year + age + age * year + sex + (years) + (years) * year + 
    status4 + status4 * year  + status4 * year * years, 
  "long4" = alkaline ~ year + age + sex + (years) + (years) * year + 
    status4 + status4 * year  + status4 * year * years,
  "long5" = SGOT ~ year + age + sex + (years) + (years) * year + 
    status4 + status4 * year  + status4 * year * years,
  "long6" = platelets ~ year + age + sex + (years) + (years) * year + 
    status4 + status4 * year  + status4 * year * years)

LongSubRandom =list(
  "long1" =  ~ year| id,   
  "long2" =  ~ year| id,    
  "long3" =  ~ year| id,    
  "long4" =  ~ year| id,    
  "long5" =  ~ year| id,    
  "long6" =  ~ year| id)

survivalVariableAll = list(
  "Tyears1",  "Tyears2", "Tyears3", "Tyears4"
)

survivalTransFunction = list(
  fun1 = function(x){log(x + 1)}, 
  fun2 = function(x){abs(x - 1)}, 
  fun3 = function(x){abs(x - 3)}, 
  fun4 = function(x){abs(x - 5)}
)

data.fit.all = list()
for(i in 1:length(LongSubFixed)){
  data.fit.all[[i]] = pbc3[pbc3$status3 == 1, ]
}
## fitting longitudinal submodel
long_fit_all = longitudinalSub(data.fit.all, LongSubFixed, LongSubRandom)

############################################################
# Risk Dynamic prediction
############################################################
i_PID = 2
data.raw.predict.1 = pbc3[pbc3$id == i_PID, ]

# prediction time: 3; horizon: 3
data.predict.all = list()
for(i in 1:length(LongSubFixed)){
  data.predict.all[[i]] = data.raw.predict.1[data.raw.predict.1$year <= 3,]
}

# predict risk probability
risk.prob = dynamicPrediction(data.predict.all, long_fit_all, survival_fit_all, 
                              prediction.time = 3, 
                              horizon = 3, time_variable = "year",
                              survivalVariableAll, survivalTransFunction,
                              bandcount1 = 10, bandcount2 = 10)

############################################################
# Biomarker Dynamic prediction
############################################################
Y_predict = dynamicPredictionBio(bio_i = 1, data.predict.all, long_fit_all, 
                                 survival_fit_all, prediction.time = time.cutoff, 
                                 horizon = prediction.horizon, time_variable = "year",
                                 survivalVariableAll, survivalTransFunction,
                                 bandcount2 = 40, bandcount3 = 400)

Y_predict[[1]]
############################################################
# Risk prediction plot
############################################################

i_PID = 2
data.raw.predict.plot = pbc3[pbc3$id == i_PID, ]
data.predict.all.pre = list(data.raw.predict.plot, data.raw.predict.plot, data.raw.predict.plot,
                            data.raw.predict.plot, data.raw.predict.plot, data.raw.predict.plot)

# plot biomarker 1 history,  predict future biomarker
predictPlot(data.predict.all.pre, long_fit_all, survival_fit_all, 
            prediction.time = 5, bio_his = 1, bio_pred = 1,
            horizon = seq(0.0, 3.0, 0.5), time_variable = "year",
            survivalVariableAll, survivalTransFunction,
            bandcount1 = 10, bandcount2 = 10, bandcount3 = 200)
