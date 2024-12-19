devtools::document()
devtools::install()
library(BJM)

#####
##### Without competing risks
#####

data(pbc3)

############################################################
# CMT plot
############################################################
pbc.cmt <- cmtPlot(data.plot.all = pbc3, condi_time2event = 5, 
    event_type_variable = NULL, event_type = NULL,
    bio_variable = "serBilir", time_variable = "year", 
    survival_variable = "years", 
    interval_time = 1/12)

############################################################
# Survival model
############################################################
data.survival.fitting =  pbc3[!duplicated(pbc3$id), ]
## survival model structures
formMarginalSurv = Surv(years, status3) ~ age + sex
formConditionalCR = NULL
 
## fitting survival submodel
survival_fit_all = survivalSub(data.survival.fitting, formMarginalSurv, 
                                formConditionalCR)

############################################################
# Longitudinal model
############################################################

## longitudinal model structure
LongSubFixed = list(
  "long1" = serBilir ~ year + age + sex +  (years) + (years) * year,  
  "long2" = prothrombin ~ year + age + sex + (years) + (years) * year,  
  "long3" = albumin ~ year + age + age * year + sex + (years) + (years) * year,  
  "long4" = alkaline ~ year + age + sex + (years) + (years) * year, 
  "long5" = SGOT ~ year + age + sex + (years) + (years) * year, 
  "long6" = platelets ~ year + age + sex + (years)  + (years) * year)

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
  fun1 = function(x){abs(x - 1)}, 
  fun2 = function(x){abs(x - 3)}, 
  fun3 = function(x){abs(x - 5)}, 
  fun4 = function(x){abs(x - 7)}
)

data.fit.all = list()
for(i in 1:length(LongSubFixed)){
  data.fit.all[[i]] = data.raw.sim.1[data.raw.sim.1$status3 == 1, ]
}
## fitting longitudinal submodel
long_fit_all = longitdinalSub(data.fit.all, LongSubFixed, LongSubRandom)

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
