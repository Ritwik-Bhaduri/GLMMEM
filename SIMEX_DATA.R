data("simGFRdata")
simGFR=simGFR[is.element(simGFR$time,c(1:12)/4) & is.element(simGFR$PID,c(1:80)*100),]

fm2=nlme::lme.formula(fixed = cfb ~ time + x1:time + trt + trt:time + trt:x1:time + 0,
                      data = simGFR, random = ~time | PID,
                      correlation = nlme::corCompSymm(0.5,form = ~time | PID, fixed = TRUE),
                      control=nlme::lmeControl(returnObject=TRUE))

start = Sys.time()
(s1 = simexlme(model=fm2, model.model=simGFR[,c("cfb","PID","time","x1","trt")],
               SIMEXvariable="x1",respvar="cfb",grpvar="PID",corform="~time | PID",
               measurement.error=res.sd,measurement.error.resp=res.sd,
               lambda = c(seq(0,2,by = 0.05), seq(1,1.5,by = 0.05)),
               B = 100, fitting.method = "linear",
               jackknife.estimation = FALSE))
end = Sys.time()
end-start
100*abs(c(fixed.time,fixed.trt,fixed.leGFR,fixed.trttime,fixed.leGFRtrt)-s1$coefficients)/abs(c(fixed.time,fixed.trt,fixed.leGFR,fixed.trttime,fixed.leGFRtrt)-fm2$coefficients$fixed)
plot(s1)
?simGFR
####################################################################################
c(fixed.time,fixed.trt,fixed.leGFR,fixed.trttime,fixed.leGFRtrt)
fm2$coefficients$fixed
s1$coefficients
# # }
# 

# fm2=nlme::lme.formula(fixed = y ~ w + z,
#                       data = data_sim, random = ~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10
#                       +X11+X12+X13+X14+X15+X16+X17+X18+X19+X20
#                       +X21+X22+X23+X24+X25+X26+X27+X28+X29+X30
#                       +X31+X32+X33+X34+X35+X36+X37+X38+X39+X40
#                       +X41+X42+X43+X44+X45+X46+X47+X48+X49+X50|cluster,
#                       correlation = nlme::corCompSymm(0.5,form = ~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10
#                                                       +X11+X12+X13+X14+X15+X16+X17+X18+X19+X20
#                                                       +X21+X22+X23+X24+X25+X26+X27+X28+X29+X30
#                                                       +X31+X32+X33+X34+X35+X36+X37+X38+X39+X40
#                                                       +X41+X42+X43+X44+X45+X46+X47+X48+X49+X50|cluster, fixed = TRUE),
#                       control=nlme::lmeControl(returnObject=TRUE))
