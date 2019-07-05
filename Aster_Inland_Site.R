library(aster)

# Each of the numbered nodes labeled Surv represent week-to-week presence/absence survival data, while the final node labeled Mass includes continuous dry biomass in grams.

vars <- c("Surv_1","Surv_2","Surv_3","Surv_4","Surv_5","Surv_6","Surv_7","Mass")
redata <- reshape(asterinland, varying = list(vars), direction = "long", timevar = "varb", times = as.factor(vars), v.names = "resp")
redata <- data.frame(redata, root=1)
massnozero <- subset(asterinland,Mass>0)
sd(massnozero$Mass)
famlist <- list(fam.bernoulli(),fam.normal.location(0.2483608))
pred <-c(0,1,2,3,4,5,6,7)
fam <- c(1,1,1,1,1,1,1,2)

# Used the following code (an alteration from one used in the Aster Tutorial) to verify that my families were appropriately assigned.

sapply(famlist,as.character)[fam]

# Now it’s time to fit the model one factor at a time and run a likelihood ratio test between each subsequent model.

aout1 <- aster(formula = resp ~ varb + Eco, pred = pred, fam = fam, varvar = varb, idvar = id, root = root, data = redata, famlist = famlist)
aout2 <- aster(formula = resp ~ varb + Eco + Treatment, pred = pred, fam = fam, varvar = varb, idvar = id, root = root, data = redata, famlist = famlist)
anova(aout1,aout2)
aout3 <- aster(formula = resp ~ varb + Eco + Treatment + Eco*Treatment, pred = pred, fam = fam, varvar = varb, idvar = id, root = root, data = redata, famlist = famlist)
anova(aout2,aout3)

# Now we want to extract the weighted biomass means, standard error, and confidence intervals from the model for each combination of our factors.

foo <- paste(as.character(redata$Eco), as.character(redata$Treatment))
bar <- sort(unique(foo))
print(bar)
baz <- redata$id[match(bar,foo)]
print(baz)
newdata <- redata[redata$id %in% baz,]
pout3 <- predict(aout3, newdata, varb, id, root, se.fit = TRUE)
qux <- data.frame(eco = newdata$Eco, treatment = newdata$Treatment, mean = pout3$fit, std_err = pout3$se.fit, low = pout3$fit - qnorm(0.975) * pout3$se.fit, high = pout3$fit + qnorm(0.975) * pout3$se.fit)
qux2 <- qux[as.character(newdata$varb) == "Mass",]
print(qux2)