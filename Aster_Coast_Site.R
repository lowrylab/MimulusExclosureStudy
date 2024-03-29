library(aster)

# Each of the numbered nodes labeled Surv represent week-to-week presence/absence survival data, while the final node labeled Mass includes continuous dry biomass in grams.

vars <- c("Surv_1","Surv_2","Surv_3","Surv_4","Surv_5","Surv_6","Surv_7","Mass")
redata <- reshape(astercoast, varying = list(vars), direction = "long", timevar = "varb", times = as.factor(vars), v.names = "resp")
redata <- data.frame(redata, root=1)
massnozero <- subset(astercoast,Mass>0)
sd(massnozero$Mass)
famlist <- list(fam.bernoulli(),fam.normal.location(3.681887))
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
aout4 <- aster(formula = resp ~ varb + Treatment, pred = pred, fam = fam, varvar = varb, idvar = id, root = root, data = redata, famlist = famlist)
anova(aout4,aout2)

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

# (FINAL) Next we need to make a bar plot of the estimated mean biomass at BOTH sites using the following ggplot code.

library(ggplot2)
plot1 <- data.frame(mean=c(1.712721, -1.256874, 3.224835, 1.526579, 0.03763703, 0.03673953, 0.18391834, 0.05572122),
                    se=c(0.2656212, 0.1921128,0.2965933, 0.2791396, 0.01104276, 0.01096556, 0.01969937, 0.01251975),
                    ecotype=c("Inland","Inland","Coastal","Coastal","Inland","Inland","Coastal","Coastal"),
                    treatment=c("Ex","Con","Ex","Con","Ex","Con","Ex","Con"),
                    site=c("Bodega Bay (coast site)","Bodega Bay (coast site)","Bodega Bay (coast site)","Bodega Bay (coast site)","Pepperwood (inland site)","Pepperwood (inland site)","Pepperwood (inland site)","Pepperwood (inland site)"))

ggplot(plot1, aes(x=ecotype, y=mean, fill=treatment))+
  geom_bar(stat="identity", position="dodge")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(.9), width=.3)+
  scale_fill_manual(values=c("#CA3542","#37AFA9"), labels=c("Control","Exclusion"), name="Treatment")+
  labs(x="Ecotype", y="Expected mean dry aboveground biomass (g)")+
  geom_hline(aes(yintercept=0), size=.3)+
  facet_grid(site~., scales="free")+
  theme_bw()+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(axis.text.x=element_text(size=16))+
  theme(axis.text.y=element_text(size=16))+
  theme(strip.text.y=element_text(size=16))+
  theme(legend.title=element_text(size=12))+
  theme(legend.text=element_text(size=12))+
  theme(legend.justification = c(0.9,1.1), legend.position=c(0.95,1), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("plot1.eps", height = 8, width = 5)

