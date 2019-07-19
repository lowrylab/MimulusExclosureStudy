# First we load up the lme4 package for coding linear, generalized linear, and nonlinear models.

library(lme4)

# Now we specify our models for both relevant response variables (final survival and dry mass) at the coast site.
# Since our survival data is essentially presence/absence data, we must fit the model around a binomial distribution.

glmminland1 <- glmer(Surv_7 ~ Treatment + Eco + Treatment:Eco + (1 | Plot) + (1 | Pop), family= "binomial", data = glmminland)

# More confusing yet is our model specification for Mass. There’s no straightforward way to analyze a zero inflated, non-normal, non-integer data set. However, upon further rumination, we may want to treat our zero values as missing data (since realistically, all replicates had some mass prior to senescence / disintegration), thus excising them from the analysis altogether. Because the remaining data is still heavily skewed to the right, this makes Mass an appropriate candidate to be modeled using a gamma distribution.
# Note: For the Mass analysis, we culled all of the samples where the plant survived to Surv_7, but for which we have no Mass data. These plant either disappeared before the final harvest or were lost during the harvest.


glmminland.adj <- glmminland[glmminland$Mass!=0,]
glmminland2 <- glmer(Mass ~ Treatment + Eco + Treatment:Eco + (1 |Plot) + (1 | Pop), family= "Gamma", data = glmminland.adj)

# Next we want to test for over dispersion and variance inflation. To do so we will be using code designed by Chad Zirbel and Christopher Warneke stored on the Brudvig Lab archives. Before testing each, we need to specify the appropriate functions.

# Code for specifying a test function for over dispersion or “overdisp_fun.”

overdisp_fun <- function(model) {
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

# Code for specifying a test function for variance inflation or “vif_glmm”.

vif_glmm <- function (fit) {
  
  v <- vcov(fit)
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}

# Now that the functions have been specified, let’s test for over dispersion and variance inflation in both models.

overdisp_fun(glmminland1)
vif_glmm(glmminland1)
overdisp_fun(glmminland2)
vif_glmm(glmminland2)

# A significant p-value in our test for over dispersion demonstrates that our model is over dispersed. The survival model comes out non-significant; however, the same cannot be said for the mass model. Regardless, a gamma GLM cannot be over dispersed, thus we aren’t to worry about the second model. For variance inflation, my understanding (at least given the Brudvig Lab paradigm) is that each output value needs to be around or below a value of 2.2 to confidently verify that the model isn’t stricken with variance inflation. Both models generally seem to pass this test.

# Next we want to test our models for significance. We would normally use the “summary” function to do so; however, this function is typically used to analyze non mixed models. Thus, we will be using the Anova Type = III function in the “car” package to analyze our models as proper mixed models. This is a more conservative approach.

library(car)
Anova(glmminland1,type="III")
Anova(glmminland2,type="III")

# To make our figures, we want to calculate total survival at the coast site for: Coastal Exclosure, Coastal Control, Inland Exclosure, Inland Control. Additionally, we will be calculated the average Biomass for each these categories as well. To begin, we will be subsetting all Exclosure and Control replicates.

inlandex<-glmminland[which(glmminland$Treatment=="Exclosure"),]
inlandcon<-glmminland[which(glmminland$Treatment=="Control"),]

# Then we want to subset those subsets further by each ecotype.

inlandcoastex<-inlandex[which(inlandex$Eco=="C"),]
inlandinlandex<-inlandex[which(inlandex$Eco=="I"),]
inlandcoastcon<-inlandcon[which(inlandcon$Eco=="C"),]
inlandinlandcon<-inlandcon[which(inlandcon$Eco=="I"),]

# We use the subsets above to calculate total survival (at the Coast) in Coastal Exclosure, Coastal Control, Inland Exclosure, and Inland Control (in that order).

sum(inlandcoastex$Surv_7)
sum(inlandcoastcon$Surv_7)
sum(inlandinlandex$Surv_7)
sum(inlandinlandcon$Surv_7)

# The same subset will be used to calculate mean biomass (at the Coast) for Coastal Exclosure, Coastal Control, Inland Exclosure, and Inland Control (in that order).

mean(inlandcoastex$Mass)
mean(inlandcoastcon$Mass)
mean(inlandinlandex$Mass)
mean(inlandinlandcon$Mass)

# Finally, we want to use the std.error function in the "plotrix" package to extract the standard error of our mass calculations (at the Coast) for Coastal Exclosure, Coastal Control, Inland Exclosure, and Inland Control (in that order). Because survival is discrete count data, no standard error can be calculated.

library(plotrix)
std.error(inlandcoastex$Mass)
std.error(inlandcoastcon$Mass)
std.error(inlandinlandex$Mass)
std.error(inlandinlandcon$Mass)

# The following will be ggplot code relevant to the graphing of my sodium data. Though it wasn't necessary, I went ahead and plugged my sodium data into a separate data frame.

saltdata <- data.frame(concentration=c(16.9278,6.4518,12.6040,6.2040,13.7008,8.2491,5.6665,5.2592,0.4503,1.1697),
                    plot=c("Coast 1","Coast 1","Coast 2","Coast 2","Coast 3","Coast 3","Inland 1","Inland 1","Na Controls","Na Controls"),
                    treatment=c("Control","Exclusion","Control","Exclusion","Control","Exclusion","Control","Exclusion","Blank","Filter"))
levels(saltdata$treatment)

# The above is just fine for graphing in a bar plot; however, the legend order will be as follows: Blank, Control, Exclusion, Filter. This is because factor levels are ordered alphabetically by R. However, this would be an unintuitive order for the legend to read. Thus, we need to go ahead and specify the order we want.

saltdata$treatment <- factor(saltdata$treatment, levels = c("Control","Exclusion","Blank","Filter"))
levels(saltdata$treatment)

# Now that the data frame is complete and our factor levels ordered, we can go ahead and plot the data in a bar plot.

ggplot(saltdata, aes(x=plot, y=concentration, fill=treatment))+
  geom_bar(stat="identity", position="dodge", width=0.7)+
  scale_fill_manual(values=c("#CA3542","#37AFA9","#BFB0B3","#49494D"), labels=c("Control","Exclusion","Blank","Filter"), name="Treatment")+
  labs(x="Plots", y="Na Concentration (ppm)")+
  geom_hline(aes(yintercept=0), size=.3)+
  theme_bw()+
  scale_y_continuous(breaks=round(seq(0, 20, by = 5)))+
  theme(axis.title.x=element_text(size=14))+
  theme(axis.title.y=element_text(size=14))+
  theme(axis.text.x=element_text(size=12))+
  theme(axis.text.y=element_text(size=12))+
  theme(strip.text.y=element_text(size=14))+
  theme(legend.title=element_text(size=12))+
  theme(legend.text=element_text(size=12))+
  theme(legend.justification = c(0.9,1.1), legend.position=c(0.95,1), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("plot7.eps", height = 5, width = 8)
