# First we load up the lme4 package for coding linear, generalized linear, and nonlinear models.

library(lme4)

# Now we specify our models for both relevant response variables (final survival and dry mass) at the coast site.
# Since our survival data is essentially presence/absence data, we must fit the model around a binomial distribution.

glmmcoast1 <- glmer(Surv_7 ~ Treatment + Eco + Treatment:Eco + (1 | Plot) + (1 | Pop), family= "binomial", data = glmmcoast)

# More confusing yet is our model specification for Mass. There’s no straightforward way to analyze a zero inflated, non-normal, non-integer data set. However, upon further rumination, we may want to treat our zero values as missing data (since realistically, all replicates had some mass prior to senescence / disintegration), thus excising them from the analysis altogether. Because the remaining data is still heavily skewed to the right, this makes Mass an appropriate candidate to be modeled using a gamma distribution.
#Note: For the Mass analysis, we culled all of the samples where the plant survived to Surv_7, but for which we have no Mass data. These plant either disappeared before the final harvest or were lost during the harvest.

glmmcoast.adj <- glmmcoast[glmmcoast$Mass!=0,]
glmmcoast2 <- glmer(Mass ~ Treatment + Eco + Treatment:Eco + (1 |Plot) + (1 | Pop), family= "Gamma", data = glmmcoast.adj)

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

overdisp_fun(glmmcoast1)
vif_glmm(glmmcoast1)
overdisp_fun(glmmcoast2)
vif_glmm(glmmcoast2)

# A significant p-value in our test for over dispersion demonstrates that our model is over dispersed. The survival model comes out non-significant; however, the same cannot be said for the mass model. Regardless, a gamma GLM cannot be over dispersed, thus we aren’t to worry about the second model. For variance inflation, my understanding (at least given the Brudvig Lab paradigm) is that each output value needs to be around or below a value of 2.2 to confidently verify that the model isn’t stricken with variance inflation. Both models generally seem to pass this test.

# Next we want to test our models for significance. We would normally use the “summary” function to do so; however, this function is typically used to analyze non mixed models. Thus, we will be using the Anova Type = III function in the “car” package to analyze our models as proper mixed models. This is a more conservative approach.

library(car)
Anova(glmmcoast1,type="III")
Anova(glmmcoast2,type="III")

# Finally, we want to calculate total survival at the coast site for: Coastal Exclosure, Coastal Control, Inland Exclosure, Inland Control. Additionally, we will be calculated the average Biomass for each these categories as well. To begin, we will be subsetting all Exclosure and Control replicates.

coastex<-glmmcoast[which(glmmcoast$Treatment=="Exclosure"),]
coastcon<-glmmcoast[which(glmmcoast$Treatment=="Control"),]

# Then we want to subset those subsets further by each ecotype.

coastcoastex<-coastex[which(coastex$Eco=="C"),]
coastinlandex<-coastex[which(coastex$Eco=="I"),]
coastcoastcon<-coastcon[which(coastcon$Eco=="C"),]
coastinlandcon<-coastcon[which(coastcon$Eco=="I"),]

# We use the subsets above to calculate total survival (at the Coast) in Coastal Exclosure, Coastal Control, Inland Exclosure, and Inland Control (in that order).

sum(coastcoastex$Surv_7)
sum(coastcoastcon$Surv_7)
sum(coastinlandex$Surv_7)
sum(coastinlandcon$Surv_7)

# (Disregard) The same subset will be used to calculate mean biomass (at the Coast) for Coastal Exclosure, Coastal Control, Inland Exclosure, and Inland Control (in that order).

mean(coastcoastex$Mass)
mean(coastcoastcon$Mass)
mean(coastinlandex$Mass)
mean(coastinlandcon$Mass)

# (Disregard) Finally, we want to use the std.error function in the "plotrix" package to extract the standard error of our mass calculations (at the Coast) for Coastal Exclosure, Coastal Control, Inland Exclosure, and Inland Control (in that order). Because survival is discrete count data, no standard error can be calculated.

library(plotrix)
std.error(coastcoastex$Mass)
std.error(coastcoastcon$Mass)
std.error(coastinlandex$Mass)
std.error(coastinlandcon$Mass)

# (FINAL) Next we need to make a bar plot of total survival at BOTH sites using the following ggplot code. First we'll build an appropriate dataframe that includes all the necessary parameters calculated on both ecotypes at both sites under both treatments (mean/total survival, standard error). Percent survival doesn't need a standard error component.

library(ggplot2)
plot2 <- data.frame(mean=c(99,87,95,11,99,94,81,81),
                    se=c(0,0,0,0,0,0,0,0),
                   ecotype=c("Coastal","Coastal","Inland","Inland","Coastal","Coastal","Inland","Inland"),
                   treatment=c("Ex","Con","Ex","Con","Ex","Con","Ex","Con"),
                   site=c("Bodega Bay (Coast Site)","Bodega Bay (Coast Site)","Bodega Bay (Coast Site)","Bodega Bay (Coast Site)","Pepperwood (Inland Site)","Pepperwood (Inland Site)","Pepperwood (Inland Site)","Pepperwood (Inland Site)"))

ggplot(plot2, aes(x=ecotype, y=mean, fill=treatment))+
  geom_bar(stat="identity", position="dodge", width=0.7)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(.9), width=0)+
  scale_fill_manual(values=c("#CA3542","#37AFA9"), labels=c("Control","Exclusion"), name="Treatment")+
  labs(x="Ecotype", y="Percent Survival (%)")+
  geom_hline(aes(yintercept=0), size=.3)+
  facet_grid(site~., scales="free")+
  theme_bw()+
  scale_y_continuous(breaks=round(seq(0, 150, by = 25)))+
  theme(strip.text.y = element_blank())+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(axis.text.x=element_text(size=16))+
  theme(axis.text.y=element_text(size=16))+
  theme(legend.title=element_text(size=12))+
  theme(legend.text=element_text(size=12))+
  theme(legend.justification = c(0.9,1.1), legend.position=c(0.68,1), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("plot2.eps", height = 8, width = 5)

# (Disregard) Next we need to make a bar plot of mean dry biomass from BOTH sites using the following ggplot code. First we'll build an appropriate dataframe that includes all the necessary parameters calculated on both ecotypes at both sites under both treatments (mean/total survival, standard error). Survival doesn't need a standard error component.

library(ggplot2)
plot3 <- data.frame(mean=c(2.737757,0.9404047,0.945464,0.02265333,0.2051113,0.06794733,0.01537,0.010746),
                    se=c(0.4111004,0.1809725,0.1423283,0.01474667,0.03436435,0.0067164,0.00193352,0.00124645),
                    ecotype=c("Coastal","Coastal","Inland","Inland","Coastal","Coastal","Inland","Inland"),
                    treatment=c("Ex","Con","Ex","Con","Ex","Con","Ex","Con"),
                    site=c("Bodega Bay (Coast Site)","Bodega Bay (Coast Site)","Bodega Bay (Coast Site)","Bodega Bay (Coast Site)","Pepperwood (Inland Site)","Pepperwood (Inland Site)","Pepperwood (Inland Site)","Pepperwood (Inland Site)"))

ggplot(plot3, aes(x=ecotype, y=mean, fill=treatment))+
  geom_bar(stat="identity", position="dodge", width=0.7)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(.7), width=.3)+
  scale_fill_manual(values=c("#CA3542","#37AFA9"), labels=c("Control","Exclusion"), name="Treatment")+
  labs(x="Ecotype", y="Mean Dry Aboveground Biomass (g)")+
  geom_hline(aes(yintercept=0), size=.3)+
  facet_grid(site~., scales="free")+
  theme_bw()+
  theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("plot3.eps", height = 8, width = 5)

# While the above code works well for making a bar plot of mean dry aboveground biomass, it may be more informative to make box plots of the data in order to compare median values rather than the mean. For example: if we use the following code to visualize our biomass data, you can note that the coastal ecotype median is comparable across treatments at the inland site. Otherwise, if you just show the mean values in a bar plot, it looks as if the treatment has a significant effect on coastal biomass. In reality, the mean is only being driven up by a handful of outliers.

hist(inlandcoastcon$Mass, main="Inland Site: Coastal Control")
hist(inlandcoastex$Mass, main="Inland Site: Coastal Exclosure")
boxplot(inlandcoastcon$Mass, main="Inland Site: Coastal Control")
boxplot(inlandcoastex$Mass, main="Inland Site: Coastal Exclosure")

# Now we want to create a formal box plot figure in ggplot2 not unlike the bar plots above. The following is the faceted box plot code for both sites.

ggplot(boxplotdata, aes(x=Treatment2, y=Mass))+
  geom_boxplot(aes(fill=Treatment))+
  scale_fill_manual(values=c("#CA3542","#37AFA9"), labels=c("Control","Exclusion"), name="Treatment")+
  labs(x="Ecotype", y="Dry Aboveground Biomass (g)")+
  facet_grid(Site~., scales="free")+
  theme_bw()+
  theme(legend.justification = c(0.9,0.9), legend.position=c(0.9,0.9), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("plot4.eps", height = 8, width = 5)

# (Disregard) As is evident in the rudimentary plot above, the data is tremendously skewed - so much that the boxes are squished in an effort to accomodate the outliers. To remedy this, we need to log transform the mass data. Because you can't log transform zero values (of which I have many in this data set) we will be adding 0.00001 to each mass value, then subsequently log transform that column. Because my knowledge of R is limited, I have already done this in the appropriate data subset (Titled: Boxplot Data (Full)). The following is the ggplot code for demonstrating this data in a boxplot.

ggplot(boxplotdata, aes(x=Treatment2, y=Logmass))+
  geom_boxplot(aes(fill=Treatment), width=0.7)+
  scale_fill_manual(values=c("#CA3542","#37AFA9"), labels=c("Control","Exclusion"), name="Treatment")+
  labs(x="Ecotype", y="Dry Aboveground Biomass (log(mg))")+
  facet_grid(Site~., scales="free")+
  theme_bw()+
  scale_y_continuous(breaks=round(seq(-2, 4, by = 0.5)))+
  theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=14))+
  theme(axis.title.y=element_text(size=14))+
  theme(axis.text.x=element_text(size=12))+
  theme(axis.text.y=element_text(size=12))+
  theme(strip.text.y=element_text(size=14))
ggsave("plot5.eps", height = 8, width = 5)

# (FINAL) So Acer believes that because I am already including the zero mass values in my aster figure, that it would be redundant to include them in my box plot. It would also make the boxplots nicer without including these zero values (no negative scale, no log transformation, no row of zero outliers evident in figure). Thus, the below code will be used to make our new box plot with no zero mass values.

boxplotdata.adj<-boxplotdata[boxplotdata$Mass!=0,]

ggplot(boxplotdata.adj, aes(x=Treatment2, y=Logmass))+
  geom_boxplot(aes(fill=Treatment), width=0.7)+
  scale_fill_manual(values=c("#CA3542","#37AFA9"), labels=c("Control","Exclusion"), name="Treatment")+
  labs(x="Ecotype", y="Dry Aboveground Biomass (g)")+
  facet_grid(Site~.)+
  theme_bw()+
  theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(axis.text.x=element_text(size=16))+
  theme(axis.text.y=element_text(size=16))+
  theme(strip.text.y=element_text(size=16))
ggsave("plot6.eps", height = 8, width = 5)
