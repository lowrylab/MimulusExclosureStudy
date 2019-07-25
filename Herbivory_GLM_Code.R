# First we load up the lme4 package for coding linear, generalized linear, and nonlinear models.

library(lme4)

# Now we specify our model for the relevant response variables (herbivory - Herb_5) at the COAST site.
# Since our herbivory data is essentially presence/absence data, we must fit the model around a binomial distribution. Because we're only interested in elucidating whether the treatment had an effect, we will be running a Generalized Linear Model without random effects.

glmcoastherb1 <- glm(Herb_5 ~ Treatment + Eco + Treatment:Eco, family= "binomial", data = glmcoastherb)

# Next, we test our model for significance.

library(car)
Anova(glmcoastherb1,type="III")

# Finally, we want to calculate the total number of herbivore damaged replicates at the end of the experiment at the coast site for: Coastal Exclosure, Coastal Control, Inland Exclosure, Inland Control.

coastherbex<-glmcoastherb[which(glmcoastherb$Treatment=="Exclosure"),]
coastherbcon<-glmcoastherb[which(glmcoastherb$Treatment=="Control"),]

# Then we want to subset those subsets further by each ecotype.

coastherbcoastex<-coastherbex[which(coastherbex$Eco=="C"),]
coastherbinlandex<-coastherbex[which(coastherbex$Eco=="I"),]
coastherbcoastcon<-coastherbcon[which(coastherbcon$Eco=="C"),]
coastherbinlandcon<-coastherbcon[which(coastherbcon$Eco=="I"),]

# We use the subsets above to calculate total survival (at the Coast) in Coastal Exclosure, Coastal Control, Inland Exclosure, and Inland Control (in that order).

sum(coastherbcoastex$Herb_5)
sum(coastherbcoastcon$Herb_5)
sum(coastherbinlandex$Herb_5)
sum(coastherbinlandcon$Herb_5)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Now we specify our model for the relevant response variables (herbivory - Herb_5) at the INLAND site.
# Since our herbivory data is essentially presence/absence data, we must fit the model around a binomial distribution. Because we're only interested in elucidating whether the treatment had an effect, we will be running a Generalized Linear Model without random effects.

glminlandherb1 <- glm(Herb_5 ~ Treatment + Eco + Treatment:Eco, family= "binomial", data = glminlandherb)

# Next, we test our model for significance.

library(car)
Anova(glminlandherb1,type="III")

# Finally, we want to calculate the total number of herbivore damaged replicates at the end of the experiment at the coast site for: Coastal Exclosure, Coastal Control, Inland Exclosure, Inland Control.

inlandherbex<-glminlandherb[which(glminlandherb$Treatment=="Exclosure"),]
inlandherbcon<-glminlandherb[which(glminlandherb$Treatment=="Control"),]

# Then we want to subset those subsets further by each ecotype.

inlandherbcoastex<-inlandherbex[which(inlandherbex$Eco=="C"),]
inlandherbinlandex<-inlandherbex[which(inlandherbex$Eco=="I"),]
inlandherbcoastcon<-inlandherbcon[which(inlandherbcon$Eco=="C"),]
inlandherbinlandcon<-inlandherbcon[which(inlandherbcon$Eco=="I"),]

# We use the subsets above to calculate total survival (at the Coast) in Coastal Exclosure, Coastal Control, Inland Exclosure, and Inland Control (in that order).

sum(inlandherbcoastex$Herb_5)
sum(inlandherbcoastcon$Herb_5)
sum(inlandherbinlandex$Herb_5)
sum(inlandherbinlandcon$Herb_5)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# (Disregard) Next we need to make a bar plot of total survival at BOTH sites using the following ggplot code. First we'll build an appropriate dataframe that includes all the necessary parameters calculated on both ecotypes at both sites under both treatments (percent herbivorized, standard error). Percent herbivorized doesn't need a standard error component.

library(ggplot2)
plot2 <- data.frame(mean=c(0,28,0,29,0,7,0.7,0.7),
                    se=c(0,0,0,0,0,0,0,0),
                    ecotype=c("Coastal","Coastal","Inland","Inland","Coastal","Coastal","Inland","Inland"),
                    treatment=c("Ex","Con","Ex","Con","Ex","Con","Ex","Con"),
                    site=c("Bodega Bay (Coast Site)","Bodega Bay (Coast Site)","Bodega Bay (Coast Site)","Bodega Bay (Coast Site)","Pepperwood (Inland Site)","Pepperwood (Inland Site)","Pepperwood (Inland Site)","Pepperwood (Inland Site)"))

ggplot(plot2, aes(x=ecotype, y=mean, fill=treatment))+
  geom_bar(stat="identity", position="dodge", width=0.7)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(.9), width=0)+
  scale_fill_manual(values=c("#CA3542","#37AFA9"), labels=c("Control","Exclusion"), name="Treatment")+
  labs(x="Ecotype", y="Percent Herbivorized (%)")+
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
ggsave("plot50.eps", height = 8, width = 5)