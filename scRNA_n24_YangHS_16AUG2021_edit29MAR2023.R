## DLPFC cellular communities (n=24): "The proportion of certain DLPFC cell subsets partially mediate the association between tau pathology and cognitive decline"
# Analysis finalized 8/16/2021. Code file summarized 3/29/2023

library(mediation)

# Analysis done with CelMod-inferred cell subtype proportions, adj. age, sex, RIN

# Dataset: d0, with 12/10/2020 phenotype data download from RADC.
# NOTE: Cell subset labels were re-arranged after these analysis. This was just re-numbering (re-labeling) cell subtypes without change in celmod proportions.
# Ast.3 --> Ast.4; End.1 --> End.2; Inh.SST --> Inh.3; Inh.1-->Inh.2; oli.1 and oli.2 annotation did not change.

# Please note that the following code used pre-update cell subset labels.

# Tau-associated cell types (FDR<0.017) that are also marginally associated with Ab (nominal p<0.05): 
# Ast.3, End.1, Inh.SST, Oli.1, Oli.2, and Inh.1

# Testing Amyloid (A) and Tau (T) in the same model: to assess 
# (1) whether A is a confounder (i.e., cell subtype no longer associated with T when A in the model), or
# (2) T is a mediator of A-cell association (T associated with cell type even after A adjustment)

celltypes<-c("Ast.3", "End.1", "Inh.SST", "oli.1", "oli.2", "Inh.1")
A<-matrix(data=NA,nrow=length(celltypes),ncol=5)
for (i in 1:length(celltypes)){
  A[i,1]<-celltypes[i]
  fmla1<-paste0(celltypes[i],"~tangles_sqrt+amyloid_sqrt+age_death+msex+RINcontinuous")
  tryCatch({
    m<-summary(lm(as.formula(fmla1),data=subset(d0)))
    A[i,2]<-coef(m)[2,1]
    A[i,3]<-coef(m)[2,2]
    A[i,4]<-coef(m)[2,3]
    A[i,5]<-coef(m)[2,4]
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
A<-as.data.frame(A)
names(A)<-c("cell_type","beta","se","t","p")
A$p<-as.numeric(as.character(A$p))
View(A)
#   cell_type                 beta                  se                 t            p
# 1     Ast.3   0.0183957245440306 0.00334573349976446   5.4982635482849 5.594518e-08
# 2     End.1   0.0110781183944231 0.00233053844454287  4.75345876415956 2.483994e-06
# 3   Inh.SST -0.00633674028286111 0.00148216219335685 -4.27533525768152 2.206804e-05
# 4     oli.1  -0.0174772871239428 0.00434367555923812   -4.023617069367 6.433344e-05
# 5     oli.2   0.0110410271043773 0.00532154998223294  2.07477654841917 3.841589e-02
# 6     Inh.1  0.00148770383137336 0.00206626562015525   0.7199964113334 4.717964e-01
# Inh.1 is mainly associated with A. No direct association with T. Difficult to resolve causal direction.
# We analyzed the rest in mediation models. (Results in Extended Data Fig. 8c, With updated cell subtype labels)

## A-->T-->Cell (Fig.7e, Ext Fig.8c)
celltypes<-c("Ast.3","End.1","Inh.SST","oli.1","oli.2")

# Note: ACME, average causal mediated effects. ADE, average direct effects. treat=independent variable.  
set.seed(1005)
m1<-lm(tangles_sqrt~amyloid_sqrt+age_death+msex+RINcontinuous,data=subset(d0,!is.na(oli.2)));summary(m1)
m2<-lm(Ast.3~amyloid_sqrt+tangles_sqrt+age_death+msex+RINcontinuous,data=subset(d0));summary(m2)
m_med_ast3<-mediate(m1,m2,sims=10000,boot=TRUE,treat="amyloid_sqrt",mediator="tangles_sqrt")
summary(m_med_ast3)
# Causal Mediation Analysis 
# Nonparametric Bootstrap Confidence Intervals with the Percentile Method
#                 Estimate 95% CI Lower 95% CI Upper p-value    
# ACME            0.008938     0.005140         0.01  <2e-16 ***
# ADE            -0.001888    -0.008727         0.00   0.598    
# Total Effect    0.007050     0.000601         0.01   0.032 *  
# Prop. Mediated  1.267722     0.471464         6.41   0.032 * 
# Sample Size Used: 631 
# Simulations: 10000 
# With 10000 simulations, min p-value that can be quantified is 1e-4 (once in 10000 simulations)

set.seed(1003)
m1<-lm(tangles_sqrt~amyloid_sqrt+age_death+msex+RINcontinuous,data=subset(d0,!is.na(oli.2)));summary(m1)
m2<-lm(End.1~amyloid_sqrt+tangles_sqrt+age_death+msex+RINcontinuous,data=subset(d0));summary(m2)
m_med_end1<-mediate(m1,m2,sims=10000,boot=TRUE,treat="amyloid_sqrt",mediator="tangles_sqrt")
summary(m_med_end1)
# Causal Mediation Analysis 
# Nonparametric Bootstrap Confidence Intervals with the Percentile Method
#                 Estimate 95% CI Lower 95% CI Upper p-value    
# ACME            0.005382     0.002948         0.01  <2e-16 ***
# ADE             0.000403    -0.004248         0.01  0.8720    
# Total Effect    0.005786     0.001399         0.01  0.0084 ** 
# Prop. Mediated  0.930324     0.422575         3.28  0.0084 ** 
# Sample Size Used: 631 
# Simulations: 10000 

m_med_end1$d0.p #to get the exact numeric value of ACME p-value
# [1] 0

set.seed(1004)
m1<-lm(tangles_sqrt~amyloid_sqrt+age_death+msex+RINcontinuous,data=subset(d0,!is.na(oli.2)));summary(m1)
m2<-lm(Inh.SST~amyloid_sqrt+tangles_sqrt+age_death+msex+RINcontinuous,data=subset(d0));summary(m2)
m_med_inhsst<-mediate(m1,m2,sims=10000,boot=TRUE,treat="amyloid_sqrt",mediator="tangles_sqrt")
summary(m_med_inhsst)
# Causal Mediation Analysis 
# Nonparametric Bootstrap Confidence Intervals with the Percentile Method
#                 Estimate 95% CI Lower 95% CI Upper p-value    
# ACME           -0.003079    -0.004621         0.00  <2e-16 ***
# ADE            -0.000339    -0.003563         0.00   0.825    
# Total Effect   -0.003417    -0.006309         0.00   0.018 *  
# Prop. Mediated  0.900934     0.347166         4.34   0.018 *  
# Sample Size Used: 631 
# Simulations: 10000 
m_med_inhsst$z0.p #to get the exact numeric value of ADE p-value
# [1] 0.8246

set.seed(1002)
m1<-lm(tangles_sqrt~amyloid_sqrt+age_death+msex+RINcontinuous,data=subset(d0,!is.na(oli.2)));summary(m1)
m2<-lm(oli.1~amyloid_sqrt+tangles_sqrt+age_death+msex+RINcontinuous,data=subset(d0));summary(m2)
m_med_oli1<-mediate(m1,m2,sims=10000,boot=TRUE,treat="amyloid_sqrt",mediator="tangles_sqrt")
summary(m_med_oli1)
# Causal Mediation Analysis 
# Nonparametric Bootstrap Confidence Intervals with the Percentile Method
#                 Estimate 95% CI Lower 95% CI Upper p-value    
# ACME           -0.00849     -0.01304         0.00  0.0002 ***
# ADE            -0.00384     -0.01274         0.01  0.3880    
# Total Effect   -0.01233     -0.02066         0.00  0.0030 ** 
# Prop. Mediated  0.68886      0.30501         2.09  0.0032 ** 
# Sample Size Used: 631 
# Simulations: 10000 
m_med_oli1$d0.p
# [1] 2e-04

set.seed(1001)
m1<-lm(tangles_sqrt~amyloid_sqrt+age_death+msex+RINcontinuous,data=subset(d0,!is.na(oli.2)));summary(m1)
m2<-lm(oli.2~amyloid_sqrt+tangles_sqrt+age_death+msex+RINcontinuous,data=subset(d0));summary(m2)
m_med_oli2<-mediate(m1,m2,sims=10000,boot=TRUE,treat="amyloid_sqrt",mediator="tangles_sqrt")
summary(m_med_oli2)
# Causal Mediation Analysis 
# Nonparametric Bootstrap Confidence Intervals with the Percentile Method
#                 Estimate 95% CI Lower 95% CI Upper p-value    
# ACME           5.36e-03     3.07e-04         0.01  0.0376 *  
# ADE            1.05e-02     1.14e-05         0.02  0.0496 *  
# Total Effect   1.59e-02     5.85e-03         0.03  0.0006 ***
# Prop. Mediated 3.37e-01     2.11e-02         1.00  0.0382 *  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Sample Size Used: 631 
# Simulations: 10000
m_med_oli2$tau.p
# [1] 6e-04

# Ast.3, End.1, Inh.SST, Oli.1: No direct A effect after adj.T. >50% of A-Cell type association mediated by T. 
# Oli.2: Mediated proportion 34%, direct A - Cell subtype association remained. Likely affected by both A and T. 
# We decided to focus on the first 4 cell subtypes with strong downstream association with tau 
#    to assess if these cell subtypes mediate tau-->cognitive decline. (Fig.7f, Ext Fig.8d)

## T-->Cell-->CogDec (Fig.7f, Ext Fig.8d)
# Note: Here, negative association with cognitive decline means more decline (as the variable is coded with the slope).
#       We reversed the sign in the manuscript to reduce confusion on the direction of effect. 

set.seed(2001)
m1<-lm(Ast.3~tangles_sqrt+age_death+msex+RINcontinuous,data=subset(d0,!is.na(cogng_demog_slope)))
m2<-lm(cogng_demog_slope~Ast.3+tangles_sqrt+age_death+msex+RINcontinuous,data=subset(d0))
m_med1_ast3<-mediate(m1,m2,sims=10000,boot=TRUE,treat="tangles_sqrt",mediator="Ast.3")
summary(m_med1_ast3)
# Causal Mediation Analysis 
# Nonparametric Bootstrap Confidence Intervals with the Percentile Method
#                Estimate 95% CI Lower 95% CI Upper p-value    
# ACME           -0.00279     -0.00469         0.00  <2e-16 ***
# ADE            -0.03716     -0.04236        -0.03  <2e-16 ***
# Total Effect   -0.03995     -0.04529        -0.03  <2e-16 ***
# Prop. Mediated  0.06977      0.03289         0.12  <2e-16 ***
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Sample Size Used: 598 
# Simulations: 10000
# Note: n=598 due to additional missingness in the cognitive decline phenotype

set.seed(2002)
m1<-lm(End.1~tangles_sqrt+age_death+msex+RINcontinuous,data=subset(d0,!is.na(cogng_demog_slope)))
m2<-lm(cogng_demog_slope~End.1+tangles_sqrt+age_death+msex+RINcontinuous,data=subset(d0))
m_med1_end1<-mediate(m1,m2,sims=10000,boot=TRUE,treat="tangles_sqrt",mediator="End.1")
summary(m_med1_end1)
# Causal Mediation Analysis 
# Nonparametric Bootstrap Confidence Intervals with the Percentile Method
#                Estimate 95% CI Lower 95% CI Upper p-value    
# ACME           -0.00249     -0.00422         0.00  <2e-16 ***
# ADE            -0.03745     -0.04252        -0.03  <2e-16 ***
# Total Effect   -0.03995     -0.04518        -0.03  <2e-16 ***
# Prop. Mediated  0.06244      0.02864         0.10  <2e-16 ***
# Sample Size Used: 598 
# Simulations: 10000

set.seed(2003)
m1<-lm(Inh.SST~tangles_sqrt+age_death+msex+RINcontinuous,data=subset(d0,!is.na(cogng_demog_slope)))
m2<-lm(cogng_demog_slope~Inh.SST+tangles_sqrt+age_death+msex+RINcontinuous,data=subset(d0))
m_med1_inhsst<-mediate(m1,m2,sims=10000,boot=TRUE,treat="tangles_sqrt",mediator="Inh.SST")
summary(m_med1_inhsst);plot(m_med1_inhsst)
# Causal Mediation Analysis 
# Nonparametric Bootstrap Confidence Intervals with the Percentile Method
#                Estimate 95% CI Lower 95% CI Upper p-value  
# ACME           -0.00141     -0.00261         0.00   8e-04 ***
# ADE            -0.03854     -0.04381        -0.03  <2e-16 ***
# Total Effect   -0.03995     -0.04522        -0.03  <2e-16 ***
# Prop. Mediated  0.03529      0.01195         0.07   8e-04 ***
# Sample Size Used: 598 
# Simulations: 10000

set.seed(2004)
m1<-lm(oli.1~tangles_sqrt+age_death+msex+RINcontinuous,data=subset(d0,!is.na(cogng_demog_slope)))
m2<-lm(cogng_demog_slope~oli.1+tangles_sqrt+age_death+msex+RINcontinuous,data=subset(d0))
m_med1_oli1<-mediate(m1,m2,sims=10000,boot=TRUE,treat="tangles_sqrt",mediator="oli.1")
summary(m_med1_oli1);plot(m_med1_oli1)
# Causal Mediation Analysis 
# Nonparametric Bootstrap Confidence Intervals with the Percentile Method
#                Estimate 95% CI Lower 95% CI Upper p-value  
# ACME           -0.00108     -0.00220         0.00   0.014 *  
# ADE            -0.03887     -0.04413        -0.03  <2e-16 ***
# Total Effect   -0.03995     -0.04515        -0.03  <2e-16 ***
# Prop. Mediated  0.02705      0.00456         0.06   0.014 *
# Sample Size Used: 598 
# Simulations: 10000

## Proportion of tau-cogdec association explained by 4 cell subtypes downstream of tau:
m<-lm(cogng_demog_slope~tangles_sqrt+age_death+msex+RINcontinuous, data=subset(d0,!is.na(Ast.3)));summary(m)
# tangles_sqrt  -0.0399487  0.0026119 -15.295  < 2e-16 ***
# Multiple R-squared:  0.2997,	Adjusted R-squared:  0.295 
m<-lm(cogng_demog_slope~Ast.3+End.1+Inh.SST+oli.1+tangles_sqrt+age_death+msex+RINcontinuous, data=subset(d0,!is.na(Ast.3)));summary(m)
# tangles_sqrt  -0.0370065  0.0026182 -14.134  < 2e-16 ***
# Multiple R-squared:  0.3359,	Adjusted R-squared:  0.3269
1-0.0370065/0.0399487 # compare the effect size (beta)
# [1] 0.07364946



