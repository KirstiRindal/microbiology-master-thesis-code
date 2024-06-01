#Power analysis 

#https://cran.r-project.org/web/packages/pwr/vignettes/pwr-vignette.html

#nstall.packages("pwr")
library(pwr)

#1: pwr.t2n.test = two-sample t-tests (unequal sample sizes)

#parameters for power analysis
effect_size <- 0.20  #effect size (f2) as medium (0.20), 
alpha <- 0.05
power <- 0.80

#sample sizes
n_HC <- 36
n_IBSM <- 28
n_IBSD <- 20
n_IBSC <- 8
n_IBS <- 56

#power analysis for each comparison

#healthy controls vs IBS-M
power_mixed <- pwr.t2n.test(n1 = n_HC, n2 = n_IBSM, d = effect_size, sig.level = alpha, power = NULL)
print(paste('Power for Healthy vs Mixed IBS:', round(power_mixed$power, 4)))

sample_size_IBSM <- pwr.t.test(d = effect_size, sig.level = alpha, power = power, type = "two.sample")$n
print(paste('Required sample size per group for Healthy vs Mixed IBS:', ceiling(sample_size_IBSM))) #64

#healthy controls vs IBS-D
power_diarrhea <- pwr.t2n.test(n1 = n_HC, n2 = n_IBSD, d = effect_size, sig.level = alpha, power = NULL)
print(paste('Power for Healthy vs Diarrhea-predominant IBS:', round(power_diarrhea$power, 4)))

sample_size_IBSD <- pwr.t.test(d = effect_size, sig.level = alpha, power = power, type = "two.sample")$n
print(paste('Required sample size per group for Healthy vs Mixed IBS:', ceiling(sample_size_IBSD))) #64

#healthy controls vs IBS-C
power_constipation <- pwr.t2n.test(n1 = n_HC, n2 = n_IBSC, d = effect_size, sig.level = alpha, power = NULL)
print(paste('Power for Healthy vs Constipation-predominant IBS:', round(power_constipation$power, 4)))

sample_size_IBSC <- pwr.t.test(d = effect_size, sig.level = alpha, power = power, type = "two.sample")$n
print(paste('Required sample size per group for Healthy vs Mixed IBS:', ceiling(sample_size_IBSC))) #64

#total power for all comparisons
total_power <- pwr.t2n.test(n1 = n_HC, n2 = n_IBSM + n_IBSD + n_IBSC, d = effect_size, sig.level = alpha, power = NULL)
print(paste('Total power for all comparisons:', round(total_power$power, 4)))

#2: pwr.f2.test = test for the general linear model.

pwr.f2.test(u = 4, v = NULL, f = 0.2, sig.level = 0.05, power = 0.8)
#v = degrees of freedom for denominator = 60 (round up). 

#what is the power of the test with 60 subjects? 
pwr.f2.test(u = 4, v = 8, f2 = 0.2, sig.level = 0.05)
#power = 0.8029916 = 80.3% 

#u=	degrees of freedom for numerator. The number of coefficients you'll have in your model (minus the intercept).
#v=	degrees of freedom for denominator. The number of error degrees of freedom: v=n-u-1, which implies n=v+u+1.
#f2=effect size. The effect size you want to be able to detect.
#power = power of the test (1 minus Type II error probability)



