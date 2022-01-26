import numpy as np
import scipy.stats as stats
from statsmodels.stats.power import NormalIndPower
from statsmodels.stats.proportion import proportion_effectsize

#ssem1: sample size estimator 1
def ssem1(p_null, p_alt, alpha, beta):
    """
    Compute the minimum number of samples needed to achieve a desired power
    level for a given effect size.
    
    Input parameters:
        p_null: base success rate under null hypothesis
        p_alt : desired success rate to be detected
        alpha : Type-I error rate
        beta  : Type-II error rate
    
    Output value:
        n : Number of samples required for each group to obtain desired power
    """
    
    # Get necessary z-scores and standard deviations (@ 1 obs per group)
    z_null = stats.norm.ppf(1 - alpha)
    z_alt  = stats.norm.ppf(beta)
    sd_null = np.sqrt(p_null * (1-p_null) + p_null * (1-p_null))
    sd_alt  = np.sqrt(p_null * (1-p_null) + p_alt  * (1-p_alt) )
    
    # Compute and return minimum sample size
    p_diff = p_alt - p_null
    n = ((z_null*sd_null - z_alt*sd_alt) / p_diff) ** 2
    return np.ceil(n)

def ssem2(p_null, p_alt, alpha, power):
    effect_size = proportion_effectsize(p_alt, p_null)
#     print(effect_size)
#     print(alpha)
#     print(power)
    n = NormalIndPower().solve_power(effect_size = effect_size, alpha = alpha, power = power,alternative = 'larger')
#     print(n)
    return np.ceil(n)


def run_sample_size_estimator():

    p_null = float(input("Baseline conversion rate? \nEnter as percentage:   "))
    p_alt = float(input("Desired (success) conversion rate? \nEnter as percentage:    "))
    alpha = float(input("What is significance level / acceptable Type I error rate (alpha)? \nEnter as percentage:       "))
    no_eval_metrics = int(input("Number of evaluation metrics used:    "))
    power = float(input("How much statistical power / Probability that we will correctly reject the null hypothesis / 100*(1 - beta) ,where beta is Type II error rate, i.e. selecting null hypothesis over alternate hypothesis \nEnter as percentage:    ")) 

    p_null = p_null/100
    p_alt = p_alt/100
    alpha = alpha/100
    power = power/100
    beta = 1 - power

    if no_eval_metrics > 1:
        print("\nBonferroni Correction applied as more than 1 evaluation metric is used.")
        alpha = alpha/no_eval_metrics

    n1 = ssem1(p_null,p_alt,alpha,beta)
    n2 = ssem2(p_null,p_alt,alpha,power)

    print(f"\nsample size estimator 1, n1 = {n1}")
    print(f"\nsample size estimator 2, n2 = {n2}")


run_sample_size_estimator()