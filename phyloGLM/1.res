[1] "((a:1,c:1):1,((b:2,d:1):3,e:4):1);"
  trait01 predictor
a       1         1
b       0         0
c       1         1
d       0         0
e       0        12

Call:
phyloglm(formula = trait01 ~ predictor, data = dat, phy = tre, 
    boot = 100)
       AIC     logLik Pen.logLik 
    11.484     -2.742     -1.238 

Method: logistic_MPLE
Mean tip height: 4
Parameter estimate(s):
alpha: 0.4325445 
      bootstrap mean: 0.8929537 (on log scale, then back transformed)
      so possible upward bias.
      bootstrap 95% CI: (0.004647496,6.878534)

Coefficients:
             Estimate    StdErr   z.value lowerbootCI upperbootCI p.value
(Intercept) -0.160692  1.223282 -0.131362   -2.143641      2.3454  0.8955
predictor   -0.057793  0.210514 -0.274534   -0.286615      0.8596  0.7837

Note: Wald-type p-values for coefficients, conditional on alpha=0.4325445
      Parametric bootstrap results based on 100 fitted replicates

(Intercept)   predictor 
-0.16069247 -0.05779306 
            (Intercept)   predictor
(Intercept)   1.4964183 -0.13019187
predictor    -0.1301919  0.04431598
