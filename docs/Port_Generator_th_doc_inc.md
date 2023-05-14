The `Port_Generator` class is a generic backtesting (out-of-sample testing)
simulator. It takes as a parameter a `ModelPipeline` object.
In fact, classes as `Port_CVaR`, `Port_Kelly`, `Port_InvVol`, etc.
(presented is the section Risk-base, Greedy and Naïve portfolio sections)
are simple (and convenient) wrappers of `Port_Generator` class instantiated
with a `ModelPipeline` constructed from a single element list containing
an object of type `CVaRAnalyzer`, `KellyEngine`, `InvVolEngine`, etc.
In other words, the call
```
pp = az.Port_CVaR(mktdata, freq=freq, fixoffset=fixoffset,
                  histoffset=hlength, pname=pname)
port = pp.set_model(alpha=alpha, coef=coef, rtype=rtype, mu0=mu0,
                    hlength=hlength)
```
is equivalent to
```
cvar = az.CVaRAnalyzer(alpha=alpha, coef=coef, rtype=rtype, mu0=mu0,
                       freq=freq, hlength=hlength)
model = az.ModelPipeline([cvar])
pp = az.Port_Generator(mktdata, freq=freq, fixoffset=fixoffset,
                       histoffset=hlength, pname=pname)
port = pp.set_model(model)
```
More examples of equivalent calls can be found in this
[script](<https://github.com/Mircea-MMXXI/azapy/blob/main/scripts/generators/Port_Generator_equivalence_test.py>).
