The `ModelPipeline` class provides a convenient mechanism to
chain together several Selector models with a Portfolio Optimization model
(*i.e.*, Risk-based, Naïve, or Greedy) to build a complex multi-stage portfolio
weights model.

Its constructor takes as a parameter a list of models,
`[S_1, ..., S_N, Opt]`, where `S_i` for `i=1, ...,N`
is a Selector Model objects and `Opt` is an Optimizer Model object.
It is imperative that the last, and only the last, element of
the list is an Optimizer Model object (e.g., instances of `az.CVaRAnalyzer`,
`az.InvVol`, `az.KellyEngine`, etc. classes).

Note the following exceptions:
* A single element list will contain only an optimizer and no selectors,
*i.e.*, `[Opt]`,
* A `None` element in the list will be ignored, therefore the following
expressions are equivalent: `[Opt]`, `[None, Opt]`,
`[NullSelector(), Opt]`. Although, in the last sequence the
`NullSelector` instance will be executed,
* In general, `Opt` is a valid instance of a Risk-based, Naïve, or
Greedy portfolio weights classes. An exception is made for Equal Weighted
Portfolio where the string `"EWP"` can be passed as a shortcut,
*e.g.*, `[S_1, ..., S_N, "EWP"]`.

Once constructed, the `ModelPipeline` object can be interrogated for the
optimal weights or it can be passed to
a `Port_Generator` object for backtesting (out-of-sample testing).
