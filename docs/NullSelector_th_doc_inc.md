This is the identity Selector. It always returns the initial market data
and the invested capital fraction set to `1`. In other words, it produces
no selection (the output is same as the input). Its only purpose is to
serve as a base class for other Selectors classes as well as to provide
consistency across the family of Selectors.

In simple applications, it can be omitted.
