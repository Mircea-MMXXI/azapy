
# NYSEgen

## Generate the NYSE (New York Stock Exchange) business days calendar.

It lists the exceptional dates: '2012-10-29', '2012-10-30' and '2001-09-11' as
non NYSE business days.

### Call:

```
NYSEgen(sdate=np.datetime64('1980-01-01'),
        edate=np.datetime64('2050-12-31'))
```

### Inputs:

* `sdate` : calendar start date. If it is not a NYSE business day it will
default to the next business day. The default is `np.datetime64('1980-01-01')`.
* `edate` : calendar end date. If it is not a NYSE business day it will
default to the previous business day. The default is
`np.datetime64('2050-12-31')`.

### Returns:
`np.busdaycalendar` object.

### Examples:

```
buscal = NYSEgen()
```
