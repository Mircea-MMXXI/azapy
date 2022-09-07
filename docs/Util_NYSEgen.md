
(Util_NYSEgen_TOP)=
# NYSEgen

## Generate the NYSE (New York Stock Exchange) business days calendar.

### Call:

```
NYSEgen(sdate='1980-01-01', edate='2050-12-31')
```

### Inputs:

* `sdate` : calendar start date. If it is not a NYSE business day it will
default to the next business day.
The default is `'1980-01-01'`.
* `edate` : calendar end date. If it is not a NYSE business day it will
default to the previous business day. The default is
`'2050-12-31'`.

### Returns:
`numpy.busdaycalendar` object.

### Examples:

```
buscal = NYSEgen()
```

[TOP](Util_NYSEgen_TOP)
