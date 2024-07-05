import azapy as az 

# New York Stock Exchange business calendar
bcal = az.calendarGen()

# Stockholm Stock Exchange (aka Nasdaq Stockholm)
bcal = az.calendarGen(name='XSTO')

# Paris Stock Exchange
bcal = az.calendarGen(name='XPAR')

# Tokyo Stock exchange
bcal = az.calendarGen(name='XTKS')

# London Stock Exchange
bcal = az.calendarGen(name='XLON')

# Singapore Stock Exchange
bcal = az.calendarGen(name='XSES')

# Shanghai Stock Exchange
bcal = az.calendarGen(name='XSHG')

# Hong Kong Stock Exchange
bcal = az.calendarGen(name='XHKG')


# The following are the valid names for exchange calendars
print(az.get_calendar_names())