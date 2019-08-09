# measurements
Extends astropy.unit.Quantity to support uncertainties associated with quantities.
The uncertainties are propagated (following standard error propagation) when doing operations with the different measurements.

```python

import measurement
from astropy import units as u

# You can simply create a new measurement (e.g. 10 +/- 2 deg) with:
ang = measurement.Measurement(value=10.0, uncertainty=2.0, unit=u.deg)

print(ang*4)
> <Measurement 40.0 +/- 8.0 deg>
```

It supports all kind of operations that one could do with astropy.unit.Quantity (e.g. change of units).


