import numpy as np
from astropy import units as u


class Measurement(u.Quantity):
    """A `Measurement` represents a number or `astropy.units.Quantity` with some
    associated uncertainties.

    Parameters
    ----------

    value : number, `numpy.ndarray`, `astropy.units.Qunatity` object (sequence)
        The numerical value of this measurement.

    """
    def __new__(cls, value, uncertainty=None, unit=None, dtype=None, copy=True, order=None,
                subok=False, ndmin=0):
        if isinstance(value, u.Quantity):
            if unit is not None:
                if unit is not value.unit:
                    value = value.to(unit).value
                else:
                    value = value.value
            else:
                value = value.value
                unit = value.unit
                
        obj = super().__new__(cls, value, unit, dtype, copy, order, subok, ndmin)
        if uncertainty is not None:
            if isinstance(uncertainty, u.Quantity):
                if obj.unit is not uncertainty.unit:
                    uncertainty = uncertainty.to(obj.unit).value
                else:
                    uncertainty = uncertainty.value
                
            obj._uncertainty = uncertainty
            # uncertainty = uncertainty
            # To implement. Check the Quantity.__new__ code!!!
        else:
            obj._uncertainty = 0.0
        return obj
 
    @property
    def uncertainty(self):
        """Returns the average value of the measurement.
        """
        return self._uncertainty

    # In Quantity you cannot set attributes... no value for example!! I SHOULD REMOVE THE SETTER!!!!!!!!
    @uncertainty.setter
    def uncertainty(self, uncertainty):
        if type(uncertainty) is u.Quantity:
            self._uncertainty = uncertainty.to(self.unit).value
        else:
            self._uncertainty = uncertainty

    @property
    def si(self):
        """
        Returns a copy of the current `Measurement` instance with SI units. The
        value of the resulting object will be scaled.
        """
        si_unit = self.unit.si
        return self._new_view(self.value*si_unit.scale, si_unit/si_unit.scale,
                             self.uncertainty*si_unit.scale)

    @property
    def cgs(self):
        """
        Returns a copy of the current `Measurement` instance with CGS units. The
        value of the resulting object will be scaled.
        """
        cgs_unit = self.unit.cgs
        return self._new_view(self.value*cgs_unit.scale, cgs_unit/cgs_unit.scale,
                              self.uncertainty*cgs_unit.scale)


    def _new_view(self, obj, unit=None, uncertainty=None):
        """Create a Measurement view of obj, and set the uncertainty
        """
        view = super()._new_view(obj, unit=unit)
        view.__array_finalize__(self)
        if uncertainty is not None:
            view.uncertainty = uncertainty
        return view

    def __quantity_subclass(self, unit):
        return Measurement, True

    def __array_finalize__(self, obj):
        if obj is None:
            return
        
        self._uncertainty = getattr(obj, '_uncertainty', 0.0)


    def __array_ufunc__(self, ufunc, method, *args, **kwargs):
        """New in version numpy 1.13 replacing __array_wrap__ and __array_prepare__
        """
        if hasattr(self, ufunc.__name__):
            # print(ufunc.__name__, *args, **kwargs)
            return getattr(self, ufunc.__name__)(*args, **kwargs)
        return super().__array_ufunc__(ufunc, method, *args, **kwargs)
        results = super().__array_ufunc__(ufunc, method, *args, **kwargs)
        if results is NotImplemented:
            return NotImplemented
        return results
        
    def __result_as_measurement(self, result, uncertainty, out):
        """Turn result into a measurement with the given uncertainty.

        If no output is given, it will take a view of the array as a measurement,
        and set the uncertainty.  If output is given, those should be measurement views
        of the result arrays, and the function will just set the uncertainty.

        Parameters
        ----------
        result : `~numpy.ndarray` or tuple of `~numpy.ndarray`
            Array(s) which need to be turned into quantity.
        uncertainty : `float` or `~astropy.units.Quantity` or None
            Uncertainty for the measurements to be returned (or `None` if the result
            should not be a measurement).  Should be tuple if result is a tuple.
        out : `~measurement.Measurement` or None
            Possible output measurement. Should be `None` or a tuple if result
            is a tuple.

        Returns
        -------
        out : `~measurement.Measurement`
           With uncertainties set.
        """
        if isinstance(result, tuple):
            if out is None:
                out = (None,) * len(result)
            return tuple(self._result_as_measurement(result_, uncertainty_, out_)
                         for (result_, uncertainty_, out_) in
                         zip(result, unit, out))

        if out is None:
            # View the result array as a Quantity with the proper unit.
            return result if uncertainty is None else self._new_view(result, uncertainty)

        # For given output, just set the unit. We know the unit is not None and
        # the output is of the correct Quantity subclass, as it was passed
        # through check_output.
        out.uncertainty = uncertainty
        return out
   
    # ufuncs prior to numpy 1.13 and other functions
    def __array_wrap__(self, obj, context=None):
        return super().__array_wrap__(self, obj, context)

    def __str__(self):
        return '{0} +/- {1}{2}'.format(self.value, self.uncertainty, self._unitstr)

    def __repr__(self):
        prefixstr = '<' + self.__class__.__name__ + ' '
        # arrstr = np.array2string(self.view(np.ndarray), separator=',', prefix=prefixstr)
        # return '{0}{1}{1:s}'.format(prefixstr, arrstr, self.__unitstr)
        return '{0}{1} +/- {2} {3:s}>'.format(prefixstr, self.value, self.uncertainty, self.unit)
    
    def _repr_latex_(self):
        """
        Generate latex representation of the quantity and its unit.
        This is used by the IPython notebook to show it all latexified.
        It only works for scalar quantities; for arrays, the standard
        reprensation is returned.

        Returns
        -------
        lstr
            LaTeX string
        """
        if not self.isscalar:
            raise NotImplementedError('Cannot represent Measurement arrays in LaTeX format')

        # Format value
        latex_value = "{0:g}".format(self.value)
        latex_uncertainty = "{0:g}".format(self.uncertainty)
        if "e" in latex_value:
            latex_value = latex_value.replace('e', '\\times 10^{') + '}'
        if "e" in latex_uncertainty:
            latex_uncertainty = latex_uncertainty.replace('e', '\\times 10^{') + '}'

        # Format unit
        # [1:-1] strips the '$' on either side needed for math mode
        latex_unit = (self.unit._repr_latex_()[1:-1]  # note this is unicode
                      if self.unit is not None
                      else "(Unit not initialised)")

        return '${0} \pm {1} \; {2}$'.format(latex_value, latex_uncertainty, latex_unit)

    def _decompose(self, allowscaledunits=False, bases=[]):
        """
        Generates a new `Measurement` with the units decomposed. Decomposed
        units have only irreducible units in them (see
        `astropy.units.UnitBase.decompose`).

        Parameters
        ----------
        allowscaledunits : bool
            If True, the resulting `Quantity` may have a scale factor
            associated with it.  If False, any scaling in the unit will
            be subsumed into the value of the resulting `Quantity`

        bases : sequence of UnitBase, optional
            The bases to decompose into.  When not provided,
            decomposes down to any irreducible units.  When provided,
            the decomposed result will only contain the given units.
            This will raises a `~astropy.units.UnitsError` if it's not possible
            to do so.

        Returns
        -------
        newq : `~measurement.Measurement`
            A new object equal to this quantity with units decomposed.

        """
        new_unit = self.unit.decompose(bases=bases)
        # Be careful here because self.value usually is a view of self;
        # be sure that the original value is not being modified.
        if not allowscaledunits and hasattr(new_unit, 'scale'):
            new_value = self.value*new_unit.scale
            new_uncertainty = self.uncertainty*new_unit.scale
            new_unit = new_unit/new_unit.scale
            return self._new_view(new_value, new_unit, new_uncertainty)
        else:
            return self._new_view(self.copy(), new_unit, self.uncertainty)

    def to(self, unit, equivalencies=[]):
        """
        Returns a new `~measurement.Measurement` object with the specified units.

        Parameters
        ----------
        unit : `~astropy.units.UnitBase` instance, str
            An object that represents the unit to convert to. Must be
            an `~astropy.units.UnitBase` object or a string parseable
            by the `~astropy.units` package.

        equivalencies : list of equivalence pairs, optional
            A list of equivalence pairs to try if the units are not
            directly convertible.  See :ref:`unit_equivalencies`.
            If not provided or ``[]``, class default equivalencies will be used
            (none for `~astropy.units.Quantity`, but may be set for subclasses)
            If `None`, no equivalencies will be applied at all, not even any
            set globally or within a context.
        """
        if equivalencies == []:
            equivalencies = self._equivalencies
        unit = u.Unit(unit)
        new_val = np.asarray(self.unit.to(unit, self.value, equivalencies=equivalencies))
        new_unc = np.asarray(self.unit.to(unit, self.uncertainty, equivalencies=equivalencies))
        return self._new_view(new_val, unit, new_unc)

    @property
    def si(self):
        si_unit = self.unit.si
        return self._new_view(self.value*si_unit.scale, si_unit/si_unit.scale, self.uncertainty*si_unit.scale)

    @property
    def cgs(self):
        cgs_unit = self.unit.cgs
        return self._new_view(self.value*cgs_unit.scale, cgs_unit/cgs_unit.scale, self.uncertainty*cgs_unit.scale)

    def __hash__(self):
        return hash(self.value)^hash(self.uncertainty)^hash(self.unit)

    def __iter__(self):
        if self.isscalar:
            raise TypeError("'{cls}' object with a scalar value is not iterable"
                            .format(cls=self.__class__.__name__))

        # Otherwise reutrn a generator
        def quantity_iter():
            for val, unc in zip(self.value, self.uncertainty):
                yield self._new_view(val, self.unit, unc)

        return quantity_iter()

    def __getitem__(self, key):
        if self.isscalar:
            raise TypeError("'{cls}' object with a scalar value does not support "
                "indexing".format(cls=self.__class__.__name__))

        out_value = self.value.__getitem__(key)
        out_uncertainty = self.uncertainty.__getitem__(key)
        return self._new_view(out_value, self.unit, out_uncertainty)

    # Arithmetic operations
    
    def __add__(self, other):
        if isinstance(other, Measurement):
            if other.unit is not self.unit:
                other = other.to(self.unit)
            uncertainty = np.sqrt(np.add(other.uncertainty**2, self.uncertainty**2))
            return self._new_view(np.add(self.value, other.value), self.unit, uncertainty)
        elif isinstance(other, u.Quantity):
            if other.unit is not self.unit:
                other = other.to(self.unit)
            return self._new_view(np.add(self.value, other.value), self.unit, self.uncertainty)
        else:
            if other == 0.0 or other == 0:
                return self
            raise UnitsError('Can only apply \'add\' function to dimensionless quantities when other'+\
                             ' argument is not a quantity (unless the latter is all zero/infinity/nan)')

    def __radd__(self, other):
        if isinstance(other, Measurement):
            if other.unit is not self.unit:
                other = other.to(self.unit)
            uncertainty = np.sqrt(np.add(other.uncertainty**2, self.uncertainty**2))
            return self._new_view(np.add(other.value, self.value), self.unit, uncertainty)
        elif isinstance(other, u.Quantity):
            if other.unit is not self.unit:
                other = other.to(self.unit)
            return self._new_view(np.add(other.value, self.value), self.unit, self.uncertainty)
        else:
            if other == 0.0 or other == 0:
                return self
            raise UnitsError('Can only apply \'add\' function to dimensionless quantities when other'+\
                             ' argument is not a quantity (unless the latter is all zero/infinity/nan)')
      
    def __iadd__(self, other):
        if isinstance(other, Measurement):
            if other.unit is not self.unit:
                other = other.to(self.unit)
            self._uncertainty = np.sqrt(np.add(other.uncertainty**2, self.uncertainty**2))
            self._value += other.value
            return self
        elif isinstance(other, u.Quantity):
            if other.unit is not self.unit:
                self._value += other.to(self.unit).value
            else:
                self._value += other.value
            return self
        else:
            if other == 0.0 or other == 0:
                return self
            raise UnitsError('Can only apply \'add\' function to dimensionless quantities when other'+\
                             ' argument is not a quantity (unless the latter is all zero/infinity/nan)')

    def __sub__(self, other):
        return self.__add__(-other)

    def __rsub__(self, other):
        return self.__radd__(-other)

    def __isub__(self, other):
        return self.__iadd__(-other)

    def __mul__(self, other):
        if isinstance(other, Measurement):
            uncertainty = np.sqrt(np.add((self.value*other.uncertainty)**2,
                                         (self.uncertainty*other.value)**2))
            return self._new_view(np.multiply(self.value, other.value), self.unit*other.unit,
                                  uncertainty).decompose()
        elif isinstance(other, u.Quantity):
            return self._new_view(np.multiply(self.value, other.value), self.unit*other.unit,
                                  np.multiply(self.uncertainty, np.abs(other.value))).decompose()
        
        elif isinstance(other, u.UnitBase):
            return self._new_view(self.value, self.unit*other, self.uncertainty).decompose()
        
        else:
            return self._new_view(np.multiply(self.value, other), self.unit,
                                  np.multiply(self.uncertainty, np.abs(other))).decompose()

    def __rmul__(self, other):
        if isinstance(other, Measurement):
            uncertainty = np.sqrt(np.add((other.uncertainty*self.value)**2,
                                         (other.value*self.uncertainty)**2))
            return self._new_view(np.multiply(other.value, self.value), other.unit*self.unit,
                                  uncertainty).decompose()
        elif isinstance(other, u.Quantity):
            return self._new_view(np.multiply(other.value, self.value), self.unit*other.unit,
                                  np.multiply(np.abs(other.value), self.uncertainty)).decompose()

        elif isinstance(other, u.UnitBase):
            return self._new_view(self.value, other*self.unit, self.uncertainty).decompose()
        
        else:
            return self._new_view(np.multiply(other, self.value), self.unit,
                                  np.multiply(np.abs(other), self.uncertainty))

    def __imul__(self, other):
        if isinstance(other, Measurement):
            self._uncertainty = np.sqrt(np.add((self.value*other.uncertainty)**2,
                                               (other.value*self.uncertainty)**2))
            self._value *= other.value
            self._unit = self.unit*other.unit
            return self
        
        elif isinstance(other, u.Quantity):
            self._uncertainty *= np.abs(other.value)
            self._value *= other.value
            self._unit = self.unit*other.unit
            return self
        
        elif isinstance(other, u.UnitBase):
            self._unit = self.unit*other
            return self
        
        else:
            self._uncertainty *= np.abs(other)
            self._value *= other
            return self

    def __truediv__(self, other):
        if isinstance(other, Measurement):
            uncertainty = np.sqrt(np.add(np.true_divide(self.uncertainty, other.value)**2,
                          np.true_divide(self.value*other.uncertainty, other.value**2)**2))
            return self._new_view(np.true_divide(self.value, other.value), self.unit/other.unit, uncertainty).decompose()
        
        elif isinstance(other, u.Quantity):
            return self._new_view(np.true_divide(self.value, other.value), self.unit/other.unit,
                                  np.true_divide(self.uncertainty, np.abs(other.value))).decompose()
        
        elif isinstance(other, u.UnitBase):
            return self._new_view(self.value, self.unit/other, self.uncertainty).decompose()
            
        else:
            return self._new_view(np.true_divide(self.value, other), self.unit,
                                  np.true_divide(self.uncertainty, np.abs(other)))

    def __rtruediv__(self, other):
        if isinstance(other, Measurement):
            uncertainty = np.sqrt(np.add((np.true_divide(other.uncertainty, self.value)**2,
                                          np.true_divide(other.value*self.uncertainty, self.value**2)**2)))
            return self._new_view(np.true_divide(other.value, self.value), other.unit/self.unit, uncertainty).decompose()
        
        elif isinstance(other, u.Quantity):
            uncertainty = np.abs(np.true_divide(other.value*self.uncertainty, self.value**2))
            return self._new_view(np.true_divide(other.value, self.value), other.unit/self.unit, uncertainty).decompose()
        
        elif isinstance(other, u.UnitBase):
            return self._new_view(1/self.value, other/self.unit, 1/self.uncertainty).decompose()
        
        else:
            return self._new_view(np.true_divide(other, self.value), 1/self.unit,
                                  np.true_divide(np.abs(other)*self.uncertainty, self.value**2))

    def __itruediv__(self, other):
        if isinstance(other, Measurement):
            self._uncertainty = np.sqrt(np.add((np.true_divide(self.uncertainty, other.value)**2,
                                          np.true_divide(self.value*other.uncertainty, other.value**2)**2)))
            self._value = np.true_divide(self.value, other.value)
            self._unit = self.unit/other.unit
            return self
        
        elif isinstance(other, u.Quantity):
            self._uncertainty = np.true_divide(self.uncertainty, np.abs(other.value))
            self._unit = self.unit/other.unit
            self._value = np.true_divide(self.value, other.value)
            return self
        
        elif isinstance(other, u.UnitBase):
            self._unit = self.unit/other
            return self
        
        else:
            self._uncertainty = np.true_divide(self.uncertainty/np.abs(other))
            self._value = np.true_divide(self.value, other)
            return self

    def __div__(self, other):
        if isinstance(other, Measurement):
            uncertainty = np.sqrt(np.add((np.divide(self.uncertainty, other.value)**2,
                                          np.divide(self.value*other.uncertainty, other.value**2)**2)))
            return self._new_view(np.divide(self.value, other.value), self.unit/other.unit, uncertainty).decompose()
        
        elif isinstance(other, u.Quantity):
            return self._new_view(np.divide(self.value, other.value), self.unit/other.unit,
                                  np.divide(self.uncertainty, np.abs(other.value))).decompose()
        
        elif isinstance(other, u.UnitBase):
            return self._new_view(self.value, self.unit/other, self.uncertainty).decompose()
        
        else:
            return self._new_view(np.divide(self.value, other), self.unit,
                                  np.divide(self.uncertainty, np.abs(other)))

    def __rdiv__(self, other):
        if isinstance(other, Measurement):
            uncertainty = np.sqrt(np.add((np.divide(other.uncertainty, self.value)**2,
                                          np.divide(other.value*self.uncertainty, self.value**2)**2)))
            return self._new_view(np.divide(other.value, self.value), other.unit/self.unit, uncertainty).decompose()
        
        elif isinstance(other, u.Quantity):
            uncertainty = np.abs(np.divide(other.value*self.uncertainty, self.value**2))
            return self._new_view(np.divide(other.value, self.value), other.unit/self.unit, uncertainty).decompose()
        
        elif isinstance(other, u.UnitBase):
            return self._new_view(self.value, self.unit/other, self.uncertainty).decompose()
        
        else:
            return self._new_view(np.divide(other, self.value), 1/self.unit,
                                  np.divide(np.abs(other)*self.uncertainty, self.value**2))

    def __idiv__(self, other):
        if isinstance(other, Measurement):
            self._uncertainty = np.sqrt(np.add((np.divide(self.uncertainty, other.value)**2,
                                          np.divide(self.value*other.uncertainty, other.value**2)**2)))
            self._value = np.divide(self.value, other.value)
            self._unit = self.unit/other.unit
            return self
        
        elif isinstance(other, u.Quantity):
            self._uncertainty = np.divide(self.uncertainty, np.abs(other.value))
            self._unit = self.unit/other.unit
            self._value = np.divide(self.value, other.value)
            return self
        
        elif isinstance(other, u.UnitBase):
            self._unit = self.unit/other
            return self
        
        else:
            self._uncertainty = np.divide(self.uncertainty/np.abs(other))
            self._value = np.divide(self.value, other)
            return self

    def __floordiv__(self, other):
        if isinstance(other, Measurement):
            uncertainty = np.floor_divide(np.sqrt(np.add((np.true_divide(self.uncertainty, other.value)**2,
                                          np.true_divide(self.value*other.uncertainty, other.value**2)**2))), 1)
            return self._new_view(np.floor_divide(self.value, other.value), self.unit/other.unit, np.uncertainty).decompose()

        elif isinstance(other, u.Quantity):
            return self._new_view(np.floor_divide(self.value, other.value), self.unit/other.unit,
                                  np.floor_divide(self.uncertainty, np.abs(other.value))).decompose()
        
        elif isinstance(other, u.UnitBase):
            return self._new_view(self.value, self.unit/other, self.uncertainty).decompose()
        
        else:
            return self._new_view(np.floor_divide(self.value, other), self.unit,
                                  np.floor_divide(self.uncertainty, np.abs(other)))

    def __rfloordiv__(self, other):
        if isinstance(other, Measurement):
            uncertainty = np.floor_divide(np.sqrt(np.add((np.true_divide(other.uncertainty, self.value)**2,
                                          np.true_divide(other.value*self.uncertainty, self.value**2)**2))), 1)
            return self._new_view(np.floor_divide(other.value, self.value), other.unit/self.unit, uncertainty).decompose()
        
        elif isinstance(other, u.Quantity):
            uncertainty = np.abs(np.floor_divide(other.value*self.uncertainty, self.value**2))
            return self._new_view(np.floor_divide(other.value, self.value), other.unit/self.unit, uncertainty).decompose()
        
        elif isinstance(other, u.UnitBase):
            return self._new_view(1/self.value, other/self.unit, 1/self.uncertainty).decompose()

        else:
            return self._new_view(np.floor_divide(other, self.value), 1/self.unit,
                                  np.floor_divide(np.abs(other)*self.uncertainty, self.value**2))

    def __ifloordiv__(self, other):
        if isinstance(other, Measurement):
            self._uncertainty = np.floor_divide(np.sqrt(np.add((np.divide(self.uncertainty, other.value)**2,
                                          np.divide(self.value*other.uncertainty, other.value**2)**2))), 1)
            self._value = np.floor_divide(self.value, other.value)
            self._unit = self.unit/other.unit
            return self
        
        elif isinstance(other, u.Quantity):
            self._uncertainty = np.floor_divide(self.uncertainty, np.abs(other.value))
            self._unit = self.unit/other.unit
            self._value = np.floor_divide(self.value, other.value)
            return self
        
        elif isinstance(other, u.UnitBase):
            self._unit = self.unit/other
            return self
       
        else:
            self._uncertainty = np.floor_divide(self.uncertainty/np.abs(other))
            self._value = np.floor_divide(self.value, other)
            return self

    # __mod__ functions are not overrrided by u.Quantity. They do not work with quantities... I think it should
    
    def __pow__(self, other):
        if isinstance(other, Measurement):
            if other.unit.physical_type != u'dimensionless':
                raise u.UnitTypeError('Can only raise something to a dimensionless quantity')

            new_value = self.value**other.value
            uncertainty = np.abs(new_value)*np.sqrt((other.value*self.uncertainty/self.value)**2 +\
                                                    (other.uncertainty*np.log(self.value))**2)
            return self._new_view(new_value, self.unit**other.value, uncertainty)
        elif isinstance(other, u.Quantity):
            if other.unit.physical_type != u'dimensionless':
                raise u.UnitTypeError('Can only raise something to a dimensionless quantity')

            new_value = self.value**other.value
            uncertainty = np.abs(new_value*other.value*self.uncertainty/self.value)
            return self._new_view(new_value, self.unit**other.value, uncertainty)
        else:
            return self._new_view(self.value**other, self.unit**other,
                                  np.abs(self.value**(other-1)*other*self.uncertainty))
    
    def __rpow__(self, other):
        if self.unit.physical_type != u'dimensionless':
            raise u.UnitTypeError('Can only raise something to a dimensionless quantity')
            
        if isinstance(other, Measurement):
            new_value = other.value**self.value
            uncertainty = np.abs(new_value)*np.sqrt((self.value*other.uncertainty/other.value)**2+
                                                    (self.uncertainty*np.log(other.value))**2)
            return self._new_view(new_value, other.unit**self.value, uncertainty)

        elif isinstance(other, u.Quantity):
            new_value = other.value**self.value
            uncertainty = np.abs(new_value*self.uncertainty*np.log(other.value))
            return self._new_view(new_value, other.unit**self.value, uncertainty)
        else:
            new_value = other**self.value
            return self._new_view(new_value, u.dimensionless_unscaled,
                                  np.abs(new_value*self.uncertainty*np.log(other)))
    
    def __ipow__(self, other):
        if isinstance(other, Measurement):
            if other.unit.physical_type != u'dimensionless':
                raise u.UnitTypeError('Can only raise something to a dimensionless quantity')

            new_value = self.value**other.value
            self._uncertainty = np.abs(new_value)*np.sqrt((other.value*self.uncertainty/self.value)**2 +\
                                                    (other.uncertainty*np.log(self.value))**2)
            self._unit = self.unit**other.value
            self._value = new_value
            return self
        elif isinstance(other, u.Quantity):
            if other.unit.physical_type != u'dimensionless':
                raise u.UnitTypeError('Can only raise something to a dimensionless quantity')

            new_value = self.value**other.value
            uncertainty = np.abs(new_value*other.value*self.uncertainty/self.value)
            return self._new_view(new_value, self.unit**other.value, uncertainty)
        else:
            return self._new_view(self.value**other, self.unit**other,
                                  np.abs(new_value*other.value*self.uncertainty/self.value))

    def __neg__(self):
        return self._new_view(-self.value, self.unit, self.uncertainty)

    def __pos__(self):
        return self._new_view(+self.value, self.unit, self.uncertainty)
    
    def __abs__(self):
        return self._new_view(np.abs(self.value), self.unit, self.uncertainty)
    

    # Functions from numpy
    def subtract(self, x1, x2, *args, **kwargs):
        # I NEED TO DO ALL THE IMPLEMENTATIONS!!!!!!!!!!!!!!!!!!  ON X1 AND X2
        if isinstance(x1, Measurement):
            if x2.unit is not x1.unit:
                x2 = x2.to(x1.unit)

            uncertainty = np.sqrt(np.add(x2.uncertainty**2, x1.uncertainty**2, **kwargs))
            return self._new_view(np.subtract(x1.value, x2.value, **kwargs), x1.unit, uncertainty).decompose()

        elif isinstance(x2, u.Quantity):
            if x2.unit is not x1.unit:
                x2 = x2.to(x1.unit)

            return self._new_view(np.subtract(x1.value, x2.value, **kwargs), x1.unit, x2.uncertainty).decompose()
        else:
            raise UnitsError('Can only apply \'add\' function to dimensionless quantities when other'+\
                             ' argument is not a quantity (unless the latter is all zero/infinity/nan)')

    def log(self, *args, **kwargs):
        if self.unit.physical_type != u'dimensionless':
            raise UnitTypeError('Can only apply \'log\' function to dimensionless quantities')
        return self._new_view(np.log(self.value, **kwargs), u.dimensionless_unscaled,
                              self.uncertainty/np.abs(self.value))

    def log10(self, *args, **kwargs):
        if self.unit.physical_type != u'dimensionless':
            raise UnitTypeError('Can only apply \'log10\' function to dimensionless quantities')
        return self._new_view(np.log10(self.value, **kwargs), u.dimensionless_unscaled,
                              self.uncertainty/(np.abs(self.value)*np.log(10)))

    def exp(self, *args, **kwargs):
        if self.unit.physical_type != u'dimensionless':
            raise UnitTypeError('Can only apply \'exp\' function to dimensionless quantities')
        new_value = np.exp(self.value, **kwargs)
        return self._new_view(new_value, u.dimensionless_unscaled, np.abs(new_value)*self.uncertainty)

    def sqrt(self, *args, **kwargs):
        return self._new_view(np.sqrt(self.value, **kwargs), self.unit**(1/2),
                              self.uncertainty/(2*np.sqrt(self.value, **kwargs)))
   
    def sin(self, *args, **kwargs):
        if self.unit.physical_type not in (u'dimensionless', u'angle'):
            raise UnitTypeError('Can only apply \'sin\' function to quantities with angle units')
        # print(self.value, self.unit, self.uncertainty)
        # print(np.sin(self.value*self.unit))
        # print(np.cos(self.value*self.unit))
        return self._new_view(np.sin(self.value*self.unit), u.dimensionless_unscaled,
                              self.uncertainty*self.unit.to(u.rad)*np.abs(np.cos(self.value*self.unit)))

    def cos(self, *args, **kwargs):
        if self.unit.physical_type not in (u'dimensionless', u'angle'):
            raise UnitTypeError('Can only apply \'sin\' function to quantities with angle units')
        return self._new_view(np.cos(self.value*self.unit), u.dimensionless_unscaled,
                              self.uncertainty*self.unit.to(u.rad)*np.abs(np.sin(self.value*self.unit)))

    def tan(self, *args, **kwargs):
        if self.unit.physical_type not in (u'dimensionless', u'angle'):
            raise UnitTypeError('Can only apply \'sin\' function to quantities with angle units')
        return self._new_view(np.tan(self.value*self.unit), u.dimensionless_unscaled,
                              self.uncertainty*self.unit.to(u.rad)/np.abs(np.cos(self.value*self.unit)**2))

    def arctan(self, *args, **kwargs):
        if self.unit.physical_type  is not u'dimensionless':
            raise UnitTypeError('Can only apply \'arctan\' function to dimensionless quantities')
        return self._new_view(np.arctan(self.value), u.rad,
                              self.uncertainty/(1+self.value**2))

    def arcsin(self, *args, **kwargs):
        if self.unit.physical_type  is not u'dimensionless':
            raise UnitTypeError('Can only apply \'arcsin\' function to dimensionless quantities')
        return self._new_view(np.arcsin(self.value), u.rad,
                              self.uncertainty/(1-self.value**2))
    
    def arccos(self, *args, **kwargs):
        if self.unit.physical_type  is not u'dimensionless':
            raise UnitTypeError('Can only apply \'arccos\' function to dimensionless quantities')
        return self._new_view(np.arccos(self.value), u.rad,
                              self.uncertainty/(1-self.value**2))

    def arctan2(self, x1, x2, *args, **kwargs):
        if isinstance(x2, Measurement):
            if x2.unit is not x1.unit:
                x2 = x2.to(x1.unit)
            uncertainty = np.sqrt((x1.uncertainty/x2.value)**2 + \
                                  (x1.value*x2.uncertainty/(x2.value**2))**2)
            uncertainty = uncertainty/(1+(x1.value/x2.value)**2)
            return self._new_view(np.arctan2(x1.value, x2.value, *args, **kwargs),
                        u.rad, uncertainty)

        elif isinstance(x2, u.Quantity):
            if x2.unit is not x1.unit:
                x2 = x2.to(x1.unit)

            uncertainty = x1.uncertainty/(1+(x1.value/x2.value)**2)
            return self._new_view(np.arctan2(x1.value, x2.value, *args, **kwargs),
                        u.rad, uncertainty)

        else: 
            return self._new_view(np.arctan2(x1.value, x2, *args, **kwargs),
                        u.rad, x1.uncertainty/(1+(x1.value/x2)**2))

    # def subtract(self, other, **kwargs):
    #     if isinstance(other, Measurement):
    #         if other.unit is not self.unit:
    #             other = other.to(self.unit)
    #
    #         uncertainty = np.sqrt(np.add(other.uncertainty**2, self.uncertainty**2, **kwargs))
    #         return self._new_view(np.subtract(self.value, other.value, **kwargs), self.unit, uncertainty)
    #
    #     elif isinstance(other, u.Quantity):
    #         if other.unit is not self.unit:
    #             other = other.to(self.unit)
        #
        #     return self._new_view(np.subtract(self.value, other.value, **kwargs), self.unit, self.uncertainty)
        # else:
        #     raise UnitsError('Can only apply \'add\' function to dimensionless quantities when other'+\
        #                      ' argument is not a quantity (unless the latter is all zero/infinity/nan)')
        #

    def sum(self, *args, **kwargs):
        """Compatible with np.sum due to the kwargs... alwasy keep or put the same structure as the
        numpy functions
        """
        return self._new_view(np.sum(self.value, **kwargs), self.unit, np.sum(self.uncertainty, **kwargs))


        
    # We use the corresponding numpy functions to evaluate the results, since
    # the methods do not always allow calling with keyword arguments.
    # For instance, np.array([0.,2.]).clip(a_min=0., a_max=1.) gives
    # TypeError: 'a_max' is an invalid keyword argument for this function.
    # def _wrap_function(self, function, *args, **kwargs):
    #     """Wrap a numpy function that processes self, returning a Quantity.
    #
    #     Parameters
    #     ----------
    #     function : callable
    #         Numpy function to wrap.
    #     args : positional arguments
    #         Any positional arguments to the function beyond the first argument
    #         (which will be set to ``self``).
    #     kwargs : keyword arguments
    #         Keyword arguments to the function.
    #
    #     If present, the following arguments are treated specially:
    #
    #     unit : `~astropy.units.Unit`
    #         Unit of the output result.  If not given, the unit of ``self``.
    #     out : `~astropy.units.Quantity`
    #         A Quantity instance in which to store the output.
    #
    #     Notes
    #     -----
    #     Output should always be assigned via a keyword argument, otherwise
    #     no proper account of the unit is taken.
    #
    #     Returns
    #     -------
    #     out : `~astropy.units.Quantity`
    #         Result of the function call, with the unit set properly.
    #     """
    #     unit = kwargs.pop('unit', self.unit)
    #     uncertainty = kwargs.pop('uncertainty', self.uncertainty)
    #     out = kwargs.pop('out', None)
    #     # Ensure we don't loop back by turning any Quantity into array views.
    #     args = (self.value,) + tuple((arg.value if isinstance(arg, Quantity)
    #                                   else arg) for arg in args)
    #     if out is not None:
    #         # If pre-allocated output is used, check it is suitable.
    #         # This also returns array view, to ensure we don't loop back.
    #         arrays = tuple(arg for arg in args if isinstance(arg, np.ndarray))
    #         kwargs['out'] = check_output(out, unit, arrays, function=function)
    #     # Apply the function and turn it back into a Quantity.
    #     result = function(*args, **kwargs)
    #     return self._result_as_quantity(result, unit, out)
    #
