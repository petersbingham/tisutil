# tisutil
Extension of matfuncutil for scattering matrices and other quantities.

## Installation

Clone the repository and install with the following commands:

    git clone https://github.com/petersbingham/tisutil.git
    cd tisutil
    python setup.py install
    
## Dependencies

Author (these will have their own dependencies):
 - [matfuncutil](https://github.com/petersbingham/matfuncutil)
 - [channelutil](https://github.com/petersbingham/channelutil)

## Overview

This package extends the containers in the [matfuncutil](https://github.com/petersbingham/matfuncutil) for the various scattering representations and provides common units and conversions. The provided discrete scattering matrix containers are `dSmat`, `dKmat`, `dTmat`, `dQmat`, `dXSmat` and `dEPhaseMat`. Also provided are the scalar containers, `dXSsca` and `dEPhaseSca`. Functions of the form `to_dSmat` etc are provided by the `dSmat`, `dKmat`, `dTmat` etc. to allow easy conversions between the representations, as well as `to_dXSsca` to calculate the total cross section. 

A continuous container is also provided, `cMatSympypolyk`, which extends `matfuncutil.cMatSympypoly` by allowing energy parametrisation using the `AsymCalc` object.

Conversions, cross sections, eigenphase sums and the Q-matrix are calculated as:<br/>

![equations](https://github.com/petersbingham/tisutil/blob/master/equations.jpg)

## Usage

A useful illustration of tisutil is given in the [matfuncutil documentation](https://github.com/petersbingham/matfuncutil). We provide an additional example here to illustate the conversion between the different scattering representations. This example plots the K-matrix and total cross section for the two channel radial well (see [twochanradialwell](https://github.com/petersbingham/twochanradialwell)).

```python
import twochanradialwell as radwell
import channelutil as chanutil
asymcalc = chanutil.AsymCalc(chanutil.hartrees, thresholds=[0.,2.])
csmat = radwell.get_Smat_fun(1., 2., 2., asymcalc, 1.)

dsmat = csmat.discretise(1.0, 8.0, 200) # Create the discrete container.

dsmat2 = dsmat.to_dSmat() # This just returns a reference to itself.
dkmat = dsmat.to_dKmat()
dkmat.plot()

dxsmat = dkmat.to_dXSmat() # or dsmat.to_dXSmat()
dtotxssca = dxsmat.to_dXSsca()
dtotxssca.plot()
```
