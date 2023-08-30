202308 Foregrounds and CMB from PanExp V1
=========================================

## Updates

## Summary

Full sky simulations for all LiteBIRD frequency channels of Galactic/Extragalactic foregrounds and CMB. The dataset includes Galactic foreground models at 3 different levels of complexity, the models have been selected by the Panexperiment Galactic Science group for common simulations across all CMB Experiments.
The instrument model assumes Gaussian beams and top-hat bandpasses.

## Instrument model

The instrument model has been extracted from `IMoV2-14June.json`, see the [notebook script for more details](utils/create_instrument_parameters.ipynb).
In particular we modified some of the center frequencies, because in the input JSON they were shifted to avoid multiple channels with the same frequency ([see the relevant issue in the private repository - requires permissions](https://github.com/litebird/litebirdXS4-private/pull/2)]:

```
M1-195 replace 196.0 GHz with 195.0 GHz
L4-140 replace 141.0 GHz with 140.0 GHz
L3-119 replace 120.0 GHz with 119.0 GHz
L4-100 replace 101.0 GHz with 100.0 GHz
L3-089 replace 90.0 GHz with 89.0 GHz
L4-078 replace 79.0 GHz with 78.0 GHz
L3-068 replace 69.0 GHz with 68.0 GHz
```

### Instrument model file format

The input and output instrument models are not public, they are available in the <https://github.com/litebird/litebirdXS4-private/> **private** repository. Ask [@giuspugl](https://github.com/giuspugl) for access.

[Text file in IPAC table format](../litebirdXS4-private/instrument_model_20230614/litebird_instrument_model.tbl), which supports units for the columns and can be read as `astropy.QTable`.

```python
from astropy import QTable
imo = QTable.read("../litebirdXS4-private/instrument_model_20230614/       litebird_instrument_model.tbl", format="ascii.ipac" )
imo.add_index("tag")
```

Example on how to access all parameters for a band, or a specific one:

```python
imo.loc["L3-089"]
imo.loc["L3-089"]["fwhm"]
```
