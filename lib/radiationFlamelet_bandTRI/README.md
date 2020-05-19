# WSGG Model

## Description


## Requirement

To use the libraries you have to modify the file

```
OpenFOAM/OpenFOAM-dev/src/thermophysicalModels/radiation/radiationModels/fvDOM/aabsorptionCoeffs/absorptionCoeffs.H 
```

Modify the line:

```
static const int nCoeffs_ = 12;

```


```
static const int nCoeffs_ = 6;

```


Then

```
wclean
wmake
```
