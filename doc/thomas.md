# Lorenz

## Abstract

GEN routine based on the Thomas attractor


## Description

``thomas`` creates a ft based on data created by a Thomas attractor. 

The Thomas attractor is a chaotic dynamical system described by a set
of nonlinear differential equations. Discovered by Pierre Thomas in
1989, it exhibits a cyclically symmetric structure and involves three
variables: x,y, and z, which typically follow a chaotic
trajectory. The Thomas attractor is known for its compact and dense
form, often resembling a distorted "butterfly" shape. Compared to
other chaotic systems like the Lorenz attractor, it features higher
symmetry and intricate cyclic interactions among its variables. Like
all chaotic systems, the Thomas attractor is highly sensitive to
initial conditions, meaning that even tiny changes in the starting
values can result in vastly different outcomes.



## Syntax


```csound

f1 0 4096 "thomas" 0 30 1600 -1 1 0.001 0 0 0.15 0.01

```
    
## Arguments
* p5 = choose axis to select values from; 0 = x, 1 = y, 2 = z
* p6 = min output; if p6 and p7 == 0 don't scale
* p7 = max output; if p6 and p7 == 0 don't scale
* p8 = normalize; 0 == normalize, -1 == don't normalize
* p9 = stepsize
* p10 = x start value
* p11 = y start value
* p12 = z start value
* [p13 = b; 0 == default -> 0.9] 
* [p14 = time delta; 0 == default -> 0.001] 


## Output


## Execution Time

* Init 

## Examples


```csound
<CsoundSynthesizer>
<CsOptions>
-odac
</CsOptions>
<CsInstruments>

sr = 44100
ksmps = 16
nchnls = 2
0dbfs = 1.0

instr 1
  iX random -1, 1
  iY random -1, 1
  iZ random -1, 1
  iStepSize = 1
  iB = 0.15
  iDt = 0.01
  iNorm = -1
  iMin = 80
  iMax = 600
  iAxis = 0
  iFreqs ftgen 0, 0, 16384, "thomas", iAxis, iMin, iMax, iNorm, iStepSize, iX,\
    iY, iZ, iB, iDt
  
  aIndex line 0,p3,1  
  aFreq table aIndex, iFreqs, 1
  aSig poscil3 0.8, aFreq  
  aEnv linseg 0,0.02,1,p3-0.02,1,0.02,0

  outs aSig, aSig
endin

</CsInstruments>
<CsScore>
i1 0 20
</CsScore>
</CsoundSynthesizer>

```


## See also

* [thomas attractor](https://en.wikipedia.org/wiki/Thomas'_cyclically_symmetric_attractor)

## Credits

Philipp von Neumann, 2024
