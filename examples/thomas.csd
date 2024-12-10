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

