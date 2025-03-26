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
  iAxis = 0
  iMin = 80
  iMax = 16000
  iNorm = -1
  iStepSize = 100
  iX = 0.1
  iY = 0
  iZ = 0
  iSigma = 0 ; default
  iRho = 0 ; default
  iBeta = 0 ; default
  iTimeDelta = 0.001 ; default
  iLorenzFreqs ftgen 0, 0, 512, "lorenz", iAxis, iMin, iMax, iNorm, iStepSize, iX,\
    iY, iZ, iSigma, iRho, iBeta, iTimeDelta
  
  aIndex line 0,p3,1  
  aFreq table3 aIndex, iLorenzFreqs, 1
  aSig poscil3 0.8, aFreq  
  aEnv linseg 0,0.02,1,p3-0.02,1,0.02,0

  outs aSig, aSig
endin

</CsInstruments>
<CsScore>
i1 0 20
</CsScore>
</CsoundSynthesizer>
