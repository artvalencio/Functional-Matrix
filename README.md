# Functional-Matrix
(beta)

Calculates the dynamical functional connectivity (DFC) of brain areas given an EEG signal.

Each data point gives the peak sliding window correlation or information-theoretical measure
between two EEG channels at a given time interval and the time lag where this peak is found. 
Such time lag could be interpreted as the time it takes for area j to receive signal
from area i.

Use DFC to calculate the dynamic functional connectivity using correlation measures;

Use DFCCAMI to calculate the dynamic functional connectivity using information-theoretical measures (CaMI,MI,TE,DI)


In construction:

(1) EEG pre-processing batch

(2) visual output of DFC results

