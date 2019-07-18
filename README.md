# Functional-Matrix
(beta)

Calculates the dynamical functional connectivity (DFC) of brain areas given an EEG signal.

Each data point gives the peak sliding window correlation between two EEG channels at
a given time interval and the time lag where this peak correlation was found. 
Such time lag could be interpreted as the time it takes for area j to receive signal
from area i.

To be included: 
(1) computation of DFC using mutual information, CaMI and transfer entropy.
(2) EEG pre-processing batch
(3) visual output of DFC results
