# Gait-phase-detection
These pages hold original data collected for experiments reported in our paper entiled
**An Online Gait Phase Detection Method Based on DWT-MEAN Templates**.
We provide the data publicly for other researchers in the same field get benefited.
# Data
The folder 'Data' includes all the experiments involved in the paper, for both training and testing data, it uses Matlab format. In testing data file, there are mainly two types of the files, one is the data file to be tested ending with the suffix .mat and another one is gait phase detection results with the suffix .fig for four states, which are the initial off-ground (IOff), the off-ground (Off), the initial on-ground (IOn), and the on-ground (On) states in a complete gait cycle. The columns of the .mat file are defined as follows:
* time, in second
* Left GCFs(The sum of the GCFs of the left heel and ball), in N
* Right GCFs(The sum of the GCFs of the right heel and ball), in N
* GCFs of the left ball,in N
* GCFs of the right ball,in N
* GCFs of the right heel,in N
* GCFs of the left heel,in N

Other useful informations are included in these .ini files:
* The detection time point of each state that are the sequential detection of the off-ground state , the initial on-ground state marked , the on-ground state , and the initial off-ground state.

# Results
The results show that the method can adaptively detect the gait phase in real time for different terrains and different payloads while walking at different speed compared with the reference methods, and gait phase detection can gain better detection results and the overall accuracy reaches about 95%.
