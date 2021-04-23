# Gait-phase-detection
These pages hold original data collected for experiments reported in our paper entiled
**Online Gait Phase Detection in Complex Environment Based on Distance and Multi-Sensors Information Fusion Using Inertial Measurement Units**.
We provide the data publicly for other researchers in the same field get benefited.
# DTWMean
The folder 'DTWMean' includes Source code also available at http://www.akt.tu-berlin.de/menue/software (Java source code)
# Data
The folder 'Data' includes all the experiments involved in the paper, for both training and testing data, it mainly uses Matlab format. In testing data file, there are mainly two types of the files, one is the data file to be tested ending with the suffix .mat and another one is gait phase detection results with the suffix .fig for four states, which are the initial off-ground (IOff), the off-ground (Off), the initial on-ground (IOn), and the on-ground (On) states in a complete gait cycle. The 1-5 at the end of the file name represent different speeds that are 2km/h,3km/h,4km/h,5km/h,and 6km/h.The columns of the .mat file are defined as follows:
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
The results show that the method can adaptively detect the gait phase in real time for different terrains and different payloads while walking at different speeds. Our proposed method can gain over 95% accuracy and time difference is 4.99±15.05ms
