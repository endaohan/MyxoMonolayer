# MyxoMonolayer
Code used in the paper Local polar order controls mechanical stress and triggers layer formation in  Myxococcus xanthus colonies. 

The codes demonstrate how the authors obtained the results reported in the paper Local polar order controls mechanical stress and triggers layer formation in Myxococcus xanthus colonies by Han et al. 
- TFM: It demonstrates how we processed the data obtained with the traction force microscopy (TFM) assay. It provides traction force, director field (cell orientation), and defect location and orientation. 
- Polarity: It demonstrates how we processed the data obtained with the polarity assay. It provides cell polarity, cell velocity, and defect tracking. 

1. System requirements
The code was developed and tested using MATLAB R2019b on a system running Windows 10 Pro, version 22H2.

2. Installation guide
The user needs to install Matlab, which could take 30 to 45 minutes.  

3. Demo
Instructions to run on data
- TFM: use the code TFM_run.m. Open the file in Matlab and click Run. The original data to analyze are in the folders "TFM/ExampleData/BFld/Images" and "TFM/ExampleData/Beads_1/Images". 
- Polarity: use the code run_polarity_analysis.m. Open the file in Matlab and click Run. The original data to analyze are in the folders" and "Polarity/ExampleData/Laser". 
- When a Matlab window jumps out saying "To run this file, you can either change the MATLAB" current folder or add its folder to the MATLAB path." click the button "Change Folder". 

Expected output
- The anticipated output is saved in the same folder as the raw data.
- The TFM code outputs displacement fields of the fluorescent beads, director (orientation) field of the cells taken with bright field imaging, detected topological defects with tracking, and reconstructed traction force. A demo video is saved with the name "Traction_hf.mp4". 
- The Polarity code outputs cell polarity, cell velocity obtained with optical flow, and defect tracking. A demo video is saved with the name Label_video.mp4". 

Expected run time for demo on a "normal" desktop computer: 
- For the TFM demo, the run time is about 45 minutes on a normal laptop computer. 
- For the Polarity demo, the run time is about 10 minutes on a normal laptop computer. 

4. Instructions for use
- To run the TFM code, the user needs to put the images of the fluorescent particles in the folder "Laser" and bright field images of the cells in the folder "BFld". The corresponding parameters in the experiments should all be saved in the file Parameters.mat, which can be read by the code. 
- To run the Polarity code, the user needs to put the images of the fluorescent particles in the folder "Laser" and bright field images of the cells in the folder "BFld". 



In this code, we have used the codes and packages written by other researchers. We appreciate all of them for sharing their work. 
- The folder track_ED contains the MATLAB Particle Tracking Code written by Daniel Blair and Eric Dufresne, which can be found here: https://site.physics.georgetown.edu/matlab/. 
- The code in the folder TFM_PIV is tweaked based on the PIV code written by Ivo Peters. We use a relatively early version. The latest version can be found here: https://ivopeters.org/code/. 
- Some functions in the code are provided by a coauthor of the paper, Katherine Copenhagen. Her GitHub page is: https://github.com/kcopenhagen. 


