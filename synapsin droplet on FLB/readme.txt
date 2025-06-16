Vesicle clamping with complexin
Last Updated: 2025.06.16

This repository contains a MATLAB script for analyzing the diffusion behavior and fluorescence intensities of particles on FLB. The workflow includes image registration, particle detection, tracking, and intensity quantification.
-------------------------
1. Main Script
-------------------------
- analysis.m
  - Loads and processes TIFF image stacks (2 channels per field).
  - Applies affine image registration using manually selected landmarks.
  - Extracts cropped regions of interest (ROI) for visualization.
  - Detects particles using local maximum, and save their coordinates.
  - Identifies colocalized particles based on the proximity of detections in both channels.
  - Measures and stores fluorescence intensity at each particle's position.
  - Generates summary plots and saves results.

-------------------------
2. Input Requirements
-------------------------
- TIFF images in `raw1_before HD\` and `raw2_after HD\` directories.
- Affine transformation information from manual image registration is stored in tform.mat.
-------------------------
3. Output Files
-------------------------
- data.mat : loaded and preprocessed image data
- Figures: visualization of colocalizaed spot, combined movie, and analysis plot
-------------------------
4. Dependencies
-------------------------
- MATLAB (R2022b or newer)
- Image Processing Toolbox, Statistics and Machine Learning toolbox
- Custom functions required in path:
  imshow3, particle detection, writeTIFrgb,  plotSpread, colocalization, find_tform, dualviewer_spliter
- bresenham.m function by Aaron Wetzler required in path (https://kr.mathworks.com/matlabcentral/fileexchange/28190-bresenham-optimized-for-matlab)

-------------------------
5. Contact
-------------------------
For questions, contact: mjshon@postech.ac.kr
