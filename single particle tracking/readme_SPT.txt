Single particle tracking on FLB
Last Updated: 2025.06.16

This repository contains a MATLAB script for analyzing the diffusion behavior and fluorescence intensities of particles on FLB. The workflow includes image registration, particle detection, tracking, and intensity quantification.
-------------------------
1. Main Script
-------------------------
- analysis.m
  - Loads and processes TIFF image stacks.
  - Extracts cropped regions of interest (ROI) for visualization.
  - Detects particles using local maximum, and save their coordinates for all frames.
  - Tracks particles across frames by setting:
	- The maximum allowed displacement per frame.
	- The minimum number of frames a particle must appear in to be considered a valid trajectory.
  - Measures and stores fluorescence intensity at each particle's position.
  - Generates summary plots and saves results.

-------------------------
2. Input Requirements
-------------------------
- TIFF images in `raw\` directories.
-------------------------
3. Output Files
-------------------------
- data.mat : loaded and preprocessed image data
- Figures: MSD plot, tracking movie, and analysis plot
-------------------------
4. Dependencies
-------------------------
- MATLAB (R2022b or newer)
- Image Processing Toolbox, Statistics and Machine Learning toolbox
- Custom functions required in path:
  imshow3, particle detection, writeTIFrgb, CalD, plotSpread, particle_tracking, trajectory
- bresenham.m function by Aaron Wetzler required in path (https://kr.mathworks.com/matlabcentral/fileexchange/28190-bresenham-optimized-for-matlab)

-------------------------
5. Contact
-------------------------
For questions, contact: mjshon@postech.ac.kr
