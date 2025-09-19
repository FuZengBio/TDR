# TDR Analysis
Working directory: First, download all code to the working directory:
WORK_DIR = 'Z:\Users\TDR\';
_____________________________________________________________________________
1. Analysis of 1s Stimulus Data in VIP and MSTd
Run code: TDR_analysis_for_VIP_MSTd_1s
Notes:
(1) If using the original TDR method (two configurations, e.g., Heading-choice in Figure 3): Set Kernel_method = 0;
To change the orthogonalization order, set reverse = 1, Default reverse = 0 corresponds to the order heading → choice
(2) If using the extended, time-varying TDR analysis (dissociating motion parameters: velocity & acceleration, and choice, e.g., Figures 4, 5, 6):
Set Kernel_method = 1;
To change the orthogonalization order, set reverse = 1 or 2
Default reverse = 0 corresponds to velocity → acceleration → choice
reverse = 1: orthogonalization order becomes choice → acceleration → velocity
reverse = 2: move acceleration to the end, i.e., velocity → choice → acceleration
_____________________________________________________________________________
2. Analysis of 2s Stimulus VIP and MSTd Data, and 1s PIVC Data (Figures 7，Supplemental Figure 1 & 2)
Run code: TDR_analysis_VIP_MSTd_2s_PIVC_1s
Order settings:
reverse = 0: velocity → acceleration → choice
reverse = 1: choice → acceleration → velocity
reverse = 2: velocity → choice → acceleration
_____________________________________________________________________________
Data Availability Note
Due to upload file size limits, the data files cannot be provided directly.
If needed, please contact the corresponding author to obtain the data and place it in the corresponding folder:
DATA_DIR = 'Z:\Users\TDR\Data\';
