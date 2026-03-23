**README: Dimensionality-reduction analyses for VIP, MSTd, and PIVC datasets**
This folder contains MATLAB scripts for running time-varying targeted dimensionality reduction (**tTDR**), **original TDR**, time-course variance significance, comparison analyses across monkeys, hemispheres as well as complementary dimensionality-reduction analyses(mTDR and dPCA).
_____________________________________________________________________________ 
**1. Time-varying targeted dimensionality reduction (tTDR)**
Notes:
(1) To obtain tTDR results for VIP and MSTd in the 1 s condition, run: tTDTR_VIP_MSTd_1s.
(2) To obtain tTDR results for VIP and MSTd in the 2 s condition, together with PIVC in the 1 s condition, run: tTDR_VIP_MSTd_2s_PIVC_1s.

To change the orthogonalization order in the tTDR scripts, set: reverse = 0; or change it to: reverse = 1; or: reverse = 2;
The meaning of each option is:
Default: reverse = 0: orthogonalization order = velocity → acceleration → choice.
reverse = 1: orthogonalization order = choice → acceleration → velocity.
reverse = 2: move acceleration to the end, i.e. velocity → choice → acceleration.

**2. Original TDR cue-heading targeted dimensionality reduction**
To obtain original TDR cue-heading targeted dimensionality reduction results, run: originalTDR_analysis_CueHeading
Regressor order for original cue-heading TDR
In originalTDR_analysis_CueHeading, the variable reverse controls the regressor order for cue and heading:
reverse = 0;
Default: reverse = 0
regressor order = cue → heading
reverse = 1
regressor order = heading → cue

**3. Original TDR heading-choice targeted dimensionality reduction**

To obtain original TDR heading-choice targeted dimensionality reduction results, run:

originalTDR_analysis_HeadingChoice

This script runs the heading-choice TDR analysis across the following datasets:

VIP 1 s
VIP 2 s
MSTd 1 s
MSTd 2 s
PIVC 1 s

PIVC is only available in the 1 s condition, so the script automatically skips the non-existing PIVC 2 s dataset.

For each dataset, the script analyzes the two stimulus modalities separately:

Vestibular
Visual
Regressor order for original heading-choice TDR

In originalTDR_analysis_HeadingChoice, the variable reverse controls the regressor order for heading and choice:

reverse = 0;
Default: reverse = 0
regressor order = heading → choice
reverse = 1
regressor order = choice → heading
4. Time-course variance, significance, and decoding
4.1 tTDR (velocity / acceleration / choice axes)

If you want to compute time-course variance, its significance, and decoding for tTDR, i.e. along the velocity, acceleration, and choice axes, run:

tTDTR_timecourse_decoding

This script computes:

time-course variance on the tTDR axes
significance of variance
heading decoding
choice decoding
4.2 Original TDR (heading / choice axes)

If you want to compute time-course variance, its significance, and decoding for original TDR, i.e. along the heading and choice axes, run:

originalTDTR_HC_timecourse_decoding

This script computes:

time-course variance on the heading and choice axes
significance of variance
heading decoding
choice decoding

Vestibular and visual conditions are analyzed separately.

5. Vestibular–visual choice-axis angle before orthogonalization in tTDR

If you want to obtain the angle between the vestibular and visual choice axes before orthogonalization in tTDR, run:

tTDR_choiceAngle_Vis_vs_Vest

This script computes the raw choice axis separately for the vestibular and visual conditions, and compares the two axes before any orthogonalization step.

It reports quantities such as:

raw angle
unsigned angle
number of shared neurons
raw vestibular choice axis
raw visual choice axis
6. Separate-monkey and left–right hemisphere analyses in tTDR
6.1 Differences between monkeys within the same brain area

If you want to examine differences between monkeys within the same brain area in tTDR, run:

For the 1 s condition:

tTDTR_VIP_MSTd_1s_separate_monkey

For the 2 s condition:

tTDTR_VIP_MSTd_2s_separate_monkey

These scripts compare tTDR results between monkeys within the same recorded area.

6.2 Differences between left and right hemispheres

If you want to examine differences between the left and right hemispheres in tTDR, run:

For the 1 s condition:

tTDTR_VIP_MSTd_1s_LR_hemisphere

For the 2 s condition:

tTDTR_VIP_MSTd_2s_LR_hemisphere

These scripts compare tTDR results between left and right recording hemispheres.

Changing the orthogonalization order

All of the scripts above use the same reverse setting:

reverse = 0;   % 0: VAC; 1: reverse parameters; 2: move acceleration to the end

The meaning of each option is:

Default: reverse = 0
order = velocity → acceleration → choice
reverse = 1
reverse the parameter order
reverse = 2
move acceleration to the end, i.e. velocity → choice → acceleration
7. Model-based targeted dimensionality reduction (mTDR) and comparison with original TDR

If you want to obtain model-based targeted dimensionality reduction (mTDR) results and compare them with original TDR in the heading–choice framework, run:

mTDR_and_compare_originalTDR_HC

For each dataset, this script:

fits an mTDR model for heading and choice
computes the corresponding original TDR projections
performs additional comparisons between mTDR and original TDR
Changing the regressor order

To change the regressor order, modify:

reverse  = 0;
orderStr = 'HC';
Default:
reverse = 0; orderStr = 'HC';
order = heading → choice
To reverse the order, change both variables accordingly so that the regression order and output labels remain consistent.
8. Demixed principal component analysis (dPCA)

If you want to perform demixed principal component analysis (dPCA), please refer to the original dPCA resources from the Machens lab.

The main analysis scripts (MATLAB) are available at:

http://github.com/machenslab/elife2016dpca

The dPCA toolkit is located at:

https://github.com/machenslab/dPCA

These resources provide the original implementation and analysis framework for dPCA, which can be adapted for the present datasets.

Data path

Please place the data files in the following folder:

DATA_DIR = 'Z:\Users\TDR\Data\';
Data Availability Note

Due to upload file size limits, the data files cannot be provided directly.
If needed, please contact the corresponding author to obtain the data and place it in the corresponding folder:

DATA_DIR = 'Z:\Users\TDR\Data\';

Working directory: First, download all code to the working directory:
WORK_DIR = 'Z:\Users\TDR\';
_____________________________________________________________________________
1. Analysis of 1s Stimulus Data in VIP and MSTd
Run code: TDR_analysis_for_VIP_MSTd_1s
Notes:
(1) If using the original TDR method (two configurations, e.g., Heading-choice in Figure 3): Set Kernel_method = 0;
To change the orthogonalization order, set reverse = 1, Default reverse = 0 corresponds to the order heading → choice;
(2) If using the extended, time-varying TDR analysis (dissociating motion parameters: velocity & acceleration, and choice, e.g., Figures 4, 5, 6):
Set Kernel_method = 1;
To change the orthogonalization order, set reverse = 1 or 2
Default reverse = 0 corresponds to velocity → acceleration → choice;
reverse = 1: orthogonalization order becomes choice → acceleration → velocity;
reverse = 2: move acceleration to the end, i.e., velocity → choice → acceleration.
_____________________________________________________________________________
2. Analysis of 2s Stimulus VIP and MSTd Data, and 1s PIVC Data (Figures 7，Supplemental Figure 1 & 2)
Run code: TDR_analysis_VIP_MSTd_2s_PIVC_1s
Order settings:
reverse = 0: velocity → acceleration → choice;
reverse = 1: choice → acceleration → velocity;
reverse = 2: velocity → choice → acceleration.
_____________________________________________________________________________
Data Availability Note
Due to upload file size limits, the data files cannot be provided directly.
If needed, please contact the corresponding author to obtain the data and place it in the corresponding folder:
DATA_DIR = 'Z:\Users\TDR\Data\';
