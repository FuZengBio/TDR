**README: Dimensionality-reduction analyses for VIP, MSTd, and PIVC datasets**

This folder contains MATLAB scripts for running time-varying targeted dimensionality reduction (**tTDR**), **original TDR**, **time-course** variance significance, comparison analyses **across monkeys**, **hemispheres** as well as complementary dimensionality-reduction analyses(**mTDR** and **dPCA**).
_____________________________________________________________________________________________________________________________________________________
**1. Time-varying targeted dimensionality reduction (tTDR)**

Notes:

(1) To obtain tTDR results for VIP and MSTd in the 1 s condition, run: **tTDTR_VIP_MSTd_1s**.

(2) To obtain tTDR results for VIP and MSTd in the 2 s condition, together with PIVC in the 1 s condition, run: **tTDR_VIP_MSTd_2s_PIVC_1s**.

To change the orthogonalization order in the tTDR scripts, set: reverse = 0; or change it to: reverse = 1; or: reverse = 2;

The meaning of each option is:

Default: reverse = 0: orthogonalization order = velocity → acceleration → choice.

reverse = 1: orthogonalization order = choice → acceleration → velocity.

reverse = 2: move acceleration to the end, i.e. velocity → choice → acceleration.

_____________________________________________________________________________________________________________________________________________________
**2. Original TDR cue-heading targeted dimensionality reduction**

To obtain original TDR cue-heading targeted dimensionality reduction results, run: **originalTDR_analysis_CueHeading**.

Regressor order for original cue-heading TDR: 

Default: reverse = 0, regressor order = cue → heading

_____________________________________________________________________________________________________________________________________________________
**3. Original TDR heading-choice targeted dimensionality reduction**

To obtain original TDR heading-choice targeted dimensionality reduction results, run: **originalTDR_analysis_HeadingChoice**.

Regressor order for heading and choice:

reverse = 0, regressor order = heading → choice

reverse = 1, regressor order = choice → heading

_____________________________________________________________________________________________________________________________________________________
**4. Time-course variance, significanceg**

(1) tTDR (i.e. along the velocity, acceleration, and choice axes), run: **tTDTR_timecourse_variance_sig**.

To change the orthogonalization order, modify **idMode** in the script.

VAC: velocity → acceleration → choice

VCA: velocity → choice → acceleration

CAV: choice → acceleration → velocity

(2) Original TDR (i.e. along the heading and choice axes), run: **originalTDTR_HC_timecourse_variance_sig**

To change the regressor order, modify both reverse and orderStr in the script.

reverse = 0; orderStr = 'HC';: heading → choice

reversed setting: choice → heading

reverse and orderStr should be changed together to keep the computation order and output labels consistent.

_____________________________________________________________________________________________________________________________________________________
**5. Vestibular–visual choice-axis angle before orthogonalization in tTDR**

To obtain the angle between the vestibular and visual choice axes before orthogonalization in tTDR, run: **tTDR_choiceAngle_Vis_vs_Vest**.

_____________________________________________________________________________________________________________________________________________________
**6. Separate-monkey and left–right hemisphere analyses in tTDR**

(1) Differences between monkeys within the same brain area

To examine differences between monkeys within the same brain area in tTDR, run:

For the 1 s condition: **tTDTR_VIP_MSTd_1s_separate_monkey**.

For the 2 s condition: **tTDTR_VIP_MSTd_2s_separate_monkey**.


(2) Differences between left and right hemispheres

To examine differences between the left and right hemispheres in tTDR, run:

For the 1 s condition: **tTDTR_VIP_MSTd_1s_LR_hemisphere**

For the 2 s condition: **tTDTR_VIP_MSTd_2s_LR_hemisphere**


Changing the orthogonalization order

All of the scripts above use the same reverse setting:

reverse = 0;   % 0: VAC; 1: reverse parameters; 2: move acceleration to the end

_____________________________________________________________________________________________________________________________________________________
**7. Model-based targeted dimensionality reduction (mTDR) and comparison with original TDR**

To obtain mTDR results and compare them with original TDR in the heading–choice framework, run: **mTDR_and_compare_originalTDR_HC**

To change the regressor order, modify both reverse and orderStr in the script.

reverse = 0; orderStr = 'HC';: heading → choice

reversed setting: choice → heading

reverse and orderStr should be changed together to keep the computation order and output labels consistent.

_____________________________________________________________________________________________________________________________________________________
**8. Demixed principal component analysis (dPCA)**

To perform dPCA, please refer to the original dPCA resources from the Machens lab.

The main analysis scripts (MATLAB) are available at: **http://github.com/machenslab/elife2016dpca**

The dPCA toolkit is located at: **https://github.com/machenslab/dPCA**

These resources provide the original implementation and analysis framework for dPCA, which can be adapted for the present datasets.
_____________________________________________________________________________________________________________________________________________________
**Data Availability Note**

Due to upload file size limits, the data files cannot be provided directly.
If needed, please contact the corresponding author to obtain the data and place it in the corresponding folder:

DATA_DIR = 'Z:\Users\TDR\Data\';

