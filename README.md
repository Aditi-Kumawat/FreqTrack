# FreqTrack
## Overview

This repository contains all the MATLAB codes used to generate the results for the frequency-dependent analytical model of ballasted rail-track systems subjected to moving loads. The model description and its findings are detailed in the paper, "Frequency-Dependent Analytical Model for Ballasted Rail-Track Systems Subjected to Moving Load."

## Model Description
The frequency-dependent analytical model proposed in this study incorporates essential substructural components, such as sleepers, rail pads, ballast, subballast, and subgrade. The model considers the effects of these components on the dynamic response of the rail-track system, including large deflections induced at the isolation frequency of the sleeper-substructure system. The model also incorporates damping effects, which significantly reduce rail-beam deflections without altering the critical velocity. The model is computationally efficient, making it suitable for engineering practice.
![Alt text](<Screenshot 2024-08-08 at 16.31.58.png>) 
Fig. 1 Definition sketch of the cross-section of railway formation 
![Alt text](<Screenshot 2024-08-08 at 16.31.46.png>)
Fig. 2 Definition sketch of the longitudinal-section of railway formation
## Validation
The proposed model's reliability was assessed by comparing its predictions with experimental data (Madshus and Kaynia, 2000; Santos et al. 2016, 2017). The comparison demonstrated that the model could accurately predict rail-beam deflections under dynamic loading conditions, making it a valuable tool for railway track design and analysis.

## Source Code
The script files given in the folder `src/` are described here briefly:
- `sc_wTimHist_*.m`: Time history plots shown in Figures 8 to 12
- `sc_contour_*.m`: generates contour plots shown in Figure 14
- `DATA_INPT/`: Directory containing input data files
- `Santos_Et_AL/`: validation codes for Santos et al. 2016, 2017 study
- `sc_w_multi_*.m`: validation codes for Madshus and Kaynia 2000

## Contributing
Contributions are welcome. Please fork the repository and create a pull request with your changes.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Contact
For any questions or issues, please contact Aditi Kumawat: aditikumawat@tum.de
