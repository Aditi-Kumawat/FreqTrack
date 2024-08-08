# FreqTrack
## Overview

This repository contains all the MATLAB codes used to generate the results for the frequency-dependent analytical model of ballasted rail-track systems subjected to moving loads. The model description and its findings are detailed in the paper, [Frequency-Dependent Analytical Model for Ballasted Rail-Track Systems Subjected to Moving Load](https://ascelibrary.org/doi/10.1061/%28ASCE%29GM.1943-5622.0001358).

## Model Description
The frequency-dependent analytical model proposed in this study incorporates essential substructural components, such as sleepers, rail pads, ballast, subballast, and subgrade. The model considers the effects of these components on the dynamic response of the rail-track system, including large deflections induced at the isolation frequency of the sleeper-substructure system. The model also incorporates damping effects, which significantly reduce rail-beam deflections without altering the critical velocity. The model is computationally efficient, making it suitable for engineering practice.
 <img width="875" alt="Screenshot 2024-08-08 at 16 31 58" src="https://github.com/user-attachments/assets/49790bbb-787d-44e1-bf2e-c9cfa4ac8cee">
Fig. 1 Definition sketch of the cross-section of railway formation 
<img width="746" alt="Screenshot 2024-08-08 at 16 31 46" src="https://github.com/user-attachments/assets/b1f0aa50-d1d3-486b-920c-78804ad012ea">
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

## References
- [Kumawat, A., Raychowdhury, P. and Chandra, S., 2019. Frequency-dependent analytical model for ballasted rail-track systems subjected to moving load. International Journal of Geomechanics, 19(4), p.04019016.](https://ascelibrary.org/doi/10.1061/%28ASCE%29GM.1943-5622.0001358)
- [Madshus, C. 5., and A. M. Kaynia. 2000. “High-speed railway lines on soft ground: Dynamic behaviour at critical train speed.” J. Sound Vib. 231 (3): 689–701.](https://doi.org/10.1006/jsvi.1999.2647)
- [Santos, N. C., A. Colaço, P. A. Costa, and R. Calçada. 2016. “Experimental analysis of track-ground vibrations on a stretch of the Portuguese railway network.” Soil Dyn. Earthquake Eng. 90: 358–380.]([https://doi.org/10.1016/j.soildyn.2017.03.004](https://doi.org/10.1016/j.soildyn.2016.09.003.)) 
- [Santos, N. C., J. Barbosa, R. Calçada, and R. Delgado. 2017. “Track-ground vibrations induced by railway traffic: Experimental validation of a 3D numerical model.” Soil Dyn. Earthquake Eng. 97: 324–344.](https://doi.org/10.1016/j.soildyn.2017.03.004) 

## Contributing
Contributions are welcome. Please fork the repository and create a pull request with your changes.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Contact
For any questions or issues, please contact Aditi Kumawat: aditikumawat@tum.de
