# KSCN in Electrolytes Compensates K Ions
### Introduction
Scripts from the manuscript W. Wang, et al. Compensating K Ions Through an Organic Salt in Electrolytes for Practical K-Ion Batteries.

### Contents
- **Data folder** includes a [molecules.xlsx](./Data/molecules.xlsx), which is a input file for hierarchical clustering.
- **Output folder** contains two tif images and two excel files, which are the data generated after running [hierarchical_clustering.py](./hierarchical_clustering.py). The excel file contains the results of the calculation of the molecular descriptors, while the tiff file is the dendrogram obtained from the hierarchical clustering analysis based on these descriptors. ''AL'' and ''SE'' in the file names indicate ''anode limit'' and ''solubility energy'', respectively.
- [hierarchical_clustering.py](./hierarchical_clustering.py), a code for hierarchical clustering.

### System Requirements
In order to run source code file in the Data folder, the following requirements need to be met:
- Windows, Mac , Linux
- Python (version 3.7.4 or newer)
- Python modules: numpy, pandas, rdkit, matplotlib, sklearn and scipy. Most version should be fine.

### Installation guide
You can download the package in zip format directly from this github site, with an estimated download time of less than 1 minute (regular office internet speed).

### Instructions for Use
- Make sure that there is a ''molecules.xlsx'' file in the Data folder, run [hierarchical_clustering.py](./hierarchical_clustering.py), and the Output folder will generate [''cluster_AL.tif''](./Output/cluster_AL.tif), [''cluster_SE.tif''](./Output/cluster_SE.tif), [''data_AL.xlsx''](./Output/data_AL.xlsx) and [''data_SE.xlsx''](./Output/data_SE.xlsx) files.

### Contributions
Y. Gao, W. Wang and G. Wu developed a workflow. G. Wu wrote the program. Y. Gao and W. Wang contributed to code debugging.

### License
This project uses the [MIT LICENSE](LICENSE).

### Disclaimer
This code is intended for educational and research purposes only. Please ensure that you comply with relevant laws and regulations as well as the terms of service of the target website when using this code. The author is not responsible for any legal liabilities or other issues arising from the use of this code.

### Contact
If you have any questions, you can contact us at: yuegao@fudan.edu.cn
