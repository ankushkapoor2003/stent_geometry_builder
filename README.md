# Stent Geometry Builders for Independent Ring (IR) and Helical Stent (HS) Designs

## Citing This Work

If you use these geometry builders in your research, please cite the paper and the software.

Research paper citations:
```
@article{Kapoor2025Surrogate,
    author = {Kapoor, Ankush and Ray, Tapabrata and Jepson, Nigel and Beier, Susann},
    title = "{A Surrogate-Assisted Multi-Concept Optimization Framework for Real-World Engineering Design}",
    journal = {Journal of Mechanical Design - Accepted 03/2025 (In Print)}
}
```

```
@article{Kapoor2024Comprehensive,
    author = {Kapoor, Ankush and Ray, Tapabrata and Jepson, Nigel and Beier, Susann},
    title = "{Comprehensive Geometric Parameterization and Computationally Efficient 3D Shape Matching Optimization of Realistic Stents}",
    journal = {Journal of Mechanical Design},
    pages = {1-24},
    year = {2024},
    month = {10},
    issn = {1050-0472},
    doi = {10.1115/1.4066961},
    url = {https://doi.org/10.1115/1.4066961},
    eprint = {https://asmedigitalcollection.asme.org/mechanicaldesign/article-pdf/doi/10.1115/1.4066961/7390188/md-24-1477.pdf},
}
```
Software citation:
```
@software{kapoor20243D,
  author       = {Kapoor, Ankush and
                  Ray, Tapabrata and
                  Jepson, Nigel and
                  Beier, Susann},
  title        = {{3D geometry builder for Independent Ring and 
                   Helical stent designs}},
  month        = oct,
  year         = 2025,
  publisher    = {Zenodo},
  version      = {0.0.2},
  doi          = {10.5281/zenodo.15117426},
  url          = {https://doi.org/10.5281/zenodo.15117426}
}
```

## Overview

This repository contains the geometry builders for the IR1, IR2, HS1, and HS2 stent designs. These tools allow for the automated generation and customization of stent geometries based on specific input parameters.

## Setup Instructions

### Prerequisites

The geometry builder requires the following software to be installed. 

- **MATLAB2023b** 
- **Anaconda**
-  **Solidworks 2024**

The builder is tested only on the software versions listed above.

### Cloning the Repository

Clone this repository to your local machine:

```bash
git clone https://github.com/yourusername/stent-geometry-builder.git
cd stent-geometry-builder
```
### Setting up the solidworks integration Python environment and corresponding modules. 

1. Open the anaconda prompt and navigate to the cloned directory location using cd command
    ```conda
    cd '\\ stent geometry builder location \\'
    ```

2. Create a new conda environent with Python 3.9.18 
    ```conda
    conda create --name Sldwrks_Integration python=3.9.18
    ```

3. Activate the created conda environemnt
    ```conda
    conda activate Sldwrks_Integration
    ```
4. In the anaconda prompt, navigate to the cloned directory location using cd command
    ```conda
    cd '\\ stent geometry builder location \\'
    ```

5. Install the sld_interface package and its dependencies in the Sldwrks_integration environment
    ```conda
    python -m pip install ".\Solidworks_library\sld_interface"
    ```

6. Retrieve the location of the python executable of the conda environement using the following command
    ```conda
    where python
    ```
    The command may show multiple python locations, copy the path within your generated conda environment

7. Open Matlab and change the matlab python environment to the newly generated one by running the following in command window
    ```matlab
    pyenv('Version', '\\ Path to the python executable retrieved from step 5 \\')
    ```

NB: 

Please note that the python library assumes that solidworks is installed in the default location i.e. "C:\Program Files\SOLIDWORKS Corp\SOLIDWORKS\SLDWORKS.exe". if the solidworks executable is at a different location, change the location in the solidworks library file located at "\stent geometry builder location\Solidworks_Library\sld_interface\sld_interface\sld.py", line 11. 

Similarly, the default solidworks template used for loading new document is assumed at the location "C:\ProgramData\SolidWorks\SOLIDWORKS 2024\templates\Part.prtdot". Kindly update line 478 in sld.py if the location is different. 

## Using Geometry builder to create stent designs

1. Navigate to the cloned geometry builder directory in windows explorer:

2. Open the corresponding matlab file for IR stent (Main_IR1) or Helical stent (Main_HS1)

3. Either keep or modify the baseline design variable values (x) in the matlab script

4. Run the matlab script

5. Three stent CAD files - Solidworks part file (.SLDPRT), neutral parasolid file (.x_t) and a stl file (.STL) will be saved in the corresponding results folder for IR stent (IR1_Results) or HS stent (HS1_Results)


