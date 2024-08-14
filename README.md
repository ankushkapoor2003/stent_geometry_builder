# Stent Geometry Builders for IR1 and HS1 Designs

## Overview

This repository contains the geometry builders for the IR1 and HS1 stent designs. These tools allow for the automated generation and customization of stent geometries based on specific input parameters.

## Setup Instructions

### Prerequisites

The geometry builder requires the following software to be installed. 

- **MATLAB2023b** 
- **Anaconda for Python integration**
-  **Solidworks 2024**

The builder is tested only on the software versions listed above.

### Cloning the Repository

Clone this repository to your local machine:

```bash
git clone https://github.com/yourusername/stent-geometry-builder.git
cd stent-geometry-builder
```
### Setting up the Python environment

1. 
### Setting Up IR1 Geometry Builder

1. Navigate to the `ir1` directory:
    ```bash
    cd ir1
    ```
2. Ensure all the required MATLAB functions are present in the `ir1` folder.

3. In MATLAB, add the `ir1` folder to your MATLAB path:
    ```matlab
    addpath(genpath('path_to_your_cloned_repo/ir1'));
    ```

4. Call the main geometry builder script:
    ```matlab
    main_ir1(input_parameters);
    ```

### Setting Up HS1 Geometry Builder

1. Navigate to the `hs1` directory:
    ```bash
    cd hs1
    ```
2. Ensure all the required MATLAB functions are present in the `hs1` folder.

3. In MATLAB, add the `hs1` folder to your MATLAB path:
    ```matlab
    addpath(genpath('path_to_your_cloned_repo/hs1'));
    ```

4. Call the main geometry builder script:
    ```matlab
    main_hs1(input_parameters);
    ```

## Running the Geometry Builders

### IR1 Geometry Builder

1. **Set up the input parameters** required for your stent design.
2. **Run the `main_ir1` script**:
    ```matlab
    main_ir1(input_parameters);
    ```
3. The output stent geometry will be saved in the `results` directory, which is automatically created in the current working directory.

### HS1 Geometry Builder

1. **Set up the input parameters** required for your stent design.
2. **Run the `main_hs1` script**:
    ```matlab
    main_hs1(input_parameters);
    ```
3. The output stent geometry will be saved in the `results` directory, which is automatically created in the current working directory.

## Citing This Work

If you use these geometry builders in your research, please cite the following paper:

```
@article{YourPaper,
    author = {Your Name and Coauthors},
    title = {Title of Your Paper},
    journal = {Journal Name},
    year = {2024},
    volume = {X},
    number = {X},
    pages = {X--X},
    doi = {DOI of the Paper}
}
```

You can also cite the software using the Zenodo DOI:

```
@software{YourSoftware,
    author = {Your Name and Coauthors},
    title = {Stent Geometry Builders for IR1 and HS1 Designs},
    year = {2024},
    doi = {DOI of the Zenodo Release}
}
```
