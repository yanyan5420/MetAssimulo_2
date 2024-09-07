# MetAssimulo 2: a web app for simulating realistic 1D & 2D Metabolomic 1H NMR spectra 

## Overview
This repository contains all the necessary code for the MetAssimulo 2 project. Please follow the instructions below to set up and run the project on your local machine.

## Prerequisites
Before you begin, ensure you have the following installed on your system:
- **Git**: Necessary for cloning the repository.
- **Anaconda or Miniconda**: Required to manage the project dependencies.

### Installing Git
If you do not have Git installed, please follow the instructions on the [Git official site](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) to install it on your operating system.

### Installing Anaconda
Anaconda is recommended for setting up and managing the project dependencies. Download and install it from the [official Anaconda website](https://www.anaconda.com/products/distribution). Alternatively, if you prefer a lighter installation, you may opt for Miniconda, available [here](https://docs.conda.io/en/latest/miniconda.html).

## Installation

### 1. Clone the Repository
```bash
git clone https://github.com/yanyan5420/MetAssimulo_2.git
cd MetAssimulo_2
```

### 2. Download the Input Data
```bash
pip install gdown
gdown https://drive.google.com/uc?id=12k-nbIcoME5GSkD8u1iYUcNzg6l-7msW
unzip Input.zip
```

### 3. Set Up the Environment
```
conda env create -f environment.yml
conda activate py37_metassimulo_2_env
```

### 4. Set Environment Variable and Run the Project
For Linux and macOS:
```
export PYTHONPATH=$HOME/MetAssimulo_2/:$PYTHONPATH
python3 apps/index.py -p Input/parameters.txt
```
For Windows:
```
set PYTHONPATH=%cd%;%PYTHONPATH%
python apps\index.py -p Input\parameters.txt
```

## Usage
After setting up the environment as described above, you can run the project by navigating to the project directory and running the application:
for Linux and macOS:
```
cd /path/to/MetAssimulo_2
conda activate py37_metassimulo_2_env
export PYTHONPATH=$HOME/MetAssimulo_2/:$PYTHONPATH
python3 apps/index.py -p Input/parameters.txt
```
for Windows:
```
cd /path/to/MetAssimulo_2
conda activate py37_metassimulo_2_env
set PYTHONPATH=%cd%;%PYTHONPATH%
python apps\index.py -p Input\parameters.txt
```
