● Preparation:  
Before establishing the scDock environment, the below system packages might be required. Please check whether these packages have been installed: build-essential, libcurl4-openssl-dev, libssl-dev, libxml2-dev, zlib1g-dev, libbz2-dev, liblzma-dev, libhdf5-dev, pkg-config, gfortran, libpcre2-dev, libgsl-dev, libglpk-dev  
For example, you can use `dpkg -l | grep build-essential` to clarify the exist of the corresponding package.  

conda is also required for scDock. Please make sure you have installed miniconda (https://www.anaconda.com/docs/getting-started/miniconda/main) or Anaconda (https://www.anaconda.com/) before operating the further process.  

Please make sure there are `conda-forge` and `bioconda` channels in your conda. The channel can be added by:  
`conda config –add channels conda-forge`  
`conda config –add channels bioconda`  
`conda config –set channel_priority strict`  

● Installation:  
To establish the virtual environment for scDock, we provide "install.sh" which can automatically create a virtual environment named "scDock" with the essential system requirment, Python (3.10.0) & R version (4.4.3), and install their modules/packages by conduct a single line command.

Please follow the instructions to install scDock:  
1. Download and place the "install.sh" to your computer.
2. Conduct `bash install.sh`

The installation process might take some time to finish.
After the installation, your computer will be ready for scDock.

We have confirm the success of scDock installation and operation in the following linux systems:  
1. Ubuntu 18.04.6
2. Ubuntu 20.04.6  
3. Ubuntu 22.04.4
4. Debian 12
