# CANUCS-Sparkles


## Purpose

Analysis of the Red Sparkler

## Installation

### Python setup

I like using a separate Python environment for each project, and encourage you to do that too. If you want to roll that way, make sure you have [miniforge](https://github.com/conda-forge/miniforge) installed. Once that is done, use `conda` to create the virtual environment and install the needed Python modules as follows:

```
conda create -n sparkly python=3.8 ipython astropy photutils jupyterlab pandas
```

To use the environment:

```
conda activate sparkly
```

When finished:

```
conda deactivate
```

### Git LFS

The version of git being used needs to have support for the Git Large File Storage (LFS) extension. We want to include some giant FITS images, and without LFS to manage the storage this isn't practical. With Git LFS, big images are no problem. I suspect that people using this repo are nerdly enough to have Git LFS already activated, but just in case you don't:

On a Mac, use Homebrew: 

```
brew install git-lfs
``` 

On Windows (using WSL) or Linux:

```
sudo apt-get update -y
sudo apt-get install -y git-lfs
```

After LFS support is added to Git, you need to activate it in the repo with:

```
git lfs install
git lfs track "*.fits"
```

If you are just cloning the repo and have already installed LFS support, I don't think you need to run any of these commands (though it won't hurt).

A useful Git LFS tutorial is here:

https://www.atlassian.com/git/tutorials/git-lfs


See information below about using the CANFAR VOS.

### FSPS

If you want to model stellar populations, you also need to install FSPS:

Firstly, make sure you have gfortran installed. (Use homebrew on a mac, or apt-get on a Linux machine
or a Windows machine using WSL).

Install the command-line version:

```
mkdir -p ~/tools
cd ~/tools
git clone https://github.com/cconroy20/fsps
echo 'export SPS_HOME="$HOME/tools/fsps/"' >> ~/.zshrc
source ~/.zshrc
cd fsps/src
make
```

Install the Python bindings using pip as follows:

```
python -m pip install fsps
```


## Images

Through the magic of Git LFS, the JWST data used by this repo are all pulled over (slowly) when you clone this repo (about 4 Gb worth of images).

If you want to pull over other images from CANFAR, you need to use the VOS tools. Hopefully, this repo already has all you need, but just in case, it's useful to know how to access additional data from our friends at the CADC.

Firstly, install the `VOS` command line tools. These are written in Python and are pip installable, so:

```
pip install vos
```

The VOS tools require a certificate. Use the CADC webpage to get the certificate (or use `getCert`). The certificate is named: `cadcproxy.pem`. The certificate should be copied to the directory: ~/.ssl  

The VOS offers commands like `vls`, `vcp`, etc. The most common use case is copying files from the CANFAR JWST storage to local storage. For example:

```
vcp "arc:projects/canucs/grizli-catalogs/SMACS0723-Mosaics/VeryRoughDraft/*_40mas.fits" .
```