# CANUCS-Sparkles


## Purpose

Analysis of the Red Sparkler

## Installation

### Python setup

Firstly, make sure you have [miniforge](https://github.com/conda-forge/miniforge) installed. Once that is done, use `conda` to create the virtual environment and install the needed Python modules as follows:

```
conda create -n sparkly python=3.8 ipython astropy photutils jupyterlab
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

The version of git being used needs to have support for git LFS as we want to access some giant FITS images and we don't want these stored in the GitHub repo, so we need to pull them in from other places.

One a Mac, use Homebrew: 

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

If you are just cloning the repo and have already installed LFS support, I don't think you need to run these commands (though it won't hurt).

A useful Git LFS tutorial is here:

https://www.atlassian.com/git/tutorials/git-lfs