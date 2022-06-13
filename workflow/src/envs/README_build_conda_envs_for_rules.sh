## Download mamba installer
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh


## Install mamba
bash Mambaforge-Linux-x86_64.sh
#NOTE: specified destination folder as /lab/work/ccrober/mambaforge


## Create environment for Snakemake
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
conda env export > workflow/envs/snakemake.yaml


## Create environment for dropkick
mamba create --name dropkick
mamba activate dropkick
mamba install python=3.9
mamba install scanpy
mamba install scikit-image
mamba install scikit-learn
cd /lab/work/ccrober/sw/dropkick
python setup.py install
cd /lab/work/ccrober/T1D_U01/workflow/src/envs/
conda env export > dropkick.yaml


## Create environment for R
mamba create --name Renv
conda env export > /lab/work/ccrober/T1D_U01/workflow/envs/Renv.yaml


## Create environment for AMULET
mamba create --name amulet
mamba activate amulet
mamba install python=3.9
mamba install numpy pandas scipy statsmodels
cd ~/bin
wget https://github.com/UcarLab/AMULET/releases/download/v1.1/AMULET-v1.1.zip
unzip AMULET-v1.1.zip
chmod +x AMULET.sh


## Allow environments to be accessed in ipython notebooks
#Following method 2 from this tutorial: https://towardsdatascience.com/get-your-conda-environment-to-show-in-jupyter-notebooks-the-easy-way-17010b76e874
mamba activate dropkick
mamba install ipykernel
ipython kernel install --name=dropkick

#Repeat above for any mamba environment that you
#want to be available as an ipython notebook kernel
