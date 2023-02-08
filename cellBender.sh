# Run aligned samples through CellBender to remove ambient RNA

# remove background gene expression
export PATH=/home/NEUROSURGERY.WISC.EDU/madrid/miniconda3/bin/:$PATH
conda create -n CellBender
source activate CellBender
conda install -c anaconda pytables
conda install pytorch torchvision -c pytorch
git clone https://github.com/broadinstitute/CellBender.git
pip install -e CellBender
cellbender remove-background --output output.h5 --expected-cells 10000 --total-droplets-included 20000 --fpr 0.01 --epochs 150 --input raw_feature_bc_matrix.h5
# must be done for each of the samples
