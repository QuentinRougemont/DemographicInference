#lines to install dadi and dependancies using conda env
mkdir dadi_test
conda create --name dadi_test
conda install numpy
conda install matplotlib=1.3.1
conda install scipy-0.13.3-cp27-cp27mu-linux_x86_64.whl
python setup.py install
python -c 'from scipy import stats'
python -c 'import matplotlib'
cd 01-scripts
python -c 'import dadi'
source activate dadi_test 
