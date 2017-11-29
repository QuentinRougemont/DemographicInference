cd $HOME
cd DemographicInference
module load gcc/4.8.5
module load python/2.7.13
mkdir temp
virtualenv temp
source temp/bin/activate
pip install numpy
pip install -Iv https://pypi.python.org/packages/2f/12/565c08132db50a0ba34a33e0901f3d1d4d72e3b432ea828e4d87be5a4991/scipy-0.13.3.tar.gz

python setup.py install 
