import os

try:
    from numpy.distutils.core import setup
except Exception as ex:
    print(ex)
    print("StochPy requires NumPy\n")
    print("See http://numpy.scipy.org/ for more information about NumPy")
    os.sys.exit(-1)

local_path = os.path.dirname(os.path.abspath(os.sys.argv[0]))		# Get the dir of setup.py
os.chdir(local_path)

from distutils.core import setup
setup(name='crn_mc',
      version='1.0',
      description='code for my research, also to learn about writing python packages',
      author='Ethan Levien',
      author_email='levien@math.utah.edu',
      url='http://www.math.utah.edu/~levien/',
      packages=['crn_mc', 'crn_mc.simulation']
      )
