=======
PyPSM
=====
Climate Proxy System Models in Python

Introduction
---------------------------
This python package does a bunch of cool stuff with climate proxies.

Dependencies
---------------------------
numpy (http://www.numpy.org/) 
scipy (http://www.scipy.org/) 
rpy2 (http://rpy.sourceforge.net/) (For BCHRON)

Optional
  matplotlib (http://matplotlib.org/) (For plotting tools)

Installation
---------------------------
Make sure the dependencies are installed, then download and unzip this package, and then:
 python setup.py install

Alternately, you can use pip:
 pip install git+https://github.com/sylvia-dee/PRYSM.git

Either method will add a module named 'psm' to your default lib/python2.7/site-packages/ directory.

If you lack root access:
 python setup.py install --user

Testing
---------------------------
From the examples/ directory, run each of the example driver scripts:
./cellulose_driver.py
./speleo_driver.py
./icecore_driver.py
./coral_driver.py

This will create output files:
XXXXX

To plot (requires matplotlib):
XXXXX

