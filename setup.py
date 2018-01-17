# The template for this setup.py came from Tim Morton,
#  who I understand took it from Dan F-M. Thanks guys!

# import our ingredients
from setuptools import setup, find_packages
import os,sys

# return the README as a string
def readme():
    with open('README.md') as f:
        return f.read()

# a little kludge to be able to get the version number from the packa
import sys
if sys.version_info[0] < 3:
    import __builtin__ as builtins
else:
    import builtins
builtins.__CRAFTROOM_SETUP__ = True
import craftroom
version = craftroom.__version__

setup(name = "craftroom",
    version = version,
    description = "Python code for astronomy (basically the equivalent of a hot glue gun, a sewing machine, and a pair of knitting needles).",
    long_description = readme(),
    author = "Zachory K. Berta-Thompson",
    author_email = "zach.bertathompson@colorado.edu",
    url = "https://github.com/zkbt/craftroom",
    packages = find_packages(),
    package_data = {'craftroom':[]},
    include_package_data=False,
    scripts = [],
    classifiers=[
      'Intended Audience :: Science/Research',
      'Programming Language :: Python',
      'Topic :: Scientific/Engineering :: Astronomy'
      ],
    install_requires=['numpy', 'astropy', 'scipy', 'matplotlib'],
    zip_safe=False,
    license='MIT',
)
