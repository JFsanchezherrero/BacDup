import os
import shutil
import sys
import glob
from setuptools import setup, find_packages, Extension

long_description_text = ""
with open("README.md", "r") as fh:
    long_description_text = fh.read()

#######
def get_require_modules():
    """
    Get main python requirements modules
    """
    with open("./BacDup/config/python_requirements.txt", 'r') as f:
        myModules = [line.strip().split(',')[0] for line in f]
    
    return myModules

#######
def get_version():
    """
    Original code: PhiSpy setup.py 
    https://github.com/linsalrob/PhiSpy/blob/master/setup.py
    """
    with open("./VERSION", 'r') as f:
        v = f.readline().strip()
    return v


setup(
    name="BacDup",
    version=get_version(),

    scripts=glob.glob('main/*'),

    ## TODO: check add several authors
    author="Jose F. Sanchez-Herrero and Alba Moya Garces",

    author_email="jfbioinformatics@gmail.com",
    description="Bacterial gene duplication analysis pipeline",

    long_description_content_type="text/markdown",
    long_description=long_description_text,
    url="https://github.com/JFsanchezherrero/BacDup",
    
    install_requires=get_require_modules(),
    packages=find_packages(),
    
    license='MIT License',

    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    include_package_data=True,
)




