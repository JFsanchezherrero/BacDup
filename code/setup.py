import setuptools
import glob

long_description_text = ""
with open("README.md", "r") as fh:
    long_description_text = fh.read()

setuptools.setup(
    name="BacDup",
    version="0.1",

    scripts=glob.glob('main/*'),
    ## TODO: check add several authors
    author="Jose F. Sanchez-Herrero and Alba Moya Garces",

    author_email="jfbioinformatics@gmail.com",
    description="Bacterial gene duplication analysis pipeline",

    long_description_content_type="text/markdown",
    long_description=long_description_text,
    url="https://github.com/JFsanchezherrero/BacDup",
    packages=setuptools.find_packages(),
    license='MIT License',

    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    include_package_data=True,

    ## TODO: check
    install_requires=[ 
        'pandas', 'patool', 'termcolor', 
	#'pysam', 'pybedtools', 'multiqc',
	'bcbio-gff', 
	'biopython', 'HCGB', 
	'ftputil', 'numpy', 'python-dateutil', 'pytz', 'six', 'wget'
    ],
)




