from setuptools import setup, find_packages

VERSION = '0.0.1'
DESCRIPTION = 'Python package for modelling stellar pulsations'

with open("README.md", "r") as f:
    long_description = f.read()

install_requires = [
    "matplotlib",
    "numpy",
    "pandas",
    "h5py"]

classifiers = [
    "Development Status :: 3 - Alpha",
    "Intended Audience :: Science/Research",
    "Programming Language :: Python :: 3",
    'Topic :: Scientific :: Astronomy']
# Setting up
setup(
        name="pypulse",
        version=VERSION,
        author="Mathias Michielsen",
        author_email="mathias.michielsen@kuleuven.be",
        description=DESCRIPTION,
        long_description=long_description,
        long_description_content_type='text/markdown',
        packages=find_packages(),
        python_requires='>=3.7',
        install_requires=install_requires,
        keywords=['python'],
        classifiers= classifiers,
        license='GPL3'
)
