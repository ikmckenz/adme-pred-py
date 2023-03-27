"""setup.py file for packaging purposes"""
import setuptools

with open("README.md", mode="r", encoding="ascii") as fh:
    long_description = fh.read()

setuptools.setup(
    name="adme_pred",
    version="0.0.2",
    author="ikmckenz",
    description="ADME/Tox Models in Python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ikmckenz/adme-pred-py",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)"
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Medical Science Apps."
    ],
)
