import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="adme-pred-py",
    version="0.0.1",
    author="ikmckenz",
    description="ADME/Tox Models in Python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ikmckenz/adme-pred-py",
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)"
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Medical Science Apps."
    ),
)