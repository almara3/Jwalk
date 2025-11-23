from setuptools import setup, find_packages
import pathlib

HERE = pathlib.Path(__file__).parent

setup(
    name="Jwalk",
    version="1.0.0",
    author="Josh Bullock",
    author_email="j.bullock@cryst.bbk.ac.uk",
    description="A tool to calculate SASDs between crosslinked residues",
    packages=find_packages("src"),
    package_dir={"": "src"},
    include_package_data=True,
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.6",
        "biopython>=1.5",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    scripts=['src/Jwalk/jwalk'],
)
