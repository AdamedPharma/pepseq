from setuptools import setup

setup(
    name="pepseq",
    version="1.0",
    description="Module to Read And Write Modified Peptide Repsesentations",
    author="MS",
    author_email="omitted_for_now",
    packages=["pepseq"],  # same as name
    install_requires=["networkx", "rdkit"],
)