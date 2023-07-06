from setuptools import setup, find_packages

setup(
    name="pepseq",
    version="1.0",
    description="Module to Read And Write Modified Peptide Repsesentations",
    author="MS",
    packages=find_packages(),
    package_data={'':['pepseq/Peptide/database/db.json','Peptide/database/db.json']},
    include_package_data=True,
    author_email="omitted_for_now",
    install_requires=["networkx", "rdkit"],

)
