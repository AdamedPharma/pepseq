from setuptools import setup, find_packages

setup(
    name="pepseq",
    version="1.0",
    description="Module to Read And Write Modified Peptide Repsesentations",
    author="MS",
    packages=(
        find_packages(where="..") + 
        find_packages(where=".")
    ),
    package_dir={"": "..", "pepseq": "./pepseq"},
    package_data={'':['pepseq/Peptide/database/db.json','Peptide/database/db.json']},
    include_package_data=True,
    author_email="omitted_for_now",
    install_requires=[
        "contourpy==1.1.1",
        "cycler==0.12.1",
        "exceptiongroup==1.1.3",
        "fonttools==4.43.1",
        "iniconfig==2.0.0",
        "kiwisolver==1.4.5",
        "matplotlib==3.8.0",
        "networkx==3.1",
        "numpy==1.26.0",
        "packaging==23.2",
        "pandas==2.1.1",
        "Pillow==10.0.1",
        "pluggy==1.3.0",
        "pycairo==1.25.0",
        "pyparsing==3.1.1",
        "pytest==7.4.2",
        "python-dateutil==2.8.2",
        "pytz==2023.3.post1",
        "rdkit==2022.3.4",
        "six==1.16.0",
        "tomli==2.0.1",
        "tqdm==4.66.1",
        "tzdata==2023.3"
        ],
)
