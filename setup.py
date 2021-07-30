from setuptools import setup, find_packages

setup(
    name="nuc_codes",
    author="Matthew Wells",
    packages=find_packages(),
    author_email="matthew.wells@canada.ca",
    entry_points={"console_scripts": ["nuc_codes=nuc_codes.nuc_codes:control_gff_finder"]},
)
