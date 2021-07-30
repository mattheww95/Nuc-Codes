from setuptools import setup

setup(
    name="nuc_codes",
    author="Matthew Wells",
    author_email="matthew.wells@canada.ca",
    entry_points={"console_scripts": ["nuc_codes=nuc_codes.nuc_codes:control_gff_finder"]},
)
