from setuptools import setup

setup(
    name="nuc-codes",
    author="Matthew Wells",
    author_email="matthew.wells@canada.ca",
    entry_points={"console_scripts": ["nuc-codes=nuc_codes.nuc-codes:control_gff_finder"]},
)
