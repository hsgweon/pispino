from setuptools import setup, find_packages
import os

__version__ = os.environ.get("VERSION", "1.0")

setup(
    name = "pispino",
    version = __version__,
    packages = ["pispino"],
    scripts = ["bin/pispino_createreadpairslist", "bin/pispino_seqprep"],
    description = "PISPINO (PIpits SPIN-Off tools): Bioinformatics toolkits for processing NGS data.",
    long_description = "An open source stand-alone suite of tools for processing NGS data. They were part of PIPITS, but subsequently became independent tools as they are extensively used for other pipelines developed by Gweon Lab.",
    author = "Hyun Soon Gweon",
    author_email = "h.s.gweon@reading.ac.uk",
	url = "https://github.com/hsgweon/pispino",
	licence = "GPA"
)
