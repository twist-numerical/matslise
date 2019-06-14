import setuptools
import sys

assert '--library' in sys.argv, "--library option has to be set"
index = sys.argv.index('--library')
sys.argv.pop(index)  # Removes the '--foo'
library = sys.argv.pop(index)
print("Using: ", library)

long_description = """
    Pyslise
"""

setuptools.setup(
    name="pyslise",
    version="${PYSLISE_VERSION}",
    author="Toon Baeyens",
    author_email="toon.baeyens@ugent.be",
    description="Python bindings for the C++ version of Matslise",
    long_description=long_description,
    url="https://github.ugent.be/tobaeyen/matslise-cpp",
    packages=[''],
    package_data={
        '': [library]
    },
    classifiers=[]
)
