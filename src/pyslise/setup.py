import sys
import setuptools
from setuptools.dist import Distribution

assert '--library' in sys.argv, "--library option has to be set"
index = sys.argv.index('--library')
sys.argv.pop(index)
library = sys.argv.pop(index)
print("Using: '%s'" % library)

move = None
if '--move' in sys.argv:
    index = sys.argv.index('--move')
    sys.argv.pop(index)
    move = sys.argv.pop(index)

long_description = """
    Pyslise
"""


class BinaryDistribution(Distribution):
    def is_pure(self):
        return False


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
    include_package_data=True,
    distclass=BinaryDistribution,
    classifiers=[]
)

if move:
    from glob import glob
    import shutil
    import os

    for w in glob('dist/*.whl'):
        to = os.path.join(move, w[5:])
#        shutil.move(w, to)
        print("Moved wheel '%s' to '%s'" % (w, to))
