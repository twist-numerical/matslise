import glob
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyslise",
    version="0.0.1",
    author="Toon Baeyens",
    author_email="toon.baeyens@ugent.be",
    description="Python bindings for the C++ version of Matslise",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.ugent.be/tobaeyen/matslise-cpp",
    packages=['pyslise'],
    package_data={
        '': [f[len('pyslise/'):] for f in glob.glob('pyslise/pyslise.*')]
    },
    classifiers=[]
)
