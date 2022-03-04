from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
import os
import shutil

assert 'PYSLISE_LIBRARY' in os.environ, "pyslise_library environment variable has to be set"
library = os.environ['PYSLISE_LIBRARY']
print("Using: '%s'" % library)

with open('description.md', 'rb') as f:
    long_description = f.read().decode('UTF-8')

class CMakeExtension(Extension):
    def __init__(self, name):
        Extension.__init__(self, name, sources=[])


class CMakeBuild(build_ext):
    def run(self):
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        path = self.get_ext_fullpath(ext.name)
        try:
            os.makedirs(os.path.dirname(path))
        except:
            pass
        shutil.copy(library, path)


setup(
    name="pyslise",
    version="${PYSLISE_VERSION}",
    author="Toon Baeyens",
    author_email="toon.baeyens@ugent.be",
    description="Python bindings for the C++ version of Matslise",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://matslise.ugent.be/",
    ext_modules=[CMakeExtension('pyslise')],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    install_requires=[] # 'numpy']
)
