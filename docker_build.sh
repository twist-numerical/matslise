#!/bin/bash

set -e

cmake_python=/opt/python/cp312-cp312
# cmake
${cmake_python}/bin/pip install cmake

# pyslise
pythons=(cp39-cp39 cp310-cp310 cp311-cp311 cp312-cp312 cp313-cp313 cp314-cp314)
for py in ${pythons[@]}; do
    mkdir /opt/matslise-build-${py}
    cd /opt/matslise-build-${py}
    /opt/python/${py}/bin/pip install setuptools
    ${cmake_python}/bin/cmake /opt/matslise/ \
        -DCMAKE_BUILD_TYPE=Release \
        -DPYTHON_EXECUTABLE=/opt/python/${py}/bin/python

    ${cmake_python}/bin/cmake --build . --target build_wheel --config Release -- -j 6

    auditwheel repair matslise/dist/*${py}*.whl -w /opt/matslise/wheelhouse
done