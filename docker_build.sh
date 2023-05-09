#!/bin/bash

# cmake
/opt/python/cp37-cp37m/bin/pip install cmake

# pyslise
pythons=(cp36-cp36m cp37-cp37m cp38-cp38 cp39-cp39 cp310-cp310)
for py in ${pythons[@]}; do
    mkdir /opt/matslise-build-${py}
    cd /opt/matslise-build-${py}
    /opt/python/cp37-cp37m/bin/cmake /opt/matslise/ \
        -DCMAKE_BUILD_TYPE=Release \
        -DPYTHON_EXECUTABLE=/opt/python/${py}/bin/python

    /opt/python/cp37-cp37m/bin/cmake --build . --target build_wheel --config Release -- -j 6

    auditwheel repair matslise/dist/*${py}*.whl --plat manylinux2010_x86_64 -w /opt/matslise/wheelhouse
done