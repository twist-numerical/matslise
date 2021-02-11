#!/bin/bash

# cmake
/opt/python/cp37-cp37m/bin/pip install cmake

# Eigen3
git clone https://gitlab.com/libeigen/eigen.git /opt/eigen
mkdir /opt/eigen/build
cd /opt/eigen/build
git checkout 3.3
/opt/python/cp37-cp37m/bin/cmake ..

# pyslise
pythons=(cp27-cp27mu cp36-cp36m cp37-cp37m cp38-cp38 cp39-cp39)
for py in ${pythons[@]}; do
    mkdir /opt/matslise-build-${py}
    cd /opt/matslise-build-${py}
    /opt/python/cp37-cp37m/bin/cmake /opt/matslise/ \
        -DCMAKE_BUILD_TYPE=Release \
        -DEigen3_DIR=/opt/eigen/build \
        -DPYTHON_EXECUTABLE=/opt/python/${py}/bin/python

    /opt/python/cp37-cp37m/bin/cmake --build . --target build_wheel --config Release -- -j 6

    auditwheel repair src/dist/*${py}*.whl --plat manylinux2010_x86_64 -w /opt/matslise/wheelhouse
done