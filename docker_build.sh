#!/bin/bash

# cmake
/opt/python/cp37-cp37m/bin/pip install cmake

# Eigen3
git clone https://github.com/eigenteam/eigen-git-mirror.git /opt/eigen
mkdir /opt/eigen/build
cd /opt/eigen/build
/opt/python/cp37-cp37m/bin/cmake ..
make

# pyslise
pythons=(cp27-cp27mu cp34-cp34m cp35-cp35m cp36-cp36m  cp37-cp37m)
for py in ${pythons[@]}; do
    mkdir /opt/matslise-build-${py}
    cd /opt/matslise-build-${py}
    /opt/python/cp37-cp37m/bin/cmake /opt/matslise/ \
        -DEigen3_DIR=/opt/eigen/build \
        -DLONG_DOUBLE=OFF \
        -DPYTHON_EXECUTABLE=/opt/python/${py}/bin/python \
        -DAUDITWHEEL_repair_plat=manylinux2010_x86_64

    /opt/python/cp37-cp37m/bin/cmake --build . --target buildWheel -- -j 6
done