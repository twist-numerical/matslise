#!/bin/bash

# cmake
/opt/python/cp37-cp37m/bin/pip install cmake

# pyslise
pythons=(cp39-cp39 cp310-cp310 cp311-cp311 cp312-cp312 cp313-cp313)
for py in ${pythons[@]}; do
    mkdir /opt/matslise-build-${py}
    cd /opt/matslise-build-${py}
    /opt/python/cp37-cp37m/bin/cmake /opt/matslise/ \
        -DCMAKE_BUILD_TYPE=Release \
        -DPYTHON_EXECUTABLE=/opt/python/${py}/bin/python

    /opt/python/cp37-cp37m/bin/cmake --build . --target build_wheel --config Release -- -j 6

    auditwheel repair matslise/dist/*${py}*.whl -w /opt/matslise/wheelhouse
done