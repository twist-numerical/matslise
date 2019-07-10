for f in $(pyenv versions --bare); do
	python=$(pyenv prefix $f)/bin/python;
	$python -m pip install wheel;
	mkdir "cmake-build-mac-$f";
	cd "cmake-build-mac-$f";
	cmake -DLONG_DOUBLE=OFF -DPYTHON_EXECUTABLE=$python \
		-DCMAKE_BUILD_TYPE=Release ..
	cmake --build . --target build_wheel --config Release -- -j 4
	cd ..
	done;
