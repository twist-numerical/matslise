# Set-ExecutionPolicy -ExecutionPolicy ByPass

$pythons = (py --list-paths)

ForEach ($line in $($pythons -split "`r`n"))
{
	if($line -match "\S") {
		$version, $pythonPath = $($Line -split "`t")
		$version = $version -replace " -",""
		$buildDir = "cmake-build-python-"+$version
		echo ("Building for: " + $version)

		iex "$($pythonPath) -m wheel version" | Out-String -OutVariable wheelVersion
		if(-Not $wheelVersion.startsWith("wheel")) {
			echo "installing wheel"
			iex "$($pythonPath) -m pip install wheel"
		}

		md -Force $buildDir
		cd $buildDir
		iex ("cmake .. "+
			"-DLONG_DOUBLE=OFF "+
			"-DPYTHON_EXECUTABLE=$($pythonPath) "+
			"-DCMAKE_BUILD_TYPE=Release ")
		cmake --build . --target buildWheel
		cd ..
	}
}
