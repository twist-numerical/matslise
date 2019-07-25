import git
import json
import os
import subprocess
from cpuinfo import get_cpu_info
from datetime import datetime

# pip install gitpython py-cpuinfo


assert 'MATSLISE_TEST_EXECUTABLE' in os.environ, "MATSLISE_TEST_EXECUTABLE environment variable has to be set"
executable = os.environ['MATSLISE_TEST_EXECUTABLE']
print("Using: '%s'" % executable)

test_output = subprocess.run([executable, '-r', 'compact', '-d', 'yes', '--use-colour', 'no'], check=False,
                             stdout=subprocess.PIPE).stdout

tests = []
next_fail = False
for line in test_output.decode('utf-8').split('\n'):
    if ' s: ' in line:
        time, name = line.split(' s: ', 1)
        try:
            time = float(time)
        except:
            next_fail = True
            continue
        tests.append({
            'name': name,
            'time': time,
            'fail': next_fail
        })
        next_fail = False
    else:
        next_fail = True

repo = git.Repo(search_parent_directories=True)
results = json.dumps({
    'timestamp': str(datetime.now()),
    'commit': repo.head.commit.hexsha,
    'changed': len(repo.head.commit.diff(None)) > 0,
    'build_type': '${CMAKE_BUILD_TYPE}',
    'mkl': '${MKL_FOUND}' == 'TRUE',
    'cpu': get_cpu_info(),
    'tests': tests
})

if 'MATSLISE_BENCHMARK_RESULTS' in os.environ:
    result_file = os.environ['MATSLISE_BENCHMARK_RESULTS']
    print("benchmark saved to: '%s'" % result_file)
    with open(result_file, 'a') as f:
        f.write(results)
