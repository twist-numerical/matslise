import pyslise
import sys
import doctest

failures, test_count = doctest.testmod(pyslise)

assert test_count > 0, "No tests were found"

sys.exit(1 if failures else 0)
