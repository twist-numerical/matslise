%module pyslise_swig

%include "std_vector.i"

%{
#include "matslise/pyslise.h"
#include "swig/python_wrapper/function.h"
%}

%feature("director") std::function<double, double>;

%template(DoubleVector) std::vector<double>;
%template(YVector) std::vector<matslise::Y>;


%include "matslise/matslise.h"
%include "swig/python_wrapper/function.h"