#include <pybind11/pybind11.h>
#include "example1.hpp"
#include "example2.hpp"

PYBIND11_MODULE(_ppc, m) {
    m.doc() = "This is _ppc's docstring.";
    m.def("add", &add, "Add up two numbers.");
    m.def("sub", &sub, "Find difference of two numbers.");
}

int main(){
    int c, d;
    
    c = add(1, 2);
    d = sub(1, 2);
}