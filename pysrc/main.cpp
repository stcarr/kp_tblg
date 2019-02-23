#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include "../src/kp_tblg_construct.h"

namespace py = pybind11;

PYBIND11_MODULE(tblg_kpy, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------
        .. currentmodule:: tblg_kpy
        .. autosummary::
           :toctree: _generate
           add
    )pbdoc";

    py::class_<Kp_tblg_construct>(m, "Kp_tblg_construct")
        .def(py::init())
        .def("setTwist", &Kp_tblg_construct::setTwist)
        .def("loadFiles", &Kp_tblg_construct::loadFiles)
        .def("setInterFac", &Kp_tblg_construct::setInterFac)
        .def("setStrainFac", &Kp_tblg_construct::setStrainFac) 
        .def("interpKP", &Kp_tblg_construct::interpKP)
	.def("prepare", &Kp_tblg_construct::prepare)
        .def("layer1Ham", &Kp_tblg_construct::layer1Ham)
        .def("layer2Ham", &Kp_tblg_construct::layer2Ham)
        .def("grapheneIntralayerTerm", &Kp_tblg_construct::grapheneIntralayerTerm)
        .def("getReciprocal", &Kp_tblg_construct::getReciprocal)
        .def("crossProd", &Kp_tblg_construct::crossProd)
        .def("getH", &Kp_tblg_construct::getH)
        .def("getSize", &Kp_tblg_construct::getSize);

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
