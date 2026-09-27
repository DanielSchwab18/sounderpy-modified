// cppimport


#include <cstring>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "include/Sounding.hpp"

namespace py = pybind11;

PYBIND11_MODULE(bindings, m) {
    py::class_<Sounding>(m, "Sounding")
        .def(py::init<std::string, bool, bool, int>(),
             py::arg("sounding_file_name"),
             py::arg("needs_interpolation"),
             py::arg("interpolate_cape"),
             py::arg("resolution"))


        //parcel buoyancy traces
        .def("get_mu_parcel_buoyancy", &Sounding::get_mu_parcel_buoyancy)
        .def("get_ml_parcel_buoyancy", &Sounding::get_ml_parcel_buoyancy)
        .def("get_sb_parcel_buoyancy", &Sounding::get_sb_parcel_buoyancy)
        .def("get_mu_entrainment_parcel_buoyancy", &Sounding::get_mu_entrainment_parcel_buoyancy)
        .def("get_ml_entrainment_parcel_buoyancy", &Sounding::get_ml_entrainment_parcel_buoyancy)
        .def("get_sb_entrainment_parcel_buoyancy", &Sounding::get_sb_entrainment_parcel_buoyancy)

        //parcel temp traces
        .def("get_mu_parcel_temp", &Sounding::get_mu_parcel_temp,
            py::arg("needs_celsius") = false)
        .def("get_ml_parcel_temp", &Sounding::get_ml_parcel_temp,
            py::arg("needs_celsius") = false)
        .def("get_sb_parcel_temp", &Sounding::get_sb_parcel_temp,
            py::arg("needs_celsius") = false)
        .def("get_mu_entrainment_parcel_temp", &Sounding::get_mu_entrainment_parcel_temp,
            py::arg("needs_celsius") = false)
        .def("get_ml_entrainment_parcel_temp", &Sounding::get_ml_entrainment_parcel_temp,
            py::arg("needs_celsius") = false)
        .def("get_sb_entrainment_parcel_temp", &Sounding::get_sb_entrainment_parcel_temp,
            py::arg("needs_celsius") = false)

        //stp
        .def("get_stp", &Sounding::compute_stp)

        //height-variant variables
        .def("get_rh", &Sounding::get_rh)
        .def("get_w", &Sounding::get_w)
        .def("get_srw", &Sounding::get_srw)
        .def("get_srh", &Sounding::get_srh)
        .def("get_shear", &Sounding::get_shear)
        .def("get_streamwise_vorticity", &Sounding::get_streamwise_vorticity)

        //capes
        .def("get_sbcape", &Sounding::get_sbcape)
        .def("get_mlcape", &Sounding::get_mlcape)
        .def("get_mucape", &Sounding::get_mucape)
        .def("get_sbecape", &Sounding::get_sbecape)
        .def("get_mlecape", &Sounding::get_mlecape)
        .def("get_muecape", &Sounding::get_muecape)

        //3capes
        .def("get_sb3cape", &Sounding::get_sb3cape)
        .def("get_ml3cape", &Sounding::get_ml3cape)
        .def("get_mu3cape", &Sounding::get_mu3cape)
        .def("get_sb3ecape", &Sounding::get_sb3ecape)
        .def("get_ml3ecape", &Sounding::get_ml3ecape)
        .def("get_mu3ecape", &Sounding::get_mu3ecape)

        //6capes
        .def("get_sb6cape", &Sounding::get_sb6cape)
        .def("get_ml6cape", &Sounding::get_ml6cape)
        .def("get_mu6cape", &Sounding::get_mu6cape)
        .def("get_sb6ecape", &Sounding::get_sb6ecape)
        .def("get_ml6ecape", &Sounding::get_ml6ecape)
        .def("get_mu6ecape", &Sounding::get_mu6ecape)

        //cins
        .def("get_sbcin", &Sounding::get_sbcin)
        .def("get_mlcin", &Sounding::get_mlcin)
        .def("get_mucin", &Sounding::get_mucin)
        .def("get_sbecin", &Sounding::get_sbecin)
        .def("get_mlecin", &Sounding::get_mlecin)
        .def("get_muecin", &Sounding::get_muecin)

        //buoyancy
        .def("get_sb_B_to_height", &Sounding::get_sb_B_to_height,
            py::arg("height_limit"))
        .def("get_ml_B_to_height", &Sounding::get_ml_B_to_height,
            py::arg("height_limit"))
        .def("get_mu_B_to_height", &Sounding::get_mu_B_to_height,
            py::arg("height_limit"))

        //lfc heights and pressures
        .def("get_sb_lfc", &Sounding::get_sb_lfc)
        .def("get_sb_lfc_p", &Sounding::get_sb_lfc_p)
        .def("get_sb_e_lfc", &Sounding::get_sb_e_lfc)
        .def("get_sb_e_lfc_p", &Sounding::get_sb_e_lfc_p)
        .def("get_ml_lfc", &Sounding::get_ml_lfc)
        .def("get_ml_lfc_p", &Sounding::get_ml_lfc_p)
        .def("get_ml_e_lfc", &Sounding::get_ml_e_lfc)
        .def("get_ml_e_lfc_p", &Sounding::get_ml_e_lfc_p)
        .def("get_mu_lfc", &Sounding::get_mu_lfc)
        .def("get_mu_lfc_p", &Sounding::get_mu_lfc_p)
        .def("get_mu_e_lfc", &Sounding::get_mu_e_lfc)
        .def("get_mu_e_lfc_p", &Sounding::get_mu_e_lfc_p)

        //lcl heights and pressures
        .def("get_sb_lcl", &Sounding::get_sb_lcl)
        .def("get_sb_lcl_p", &Sounding::get_sb_lcl_p)
        .def("get_sb_e_lcl", &Sounding::get_sb_e_lcl)
        .def("get_sb_e_lcl_p", &Sounding::get_sb_e_lcl_p)
        .def("get_ml_lcl", &Sounding::get_ml_lcl)
        .def("get_ml_lcl_p", &Sounding::get_ml_lcl_p)
        .def("get_ml_e_lcl", &Sounding::get_ml_e_lcl)
        .def("get_ml_e_lcl_p", &Sounding::get_ml_e_lcl_p)
        .def("get_mu_lcl", &Sounding::get_mu_lcl)
        .def("get_mu_lcl_p", &Sounding::get_mu_lcl_p)
        .def("get_mu_e_lcl", &Sounding::get_mu_e_lcl)
        .def("get_mu_e_lcl_p", &Sounding::get_mu_e_lcl_p)

        //el heights, pressures, and temps
        .def("get_sb_el", &Sounding::get_sb_el)
        .def("get_sb_el_p", &Sounding::get_sb_el_p)
        .def("get_sb_el_t", &Sounding::get_sb_el_t)
        .def("get_sb_e_el", &Sounding::get_sb_e_el)
        .def("get_sb_e_el_p", &Sounding::get_sb_e_el_p)
        .def("get_sb_e_el_t", &Sounding::get_sb_e_el_t)
        .def("get_ml_el", &Sounding::get_ml_el)
        .def("get_ml_el_p", &Sounding::get_ml_el_p)
        .def("get_ml_el_t", &Sounding::get_ml_el_t)
        .def("get_ml_e_el", &Sounding::get_ml_e_el)
        .def("get_ml_e_el_p", &Sounding::get_ml_e_el_p)
        .def("get_ml_e_el_t", &Sounding::get_ml_e_el_t)
        .def("get_mu_el", &Sounding::get_mu_el)
        .def("get_mu_el_p", &Sounding::get_mu_el_p)
        .def("get_mu_el_t", &Sounding::get_mu_el_t)
        .def("get_mu_e_el", &Sounding::get_mu_e_el)
        .def("get_mu_e_el_p", &Sounding::get_mu_e_el_p)
        .def("get_mu_e_el_t", &Sounding::get_mu_e_el_t)
        

        //Temps and stuff
        .def("get_temperature", &Sounding::get_temperature)
        .def("get_mixing_ratio", &Sounding::get_mixing_ratio,
            py::arg("grams_per_kilogram") = false)
        .def("get_specific_humidity", &Sounding::get_specific_humidity,
            py::arg("grams_per_kilogram") = false)
        .def("get_relative_humidity", &Sounding::get_relative_humidity)
        .def("get_virtual_temperature", &Sounding::get_virtual_temperature,
            py::arg("needs_celsius") = false)
        .def("get_potential_temperature", &Sounding::get_potential_temperature)
        .def("get_sfc_pressure", &Sounding::get_sfc_pressure)


        .def("get_height", &Sounding::get_height)
        .def("get_pressure", &Sounding::get_pressure);
}