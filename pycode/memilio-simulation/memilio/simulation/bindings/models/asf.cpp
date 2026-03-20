/* 
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Kilian Volmer
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#include "pybind_util.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

// Macro to skip main() in jolly_bindings.cpp
#define ASF_BINDINGS_SKIP_MAIN
#include "asf.cpp"

namespace py = pybind11;

PYBIND11_MODULE(_simulation_asf, m)
{
    m.def("simulate", &simulate, "Simulates the asf model", py::arg("epsilon") = 40, py::arg("beta") = 0.0008,
          py::arg("beta2") = 0.0002);

    m.attr("__version__") = "dev";
}
