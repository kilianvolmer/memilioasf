/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Julia Bicker, René Schmieding, Kilian Volmer
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

#include "memilio/config.h"
#include "memilio/epidemiology/age_group.h"
#include "smm/simulation.h"
#include "smm/model.h"
#include "smm/parameters.h"
#include "memilio/data/analyze_result.h"

enum class InfectionState
{
    S,
    E,
    Iu,
    Id,
    R,
    Cu,
    Cd,
    D,
    B,
    Count
};

mio::TimeSeries<ScalarType> simulate(ScalarType epsilon = 0.2, ScalarType beta = 0.00008, ScalarType beta2 = 0.00002)
{
    using Status = mio::Index<InfectionState>;
    using enum InfectionState;
    using mio::regions::Region;

    /* Example how to run the stochastic metapopulation models with four regions. Within each region we differentiate by
       age groups, species and infection states. The infection states are S, E, C, I, R, D. For the number of age groups
       and species we choose: */
    const size_t num_regions = 1;
    using Model              = mio::smm::Model<ScalarType, InfectionState, Status>;

    /* Define the parameters that are used in the equations */

    // ScalarType beta    = 0.00008;
    // ScalarType beta2   = 0.00002;
    ScalarType lambda = 1. / 4.;
    ScalarType gamma  = 1. / 5.;
    // ScalarType epsilon = 0.2;
    ScalarType sigma = 1. / (38 * 7);
    ScalarType mu    = 1. / (2 * 365);
    ScalarType theta = 0.95;
    ScalarType d     = 1. / 90.;
    ScalarType r     = 1. / 30.;

    Model model(Status{Count}, Region(num_regions));
    for (size_t reg = 0; reg < num_regions; ++reg) {
        model.populations[{Region(reg), S}]  = 10000;
        model.populations[{Region(reg), E}]  = 100;
        model.populations[{Region(reg), Iu}] = 50;
        model.populations[{Region(reg), Id}] = 0;
        model.populations[{Region(reg), R}]  = 0;
        model.populations[{Region(reg), Cu}] = 0;
        model.populations[{Region(reg), Cd}] = 0;
        model.populations[{Region(reg), D}]  = 0;
        model.populations[{Region(reg), B}]  = 10000;
    }

    using AR = mio::smm::AdoptionRates<ScalarType, Status, Region>;

    //Set infection state adoption rates. Adoptions only happen within a region.
    AR::Type adoption_rates;
    for (size_t reg = 0; reg < num_regions; ++reg) {
        adoption_rates.push_back({S, E, Region(reg), 1., {{Id, beta}, {Iu, beta}}});
        adoption_rates.push_back({S, E, Region(reg), 1., {{Cd, beta2}, {Cu, beta2}}});
        adoption_rates.push_back({S, D, Region(reg), mu, {}});
        adoption_rates.push_back({E, Iu, Region(reg), lambda, {}});
        adoption_rates.push_back({Iu, Cu, Region(reg), gamma * theta, {}});
        adoption_rates.push_back({Iu, Id, Region(reg), epsilon, {}});
        adoption_rates.push_back({Cu, D, Region(reg), d, {}});
        adoption_rates.push_back({Id, Cd, Region(reg), gamma * theta, {}});
        adoption_rates.push_back({Cd, D, Region(reg), r, {}});
        adoption_rates.push_back({Iu, R, Region(reg), gamma * (1. - theta), {}});
        adoption_rates.push_back({Id, R, Region(reg), gamma * (1. - theta), {}});
        adoption_rates.push_back({R, S, Region(reg), sigma, {}});
        adoption_rates.push_back({R, D, Region(reg), mu, {}});
        adoption_rates.push_back({B, S, Region(reg), 1, {{{S}, mu}}});
    }

    model.parameters.get<AR>() = adoption_rates;

    ScalarType dt   = 1.;
    ScalarType tmax = 300.0;

    auto sim = mio::smm::Simulation(model, 0.0, dt);
    sim.advance(tmax);

    auto interpolated_results = mio::interpolate_simulation_result(sim.get_result());
    interpolated_results.print_table({"S", "E", "C", "I", "R", "D "});
    return interpolated_results;
}

#ifndef ASF_BINDINGS_SKIP_MAIN
int main()
{
    auto result = simulate();
    mio::log_debug("The result has simulated until {}.", result.get_last_value[0]);
    return 0;
}
#endif