#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <unordered_set>

#include "karst_core.hpp"
#include "rng.hpp"
#include "models/karst_models.hpp"

namespace py = pybind11;

namespace {
py::array_t<int> parse_karyotypes(const std::vector<std::string>& labels, bool unique) {
    if (labels.empty()) throw std::invalid_argument("At least one karyotype label is required.");
    const auto first = karst::parse_label(labels.front());
    const auto columns = static_cast<py::ssize_t>(first.size());
    py::array_t<int> output({static_cast<py::ssize_t>(labels.size()), columns});
    auto values = output.mutable_unchecked<2>();
    std::unordered_set<std::string> seen;
    for (size_t row = 0; row < labels.size(); ++row) {
        const auto result = row == 0 ? first : karst::parse_label(labels[row]);
        if (result.size() != first.size())
            throw std::invalid_argument("All karyotypes must have the same chromosome dimension.");
        if (unique) {
            const auto key = karst::make_key(result);
            if (!seen.insert(key).second) throw std::invalid_argument("Duplicate karyotype label: " + key);
        }
        for (py::ssize_t column = 0; column < columns; ++column)
            values(row, column) = result[column];
    }
    return output;
}

std::vector<std::string> format_karyotypes(const py::array_t<int>& matrix) {
    const auto values = matrix.unchecked<2>();
    std::vector<std::string> labels;
    labels.reserve(values.shape(0));
    std::vector<int> row(values.shape(1));
    for (py::ssize_t i = 0; i < values.shape(0); ++i) {
        for (py::ssize_t j = 0; j < values.shape(1); ++j) row[j] = values(i, j);
        labels.push_back(karst::make_key(row));
    }
    return labels;
}
}

PYBIND11_MODULE(_core, module) {
    module.doc() = "Native extension for karstpy.";
    bind_rng(module);
    bind_models(module);
    module.def("parse_karyotypes", &parse_karyotypes, py::arg("labels"), py::arg("unique") = true);
    module.def("format_karyotypes", &format_karyotypes);
    module.def("pij", [](int parent, int child, double beta) {
        if (!(beta >= 0.0 && beta <= 1.0)) throw std::invalid_argument("`beta` must be in [0, 1].");
        return karst::chromosome_transition_probability(parent, child, beta);
    });
}
