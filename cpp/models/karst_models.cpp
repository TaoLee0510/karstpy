#include "../tmb_ad.hpp"
#include "karst_models.hpp"

#include <array>
#include <map>
#include <string>
#include <utility>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <Eigen/SparseCholesky>
#include <simplicial_inverse_subset.hpp>

namespace py = pybind11;

namespace {
template<class Type> using vector = Eigen::Array<Type, Eigen::Dynamic, 1>;
template<class Type> using matrix = Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>;
using ad = TMBad::ad_aug;

template<class Type>
vector<Type> operator*(const matrix<Type>& mat, const vector<Type>& values) {
    return (mat * values.matrix()).array();
}

template<class Type>
Type dnorm(const Type& x, const Type& mean, const Type& sd, bool) {
    const Type z = (x - mean) / sd;
    return Type(-0.91893853320467274178) - log(sd) - Type(0.5) * z * z;
}

template<class Type>
Type dpois(const Type& count, const Type& mean, bool) {
    const auto y = TMBad::Value(count);
    if (y == 0.0) return -mean;
    return count * log(mean) - mean - Type(std::lgamma(y + 1));
}

template<class Type>
Type logspace_add(const Type& x, const Type& y) {
    using std::exp;
    using std::log1p;
    const Type top = TMBad::CondExpGt(x, y, x, y);
    const Type bottom = TMBad::CondExpGt(x, y, y, x);
    return top + log1p(exp(bottom - top));
}

// Taylor coefficients of log(sinh(sqrt(s) / 2) / (sqrt(s) / 2)) in s.
// At |gap| <= 0.5, truncating at s^8 also preserves the first three derivatives.
template<int order>
double log_sinhc_series(double s) {
    constexpr std::array<double, 9> coefficients{
        0.0, 1.0 / 24, -1.0 / 2880, 5.511463844797177e-6,
        -1.0333994708994705e-7, 2.0876756987868084e-9, -4.403491782239573e-11,
        9.558954664774758e-13, -2.1185501852016107e-14
    };
    double value = 0.0;
    for (int degree = static_cast<int>(coefficients.size()) - 1; degree >= order; --degree) {
        double coefficient = coefficients[degree];
        for (int k = 0; k < order; ++k) coefficient *= degree - k;
        value = value * s + coefficient;
    }
    return value;
}

template<int order>
double smoothed_log_sinhc(double gap) {
    static_assert(order >= 0 && order <= 3);
    constexpr double smoothing = 1e-12;
    if (std::abs(gap) < 0.5) {
        const double gap_squared = gap * gap;
        const double s = gap_squared + smoothing;
        if constexpr (order == 0) return log_sinhc_series<0>(s);
        if constexpr (order == 1) return 2 * gap * log_sinhc_series<1>(s);
        if constexpr (order == 2) return 2 * log_sinhc_series<1>(s) + 4 * gap_squared * log_sinhc_series<2>(s);
        if constexpr (order == 3) return gap * (12 * log_sinhc_series<2>(s) + 8 * gap_squared * log_sinhc_series<3>(s));
    }
    const double x = std::hypot(gap, std::sqrt(smoothing));
    if constexpr (order == 0) return 0.5 * x + std::log(-std::expm1(-x) / x);

    // Differentiate in x, then apply the chain rule for x = sqrt(gap^2 + smoothing).
    const double inverse_x = 1 / x;
    const double dx = gap / x;
    const double tail_exp = std::exp(-x);
    const double tail = tail_exp / (1 - tail_exp);
    const double first = 0.5 + tail - inverse_x;
    if constexpr (order == 1) return first * dx;
    const double d2x = smoothing * inverse_x * inverse_x * inverse_x;
    const double tail_product = tail * (1 + tail);
    const double second = inverse_x * inverse_x - tail_product;
    if constexpr (order == 2) return second * dx * dx + first * d2x;
    const double d3x = -3 * d2x * dx * inverse_x;
    const double third = tail_product * (1 + 2 * tail) - 2 * inverse_x * inverse_x * inverse_x;
    return third * dx * dx * dx + 3 * second * dx * d2x + first * d3x;
}

// Keep each derivative as a unary operator during taping. The numerical
// branch is selected on every sweep. The Hessian's reverse sweep uses order 3.
template<int order>
struct LogSinhcOp : TMBad::global::UnaryOperator {
    static constexpr bool have_eval = true;
    static constexpr int max_fuse_depth = 0;

    double eval(double gap) { return smoothed_log_sinhc<order>(gap); }
    ad eval(const ad& gap) {
        if (gap.constant()) return eval(gap.Value());
        return TMBad::get_glob()->add_to_stack<LogSinhcOp<order>>(TMBad::ad_plain(gap));
    }
    TMBad::Writer eval(const TMBad::Writer&) {
        throw std::logic_error("LogSinhcOp does not support code generation");
    }
    void reverse(TMBad::ReverseArgs<double>& args) {
        if (args.dy(0) != 0) args.dx(0) += args.dy(0) * smoothed_log_sinhc<order + 1>(args.x(0));
    }
    void reverse(TMBad::ReverseArgs<ad>& args) {
        if constexpr (order < 2) {
            args.dx(0) += args.dy(0) * LogSinhcOp<order + 1>().eval(args.x(0));
        } else {
            throw std::logic_error("LogSinhcOp supports up to three numerical derivatives");
        }
    }
    void reverse(TMBad::ReverseArgs<TMBad::Writer>&) {
        throw std::logic_error("LogSinhcOp does not support code generation");
    }
    const char* op_name() { return "LogSinhcOp"; }
};

using Data = std::map<std::string, Eigen::MatrixXd>;
struct Shape {
    Eigen::Index rows;
    Eigen::Index columns;
    size_t offset;
};

template<class Type>
struct objective_function {
    const Data& data;
    const std::map<std::string, Shape>& shapes;
    const std::vector<Type>& parameters;
    py::dict* reports;

    template<class Output>
    Output matrix_data(const std::string& name) const {
        const auto& item = data.at(name);
        return Eigen::Map<const Eigen::MatrixXd>(item.data(), item.rows(), item.cols())
            .template cast<typename Output::Scalar>();
    }

    template<class Output>
    Output vector_data(const std::string& name) const {
        const auto& item = data.at(name);
        return Eigen::Map<const Eigen::VectorXd>(item.data(), item.size())
            .template cast<typename Output::Scalar>().array();
    }

    template<class Output>
    Output scalar_data(const std::string& name) const { return Output(data.at(name)(0)); }

    matrix<Type> parameter_matrix(const std::string& name) const {
        const auto& shape = shapes.at(name);
        return Eigen::Map<const matrix<Type>>(parameters.data() + shape.offset, shape.rows, shape.columns);
    }

    vector<Type> parameter_vector(const std::string& name) const {
        const auto& shape = shapes.at(name);
        return Eigen::Map<const vector<Type>>(parameters.data() + shape.offset, shape.rows * shape.columns);
    }

    Type parameter_scalar(const std::string& name) const { return parameters.at(shapes.at(name).offset); }

    template<class Output>
    void report(const char* name, const Output& value) const {
        if constexpr (std::is_same_v<Type, double>) {
            if (reports) (*reports)[name] = py::cast(value);
        }
    }
};

#define DATA_VECTOR(name) auto name = obj->template vector_data<vector<Type>>(#name)
#define DATA_IVECTOR(name) auto name = obj->template vector_data<vector<int>>(#name)
#define DATA_MATRIX(name) auto name = obj->template matrix_data<matrix<Type>>(#name)
#define DATA_SCALAR(name) auto name = obj->template scalar_data<Type>(#name)
#define DATA_INTEGER(name) auto name = obj->template scalar_data<int>(#name)
#define PARAMETER_VECTOR(name) auto name = obj->parameter_vector(#name)
#define PARAMETER_MATRIX(name) auto name = obj->parameter_matrix(#name)
#define PARAMETER(name) auto name = obj->parameter_scalar(#name)
#define REPORT(name) obj->report(#name, name)

#include "sparse_joint_model.h"
#include "sparse_joint_series_model.h"
#include "sparse_joint_pairwise_model.h"
#include "sparse_joint_state_space_model.h"

class ModelTape {
public:
    ModelTape(const std::string& name, const py::dict& input, const py::dict& parameters)
        : name(name) {
        for (const auto& entry : input)
            data.emplace(entry.first.cast<std::string>(), entry.second.cast<Eigen::MatrixXd>());
        for (const auto& entry : parameters) {
            const auto values = entry.second.cast<Eigen::MatrixXd>();
            shapes.emplace(entry.first.cast<std::string>(), Shape{values.rows(), values.cols(), start.size()});
            start.insert(start.end(), values.data(), values.data() + values.size());
        }
        auto func = [this](const std::vector<ad>& x) { return std::vector<ad>{evaluate(x, nullptr)}; };
        objective = TMBad::ADFun<>(func, start);
        objective.optimize();
        gradient = objective.JacFun();
        gradient.optimize();
        hessian = gradient.SpJacFun();
        hessian.optimize();
    }

    double value(const std::vector<double>& x) { return objective(x).at(0); }
    Eigen::VectorXd grad(const std::vector<double>& x) {
        const auto values = gradient(x);
        return Eigen::Map<const Eigen::VectorXd>(values.data(), values.size());
    }

    Eigen::VectorXd hessian_values(const std::vector<double>& x) {
        const auto values = hessian(x);
        return Eigen::Map<const Eigen::VectorXd>(values.data(), values.size());
    }

    Eigen::VectorXd hessian_reverse(const std::vector<double>& x, const std::vector<double>& weights) {
        hessian(x);
        const auto values = hessian.reverse(weights);
        return Eigen::Map<const Eigen::VectorXd>(values.data(), values.size());
    }

    py::dict report(const std::vector<double>& x) {
        py::dict result;
        evaluate(x, &result);
        return result;
    }

    auto pattern() const { return std::make_pair(hessian.i, hessian.j); }
    std::vector<double> start;

private:
    template<class Type>
    Type evaluate(const std::vector<Type>& x, py::dict* report) const {
        objective_function<Type> obj{data, shapes, x, report};
        if (name == "sparse_joint") return sparse_joint_model(&obj);
        if (name == "sparse_joint_series") return sparse_joint_series_model(&obj);
        if (name == "sparse_joint_pairwise") return sparse_joint_pairwise_model(&obj);
        if (name == "sparse_joint_state_space") return sparse_joint_state_space_model(&obj);
        throw std::invalid_argument("Unknown sparse joint model: " + name);
    }

    std::string name;
    Data data;
    std::map<std::string, Shape> shapes;
    TMBad::ADFun<> objective, gradient;
    TMBad::Sparse<TMBad::ADFun<>> hessian;
};

using SparseMatrix = Eigen::SparseMatrix<double>;
class Cholesky {
public:
    explicit Cholesky(const SparseMatrix& matrix) : inverse(nullptr) {
        inverse.factor = std::make_shared<Eigen::SimplicialLLT<SparseMatrix>>();
        inverse.factor->analyzePattern(matrix);
    }

    bool factorize(const SparseMatrix& matrix, double shift) {
        inverse.factor->setShift(shift);
        inverse.factor->factorize(matrix);
        return inverse.factor->info() == Eigen::Success;
    }

    Eigen::MatrixXd solve(const Eigen::MatrixXd& x) const { return inverse.factor->solve(x); }

    double log_det() const {
        const SparseMatrix lower = inverse.factor->matrixL();
        return 2 * lower.diagonal().array().log().sum();
    }

    Eigen::MatrixXd draw(const Eigen::MatrixXd& normal) const {
        return inverse.factor->permutationPinv() * inverse.factor->matrixU().solve(normal);
    }

    Eigen::VectorXd inverse_entries() {
        const auto selected = inverse.chol2inv();
        const auto perm = inverse.factor->permutationP().indices();
        Eigen::VectorXd diagonal(perm.size());
        for (int i = 0; i < perm.size(); ++i) diagonal[i] = selected.coeff(perm[i], perm[i]);
        return diagonal;
    }

    Eigen::VectorXd inverse_weights(Eigen::VectorXi row, Eigen::VectorXi col) {
        const auto selected = inverse.chol2inv();
        const auto perm = inverse.factor->permutationP().indices();
        Eigen::VectorXd result(row.size());
        for (int i = 0; i < row.size(); ++i) {
            const auto r = perm[row[i]], c = perm[col[i]];
            result[i] = selected.coeff(std::max(r, c), std::min(r, c));
        }
        return result;
    }

private:
    Eigen::SimplicialInverseSubset<double> inverse;
};
}

void bind_models(py::module_& module) {
    py::class_<ModelTape>(module, "ModelTape")
        .def(py::init<const std::string&, const py::dict&, const py::dict&>())
        .def_readonly("start", &ModelTape::start)
        .def_property_readonly("pattern", &ModelTape::pattern)
        .def("value", &ModelTape::value)
        .def("gradient", &ModelTape::grad)
        .def("hessian_values", &ModelTape::hessian_values)
        .def("hessian_reverse", &ModelTape::hessian_reverse)
        .def("report", &ModelTape::report);
    py::class_<Cholesky>(module, "Cholesky")
        .def(py::init<const SparseMatrix&>())
        .def("factorize", &Cholesky::factorize, py::arg("matrix"), py::arg("shift") = 0.0)
        .def("solve", &Cholesky::solve)
        .def("log_det", &Cholesky::log_det)
        .def("draw", &Cholesky::draw)
        .def("inverse_diagonal", &Cholesky::inverse_entries)
        .def("inverse_weights", &Cholesky::inverse_weights);
}
