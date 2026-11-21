// The seed, twist and rejection algorithms below are adapted from R 4.4.2
// src/main/RNG.c, Copyright (C) 1997-2023 The R Core Team, GPL-2.0-or-later.
// Mersenne Twister: Copyright (C) 1997, 1999 Makoto Matsumoto and Takuji Nishimura.
#include "rng.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>
#include <vector>
#include <pybind11/stl.h>
#include <Rmath.h>

namespace {
class RandomStream;
thread_local RandomStream* active_stream = nullptr;

class RandomStream {
public:
    explicit RandomStream(std::uint32_t seed) {
        // R scrambles for 50 iterations and saves one more value in the
        // position slot before filling the 624-word MT19937 state array.
        for (int i = 0; i < 51; ++i) seed = 69069u * seed + 1u;
        for (auto& word : state) {
            seed = 69069u * seed + 1u;
            word = seed;
        }
    }

    double uniform() {
        if (position == state.size()) twist();
        auto word = state[position++];
        word ^= word >> 11;
        word ^= (word << 7) & 0x9d2c5680u;
        word ^= (word << 15) & 0xefc60000u;
        word ^= word >> 18;
        return word == 0 ? 0.5 * 2.328306437080797e-10 : word * 0x1p-32;
    }

    int index(int n) {
        const int bits = static_cast<int>(std::ceil(std::log2(n)));
        const auto mask = (std::uint64_t{1} << bits) - 1u;
        std::uint64_t value;
        do {
            value = 0;
            for (int bit = 0; bit <= bits; bit += 16)
                value = 65536u * value + static_cast<unsigned>(uniform() * 65536);
            value &= mask;
        } while (value >= static_cast<unsigned>(n));
        return static_cast<int>(value);
    }

    std::vector<int> sample(int n, int size, bool replace) {
        if (n < 1 || size < 0 || (!replace && size > n))
            throw std::invalid_argument("Invalid sampling population or size.");
        std::vector<int> output;
        output.reserve(size);
        if (replace) {
            for (int i = 0; i < size; ++i) output.push_back(index(n));
        } else {
            // Match R's swap with the last available item, storing only the
            // positions that change over the course of the draw.
            std::unordered_map<int, int> moved;
            const auto element = [&moved](int i) {
                const auto found = moved.find(i);
                return found == moved.end() ? i : found->second;
            };
            for (int i = 0; i < size; ++i) {
                const auto position = index(n);
                output.push_back(element(position));
                moved[position] = element(--n);
            }
        }
        return output;
    }

    template<class Function>
    auto draw(Function function) {
        struct Scope {
            RandomStream* previous;
            ~Scope() { active_stream = previous; }
        } scope{active_stream};
        active_stream = this;
        return function();
    }

private:
    void twist() {
        for (size_t i = 0; i < state.size(); ++i) {
            const auto word = (state[i] & 0x80000000u) | (state[(i + 1) % 624] & 0x7fffffffu);
            state[i] = state[(i + 397) % 624] ^ (word >> 1) ^ ((word & 1u) * 0x9908b0dfu);
        }
        position = 0;
    }

    std::array<std::uint32_t, 624> state;
    size_t position = state.size();
};
}

extern "C" double unif_rand() {
    return active_stream->uniform();
}

void bind_rng(pybind11::module_& module) {
    namespace py = pybind11;
    py::class_<RandomStream>(module, "RandomStream")
        .def(py::init<std::uint32_t>())
        .def("uniform", &RandomStream::uniform)
        .def("sample", &RandomStream::sample)
        .def("normal", [](RandomStream& stream) { return stream.draw([] { return norm_rand(); }); })
        .def("binomial", [](RandomStream& stream, int n, double p) {
            return stream.draw([&] { return rbinom(n, p); });
        })
        .def("poisson", [](RandomStream& stream, double mean) {
            return stream.draw([&] { return rpois(mean); });
        })
        .def("gamma", [](RandomStream& stream, double shape) {
            return stream.draw([&] { return rgamma(shape, 1.0); });
        })
        .def("multinomial", [](RandomStream& stream, int size, std::vector<double> probs) {
            std::vector<int> result(probs.size());
            stream.draw([&] { rmultinom(size, probs.data(), static_cast<int>(probs.size()), result.data()); });
            return result;
        });
}
