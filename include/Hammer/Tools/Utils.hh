///
/// @file  Utils.hh
/// @brief Hammer utility functions
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_TOOLS_UTILS
#define HAMMER_TOOLS_UTILS

#include <exception>
#include <map>
#include <unordered_map>
#include <string>
#include <algorithm>
#include <iterator>
#include <vector>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <optional>

#include <boost/functional/hash.hpp>

#include "Hammer/Config/HammerConfig.hh"

#define UNUSED(x) ((void) (x))

namespace Hammer {

    using MaybeBool = std::optional<bool>;

    /// @brief
    /// @return
    inline std::string version() {
        return HAMMER_NAME " ver. " HAMMER_VERSION;
    }

    template <typename KeyType, typename ValueType>
    ValueType getOrDefault(const std::map<KeyType, ValueType>& data, const KeyType& key, const ValueType& fallback) {
        auto it = data.find(key);
        return (it == data.end()) ? fallback : it->second;
    }

    template <typename KeyType, typename ValueType>
    ValueType getOrDefault(const std::unordered_map<KeyType, ValueType>& data, const KeyType& key,
                           const ValueType& fallback) {
        auto it = data.find(key);
        return (it == data.end()) ? fallback : it->second;
    }

    template <typename KeyType, typename ValueType>
    auto const& getOrThrow(const std::map<KeyType, ValueType>& data, const KeyType& key, const std::exception& error) {
        auto it = data.find(key);
        if (it == data.end()) {
            throw error;
        }
        return it->second;
    }

    template <typename KeyType, typename ValueType>
    auto& getOrThrow(std::map<KeyType, ValueType>& data, const KeyType& key, const std::exception& error) {
        auto it = data.find(key);
        if (it == data.end()) {
            throw error;
        }
        return it->second;
    }

    template <typename KeyType, typename ValueType>
    auto const& getOrThrow(const std::unordered_map<KeyType, ValueType>& data, const KeyType& key,
                           const std::exception& error) {
        auto it = data.find(key);
        if (it == data.end()) {
            throw error;
        }
        return it->second;
    }

    template <typename KeyType, typename ValueType>
    auto& getOrThrow(std::unordered_map<KeyType, ValueType>& data, const KeyType& key, const std::exception& error) {
        auto it = data.find(key);
        if (it == data.end()) {
            throw error;
        }
        return it->second;
    }

    template <typename ContainerType, typename KeyType>
    bool checkExistence(ContainerType& container, KeyType const& key) {
        return (container.find(key) != container.end());
    }

    template <typename ContainerType, typename KeyType, typename... ResidualKeyTypes>
    bool checkExistence(ContainerType& container, KeyType const& key, ResidualKeyTypes const&... rks) {
        auto it = container.find(key);
        if (it == container.end()) {
            return false;
        }
        return checkExistence(it->second, rks...);
    }

    template <typename InputIterator, typename OutputIterator, typename UnaryOperation>
    OutputIterator transform_n(InputIterator _first, size_t _n, OutputIterator _result, UnaryOperation _op) {
        return std::generate_n(_result, _n, [&_first, &_op]() -> decltype(auto) { return _op(*_first++); });
    }

    template <typename T>
    struct reversion_wrapper {
        T& iterable;
    };

    template <typename T>
    auto begin(reversion_wrapper<T> w) {
        return std::rbegin(w.iterable);
    }

    template <typename T>
    auto end(reversion_wrapper<T> w) {
        return std::rend(w.iterable);
    }

    template <typename T>
    reversion_wrapper<T> reverse_range(T&& iterable) {
        return {iterable};
    }

    template <typename T>
    std::ostream& operator<<(std::ostream& out, const std::vector<T>& v) {
        if (!v.empty()) {
            out << '[';
            std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
            out << "\b\b]";
        }
        return out;
    }


    template <typename T>
    inline void combine_hash(std::size_t& /*unused*/, T const& /*unused*/) {
        throw std::bad_typeid{};
    }

    template <>
    inline void combine_hash(std::size_t& seed, size_t const& value) {
        seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }

    template <>
    inline void combine_hash(std::size_t& seed, int const& value) {
        combine_hash(seed, static_cast<size_t const&>(value));
    }

    template <>
    inline void combine_hash(std::size_t& seed, uint32_t const& value) {
        combine_hash(seed, static_cast<size_t const&>(value));
    }

    template <typename K, typename V>
    using UMap = std::unordered_map<K, V, boost::hash<K>>;


    namespace Hash {

        uint32_t murmur3_32(const uint8_t* key, size_t len, uint32_t seed);

    } // namespace Hash


    template <typename T>
    class DisjointSet {

        static_assert(std::is_integral_v<T>, "Integer-like type is required for DisjointSet");

    private:

        std::unordered_map<T, T> _data = {};

    public:

        T find(T x) {
            if (!_data.count(x)) {
                _data[x] = x;
            }
            if (_data[x] != x) {
                _data[x] = find(_data[x]);
            }
            return _data[x];
        }

        void unite(T x, T y) {
            T px = find(x);
            T py = find(y);
            if (px != py) {
                _data[py] = px;
            }
        }

        void setParent(T x, T root) {
            _data[x] = root;
        }

        // T parent(T x) {
        //     return _data[x];
        // }
    };

    struct hash_pair {
        template <class T1, class T2>
        std::size_t operator()(const std::pair<T1, T2>& p) const {
            auto h1 = std::hash<T1>{}(p.first);
            auto h2 = std::hash<T2>{}(p.second);

            std::size_t seed = 0;
            boost::hash_combine(seed, h1);
            boost::hash_combine(seed, h2);
            return seed;
        }
    };


    namespace Spec {
        inline std::string none() { // readability wrapper, cannot be changed from "" to preserve general index labels
            return "";
        }
    } // namespace Spec


    class ParseOutput final {

    public:

        ParseOutput() = default;
        ParseOutput(const ParseOutput&) = delete;
        ParseOutput& operator=(const ParseOutput&) = delete;
        ParseOutput(ParseOutput&&) noexcept = default;
        ParseOutput& operator=(ParseOutput&&) noexcept = default;
        ~ParseOutput() noexcept {
            _oldBuffer = nullptr;
        }


        void divert() {
            _oldBuffer = std::cout.rdbuf();
            std::cout.rdbuf(_buffer.rdbuf());
        }

        const std::stringstream& getStream() const {
            return _buffer;
        }

        void restore() {
            std::cout.rdbuf(_oldBuffer);
        }

    private:

        std::streambuf* _oldBuffer;
        std::stringstream _buffer;
    };

} // namespace Hammer

#endif
