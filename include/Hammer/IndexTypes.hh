///
/// @file  IndexTypes.hh
/// @brief Hammer data types declarations
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_IndexTypes_HH
#define HAMMER_IndexTypes_HH

#include <cstdint>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <boost/algorithm/string.hpp>

#include "Hammer/Exceptions.hh"
#include "Hammer/Tools/Utils.hh"

namespace Hammer {

    enum class WTerm : std::uint8_t { COMMON, NUMERATOR, DENOMINATOR };

    using ParticleIndex = size_t;
    using ParticleIndices = std::vector<ParticleIndex>;
    using UniqueParticleIndices = std::set<ParticleIndex>;

    using HashId = size_t;

    using AmplitudeUID = HashId;
    template <typename T>
    using AmplitudeIdDict = std::map<AmplitudeUID, T>;
    using HadronicUID = HashId;
    template <typename T>
    using HadronicIdDict = std::map<HadronicUID, T>;

    using VertexUID = HashId;
    template <typename T>
    using VertexIdDict = std::map<VertexUID, T>;
    using VertexUIDSet = std::set<VertexUID>;
    template <typename T>
    using VertexDict = std::map<ParticleIndex, T>;
    using VertexName = std::string;


    using ProcessUID = HashId;
    template <typename T>
    using ProcIdDict = std::map<ProcessUID, T>;

    using EventUID = std::set<ProcessUID>;
    template <typename T>
    using EventIdDict = UMap<EventUID, T>;
    using EventUIDGroup = std::set<EventUID>;
    template <typename T>
    using EventIdGroupDict = UMap<EventUIDGroup, T>;


    using SchemeName = std::string;
    using SchemeNameList = std::vector<SchemeName>;
    template <typename T>
    using SchemeDict = std::map<SchemeName, T>;

    using SpecializationName = std::string;
    using SpecializationNameList = std::vector<SpecializationName>;
    template <typename T>
    using SpecializationDict = std::unordered_map<SpecializationName, T>;

    template <typename T>
    using WCPrefixDict = std::map<std::string, T>;

    using FFIndex = size_t;
    template <typename T>
    using FFIndexDict = std::map<FFIndex, T>;

    struct FFPrefixGroup {
        std::string prefix;
        std::string group;
        [[nodiscard]] std::string get() const {
            return prefix + group;
        }
        bool operator<(const FFPrefixGroup& r) const {
            return (prefix < r.prefix) || ((prefix == r.prefix) && (group < r.group));
        }
        bool operator==(const FFPrefixGroup& r) const {
            return (prefix == r.prefix) && (group == r.group);
        }
    };
    template <typename T>
    using FFPrefixGroupDict = std::map<FFPrefixGroup, T>;

    struct SpecPrefixId {

        std::string prefix;
        std::string id;

        SpecPrefixId() = default;

        explicit SpecPrefixId(const std::pair<std::string, std::string>& values)
            : prefix{values.first}, id{values.second} {
        }

        SpecPrefixId(const std::string& otherPrefix, const std::string& otherId) : prefix{otherPrefix}, id{otherId} {
        }

        SpecPrefixId(std::string&& otherPrefix, std::string&& otherId) : prefix{otherPrefix}, id{otherId} {
        }

        static SpecPrefixId toSpecPrefixId(const std::string& fullId) {
            std::vector<std::string> chunks;
            boost::algorithm::split(chunks, fullId, [](char c) { return c == '@'; });
            auto specId = (chunks.size() > 1) ? chunks[1] : Spec::none();
            auto basePrefix = chunks[0];
            return SpecPrefixId{basePrefix, specId};
        }

        SpecPrefixId& operator=(const std::pair<std::string, std::string>& values) {
            prefix = values.first;
            id = values.second;
            return *this;
        }

        [[nodiscard]] bool hasSpecId() const {
            return id != Spec::none();
        }

        [[nodiscard]] std::string get() const {
            if (!id.empty()) {
                return prefix + "@" + id;
            }
            return prefix;
        }

        bool operator<(const SpecPrefixId& r) const {
            return (prefix < r.prefix) || ((prefix == r.prefix) && (id < r.id));
        }

        bool operator==(const SpecPrefixId& r) const {
            return (prefix == r.prefix) && (id == r.id);
        }
    };

    template <typename T>
    using SpecPrefixIdDict = std::map<SpecPrefixId, T>;

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wswitch-default"
    template <typename T>
    struct NumDenPair {

        NumDenPair() = default;

        T numerator;
        T denominator;

        T& get(WTerm what) {
            switch (what) {
            case WTerm::NUMERATOR:
                return numerator;
            case WTerm::DENOMINATOR:
                return denominator;
            case WTerm::COMMON:
                throw Error("Invalid option");
            }
            throw Error("Invalid option");
        }

        [[nodiscard]] const T& get(WTerm what) const {
            switch (what) {
            case WTerm::NUMERATOR:
                return numerator;
            case WTerm::DENOMINATOR:
                return denominator;
            case WTerm::COMMON:
                throw Error("Invalid option");
            }
            throw Error("Invalid option");
        }
    };

#pragma clang diagnostic pop

} // namespace Hammer

#endif
