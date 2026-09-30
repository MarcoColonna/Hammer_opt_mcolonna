///
/// @file  HammerYaml.hh
/// @brief Hammer YaML utility functions
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_TOOLS_HAMMERYAML
#define HAMMER_TOOLS_HAMMERYAML

#include <map>
#include <string>
#include <complex>
#include <vector>

#include "yaml-cpp/yaml.h"

namespace Hammer {

    template <typename T>
    void writeEntry(YAML::Emitter& emitter, const std::pair<std::string, const T*>& elem) {
        emitter << YAML::Key << elem.first;
        emitter << YAML::Value << *elem.second;
    }

    template <typename T>
    void writeDict(YAML::Emitter& emitter, const std::map<std::string, const T*>& dict) {
        for (auto& elem : dict) {
            writeEntry<T>(emitter, elem);
        }
    }

    template <typename T>
    void writeDict2(YAML::Emitter& emitter, const std::map<std::string, std::map<std::string, const T*>>& dict2) {
        for (const auto& elem : dict2) {
            emitter << YAML::Key << elem.first;
            emitter << YAML::Value << YAML::BeginMap;
            writeDict<T>(emitter, elem.second);
            emitter << YAML::EndMap;
        }
    }


} // namespace Hammer

namespace YAML {

    template <>
    struct convert<std::complex<double>> {

        static Node encode(const std::complex<double>& value) {
            Node node;
            node.push_back(value.real());
            node.push_back(value.imag());
            node.SetStyle(EmitterStyle::Flow);
            return node;
        }

        static bool decode(const Node& node, std::complex<double>& value) {
            try {
                auto tmp = node.as<std::vector<double>>();
                if (tmp.size() == 2) {
                    value = std::complex<double>{tmp[0], tmp[1]};
                } else {
                    return false;
                }
            } catch (YAML::Exception&) {
                return false;
            }
            return true;
        }
    };

} // namespace YAML

#endif
