///
/// @file  SettingVisitors.hh
/// @brief Various operations on Setting class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_SettingVisitors_HH
#define HAMMER_SettingVisitors_HH

#include "Hammer/Tools/Setting.hh"

namespace Hammer {

    struct SettingChecker : public boost::static_visitor<bool> {
        bool operator()(boost::blank a, boost::blank b) const;
        bool operator()(bool a, bool b) const;
        bool operator()(int a, int b) const;
        bool operator()(double a, double b) const;
        bool operator()(std::complex<double> a, std::complex<double> b) const;
        bool operator()(const std::string& a, const std::string& b) const;
        bool operator()(const std::vector<double>& a, const std::vector<double>& b) const;
        bool operator()(const std::vector<std::string>& a, const std::vector<std::string>& b) const;
        bool operator()(const Setting::MatrixType& a, const Setting::MatrixType& b) const;
        template <typename T, typename S>
        typename std::enable_if_t<!std::is_same_v<T, S>, bool> operator()(const T& /*unused*/,
                                                                          const S& /*unused*/) const {
            return false;
        }
    };

    class SettingWriter : public boost::static_visitor<Setting::WrittenSettingType> {
    public:

        SettingWriter(flatbuffers::FlatBufferBuilder* builder);

        Setting::WrittenSettingType operator()(boost::blank a) const;
        Setting::WrittenSettingType operator()(bool a) const;
        Setting::WrittenSettingType operator()(int a) const;
        Setting::WrittenSettingType operator()(double a) const;
        Setting::WrittenSettingType operator()(std::complex<double> a) const;
        Setting::WrittenSettingType operator()(const std::string& a) const;
        Setting::WrittenSettingType operator()(const std::vector<double>& a) const;
        Setting::WrittenSettingType operator()(const std::vector<std::string>& a) const;
        Setting::WrittenSettingType operator()(const Setting::MatrixType& a) const;

    private:

        flatbuffers::FlatBufferBuilder* _builder;
    };


    struct SettingEncoder : public boost::static_visitor<YAML::Node> {
        YAML::Node operator()(boost::blank a) const;
        YAML::Node operator()(bool a) const;
        YAML::Node operator()(int a) const;
        YAML::Node operator()(double a) const;
        YAML::Node operator()(std::complex<double> a) const;
        YAML::Node operator()(const std::string& a) const;
        YAML::Node operator()(const std::vector<double>& a) const;
        YAML::Node operator()(const std::vector<std::string>& a) const;
        YAML::Node operator()(const Setting::MatrixType& a) const;
    };

    struct SettingStringConverter : public boost::static_visitor<std::string> {
        std::string operator()(boost::blank a) const;
        std::string operator()(bool a) const;
        std::string operator()(int a) const;
        std::string operator()(double a) const;
        std::string operator()(std::complex<double> a) const;
        std::string operator()(const std::string& a) const;
        std::string operator()(const std::vector<double>& a) const;
        std::string operator()(const std::vector<std::string>& a) const;
        std::string operator()(const Setting::MatrixType& a) const;
    };

} // namespace Hammer


#endif
