///
/// @file  Setting.cc
/// @brief Hammer setting class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/Tools/Setting.hh"
#include "Hammer/Tools/SettingVisitors.hh"
#include "Hammer/Tools/HammerSerial.hh"
#include "Hammer/Tools/Logging.hh"
#include "Hammer/Math/Utils.hh"

using namespace std;

namespace Hammer {

    bool Setting::_encodeUseDefault = false;

    bool SettingChecker::operator()(boost::blank /*unused*/, boost::blank /*unused*/) const {
        return true;
    }

    bool SettingChecker::operator()(bool a, bool b) const {
        return a == b;
    }

    bool SettingChecker::operator()(int a, int b) const {
        return a == b;
    }

    bool SettingChecker::operator()(double a, double b) const {
        return isZero(a - b);
    }

    bool SettingChecker::operator()(complex<double> a, complex<double> b) const {
        return isZero(a - b);
    }

    bool SettingChecker::operator()(const string& a, const string& b) const {
        return a == b;
    }

    bool SettingChecker::operator()(const vector<double>& a, const vector<double>& b) const {
        if (a.size() != b.size()) {
            return false;
        }
        return equal(a.begin(), a.end(), b.begin(), [](double n1, double n2) -> bool { return isZero(n1 - n2); });
    }

    bool SettingChecker::operator()(const vector<string>& a, const vector<string>& b) const {
        return a == b;
    }

    bool SettingChecker::operator()(const Setting::MatrixType& a, const Setting::MatrixType& b) const {
        return a == b;
    }

    SettingWriter::SettingWriter(flatbuffers::FlatBufferBuilder* builder) : _builder{builder} {
    }

    Setting::WrittenSettingType SettingWriter::operator()(boost::blank /*unused*/) const {
        auto resc = Serial::FBSBool{false};
        auto out = _builder->CreateStruct(resc);
        return make_pair(out.Union(), Serial::FBSettingTypes::FBSBool);
    }

    Setting::WrittenSettingType SettingWriter::operator()(bool a) const {
        auto resc = Serial::FBSBool{a};
        auto out = _builder->CreateStruct(resc);
        return make_pair(out.Union(), Serial::FBSettingTypes::FBSBool);
    }

    Setting::WrittenSettingType SettingWriter::operator()(int a) const {
        auto resc = Serial::FBSInt{a};
        auto out = _builder->CreateStruct(resc);
        return make_pair(out.Union(), Serial::FBSettingTypes::FBSInt);
    }

    Setting::WrittenSettingType SettingWriter::operator()(double a) const {
        auto resc = Serial::FBSDouble{a};
        auto out = _builder->CreateStruct(resc);
        return make_pair(out.Union(), Serial::FBSettingTypes::FBSDouble);
    }

    Setting::WrittenSettingType SettingWriter::operator()(std::complex<double> a) const {
        auto resc = Serial::FBComplex{a.real(), a.imag()};
        auto out = _builder->CreateStruct(resc);
        return make_pair(out.Union(), Serial::FBSettingTypes::FBComplex);
    }

    Setting::WrittenSettingType SettingWriter::operator()(const std::string& a) const {
        auto serialval = _builder->CreateString(a);
        Serial::FBSStringBuilder serialBuilder{*_builder};
        serialBuilder.add_value(serialval);
        auto out = serialBuilder.Finish();
        return make_pair(out.Union(), Serial::FBSettingTypes::FBSString);
    }

    Setting::WrittenSettingType SettingWriter::operator()(const std::vector<double>& a) const {
        auto serialvals = _builder->CreateVector(a);
        Serial::FBVecDoubleBuilder serialBuilder{*_builder};
        serialBuilder.add_value(serialvals);
        auto out = serialBuilder.Finish();
        return make_pair(out.Union(), Serial::FBSettingTypes::FBVecDouble);
    }

    Setting::WrittenSettingType SettingWriter::operator()(const std::vector<std::string>& a) const {
        auto serialvals = _builder->CreateVectorOfStrings(a);
        Serial::FBVecStringBuilder serialBuilder{*_builder};
        serialBuilder.add_value(serialvals);
        auto out = serialBuilder.Finish();
        return make_pair(out.Union(), Serial::FBSettingTypes::FBVecString);
    }

    Setting::WrittenSettingType SettingWriter::operator()(const Setting::MatrixType& a) const {
        vector<flatbuffers::Offset<Serial::FBVecDouble>> vecs;
        for (const auto& elem : a) {
            auto serialvals = _builder->CreateVector(elem);
            Serial::FBVecDoubleBuilder serialBuilder{*_builder};
            serialBuilder.add_value(serialvals);
            vecs.push_back(serialBuilder.Finish());
        }
        auto serialvecs = _builder->CreateVector(vecs);
        Serial::FBSVecVecDoubleBuilder serialBuilder{*_builder};
        serialBuilder.add_value(serialvecs);
        auto out = serialBuilder.Finish();
        return make_pair(out.Union(), Serial::FBSettingTypes::FBSVecVecDouble);
    }

    YAML::Node SettingEncoder::operator()(boost::blank /*unused*/) const {
        YAML::Node node;
        return node;
    }

    YAML::Node SettingEncoder::operator()(bool a) const {
        YAML::Node node;
        node = a;
        return node;
    }

    YAML::Node SettingEncoder::operator()(int a) const {
        YAML::Node node;
        node = a;
        return node;
    }

    YAML::Node SettingEncoder::operator()(double a) const {
        YAML::Node node;
        node = a;
        return node;
    }

    YAML::Node SettingEncoder::operator()(std::complex<double> a) const {
        YAML::Node node;
        node.push_back(a.real());
        node.push_back(a.imag());
        node.SetStyle(YAML::EmitterStyle::Flow);
        return node;
    }

    YAML::Node SettingEncoder::operator()(const std::string& a) const {
        YAML::Node node;
        node = a;
        return node;
    }

    YAML::Node SettingEncoder::operator()(const std::vector<double>& a) const {
        YAML::Node node;
        for (auto elem : a) {
            node.push_back(elem);
        }
        node.SetStyle(YAML::EmitterStyle::Flow);
        return node;
    }

    YAML::Node SettingEncoder::operator()(const std::vector<std::string>& a) const {
        YAML::Node node;
        for (const auto& elem : a) {
            node.push_back(elem);
        }
        node.SetStyle(YAML::EmitterStyle::Flow);
        return node;
    }

    YAML::Node SettingEncoder::operator()(const Setting::MatrixType& a) const {
        YAML::Node node;
        for (const auto& elem : a) {
            YAML::Node tmpNode;
            for (auto elem2 : elem) {
                tmpNode.push_back(elem2);
            }
            tmpNode.SetStyle(YAML::EmitterStyle::Flow);
            node.push_back(tmpNode);
        }
        node.SetStyle(YAML::EmitterStyle::Flow);
        return node;
    }

    std::string SettingStringConverter::operator()(boost::blank /*unused*/) const {
        return "BLANK";
    }

    std::string SettingStringConverter::operator()(bool a) const {
        return (a ? "On" : "Off");
    }

    std::string SettingStringConverter::operator()(int a) const {
        return to_string(a);
    }

    std::string SettingStringConverter::operator()(double a) const {
        return to_string(a);
    }

    std::string SettingStringConverter::operator()(std::complex<double> a) const {
        string tmp = "[ ";
        tmp += to_string(a.real());
        tmp += ", ";
        tmp += to_string(a.imag());
        tmp += "]";
        return tmp;
    }

    std::string SettingStringConverter::operator()(const std::string& a) const {
        return a;
    }

    std::string SettingStringConverter::operator()(const std::vector<double>& a) const {
        string tmp = "[ ";
        for (auto elem : a) {
            tmp += to_string(elem) + ", ";
        }
        tmp += "]";
        return tmp;
    }

    std::string SettingStringConverter::operator()(const std::vector<std::string>& a) const {
        string tmp = "[ ";
        for (const auto& elem : a) {
            tmp += elem + ", ";
        }
        tmp += "]";
        return tmp;
    }

    std::string SettingStringConverter::operator()(const Setting::MatrixType& a) const {
        string tmp = "[ ";
        for (const auto& elem : a) {
            tmp += "[ ";
            for (const auto& elem2 : elem) {
                tmp += to_string(elem2) + ", ";
            }
            tmp += "], ";
        }
        tmp += "]";
        return tmp;
    }

    Setting::Setting(const Serial::FBSetting* msgreader) {
        read(msgreader);
    }

    Setting& Setting::operator=(const Setting& other) {
        if (_value.which() != 0 && _value.which() != other._value.which()) {
            throw InitializationError("Setting types do not match!");
        }
        _value = other._value;
        if (_default.which() == 0) {
            _default = other._default;
        }
        return *this;
    }

    void Setting::reset() {
        _value = _default;
    }

    void Setting::setEncodeUseDefault(bool useDefault) {
        _encodeUseDefault = useDefault;
    }

    bool Setting::wasChanged() const {
        SettingChecker checker;
        return !(boost::apply_visitor(checker, _value, _default));
    }

    bool Setting::isSame(const Setting& other) const {
        SettingChecker checker;
        return boost::apply_visitor(checker, _value, other._value);
    }

    void Setting::setDefault() {
        _default = _value;
    }

    void Setting::update(const Setting& other) {
        if (_value.which() != 0 && _value.which() != other._value.which()) {
            // deal with complex<double>
            if (_value.which() == 5 && other._value.which() == 7) {
                const vector<double>* tmp = boost::get<vector<double>>(&other._value);
                if (tmp != nullptr && tmp->size() == 2) {
                    _value = complex<double>((*tmp)[0], (*tmp)[1]);
                    return;
                }
            }
            throw InitializationError("Setting types do not match!");
        }
        bool updateDefault = (_value.which() == 0);
        _value = other._value;
        if (updateDefault) {
            setDefault();
        }
    }

    Setting::WrittenSettingType Setting::write(flatbuffers::FlatBufferBuilder* msgwriter) const {
        SettingWriter writer{msgwriter};
        return boost::apply_visitor(writer, _value);
    }

    void Setting::read(const Serial::FBSetting* msgreader) {
        switch (msgreader->value_type()) {
        case Serial::FBSettingTypes::FBSBool:
            setValue<bool>(msgreader->value_as_FBSBool()->value());
            break;
        case Serial::FBSettingTypes::FBSInt:
            setValue<int>(msgreader->value_as_FBSInt()->value());
            break;
        case Serial::FBSettingTypes::FBSDouble:
            setValue<double>(msgreader->value_as_FBSDouble()->value());
            break;
        case Serial::FBSettingTypes::FBSString: {
            const auto* strField = msgreader->value_as_FBSString()->value();
            if (strField == nullptr) {
                Log::getLog("Hammer.Setting") << Log::ERROR << "Missing value in read: call Peter, Paul and Mary!\n";
                return;
            }
            setValue(string(strField->c_str()));
            break;
        }
        case Serial::FBSettingTypes::FBComplex: {
            const auto* serialcompl = msgreader->value_as_FBComplex();
            complex<double> val{serialcompl->re(), serialcompl->im()};
            setValue(val);
            break;
        }
        case Serial::FBSettingTypes::FBVecDouble: {
            const auto* serialvecd = msgreader->value_as_FBVecDouble()->value();
            if (serialvecd == nullptr) {
                Log::getLog("Hammer.Setting") << Log::ERROR << "Missing value in read: call Peter, Paul and Mary!\n";
                return;
            }
            vector<double> vecd;
            vecd.reserve(serialvecd->size());
            for (unsigned int j = 0; j < serialvecd->size(); ++j) {
                vecd.push_back(serialvecd->Get(j));
            }
            setValue(vecd);
            break;
        }
        case Serial::FBSettingTypes::FBVecString: {
            const auto* serialvecs = msgreader->value_as_FBVecString()->value();
            if (serialvecs == nullptr) {
                Log::getLog("Hammer.Setting") << Log::ERROR << "Missing value in read: call Peter, Paul and Mary!\n";
                return;
            }
            vector<string> vecs;
            vecs.reserve(serialvecs->size());
            for (unsigned int j = 0; j < serialvecs->size(); ++j) {
                const auto* sstr = serialvecs->Get(j);
                if (sstr == nullptr) {
                    continue;
                }
                vecs.emplace_back(sstr->c_str());
            }
            setValue(vecs);
            break;
        }
        case Serial::FBSettingTypes::FBSVecVecDouble: {
            const auto* serialvecvecd = msgreader->value_as_FBSVecVecDouble()->value();
            if (serialvecvecd == nullptr) {
                Log::getLog("Hammer.Setting") << Log::ERROR << "Missing value in read: call Peter, Paul and Mary!\n";
                return;
            }
            vector<vector<double>> vecvecd;
            vecvecd.reserve(serialvecvecd->size());
            for (unsigned int i = 0; i < serialvecvecd->size(); ++i) {
                const auto* innerElem = serialvecvecd->Get(i);
                if (innerElem == nullptr || innerElem->value() == nullptr) {
                    continue;
                }
                const auto* serialvecd = innerElem->value();
                vector<double> vecd;
                vecd.reserve(serialvecd->size());
                for (unsigned int j = 0; j < serialvecd->size(); ++j) {
                    vecd.push_back(serialvecd->Get(j));
                }
                vecvecd.push_back(vecd);
            }
            setValue(vecvecd);
            break;
        }
        case Serial::FBSettingTypes::NONE:
            break;
        }
        if (_default.empty()) {
            _default = _value;
        }
    }

    string Setting::toString(const bool useDefault) const {
        SettingStringConverter strConverter;
        return boost::apply_visitor(strConverter, (useDefault ? _default : _value));
    }


    YAML::Emitter& operator<<(YAML::Emitter& out, const Setting& s) {
        out << YAML::convert<Setting>::encode(s);
        return out;
    }

} // namespace Hammer

namespace YAML {

    Node convert<::Hammer::Setting>::encode(const ::Hammer::Setting& value) {
        ::Hammer::SettingEncoder encoder;
        if (::Hammer::Setting::_encodeUseDefault) {
            return boost::apply_visitor(encoder, value._default);
        }
        return boost::apply_visitor(encoder, value._value);
    }

    bool convert<::Hammer::Setting>::decode(const Node& node, ::Hammer::Setting& value) {
        if (node.IsScalar()) {
            try {
                bool tmp = node.as<bool>();
                value._value = tmp;
            } catch (YAML::Exception&) {
                try {
                    int tmp = node.as<int>();
                    value._value = tmp;
                } catch (YAML::Exception&) {
                    try {
                        auto tmp = node.as<double>();
                        value._value = tmp;
                    } catch (YAML::Exception&) {
                        auto tmp = node.as<string>();
                        value._value = tmp;
                    }
                }
            }
            return true;
        }
        if (node.IsSequence()) {
            if (node.size() > 0 && node.begin()->Type() == YAML::NodeType::Sequence &&
                node.begin()->begin()->Type() == YAML::NodeType::Scalar) {
                ::Hammer::Setting::MatrixType tmp;
                try {
                    for (YAML::const_iterator itvv = node.begin(); itvv != node.end(); ++itvv) {
                        auto tmpvec = itvv->as<vector<double>>();
                        tmp.push_back(tmpvec);
                    }
                    value._value = tmp;
                } catch (YAML::Exception&) {
                    return false;
                }
                return true;
            }
            if (node.size() > 0 && node.begin()->Type() != YAML::NodeType::Scalar) {
                return false;
            }
            try {
                auto tmp = node.as<vector<double>>();
                if (tmp.size() == 2 && value.getDefault<complex<double>>() != nullptr) {
                    value._value = complex<double>{tmp[0], tmp[1]};
                } else {
                    value._value = tmp;
                }
            } catch (YAML::Exception&) {
                try {
                    auto tmp = node.as<vector<string>>();
                    value._value = tmp;
                } catch (YAML::Exception&) {
                    return false;
                }
            }
            return true;
        }
        return false;
    }

} // namespace YAML
