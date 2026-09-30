///
/// @file  TestSetting.cc
/// @brief Unit tests for Hammer::Setting
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "gtest/gtest.h"
#include "Hammer/Tools/Setting.hh"
#include "Hammer/Tools/HammerSerial.hh"
#include "Hammer/Exceptions.hh"

#include <complex>
#include <string>
#include <vector>

using namespace std;

namespace Hammer {

    TEST(SettingTest, DefaultConstructedIsBlank) {
        Setting s;
        EXPECT_EQ(s.toString(), "BLANK");
    }

    TEST(SettingTest, SetAndGetBool) {
        Setting s;
        s.setValue<bool>(true);
        auto* p = s.getValue<bool>();
        ASSERT_NE(p, nullptr);
        EXPECT_TRUE(*p);
        EXPECT_EQ(s.toString(), "On");

        s.setValue<bool>(false);
        EXPECT_EQ(s.toString(), "Off");
    }

    TEST(SettingTest, SetAndGetInt) {
        Setting s;
        s.setValue<int>(42);
        auto* p = s.getValue<int>();
        ASSERT_NE(p, nullptr);
        EXPECT_EQ(*p, 42);
        EXPECT_EQ(s.toString(), "42");
    }

    TEST(SettingTest, SetAndGetDouble) {
        Setting s;
        s.setValue<double>(3.14);
        auto* p = s.getValue<double>();
        ASSERT_NE(p, nullptr);
        EXPECT_NEAR(*p, 3.14, 1e-10);
        EXPECT_NE(s.toString(), "");
    }

    TEST(SettingTest, SetAndGetString) {
        Setting s;
        s.setValue<string>(string("hello"));
        auto* p = s.getValue<string>();
        ASSERT_NE(p, nullptr);
        EXPECT_EQ(*p, "hello");
        EXPECT_EQ(s.toString(), "hello");
    }

    TEST(SettingTest, SetAndGetComplex) {
        Setting s;
        complex<double> val{1.0, 2.0};
        s.setValue(val);
        auto* p = s.getValue<complex<double>>();
        ASSERT_NE(p, nullptr);
        EXPECT_NEAR(p->real(), 1.0, 1e-10);
        EXPECT_NEAR(p->imag(), 2.0, 1e-10);
        string str = s.toString();
        EXPECT_NE(str.find('1'), string::npos);
    }

    TEST(SettingTest, SetAndGetVectorDouble) {
        Setting s;
        vector<double> vals{1.0, 2.0, 3.0};
        s.setValue(vals);
        auto* p = s.getValue<vector<double>>();
        ASSERT_NE(p, nullptr);
        EXPECT_EQ(p->size(), 3u);
        EXPECT_NEAR((*p)[0], 1.0, 1e-10);
        EXPECT_NEAR((*p)[2], 3.0, 1e-10);
        EXPECT_NE(s.toString(), "");
    }

    TEST(SettingTest, SetAndGetVectorString) {
        Setting s;
        vector<string> vals{"a", "b", "c"};
        s.setValue(vals);
        auto* p = s.getValue<vector<string>>();
        ASSERT_NE(p, nullptr);
        EXPECT_EQ(p->size(), 3u);
        EXPECT_EQ((*p)[1], "b");
        string str = s.toString();
        EXPECT_NE(str.find('b'), string::npos);
    }

    TEST(SettingTest, SetAndGetMatrix) {
        Setting s;
        vector<vector<double>> mat{{1.0, 2.0}, {3.0, 4.0}};
        s.setValue(mat);
        auto* p = s.getValue<vector<vector<double>>>();
        ASSERT_NE(p, nullptr);
        EXPECT_EQ(p->size(), 2u);
        EXPECT_NEAR((*p)[0][0], 1.0, 1e-10);
        EXPECT_NEAR((*p)[1][1], 4.0, 1e-10);
        EXPECT_NE(s.toString(), "");
    }

    TEST(SettingTest, GetWrongTypeReturnsNull) {
        Setting s;
        s.setValue<int>(5);
        EXPECT_EQ(s.getValue<bool>(), nullptr);
        EXPECT_EQ(s.getValue<double>(), nullptr);
        EXPECT_EQ(s.getValue<string>(), nullptr);
    }

    TEST(SettingTest, SetWrongTypeThrows) {
        Setting s;
        s.setValue<int>(5);
        EXPECT_THROW(s.setValue<bool>(true), InitializationError);
        EXPECT_THROW(s.setValue<double>(1.0), InitializationError);
    }

    TEST(SettingTest, DefaultValueAndReset) {
        Setting s;
        s.setValue<int>(10);
        s.setDefault();
        s.setValue<int>(99);
        EXPECT_EQ(*s.getValue<int>(), 99);
        s.reset();
        EXPECT_EQ(*s.getValue<int>(), 10);
    }

    TEST(SettingTest, SetDefaultTemplate) {
        Setting s;
        s.setValue<int>(5);
        s.setDefault<int>(7);
        EXPECT_EQ(*s.getDefault<int>(), 7);
        s.reset();
        EXPECT_EQ(*s.getValue<int>(), 7);
    }

    TEST(SettingTest, WasChanged) {
        Setting s;
        s.setValue<int>(5);
        s.setDefault();
        EXPECT_FALSE(s.wasChanged());
        s.setValue<int>(6);
        EXPECT_TRUE(s.wasChanged());
        s.reset();
        EXPECT_FALSE(s.wasChanged());
    }

    TEST(SettingTest, IsSame) {
        Setting a;
        Setting b;
        a.setValue<int>(5);
        b.setValue<int>(5);
        EXPECT_TRUE(a.isSame(b));
        b.setValue<int>(6);
        EXPECT_FALSE(a.isSame(b));
    }

    TEST(SettingTest, IsSameDouble) {
        Setting a;
        Setting b;
        a.setValue<double>(1.0);
        b.setValue<double>(1.0);
        EXPECT_TRUE(a.isSame(b));
        b.setValue<double>(2.0);
        EXPECT_FALSE(a.isSame(b));
    }

    TEST(SettingTest, IsSameComplex) {
        Setting a;
        Setting b;
        a.setValue(complex<double>{1.0, 2.0});
        b.setValue(complex<double>{1.0, 2.0});
        EXPECT_TRUE(a.isSame(b));
        b.setValue(complex<double>{1.0, 3.0});
        EXPECT_FALSE(a.isSame(b));
    }

    TEST(SettingTest, IsSameVectorDouble) {
        Setting a;
        Setting b;
        a.setValue(vector<double>{1.0, 2.0});
        b.setValue(vector<double>{1.0, 2.0});
        EXPECT_TRUE(a.isSame(b));

        Setting c;
        c.setValue(vector<double>{1.0, 3.0});
        EXPECT_FALSE(a.isSame(c));

        Setting d;
        d.setValue(vector<double>{1.0});
        EXPECT_FALSE(a.isSame(d));
    }

    TEST(SettingTest, IsSameVectorString) {
        Setting a;
        Setting b;
        a.setValue(vector<string>{"x", "y"});
        b.setValue(vector<string>{"x", "y"});
        EXPECT_TRUE(a.isSame(b));
        b.setValue(vector<string>{"x", "z"});
        EXPECT_FALSE(a.isSame(b));
    }

    TEST(SettingTest, IsSameMatrix) {
        Setting a;
        Setting b;
        a.setValue(vector<vector<double>>{{1.0, 2.0}});
        b.setValue(vector<vector<double>>{{1.0, 2.0}});
        EXPECT_TRUE(a.isSame(b));
        b.setValue(vector<vector<double>>{{1.0, 3.0}});
        EXPECT_FALSE(a.isSame(b));
    }

    TEST(SettingTest, IsSameBool) {
        Setting a;
        Setting b;
        a.setValue<bool>(true);
        b.setValue<bool>(true);
        EXPECT_TRUE(a.isSame(b));
        b.setValue<bool>(false);
        EXPECT_FALSE(a.isSame(b));
    }

    TEST(SettingTest, AssignmentSameType) {
        Setting a;
        Setting b;
        a.setValue<int>(10);
        b.setValue<int>(20);
        a = b;
        EXPECT_EQ(*a.getValue<int>(), 20);
    }

    TEST(SettingTest, AssignmentBlankReceivesAny) {
        Setting a;
        Setting b;
        b.setValue<int>(7);
        a = b;
        EXPECT_EQ(*a.getValue<int>(), 7);
    }

    TEST(SettingTest, AssignmentTypeMismatchThrows) {
        Setting a;
        Setting b;
        a.setValue<int>(1);
        b.setValue<double>(2.0);
        EXPECT_THROW(a = b, InitializationError);
    }

    TEST(SettingTest, UpdateBlankAdoptsValueAndDefault) {
        Setting a;
        Setting b;
        b.setValue<int>(5);
        a.update(b);
        EXPECT_EQ(*a.getValue<int>(), 5);
    }

    TEST(SettingTest, UpdateSameTypeWorks) {
        Setting a;
        Setting b;
        a.setValue<int>(1);
        a.setDefault();
        b.setValue<int>(2);
        a.update(b);
        EXPECT_EQ(*a.getValue<int>(), 2);
    }

    TEST(SettingTest, UpdateTypeMismatchThrows) {
        Setting a;
        Setting b;
        a.setValue<int>(1);
        b.setValue<double>(1.0);
        EXPECT_THROW(a.update(b), InitializationError);
    }

    TEST(SettingTest, UpdateComplexFromVecDouble) {
        Setting a;
        Setting b;
        a.setValue(complex<double>{0.0, 0.0});
        b.setValue(vector<double>{3.0, 4.0});
        a.update(b);
        auto* p = a.getValue<complex<double>>();
        ASSERT_NE(p, nullptr);
        EXPECT_NEAR(p->real(), 3.0, 1e-10);
        EXPECT_NEAR(p->imag(), 4.0, 1e-10);
    }

    TEST(SettingTest, UpdateComplexFromVecDoubleBadSizeThrows) {
        Setting a;
        Setting b;
        a.setValue(complex<double>{0.0, 0.0});
        b.setValue(vector<double>{3.0, 4.0, 5.0});
        EXPECT_THROW(a.update(b), InitializationError);
    }

    TEST(SettingTest, CopyConstructor) {
        Setting a;
        a.setValue<int>(55);
        Setting b{a};
        EXPECT_EQ(*b.getValue<int>(), 55);
    }

    TEST(SettingTest, MoveAssignment) {
        Setting a;
        a.setValue<int>(77);
        Setting b;
        b = std::move(a);
        EXPECT_EQ(*b.getValue<int>(), 77);
    }

    TEST(SettingTest, ToStringDefault) {
        Setting s;
        s.setValue<int>(100);
        s.setDefault();
        s.setValue<int>(200);
        EXPECT_EQ(s.toString(false), "200");
        EXPECT_EQ(s.toString(true), "100");
    }

    TEST(SettingTest, EncodeUseDefault) {
        Setting::setEncodeUseDefault(false);
        Setting s;
        s.setValue<int>(10);
        Setting::setEncodeUseDefault(true);
        Setting::setEncodeUseDefault(false);
    }

    TEST(SettingTest, YamlEncodeDecodeBool) {
        Setting s;
        s.setValue<bool>(true);
        YAML::Node node = YAML::convert<Setting>::encode(s);
        Setting s2;
        EXPECT_TRUE(YAML::convert<Setting>::decode(node, s2));
        auto* p = s2.getValue<bool>();
        ASSERT_NE(p, nullptr);
        EXPECT_TRUE(*p);
    }

    TEST(SettingTest, YamlEncodeDecodeInt) {
        Setting s;
        s.setValue<int>(42);
        YAML::Node node = YAML::convert<Setting>::encode(s);
        Setting s2;
        EXPECT_TRUE(YAML::convert<Setting>::decode(node, s2));
        auto* p = s2.getValue<int>();
        ASSERT_NE(p, nullptr);
        EXPECT_EQ(*p, 42);
    }

    TEST(SettingTest, YamlEncodeDecodeDouble) {
        Setting s;
        s.setValue<double>(2.718);
        YAML::Node node = YAML::convert<Setting>::encode(s);
        Setting s2;
        EXPECT_TRUE(YAML::convert<Setting>::decode(node, s2));
        auto* p = s2.getValue<double>();
        ASSERT_NE(p, nullptr);
        EXPECT_NEAR(*p, 2.718, 1e-6);
    }

    TEST(SettingTest, YamlEncodeDecodeString) {
        Setting s;
        s.setValue<string>(string("world"));
        YAML::Node node = YAML::convert<Setting>::encode(s);
        Setting s2;
        EXPECT_TRUE(YAML::convert<Setting>::decode(node, s2));
        auto* p = s2.getValue<string>();
        ASSERT_NE(p, nullptr);
        EXPECT_EQ(*p, "world");
    }

    TEST(SettingTest, YamlEncodeDecodeVecDouble) {
        Setting s;
        s.setValue(vector<double>{1.0, 2.0, 3.0});
        YAML::Node node = YAML::convert<Setting>::encode(s);
        Setting s2;
        EXPECT_TRUE(YAML::convert<Setting>::decode(node, s2));
        auto* p = s2.getValue<vector<double>>();
        ASSERT_NE(p, nullptr);
        EXPECT_EQ(p->size(), 3u);
    }

    TEST(SettingTest, YamlEncodeDecodeVecString) {
        Setting s;
        s.setValue(vector<string>{"foo", "bar"});
        YAML::Node node = YAML::convert<Setting>::encode(s);
        Setting s2;
        EXPECT_TRUE(YAML::convert<Setting>::decode(node, s2));
        auto* p = s2.getValue<vector<string>>();
        ASSERT_NE(p, nullptr);
        EXPECT_EQ(p->size(), 2u);
        EXPECT_EQ((*p)[0], "foo");
    }

    TEST(SettingTest, YamlEncodeDecodeMatrix) {
        Setting s;
        s.setValue(vector<vector<double>>{{1.0, 2.0}, {3.0, 4.0}});
        YAML::Node node = YAML::convert<Setting>::encode(s);
        Setting s2;
        EXPECT_TRUE(YAML::convert<Setting>::decode(node, s2));
        auto* p = s2.getValue<vector<vector<double>>>();
        ASSERT_NE(p, nullptr);
        EXPECT_EQ(p->size(), 2u);
    }

    TEST(SettingTest, YamlDecodeComplexFromVecDouble) {
        Setting s;
        s.setValue(complex<double>{0.0, 0.0});
        s.setDefault();
        YAML::Node node;
        node.push_back(5.0);
        node.push_back(6.0);
        EXPECT_TRUE(YAML::convert<Setting>::decode(node, s));
        auto* p = s.getValue<complex<double>>();
        ASSERT_NE(p, nullptr);
        EXPECT_NEAR(p->real(), 5.0, 1e-10);
        EXPECT_NEAR(p->imag(), 6.0, 1e-10);
    }

    TEST(SettingTest, YamlDecodeMapNodeReturnsFalse) {
        Setting s;
        YAML::Node mapNode;
        mapNode["key"] = "value";
        EXPECT_FALSE(YAML::convert<Setting>::decode(mapNode, s));
    }

    TEST(SettingTest, YamlEmitterOperator) {
        Setting s;
        s.setValue<int>(99);
        YAML::Emitter emitter;
        emitter << s;
        EXPECT_NE(string(emitter.c_str()).find("99"), string::npos);
    }

    TEST(SettingTest, YamlEncodeDecodeBlank) {
        Setting s;
        YAML::Node node = YAML::convert<Setting>::encode(s);
        EXPECT_TRUE(node.IsNull() || !node.IsDefined() || node.IsScalar());
    }

    TEST(SettingTest, YamlEncodeDecodeComplex) {
        Setting s;
        s.setValue(complex<double>{1.5, -2.5});
        YAML::Node node = YAML::convert<Setting>::encode(s);
        EXPECT_TRUE(node.IsSequence());
        EXPECT_EQ(node.size(), 2u);
    }

    TEST(SettingTest, GetDefaultNullWhenBlank) {
        Setting s;
        EXPECT_EQ(s.getDefault<int>(), nullptr);
        EXPECT_EQ(s.getDefault<bool>(), nullptr);
    }

    TEST(SettingTest, SetDefaultTypeMismatchThrows) {
        Setting s;
        s.setValue<int>(1);
        s.setDefault<int>(1);
        EXPECT_THROW(s.setDefault<double>(2.0), InitializationError);
    }

    TEST(SettingTest, SetValueTypeMismatchThrows) {
        Setting s;
        s.setValue<bool>(true);
        EXPECT_THROW(s.setValue<int>(1), InitializationError);
    }

    TEST(SettingTest, WasChangedBlankDefault) {
        Setting s;
        EXPECT_FALSE(s.wasChanged());
    }

    namespace {
        Setting fbRoundTrip(Setting& orig) {
            flatbuffers::FlatBufferBuilder builder;
            auto sval = orig.write(&builder);
            Serial::FBSettingBuilder sb{builder};
            sb.add_value_type(sval.second);
            sb.add_value(sval.first);
            builder.Finish(sb.Finish());
            const auto* fbset = flatbuffers::GetRoot<Serial::FBSetting>(builder.GetBufferPointer());
            return Setting(fbset);
        }
    } // namespace

    TEST(SettingTest, FlatBuffersWriteReadBlank) {
        Setting orig;
        auto s = fbRoundTrip(orig);
        ASSERT_NE(s.getValue<bool>(), nullptr);
        EXPECT_FALSE(*s.getValue<bool>());
    }

    TEST(SettingTest, FlatBuffersWriteReadBool) {
        Setting orig;
        orig.setValue<bool>(true);
        auto s = fbRoundTrip(orig);
        ASSERT_NE(s.getValue<bool>(), nullptr);
        EXPECT_TRUE(*s.getValue<bool>());
    }

    TEST(SettingTest, FlatBuffersWriteReadInt) {
        Setting orig;
        orig.setValue<int>(42);
        auto s = fbRoundTrip(orig);
        ASSERT_NE(s.getValue<int>(), nullptr);
        EXPECT_EQ(*s.getValue<int>(), 42);
    }

    TEST(SettingTest, FlatBuffersWriteReadDouble) {
        Setting orig;
        orig.setValue<double>(2.718);
        auto s = fbRoundTrip(orig);
        ASSERT_NE(s.getValue<double>(), nullptr);
        EXPECT_NEAR(*s.getValue<double>(), 2.718, 1e-10);
    }

    TEST(SettingTest, FlatBuffersWriteReadString) {
        Setting orig;
        orig.setValue(string("hello"));
        auto s = fbRoundTrip(orig);
        ASSERT_NE(s.getValue<string>(), nullptr);
        EXPECT_EQ(*s.getValue<string>(), "hello");
    }

    TEST(SettingTest, FlatBuffersWriteReadComplex) {
        Setting orig;
        orig.setValue(complex<double>{1.5, -2.5});
        auto s = fbRoundTrip(orig);
        ASSERT_NE(s.getValue<complex<double>>(), nullptr);
        EXPECT_NEAR(s.getValue<complex<double>>()->real(), 1.5, 1e-10);
        EXPECT_NEAR(s.getValue<complex<double>>()->imag(), -2.5, 1e-10);
    }

    TEST(SettingTest, FlatBuffersWriteReadVecDouble) {
        Setting orig;
        orig.setValue(vector<double>{1.0, 2.0, 3.0});
        auto s = fbRoundTrip(orig);
        ASSERT_NE(s.getValue<vector<double>>(), nullptr);
        EXPECT_EQ(s.getValue<vector<double>>()->size(), 3u);
        EXPECT_NEAR((*s.getValue<vector<double>>())[1], 2.0, 1e-10);
    }

    TEST(SettingTest, FlatBuffersWriteReadVecString) {
        Setting orig;
        orig.setValue(vector<string>{"foo", "bar"});
        auto s = fbRoundTrip(orig);
        ASSERT_NE(s.getValue<vector<string>>(), nullptr);
        EXPECT_EQ((*s.getValue<vector<string>>())[0], "foo");
    }

    TEST(SettingTest, FlatBuffersWriteReadMatrix) {
        Setting orig;
        orig.setValue(vector<vector<double>>{{1.0, 2.0}, {3.0, 4.0}});
        auto s = fbRoundTrip(orig);
        auto* p = s.getValue<vector<vector<double>>>();
        ASSERT_NE(p, nullptr);
        EXPECT_EQ(p->size(), 2u);
        EXPECT_NEAR((*p)[1][1], 4.0, 1e-10);
    }

    TEST(SettingTest, YamlDecodeMatrixNonNumericInner) {
        YAML::Node node;
        YAML::Node inner;
        inner.push_back("not_a_number");
        inner.push_back("also_not");
        node.push_back(inner);
        Setting s;
        EXPECT_FALSE(YAML::convert<Setting>::decode(node, s));
    }

    TEST(SettingTest, YamlDecodeSequenceWithMapElement) {
        YAML::Node node;
        YAML::Node mapElem;
        mapElem["key"] = "val";
        node.push_back(mapElem);
        Setting s;
        EXPECT_FALSE(YAML::convert<Setting>::decode(node, s));
    }

    TEST(SettingTest, YamlDecodeSequenceVecStringFails) {
        YAML::Node node;
        node.push_back("hello");
        YAML::Node mapElem;
        mapElem["key"] = "val";
        node.push_back(mapElem);
        Setting s;
        EXPECT_FALSE(YAML::convert<Setting>::decode(node, s));
    }

    TEST(SettingTest, FlatBuffersReadNone) {
        flatbuffers::FlatBufferBuilder builder;
        Serial::FBSettingBuilder sb{builder};
        builder.Finish(sb.Finish());
        const auto* fbset = flatbuffers::GetRoot<Serial::FBSetting>(builder.GetBufferPointer());
        Setting s(fbset);
        EXPECT_EQ(s.getValue<bool>(), nullptr);
        EXPECT_EQ(s.getValue<int>(), nullptr);
    }

} // namespace Hammer
