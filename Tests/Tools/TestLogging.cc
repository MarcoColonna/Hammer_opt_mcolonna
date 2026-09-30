///
/// @file  TestLogging.cc
/// @brief Tests for Logging
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "gtest/gtest.h"
#include "Hammer/Tools/Logging.hh"
#include "Hammer/Exceptions.hh"

using namespace std;

namespace Hammer {

    TEST(LoggingTest, ConstructionDestruction) {
        Log& lg = Log::getLog("TestLogging");
        EXPECT_TRUE(lg.isActive(Log::WARN));
        EXPECT_FALSE(lg.isActive(Log::DEBUG));
    }

    TEST(LoggingTest, Messages) {
        Log& lg = Log::getLog("TestLogging");
        lg.info("Info Message");
        lg.warn("Warn Message");
        lg.error("Error Message");
    }

    TEST(LoggingTest, ConstructorWithLevel) {
        // Exercise the Log(name, level) constructor path via getLog after setLevel
        Log::setLevel("TestLogging.WithLevel", Log::WARNING);
        Log& lg = Log::getLog("TestLogging.WithLevel");
        EXPECT_TRUE(lg.isActive(Log::WARNING));
        EXPECT_FALSE(lg.isActive(Log::DEBUG));
        EXPECT_EQ(lg.getLevel(), Log::WARNING);
    }

    TEST(LoggingTest, SetLevelsMap) {
        Log& lgA = Log::getLog("TestLogging.MapA");
        Log& lgB = Log::getLog("TestLogging.MapB");
        Log::LevelMap lm;
        lm["TestLogging.MapA"] = Log::ERROR;
        lm["TestLogging.MapB"] = Log::DEBUG;
        Log::setLevels(lm);
        EXPECT_EQ(lgA.getLevel(), Log::ERROR);
        EXPECT_EQ(lgB.getLevel(), Log::DEBUG);
    }

    TEST(LoggingTest, UpdateCounters) {
        Log& lg = Log::getLog("TestLogging.Counter");
        Log::setWarningMaxCount("TestLogging.Counter", 5);
        Log::resetWarningCounters();
        EXPECT_TRUE(lg.isActive(Log::WARN));
    }

    TEST(LoggingTest, WarnCounterOverflow) {
        Log::setWarningMaxCount("TestLogging.Overflow", 2);
        Log& lg = Log::getLog("TestLogging.Overflow");
        for (int i = 0; i < 5; ++i) {
            lg.warn("overflow warn " + to_string(i));
        }
        EXPECT_TRUE(lg.isActive(Log::WARN));
        Log::resetWarningCounters();
    }

    TEST(LoggingTest, StreamingOperatorWarnOverflow) {
        Log::setWarningMaxCount("TestLogging.StreamOverflow", 2);
        Log& lg = Log::getLog("TestLogging.StreamOverflow");
        for (int i = 0; i < 5; ++i) {
            lg << Log::WARN << "stream overflow " << i << '\n';
        }
        EXPECT_TRUE(lg.isActive(Log::WARN));
        Log::resetWarningCounters();
    }

    TEST(LoggingTest, SetMaxWarningsWithExistingLog) {
        Log& lgX = Log::getLog("Hammer.TestX");
        Log::setWarningMaxCount("Hammer", 3);
        Log::resetWarningCounters();
        EXPECT_TRUE(lgX.isActive(Log::WARN));
    }

    TEST(LoggingTest, SetShowFlags) {
        Log::setShowTimestamp(true);
        Log::setShowLevel(false);
        Log::setShowLoggerName(false);
        Log::setUseColors(false);
        Log& lg = Log::getLog("TestLogging.Flags");
        lg.info("flags test message");
        Log::setShowTimestamp(false);
        Log::setShowLevel(true);
        Log::setShowLoggerName(true);
        Log::setUseColors(true);
    }

    TEST(LoggingTest, GetLevelFromName) {
        Log::setLevel("TestLogging.Level", Log::TRACE);
        EXPECT_EQ(Log::getLog("TestLogging.Level").getLevel(), Log::TRACE);
        Log::setLevel("TestLogging.Level", Log::DEBUG);
        EXPECT_EQ(Log::getLog("TestLogging.Level").getLevel(), Log::DEBUG);
        Log::setLevel("TestLogging.Level", Log::INFO);
        EXPECT_EQ(Log::getLog("TestLogging.Level").getLevel(), Log::INFO);
        Log::setLevel("TestLogging.Level", Log::WARN);
        EXPECT_EQ(Log::getLog("TestLogging.Level").getLevel(), Log::WARN);
        Log::setLevel("TestLogging.Level", Log::ERROR);
        EXPECT_EQ(Log::getLog("TestLogging.Level").getLevel(), Log::ERROR);
    }

    TEST(LoggingTest, StreamingOperatorInactive) {
        Log::setLevel("TestLogging.StreamInactive", Log::ERROR);
        Log& lg = Log::getLog("TestLogging.StreamInactive");
        lg << Log::DEBUG << "should be discarded" << '\n';
        EXPECT_FALSE(lg.isActive(Log::DEBUG));
    }

    TEST(LoggingTest, LevelFilteringSuppressesLower) {
        Log::setLevel("TestLogging.FilterTest", Log::ERROR);
        Log& lg = Log::getLog("TestLogging.FilterTest");
        EXPECT_FALSE(lg.isActive(Log::WARN));
        EXPECT_FALSE(lg.isActive(Log::INFO));
        EXPECT_FALSE(lg.isActive(Log::DEBUG));
        EXPECT_FALSE(lg.isActive(Log::TRACE));
        EXPECT_TRUE(lg.isActive(Log::ERROR));
        lg.warn("this should be suppressed");
        lg.info("this should be suppressed");
        lg.debug("this should be suppressed");
    }

    TEST(LoggingTest, NamedLoggerDistinctLevels) {
        Log::setLevel("TestLogging.Alpha", Log::TRACE);
        Log::setLevel("TestLogging.Beta", Log::ERROR);
        Log& lgA = Log::getLog("TestLogging.Alpha");
        Log& lgB = Log::getLog("TestLogging.Beta");
        EXPECT_EQ(lgA.getLevel(), Log::TRACE);
        EXPECT_EQ(lgB.getLevel(), Log::ERROR);
        EXPECT_TRUE(lgA.isActive(Log::TRACE));
        EXPECT_FALSE(lgB.isActive(Log::WARN));
    }

    TEST(LoggingTest, AllLevelMessages) {
        Log::setLevel("TestLogging.AllLevels", Log::TRACE);
        Log& lg = Log::getLog("TestLogging.AllLevels");
        lg.trace("trace message");
        lg.debug("debug message");
        lg.info("info message");
        lg.warn("warn message");
        lg.error("error message");
        EXPECT_TRUE(lg.isActive(Log::TRACE));
    }

    TEST(LoggingTest, GlobalLevelOverrideViaSetLevel) {
        Log::setLevel("TestLogging.Global", Log::WARNING);
        Log& lg = Log::getLog("TestLogging.Global");
        EXPECT_EQ(lg.getLevel(), Log::WARNING);
        EXPECT_TRUE(lg.isActive(Log::WARNING));
        EXPECT_TRUE(lg.isActive(Log::ERROR));
        EXPECT_FALSE(lg.isActive(Log::INFO));
        EXPECT_FALSE(lg.isActive(Log::DEBUG));

        Log::setLevel("TestLogging.Global", Log::INFO);
        EXPECT_EQ(lg.getLevel(), Log::INFO);
        EXPECT_TRUE(lg.isActive(Log::INFO));
        EXPECT_FALSE(lg.isActive(Log::DEBUG));
    }

    TEST(LoggingTest, FormatMessageWithAllFlags) {
        Log::setShowTimestamp(false);
        Log::setShowLevel(true);
        Log::setShowLoggerName(true);
        Log::setUseColors(false);
        Log& lg = Log::getLog("TestLogging.Format");
        lg.info("format message");
        Log::setShowTimestamp(true);
        lg.info("with timestamp");
        Log::setShowTimestamp(false);
        Log::setUseColors(true);
    }

    TEST(LoggingTest, StreamingAtActiveLevel) {
        Log::setLevel("TestLogging.Stream", Log::INFO);
        Log& lg = Log::getLog("TestLogging.Stream");
        lg << Log::INFO << "streamed info message" << '\n';
        lg << Log::WARN << "streamed warn message" << '\n';
        lg << Log::ERROR << "streamed error message" << '\n';
        EXPECT_TRUE(lg.isActive(Log::INFO));
    }

    TEST(LoggingTest, ParentLoggerLevelInheritance) {
        Log::setLevel("ParentLogger", Log::ERROR);
        Log& child = Log::getLog("ParentLogger.Child");
        EXPECT_EQ(child.getLevel(), Log::ERROR);
    }

} // namespace Hammer

namespace AssertTest {

    TEST(AssertHandleAssertTest, WritesToCerrAndReturns) {
        testing::internal::CaptureStderr();
        Assert::HandleAssert("test message", "x == 0", "TestLogging.cc", 42);
        string output = testing::internal::GetCapturedStderr();
        EXPECT_NE(output.find("test message"), string::npos);
        EXPECT_NE(output.find("x == 0"), string::npos);
        EXPECT_NE(output.find("TestLogging.cc"), string::npos);
        EXPECT_NE(output.find("42"), string::npos);
    }

    TEST(AssertHandleAssertTest, WritesTerminatingMessageToCerr) {
        testing::internal::CaptureStderr();
        Assert::HandleAssert("another assertion", "ptr != nullptr", "SomeFile.cc", 100);
        string output = testing::internal::GetCapturedStderr();
        EXPECT_NE(output.find("Application now terminating"), string::npos);
    }

} // namespace AssertTest
