///
/// @file  TestIOTypes.cc
/// @brief Unit tests for Hammer::IOBuffer and Hammer::IOBuffers
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "gtest/gtest.h"
#include "Hammer/Tools/IOTypes.hh"

#include <sstream>
#include <cstring>

using namespace std;

namespace Hammer {

    TEST(IOBufferTest, DefaultConstruction) {
        IOBuffer buf;
        EXPECT_EQ(buf.kind, static_cast<char>(RecordType::UNDEFINED));
        EXPECT_EQ(buf.length, 0u);
        EXPECT_EQ(buf.start, nullptr);
    }

    TEST(IOBufferTest, Init) {
        IOBuffer buf;
        buf.init(16);
        EXPECT_EQ(buf.length, 16u);
        EXPECT_NE(buf.start, nullptr);
        EXPECT_EQ(buf.kind, static_cast<char>(RecordType::UNDEFINED));
        delete[] buf.start;
        buf.start = nullptr;
    }

    TEST(IOBufferTest, InitReplacesExisting) {
        IOBuffer buf;
        buf.init(8);
        buf.init(32);
        EXPECT_NE(buf.start, nullptr);
        EXPECT_EQ(buf.length, 32u);
        delete[] buf.start;
        buf.start = nullptr;
    }

    TEST(IOBufferTest, SaveAndLoad) {
        IOBuffer writer;
        writer.init(4);
        writer.kind = RecordType::EVENT;
        writer.start[0] = 0xDE;
        writer.start[1] = 0xAD;
        writer.start[2] = 0xBE;
        writer.start[3] = 0xEF;

        ostringstream oss;
        writer.save(oss);

        istringstream iss(oss.str());
        IOBuffer reader;
        reader.load(iss);

        EXPECT_EQ(reader.kind, RecordType::EVENT);
        EXPECT_EQ(reader.length, 4u);
        ASSERT_NE(reader.start, nullptr);
        EXPECT_EQ(reader.start[0], 0xDE);
        EXPECT_EQ(reader.start[1], 0xAD);
        EXPECT_EQ(reader.start[2], 0xBE);
        EXPECT_EQ(reader.start[3], 0xEF);

        delete[] writer.start;
        writer.start = nullptr;
        delete[] reader.start;
        reader.start = nullptr;
    }

    TEST(IOBufferTest, StreamOutAndIn) {
        IOBuffer writer;
        writer.init(3);
        writer.kind = RecordType::HISTOGRAM;
        writer.start[0] = 1;
        writer.start[1] = 2;
        writer.start[2] = 3;

        ostringstream oss;
        oss << writer;

        istringstream iss(oss.str());
        IOBuffer reader;
        iss >> reader;

        EXPECT_EQ(reader.kind, RecordType::HISTOGRAM);
        EXPECT_EQ(reader.length, 3u);
        ASSERT_NE(reader.start, nullptr);
        EXPECT_EQ(reader.start[0], 1u);
        EXPECT_EQ(reader.start[1], 2u);
        EXPECT_EQ(reader.start[2], 3u);

        delete[] writer.start;
        writer.start = nullptr;
        delete[] reader.start;
        reader.start = nullptr;
    }

    TEST(IOBufferTest, LoadReusesSmallerBuffer) {
        IOBuffer writer;
        writer.init(4);
        writer.kind = RecordType::RATE;
        writer.start[0] = 10;
        writer.start[1] = 20;
        writer.start[2] = 30;
        writer.start[3] = 40;

        ostringstream oss;
        oss << writer;

        IOBuffer reader;
        reader.init(8);
        uint8_t* origStart = reader.start;

        istringstream iss(oss.str());
        iss >> reader;

        EXPECT_EQ(reader.start, origStart);
        EXPECT_EQ(reader.kind, RecordType::RATE);
        EXPECT_EQ(reader.start[0], 10u);

        delete[] writer.start;
        writer.start = nullptr;
        delete[] reader.start;
        reader.start = nullptr;
    }

    TEST(IOBufferTest, LoadFromEmptyStreamSetsUndefined) {
        IOBuffer reader;
        istringstream iss("");
        iss >> reader;
        EXPECT_EQ(reader.kind, static_cast<char>(RecordType::UNDEFINED));
    }

    TEST(IOBufferTest, RecordTypeValues) {
        EXPECT_EQ(RecordType::UNDEFINED, 'u');
        EXPECT_EQ(RecordType::HEADER, 'b');
        EXPECT_EQ(RecordType::EVENT, 'e');
        EXPECT_EQ(RecordType::HISTOGRAM, 'h');
        EXPECT_EQ(RecordType::RATE, 'r');
        EXPECT_EQ(RecordType::HISTOGRAM_DEFINITION, 'd');
    }

    TEST(IOBuffersTest, DefaultConstruction) {
        IOBuffers bufs;
        EXPECT_TRUE(bufs.empty());
        EXPECT_EQ(bufs.size(), 0u);
    }

    TEST(IOBuffersTest, MoveConstruction) {
        IOBuffers a;
        IOBuffers b(std::move(a));
        EXPECT_TRUE(b.empty());
    }

    TEST(IOBuffersTest, MoveAssignment) {
        IOBuffers a;
        IOBuffers b;
        b = std::move(a);
        EXPECT_TRUE(b.empty());
    }

    TEST(IOBuffersTest, SaveEmpty) {
        IOBuffers bufs;
        ostringstream oss;
        bufs.save(oss);
        EXPECT_EQ(oss.str().size(), 0u);
    }

    TEST(IOBuffersTest, StreamOutEmpty) {
        IOBuffers bufs;
        ostringstream oss;
        oss << bufs;
        EXPECT_EQ(oss.str().size(), 0u);
    }

    TEST(IOBuffersTest, IteratorsOnEmpty) {
        IOBuffers bufs;
        EXPECT_EQ(bufs.begin(), bufs.end());
        EXPECT_EQ(bufs.cbegin(), bufs.cend());
        EXPECT_EQ(bufs.rbegin(), bufs.rend());
        EXPECT_EQ(bufs.crbegin(), bufs.crend());
    }

    TEST(IOBuffersTest, ConstIteratorsOnEmpty) {
        const IOBuffers bufs;
        EXPECT_EQ(bufs.begin(), bufs.end());
        EXPECT_EQ(bufs.cbegin(), bufs.cend());
        EXPECT_EQ(bufs.rbegin(), bufs.rend());
        EXPECT_EQ(bufs.crbegin(), bufs.crend());
    }

    TEST(IOBuffersTest, ClearEmpty) {
        IOBuffers bufs;
        bufs.clear();
        EXPECT_TRUE(bufs.empty());
    }

    TEST(HistoInfoTest, DefaultFields) {
        HistoInfo hi;
        EXPECT_EQ(hi.name, "");
        EXPECT_EQ(hi.scheme, "");
        EXPECT_EQ(hi.specialization, "");
        EXPECT_TRUE(hi.eventGroupId.empty());
    }

    TEST(BinContentsTest, Fields) {
        BinContents bc{1.0, 2.0, 3u};
        EXPECT_NEAR(bc.sumWi, 1.0, 1e-10);
        EXPECT_NEAR(bc.sumWi2, 2.0, 1e-10);
        EXPECT_EQ(bc.n, 3u);
    }

    TEST(IOBufferTest, LoadFromEofStreamDeletesStart) {
        IOBuffer buf;
        buf.init(4);
        istringstream iss("");
        iss.setstate(ios::eofbit);
        iss >> buf;
        EXPECT_EQ(buf.kind, static_cast<char>(RecordType::UNDEFINED));
        EXPECT_EQ(buf.length, 0u);
        EXPECT_EQ(buf.start, nullptr);
    }

} // namespace Hammer
