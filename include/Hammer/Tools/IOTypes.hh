///
/// @file  IOTypes.hh
/// @brief Declarations for Hammer IO structs
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_TOOLS_IOTYPES_HH
#define HAMMER_TOOLS_IOTYPES_HH

#include <vector>
#include <memory>
#include <iostream>
#include <set>
#include <string>
#include <cstdint>

#include "Hammer/Config/HammerConfig.hh"


namespace Hammer {

    enum RecordType : char {
        UNDEFINED = 'u',
        HEADER = 'b',
        EVENT = 'e',
        HISTOGRAM = 'h',
        RATE = 'r',
        HISTOGRAM_DEFINITION = 'd'
    };

    class IOBuffer {
    public:

        char kind = RecordType::UNDEFINED;
        uint32_t length = 0ul;
        uint8_t* start = nullptr;

        void load(std::istream& is);
        void save(std::ostream& os) const;
        void init(size_t dataSize);
    };

    std::ostream& operator<<(std::ostream& os, const IOBuffer& buf);
    std::istream& operator>>(std::istream& is, IOBuffer& buf);

    namespace Serial {

        class DetachedBuffers;

    }

    class IOBuffers {
    public:

        IOBuffers();
        explicit IOBuffers(std::unique_ptr<Serial::DetachedBuffers>&& data);
        IOBuffers(const IOBuffers& /*unused*/) = delete;
        IOBuffers& operator=(const IOBuffers& /*unused*/) = delete;
        IOBuffers(IOBuffers&& /*unused*/) noexcept;
        IOBuffers& operator=(IOBuffers&& /*unused*/) noexcept;
        ~IOBuffers();


        using iterator = std::vector<IOBuffer>::iterator;
        using const_iterator = std::vector<IOBuffer>::const_iterator;
        using reverse_iterator = std::vector<IOBuffer>::reverse_iterator;
        using const_reverse_iterator = std::vector<IOBuffer>::const_reverse_iterator;

        IOBuffer& at(size_t pos);
        [[nodiscard]] const IOBuffer& at(size_t pos) const;

        IOBuffer& operator[](size_t pos);
        const IOBuffer& operator[](size_t pos) const;

        IOBuffer& front();
        [[nodiscard]] const IOBuffer& front() const;

        IOBuffer& back();
        [[nodiscard]] const IOBuffer& back() const;

        [[nodiscard]] size_t size() const;
        [[nodiscard]] bool empty() const;

        iterator begin() noexcept;
        [[nodiscard]] const_iterator begin() const noexcept;
        [[nodiscard]] const_iterator cbegin() const noexcept;

        iterator end() noexcept;
        [[nodiscard]] const_iterator end() const noexcept;
        [[nodiscard]] const_iterator cend() const noexcept;

        reverse_iterator rbegin() noexcept;
        [[nodiscard]] const_reverse_iterator rbegin() const noexcept;
        [[nodiscard]] const_reverse_iterator crbegin() const noexcept;

        reverse_iterator rend() noexcept;
        [[nodiscard]] const_reverse_iterator rend() const noexcept;
        [[nodiscard]] const_reverse_iterator crend() const noexcept;

        void clear();

        void save(std::ostream& os) const;

    private:

        void init();


        std::vector<IOBuffer> _buffers;
        std::unique_ptr<Serial::DetachedBuffers> _pOwner;
    };

    std::ostream& operator<<(std::ostream& os, const IOBuffers& buf);

#ifdef HAVE_ROOT

    class RootIOBuffer {
    public:

        std::uint8_t kind = RecordType::UNDEFINED;
        std::uint32_t length = 0ul;
        std::uint32_t maxLength = 0ul;
        std::uint8_t* start = nullptr;

        void init(size_t maxSize);
        RootIOBuffer& operator=(const IOBuffer& other);
    };

#endif

    ///@brief  contents of a histogram bin after full contraction (real weights)
    ///        to be used to export the histogram outside Hammer
    struct BinContents {
        double sumWi;  ///< sum of weights
        double sumWi2; ///< sum of squared weights
        size_t n;      ////< number of entries in the bin
    };

    using IOHistogram = std::vector<BinContents>;

    struct HistoInfo {

        HistoInfo() = default;

        std::string name;
        std::string scheme;
        std::string specialization;
        std::set<std::set<size_t>> eventGroupId;
    };

} // namespace Hammer

#endif
