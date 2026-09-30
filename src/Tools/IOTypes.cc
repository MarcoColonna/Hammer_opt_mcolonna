///
/// @file  IOTypes.cc
/// @brief Declarations for Hammer IO structs
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <algorithm>
#include <cstdint>

#include "Hammer/Tools/IOTypes.hh"
#include "Hammer/Tools/HammerSerial.hh"

using namespace std;

namespace Hammer {

    ostream& operator<<(ostream& os, const IOBuffer& buf) {
        os.write((&buf.kind), sizeof(buf.kind));
        os.write(reinterpret_cast<const char*>(&buf.length), sizeof(buf.length));
        os.write(reinterpret_cast<const char*>(buf.start), buf.length);
        return os;
    }

    istream& operator>>(istream& is, IOBuffer& buf) {
        if (is.eof()) {
            buf.kind = RecordType::UNDEFINED;
            buf.length = 0;
            if (buf.start != nullptr) {
                delete[] buf.start;
                buf.start = nullptr;
            }
            return is;
        }
        uint32_t tmplength = 0;
        try {
            is.read((&buf.kind), sizeof(buf.kind));
            is.read(reinterpret_cast<char*>(&tmplength), sizeof(tmplength));
            if (tmplength > buf.length) {
                delete[] buf.start;
                buf.start = new uint8_t[tmplength];
                if (buf.start != nullptr) {
                    buf.length = tmplength;
                } else {
                    buf.length = 0;
                }
            }
            is.read(reinterpret_cast<char*>(buf.start), tmplength);
        } catch (ios_base::failure&) {
            tmplength = 0;
        }
        if (tmplength == 0 && buf.start != nullptr) {
            delete[] buf.start;
            buf.start = nullptr;
            buf.kind = RecordType::UNDEFINED;
            buf.length = 0;
        }
        return is;
    }

    void IOBuffer::load(std::istream& is) {
        is >> (*this);
    }

    void IOBuffer::save(std::ostream& os) const {
        os << (*this);
    }

    void IOBuffer::init(size_t dataSize) {
        if (start != nullptr) {
            delete[] start;
            start = nullptr;
            length = 0ul;
        }
        start = new uint8_t[dataSize];
        if (start != nullptr) {
            length = static_cast<uint32_t>(dataSize);
        } else {
            length = 0;
        }
        kind = RecordType::UNDEFINED;
    }

    IOBuffers::IOBuffers() : _pOwner{nullptr} {
    }

    IOBuffers::IOBuffers(unique_ptr<Serial::DetachedBuffers>&& data) : _pOwner{std::move(data)} {
        init();
    }

    IOBuffers::IOBuffers(IOBuffers&& other) noexcept
        : _buffers{std::move(other._buffers)}, _pOwner{std::move(other._pOwner)} {
    }

    IOBuffers& IOBuffers::operator=(IOBuffers&& other) noexcept {
        _pOwner = std::move(other._pOwner);
        _buffers.swap(other._buffers);
        return *this;
    }

    IOBuffers::~IOBuffers() {
        clear();
    }

    IOBuffer& IOBuffers::at(size_t pos) {
        return _buffers.at(pos);
    }

    const IOBuffer& IOBuffers::at(size_t pos) const {
        return _buffers.at(pos);
    }

    IOBuffer& IOBuffers::operator[](size_t pos) {
        return _buffers[pos];
    }

    const IOBuffer& IOBuffers::operator[](size_t pos) const {
        return _buffers[pos];
    }

    IOBuffer& IOBuffers::front() {
        return _buffers.front();
    }

    const IOBuffer& IOBuffers::front() const {
        return _buffers.front();
    }

    IOBuffer& IOBuffers::back() {
        return _buffers.back();
    }

    const IOBuffer& IOBuffers::back() const {
        return _buffers.back();
    }

    size_t IOBuffers::size() const {
        return _buffers.size();
    }

    bool IOBuffers::empty() const {
        return _buffers.empty();
    }

    IOBuffers::iterator IOBuffers::begin() noexcept {
        return _buffers.begin();
    }

    IOBuffers::const_iterator IOBuffers::begin() const noexcept {
        return _buffers.begin();
    }

    IOBuffers::const_iterator IOBuffers::cbegin() const noexcept {
        return _buffers.cbegin();
    }

    IOBuffers::iterator IOBuffers::end() noexcept {
        return _buffers.end();
    }

    IOBuffers::const_iterator IOBuffers::end() const noexcept {
        return _buffers.end();
    }

    IOBuffers::const_iterator IOBuffers::cend() const noexcept {
        return _buffers.cend();
    }

    IOBuffers::reverse_iterator IOBuffers::rbegin() noexcept {
        return _buffers.rbegin();
    }

    IOBuffers::const_reverse_iterator IOBuffers::rbegin() const noexcept {
        return _buffers.rbegin();
    }

    IOBuffers::const_reverse_iterator IOBuffers::crbegin() const noexcept {
        return _buffers.crbegin();
    }

    IOBuffers::reverse_iterator IOBuffers::rend() noexcept {
        return _buffers.rend();
    }

    IOBuffers::const_reverse_iterator IOBuffers::rend() const noexcept {
        return _buffers.rend();
    }

    IOBuffers::const_reverse_iterator IOBuffers::crend() const noexcept {
        return _buffers.crend();
    }

    void IOBuffers::clear() {
        _buffers.clear();
        _pOwner.reset();
    }

    void IOBuffers::init() {
        _buffers.clear();
        _buffers.reserve(_pOwner->size());
        for (size_t i = 0; i < _pOwner->size(); ++i) {
            IOBuffer buf;
            buf.kind = static_cast<RecordType>(_pOwner->bufferType(i));
            buf.start = const_cast<uint8_t*>(_pOwner->buffer(i).data());
            buf.length = static_cast<uint32_t>(_pOwner->buffer(i).size());
            _buffers.push_back(buf);
        }
    }

    ostream& operator<<(ostream& os, const IOBuffers& buf) {
        for (const auto& elem : buf) {
            os << elem;
        }
        return os;
    }

    void IOBuffers::save(std::ostream& os) const {
        os << (*this);
    }

#ifdef HAVE_ROOT

    RootIOBuffer& RootIOBuffer::operator=(const IOBuffer& other) {
        kind = static_cast<uint8_t>(other.kind);
        if (other.length > static_cast<uint32_t>(maxLength)) {
            throw out_of_range("Need reallocating: " + to_string(other.length) + " vs " + to_string(maxLength));
        }
        length = static_cast<uint32_t>(other.length);
        copy(other.start, other.start + other.length, start);
        return *this;
    }

    void RootIOBuffer::init(size_t maxSize) {
        if (start != nullptr) {
            delete[] start;
            start = nullptr;
            length = 0ul;
            maxLength = 0ul;
        }
        length = 0ul;
        if (maxSize > 0) {
            start = new uint8_t[maxSize];
            if (start != nullptr) {
                maxLength = static_cast<uint32_t>(maxSize);
            } else {
                maxLength = 0ul;
            }
        }
        kind = RecordType::UNDEFINED;
    }

#endif

} // namespace Hammer
