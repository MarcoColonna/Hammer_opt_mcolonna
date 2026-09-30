///
/// @file  HistogramSet.hh
/// @brief Container class for histograms belonging to different event types
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_HISTOGRAMSET_HH
#define HAMMER_HISTOGRAMSET_HH

#include <cstdint>
#include <map>
#include <vector>
#include <memory>

#include "Hammer/Tools/HammerSerial.fhh"
#include "Hammer/Tools/Utils.hh"
#include "Hammer/Tools/IOTypes.hh"
#include "Hammer/IndexTypes.hh"

#include "Hammer/Math/Tensor.hh"

namespace Hammer {

    class Histogram;
    class HistogramDefinition;

    class HistogramSet {

    public:

        explicit HistogramSet(bool compressed = false);
        HistogramSet(const HistogramSet& other, const std::string& newName, const HistogramDefinition& newDef,
                     const std::set<uint16_t>& collapsedIndexPositions);

        Histogram* getHistogram(const EventUID& id);
        Histogram* getHistogram(size_t id);
        size_t addEventId(const EventUID& id, const LabelsList& labels);
        Histogram* addHistogram(size_t id, std::unique_ptr<Histogram> hist);

        EventUIDGroup read(const Serial::FBHistogram* msgreader, const HistogramDefinition& def, bool merge);
        std::unique_ptr<Serial::FBHistogramBuilder> write(flatbuffers::FlatBufferBuilder* msgwriter,
                                                          const EventUID& id) const;

        [[nodiscard]] EventIdGroupDict<IOHistogram>
        specializeEventHistograms(const std::vector<Tensor>& externalData) const;
        [[nodiscard]] IOHistogram specializeSumHistogram(const std::vector<Tensor>& externalData) const;

        [[nodiscard]] EventUIDGroup getEventIdsInHistogram() const;
        [[nodiscard]] std::vector<LabelsList> getHistogramLabels() const;
        [[nodiscard]] std::vector<EventUID> getEventUIDRepresentatives() const;

        void clear();

    private:

        EventIdDict<size_t> _compressionDict;
        std::map<size_t, EventUIDGroup> _compressionReverseDict;

        std::map<size_t, LabelsList> _labelsDict;
        UMap<LabelsList, size_t> _labelsReverseDict;

        std::vector<std::unique_ptr<Histogram>> _data;

        bool _isCompressed;
    };

} // namespace Hammer

#endif
