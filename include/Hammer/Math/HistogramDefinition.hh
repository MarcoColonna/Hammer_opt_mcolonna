///
/// @file  HistogramDefinition.hh
/// @brief Histogram definition class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_MATH_HISTOGRAMDEFINITION_HH
#define HAMMER_MATH_HISTOGRAMDEFINITION_HH

#include "Hammer/Tools/HammerSerial.fhh"
#include "Hammer/Tools/Utils.hh"
#include "Hammer/IndexTypes.hh"
#include "Hammer/Math/Tensor.hh"
#include "Hammer/Math/MultiDim/BinnedIndexing.hh"
#include "Hammer/Math/MultiDim/SequentialIndexing.hh"

namespace Hammer {

    namespace MD = MultiDimensional;

    class Log;

    class HistogramDefinition {

    public:

        HistogramDefinition(bool hasUnderOverFlow, const IndexList& numBins, const MD::BinEdgeList& binEdges,
                            bool shouldCompress, bool keepErrors);
        HistogramDefinition(bool hasUnderOverFlow, const IndexList& numBins, const MD::BinRangeList& ranges,
                            bool shouldCompress, bool keepErrors);
        HistogramDefinition(bool hasUnderOverFlow, const MD::BinEdgeList& binEdges, bool shouldCompress,
                            bool keepErrors);

        explicit HistogramDefinition(const Serial::FBHistoDefinition* msgreader);

        void addFixedData(const MultiDimensional::SharedTensorData& data);
        void resetFixedData();
        const Tensor& getFixedData(const EventUID& id, const std::string& specializationName,
                                   const std::string& schemeName, const LabelsList& tensorLabels) const;

        const MD::BinnedIndexing<MD::SequentialIndexing>& getIndexing() const;
        bool shouldCompress() const;
        void setCompression(bool value);
        bool shouldKeepErrors() const;
        void setKeepErrors(bool value);

        void write(flatbuffers::FlatBufferBuilder* msgwriter, const std::string& name) const;

        bool checkDefinition(const Serial::FBHistoDefinition* msgreader) const;

    private:

        /// @brief    logging facility
        /// @return   stream to be used for logging
        static Log& getLog();

        void read(const Serial::FBHistoDefinition* msgreader);

        MD::BinnedIndexing<MD::SequentialIndexing> _indexing;
        bool _compressHistograms;
        bool _keepErrors;

        using FixedKey = std::tuple<EventUID, SchemeName, SpecializationName>;
        mutable UMap<FixedKey, Tensor> _fixedData;

        std::map<IndexLabel, MultiDimensional::SharedTensorData> _fixedDataPool;
    };


} // namespace Hammer

#endif
