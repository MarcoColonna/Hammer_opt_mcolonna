///
/// @file  HistogramDefinition.cc
/// @brief Histogram definition class
///


#include "Hammer/Math/HistogramDefinition.hh"
#include "Hammer/Exceptions.hh"
#include "Hammer/Tools/HammerSerial.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Tools/Logging.hh"
#include "Hammer/Math/MultiDim/VectorContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    HistogramDefinition::HistogramDefinition(bool hasUnderOverFlow, const IndexList& numBins,
                                             const MD::BinEdgeList& binEdges, bool shouldCompress, bool keepErrors)
        : _indexing{numBins, binEdges, hasUnderOverFlow}, _compressHistograms{shouldCompress}, _keepErrors{keepErrors} {
    }

    HistogramDefinition::HistogramDefinition(bool hasUnderOverFlow, const IndexList& numBins,
                                             const MD::BinRangeList& ranges, bool shouldCompress, bool keepErrors)
        : _indexing{numBins, ranges, hasUnderOverFlow}, _compressHistograms{shouldCompress}, _keepErrors{keepErrors} {
    }

    HistogramDefinition::HistogramDefinition(bool hasUnderOverFlow, const MD::BinEdgeList& binEdges,
                                             bool shouldCompress, bool keepErrors)
        : _indexing{binEdges, hasUnderOverFlow}, _compressHistograms{shouldCompress}, _keepErrors{keepErrors} {
    }

    HistogramDefinition::HistogramDefinition(const Serial::FBHistoDefinition* msgreader) {
        read(msgreader);
    }

    Log& HistogramDefinition::getLog() {
        return Log::getLog("Hammer.HistogramDefinition");
    }

    void HistogramDefinition::addFixedData(const MD::SharedTensorData& data) {
        ASSERT(data->rank() == 1);
        IndexLabel lab = data->labels()[0];
        if (_fixedDataPool.find(lab) != _fixedDataPool.end()) {
            MSG_WARNING("Fixed data already present for label " + to_string(lab) + ". Overwriting. ");
        }
        _fixedDataPool[lab] = data;
    }

    void HistogramDefinition::resetFixedData() {
        _fixedDataPool.clear();
    }

    const Tensor& HistogramDefinition::getFixedData(const EventUID& id, const SchemeName& specializationName,
                                                    const SchemeName& schemeName,
                                                    const LabelsList& tensorLabels) const {
        auto it = _fixedData.find({id, schemeName, specializationName});
        if (it != _fixedData.end()) {
            return it->second;
        }
        vector<pair<MD::SharedTensorData, bool>> pTlist;
        for (const auto& elem : tensorLabels) {
            auto itt = _fixedDataPool.find(elem);
            if (itt != _fixedDataPool.end()) {
                pTlist.emplace_back(itt->second, false);
                continue;
            }
            itt = _fixedDataPool.find(static_cast<IndexLabel>(-elem));
            if (itt != _fixedDataPool.end()) {
                pTlist.emplace_back(itt->second, true);
            }
        }
        Tensor t{"FixedData", std::move(pTlist)};
        auto res = _fixedData.insert({{id, schemeName, specializationName}, t});
        if (res.second) {
            return res.first->second;
        }
        throw Error("Unable to build fixed Data for histogram");
    }

    const MD::BinnedIndexing<MD::SequentialIndexing>& HistogramDefinition::getIndexing() const {
        return _indexing;
    }

    bool HistogramDefinition::shouldCompress() const {
        return _compressHistograms;
    }

    void HistogramDefinition::setCompression(bool value) {
        _compressHistograms = value;
    }

    bool HistogramDefinition::shouldKeepErrors() const {
        return _keepErrors;
    }

    void HistogramDefinition::setKeepErrors(bool value) {
        _keepErrors = value;
    }


    void HistogramDefinition::write(flatbuffers::FlatBufferBuilder* msgwriter, const string& name) const {
        IndexList dims;
        dims.reserve(_indexing.rank());
        const auto& fullDimensions = _indexing.dims();
        copy(fullDimensions.begin(), fullDimensions.end(), back_inserter(dims));
        vector<flatbuffers::Offset<Serial::FBVecDouble>> vecsedges;
        auto edges = _indexing.edges();
        for (auto& elem : edges) {
            auto serialedges1d = msgwriter->CreateVector(elem);
            Serial::FBVecDoubleBuilder serialBuilder{*msgwriter};
            serialBuilder.add_value(serialedges1d);
            vecsedges.push_back(serialBuilder.Finish());
        }
        vector<int16_t> labels;
        labels.reserve(_fixedDataPool.size());
        vector<flatbuffers::Offset<Serial::FBSingleTensor>> vectensor;
        vectensor.reserve(_fixedDataPool.size());
        for (const auto& elem : _fixedDataPool) {
            labels.push_back(static_cast<int16_t>(elem.first));
            auto res = elem.second->write(msgwriter);
            vectensor.emplace_back(res.first.o);
        }
        auto serialedges = msgwriter->CreateVector(vecsedges);
        auto serialdims = msgwriter->CreateVector(dims);
        auto seriallabels = msgwriter->CreateVector(labels);
        auto serialtensors = msgwriter->CreateVector(vectensor);
        auto serialname = msgwriter->CreateString(name);
        Serial::FBHistoDefinitionBuilder serialDefHisto{*msgwriter};
        serialDefHisto.add_name(serialname);
        serialDefHisto.add_dims(serialdims);
        serialDefHisto.add_uoverflow(_indexing.hasUnderOverFlow());
        serialDefHisto.add_edges(serialedges);
        serialDefHisto.add_compress(_compressHistograms);
        serialDefHisto.add_keeperrors(_keepErrors);
        serialDefHisto.add_fixedlabels(seriallabels);
        serialDefHisto.add_fixeddata(serialtensors);
        auto out = serialDefHisto.Finish();
        msgwriter->Finish(out);
    }


    void HistogramDefinition::read(const Serial::FBHistoDefinition* msgreader) {
        if (msgreader != nullptr) {
            IndexList dims;
            const auto* sDims = msgreader->dims();
            if (sDims == nullptr) {
                MSG_ERROR("Missing dims on read. Call Peter, Paul and Mary.");
                return;
            }
            for (unsigned int j = 0; j < sDims->size(); ++j) {
                dims.push_back(sDims->Get(j));
            }
            MD::BinEdgeList edges;
            const auto* sEdges = msgreader->edges();
            if (sEdges == nullptr) {
                MSG_ERROR("Missing edges on read. Call Peter, Paul and Mary.");
                return;
            }
            for (unsigned int j = 0; j < sEdges->size(); ++j) {
                MD::BinValue edge;
                const auto* edgeElem = sEdges->Get(j);
                if (edgeElem == nullptr || edgeElem->value() == nullptr) {
                    continue;
                }
                const auto* sEdge = edgeElem->value();
                for (unsigned int k = 0; k < sEdge->size(); ++k) {
                    edge.push_back(sEdge->Get(k));
                }
                edges.push_back(edge);
            }
            _indexing = MD::BinnedIndexing<MD::SequentialIndexing>{dims, edges, msgreader->uoverflow()};
            _compressHistograms = msgreader->compress();
            _keepErrors = msgreader->keeperrors();
            const auto* sLabs = msgreader->fixedlabels();
            const auto* sTens = msgreader->fixeddata();
            if (sLabs == nullptr || sTens == nullptr) {
                MSG_ERROR("Missing fixedlabels or fixeddata on read. Call Peter, Paul and Mary.");
                return;
            }
            for (unsigned int j = 0; j < sLabs->size(); ++j) {
                _fixedDataPool.insert({static_cast<IndexLabel>(sLabs->Get(j)),
                                       MD::SharedTensorData{new MD::VectorContainer{sTens->Get(j)}}});
            }
        } else {
            throw InitializationError(
                "Invalid record for histogram definition brought to you by the Ministry of Silly Walks");
        }
    }

    bool HistogramDefinition::checkDefinition(const Serial::FBHistoDefinition* msgreader) const {
        if (msgreader != nullptr) {
            const auto* sDims = msgreader->dims();
            if (sDims == nullptr) {
                MSG_ERROR("Missing dims on read. Call Peter, Paul and Mary.");
                return false;
            }
            bool result = true;
            result &= _compressHistograms == msgreader->compress();
            result &= _keepErrors == msgreader->keeperrors();
            for (unsigned int j = 0; j < sDims->size() && result; ++j) {
                result &= _indexing.dims()[j] == sDims->Get(j);
            }
            const auto* sEdges = msgreader->edges();
            if (sEdges == nullptr) {
                MSG_ERROR("Missing edges on read. Call Peter, Paul and Mary.");
                return false;
            }
            for (IndexType j = 0; j < sEdges->size() && result; ++j) {
                const auto& edge = _indexing.edge(j);
                const auto* edgeElem = sEdges->Get(j);
                if (edgeElem == nullptr || edgeElem->value() == nullptr) {
                    continue;
                }
                const auto* sEdge = edgeElem->value();
                for (unsigned int k = 0; k < sEdge->size() && result; ++k) {
                    result &= isZero(fabs(edge[k] - sEdge->Get(k)));
                }
            }
            const auto* sLabs = msgreader->fixedlabels();
            const auto* sTens = msgreader->fixeddata();
            if (sLabs == nullptr || sTens == nullptr) {
                MSG_ERROR("Missing fixedlabels or fixeddata on read. Call Peter, Paul and Mary.");
                return false;
            }
            result &= sLabs->size() == _fixedDataPool.size();
            for (unsigned int j = 0; j < sLabs->size() && result; ++j) {
                auto it = _fixedDataPool.find(static_cast<IndexLabel>(sLabs->Get(j)));
                result &= it != _fixedDataPool.end();
                if (result) {
                    MD::VectorContainer t{sTens->Get(j)};
                    result &= t.compare(*(it->second));
                }
            }
            return result;
        }
        return true;
    }

} // namespace Hammer
