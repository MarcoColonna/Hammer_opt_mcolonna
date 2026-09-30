///
/// @file  AmplBToQLepNuBase.cc
/// @brief \f$ b -> c \tau\nu \f$ base amplitude
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/Amplitudes/AmplBToQLepNuBase.hh"
#include "Hammer/Math/MultiDim/IContainer.hh"
#include "Hammer/IndexLabels.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;

    AmplBToQLepNuBase::AmplBToQLepNuBase()
        : _perms{{6, 9, 5, 4, 7, 3, 2}}, _flips{{1., -1., -1., -1., -1., 1., 1., 1., 1., -1., -1.}} {
    }

    void AmplBToQLepNuBase::updateWilsonCoeffLabelPrefix() {
        _mWCLabel = NONE;
        _mWCPrefix = "None";
        _multiplicity = 1ul;
        for (auto elem : getTensor().labels()) {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wswitch-enum"
            switch (elem) {
            case WILSON_BCTAUNU: {
                _mWCLabel = WILSON_BCTAUNU;
                _mWCPrefix = "BtoCTauNu";
                break;
            }
            case WILSON_BCMUNU: {
                _mWCLabel = WILSON_BCMUNU;
                _mWCPrefix = "BtoCMuNu";
                break;
            }
            case WILSON_BCENU: {
                _mWCLabel = WILSON_BCENU;
                _mWCPrefix = "BtoCENu";
                break;
            }
            case WILSON_BUTAUNU: {
                _mWCLabel = WILSON_BUTAUNU;
                _mWCPrefix = "BtoUTauNu";
                break;
            }
            case WILSON_BUMUNU: {
                _mWCLabel = WILSON_BUMUNU;
                _mWCPrefix = "BtoUMuNu";
                break;
            }
            case WILSON_BUENU: {
                _mWCLabel = WILSON_BUENU;
                _mWCPrefix = "BtoUENu";
                break;
            }
            default:
                break;
            }
#pragma clang diagnostic pop
        }
    }


    void AmplBToQLepNuBase::defineSettings() {
        //_mWCNames = {"SM",  "S_aRbL", "S_aRbR", "S_aLbL", "S_aLbR", "V_aRbL", "V_aRbR", "V_aLbL", "V_aLbR", "T_aRbL",
        //"T_aLbR"};
        _mWCNames = {"SM",     "S_qLlL", "S_qRlL", "V_qLlL", "V_qRlL", "T_qLlL",
                     "S_qLlR", "S_qRlR", "V_qLlR", "V_qRlR", "T_qRlR"};
        setPath(_mWCPrefix);
        addSetting<complex<double>>("SM", 1.0);
        for (const auto& elem : _mWCNames) {
            if (elem != "SM") {
                addSetting<complex<double>>(elem, 0.);
            }
        }
    }

    void AmplBToQLepNuBase::preProcessWCValues(vector<complex<double>>& data, bool reverse) const {
        // getWC order: conj( {"SM", "S_qLlL", "S_qRlL", "V_qLlL", "V_qRlL", "T_qLlL", "S_qLlR", "S_qRlR", "V_qLlR",
        // "V_qRlR", "T_qRlR"} ) order for amplitudes: {SM, aSR bSL, aSR bSR, aSL bSL, aSL bSR, aVR bVL, aVR bVR, aVL
        // bVL,  aVL bVR, aTR bTL, aTL bTR}
        if (!reverse) {
            auto it = _perms.begin();
            MD::IContainer::ElementType tmp = data[*it];
            IndexType previous = *it;
            ++it;
            for (; it != _perms.end(); ++it) {
                data[previous] = data[*it];
                previous = *it;
            }
            data[previous] = tmp;
            for (IndexType i = 0; i < static_cast<IndexType>(_flips.size()); ++i) {
                data[i] = _flips[i] * conj(data[i]);
            }
            // auto temp = ext.data();
            // ext.data() = {conj(temp[0]), -conj(temp[1]), -conj(temp[6]), -conj(temp[2]), -conj(temp[7]),
            // conj(temp[4]),
            //               conj(temp[9]), conj(temp[3]),  conj(temp[8]),  -conj(temp[5]), -conj(temp[10])};

        } else {
            for (IndexType i = 0; i < static_cast<IndexType>(_flips.size()); ++i) {
                data[i] = _flips[i] * conj(data[i]);
            }
            auto it = _perms.rbegin();
            MD::IContainer::ElementType tmp = data[*it];
            IndexType previous = *it;
            ++it;
            for (; it != _perms.rend(); ++it) {
                data[previous] = data[*it];
                previous = *it;
            }
            data[previous] = tmp;
        }
    }


} // namespace Hammer
