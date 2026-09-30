///
/// @file  FFEval.cc
/// @brief Hammer FFEval class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include "Hammer/FFEval.hh"
#include "Hammer/Tools/HammerSerial.hh"
#include "Hammer/Tools/Logging.hh"
#include "Hammer/Math/Tensor.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Tools/Utils.hh"
#include "Hammer/DictionaryManager.hh"
#include "Hammer/SchemeDefinitions.hh"
#include "Hammer/ProvidersRepo.hh"
#include "Hammer/FormFactorBase.hh"
#include "Hammer/ExternalData.hh"
#include "Hammer/Math/MultiDim/IContainer.hh"

using namespace std;

namespace Hammer {

    namespace MD = MultiDimensional;
    
    FFEval::FFEval(DictionaryManager* dict) : _dictionaries{dict} {
       getLog().setWarningMaxCount("Hammer.FFEval", 1);
    }

    FFEval::~FFEval() noexcept {
        _dictionaries = nullptr;    
    }
    
    vector<double> FFEval::evalFormFactors(const string& scheme, const HashId& process, const vector<double>& point,
                                                  const vector<double>& masses) const {
        MSG_WARNING("Computing elements of internal form factor tensor. " 
                    "Consult the relevant amplitude and/or form factor classes for the definition of the form factor basis.");
        vector<double> result{};
        auto ffs = _dictionaries->providers().getFormFactor(process);
        auto defs = _dictionaries->schemeDefs().getFFSchemesForProcess(process);
        auto it = defs.find(scheme);
        if(it != defs.end()){
            auto ff = ffs[it->second];
            ff->evalAtPSPoint(point, masses);
            auto fften = ff->getTensor();
            if(fften.hasFFVarLabels()){
	           auto extptr = _dictionaries->externalData().getFFEigenVectors(ff, scheme); 
               if (extptr && extptr->rank() > 0) {
                   Tensor ext{"ext", {{extptr, false}}}; 
                   fften.dot(ext);
               }  	
            }
            if(fften.rank() > 1) {
		          MSG_ERROR("Form factor tensor is not rank 1. There is no song to sing in the evening.");
                  return result;    
            } else if (fften.rank() == 1){
                result.reserve(fften.dims()[0]);
                for(IndexType idx = 0; idx < fften.dims()[0]; ++idx){
                    result.push_back(fften.element({idx}).real());
                }
            }   
        } else {
            MSG_ERROR("Scheme '" + scheme + "' not found for requested process. There is no song to sing in the morning.");
        }
        return result;
    }
    
    Log& FFEval::getLog() {
        return Log::getLog("Hammer.FFEval");
    }

    void FFEval::defineSettings() {
        setPath("Hammer.FFEval");
    }

} // namespace Hammer
