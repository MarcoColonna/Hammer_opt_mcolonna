///
/// @file  FFEval.hh
/// @brief Hammer FFEval class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_FFEVAL_HH
#define HAMMER_FFEVAL_HH

#include <string>
#include <vector>

#include "Hammer/IndexTypes.hh"
#include "Hammer/Tools/HammerSerial.fhh"
#include "Hammer/Tools/SettingsConsumer.hh"

namespace Hammer {

    class Log;
    class DictionaryManager;

    /// @brief FFEval container class
    ///
    /// Computes the set of form factor values in a FF class at particular phase space point
    ///
    /// @ingroup Core
    class FFEval : public SettingsConsumer {

    public:

        /// @brief
        /// @param[in] histograms
        /// @param[in] dict
        explicit FFEval(DictionaryManager* dict = nullptr);

        FFEval(const FFEval& other) = delete;
        FFEval& operator=(const FFEval& other) = delete;
        FFEval(FFEval&& other) = delete;
        FFEval& operator=(FFEval&& other) = delete;

        ~FFEval() noexcept override;
        
        /// @brief
        /// @param[in] scheme
        /// @param[in] process
        /// @return
        [[nodiscard]] std::vector<double> evalFormFactors(const std::string& scheme, const HashId& process, const std::vector<double>& point, const std::vector<double>& masses = {}) const;

    protected:

        /// @brief logging facility
        /// @return   stream to be used for logging
        static Log& getLog();
        
        /// @brief
        void defineSettings() override;
        
    private:

        DictionaryManager* _dictionaries;

    };

} // namespace Hammer

#endif
