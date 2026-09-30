///
/// @file  Logging.cc
/// @brief Message logging routines
/// @brief based on Logging from Rivet v.2.1.1
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#include <ctime>
#include <unistd.h>

#include "Hammer/Tools/Logging.hh"
#include "Hammer/Exceptions.hh"

using namespace std;

namespace Hammer {

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wexit-time-destructors"

    Log::LogMap& Log::existingLogs() {
        static Log::LogMap x;
        return x;
    }

    Log::LevelMap& Log::defaultLevels() {
        static Log::LevelMap x;
        return x;
    }

    Log::WarningCountMap& Log::defaultMaxWarnings() {
        static Log::WarningCountMap x;
        return x;
    }

    Log::ColorCodes& Log::colorCodes() {
        static Log::ColorCodes x;
        return x;
    }

    string& Log::endColorCode() {
        static string x;
        return x;
    }

    mutex& Log::lock() {
        static mutex x;
        return x;
    }

    mutex& Log::displayLock() {
        static mutex x;
        return x;
    }

#pragma clang diagnostic pop

    bool Log::showTimestamp = false;
    bool Log::showLogLevel = true;
    bool Log::showLoggerName = true;
    bool Log::useShellColors = true;

    void Log::updateLevels() {
        for (auto& lev : defaultLevels()) {
            for (auto& logIt : existingLogs()) {
                if (logIt.first.find(lev.first) == 0) {
                    logIt.second->setLevel(lev.second);
                }
            }
        }
    }

    void Log::updateCounters() {
        if (defaultMaxWarnings().empty()) {
            defaultMaxWarnings()["Default"] = 10;
        }
        for (auto& lev : defaultMaxWarnings()) {
            if (lev.first.find("Default") == 0) {
                continue;
            }
            for (auto& logIt : existingLogs()) {
                if (logIt.first.find(lev.first) == 0) {
                    logIt.second->setWarnCount(lev.second);
                } else {
                    logIt.second->setWarnCount(10);
                }
            }
        }
    }

    Log::Log(string name)
        : _nostream(new ostream(nullptr)), _level(INFO), _warnCounter(0), _maxWarning(10), _name(std::move(name)) {
    }


    Log::Log(string name, int level)
        : _nostream(new ostream(nullptr)), _level(level), _warnCounter(0), _maxWarning(10), _name(std::move(name)) {
    }

    Log::Log(string name, int level, int maxCount)
        : _nostream(new ostream(nullptr)), _level(level), _warnCounter(0), _maxWarning(maxCount),
          _name(std::move(name)) {
    }

    void Log::setLevel(const string& name, int level) {
        {
            lock_guard<mutex> guard(lock());
            defaultLevels()[name] = level;
        }
        updateLevels();
    }

    void Log::setLevels(const LevelMap& logLevels) {
        {
            lock_guard<mutex> guard(lock());
            for (const auto& logLevel : logLevels) {
                defaultLevels()[logLevel.first] = logLevel.second;
            }
        }
        updateLevels();
    }


    void Log::setWarningMaxCount(const std::string& name, int maxCount) {
        {
            lock_guard<mutex> guard(lock());
            defaultMaxWarnings()[name] = maxCount;
        }
        updateCounters();
    }

    void Log::resetWarningCounters() {
        lock_guard<mutex> guard(lock());
        for (auto& logIt : existingLogs()) {
            logIt.second->_warnCounter = 0;
        }
    }

    Log& Log::getLog(const string& name) {
        if (existingLogs().find(name) == existingLogs().end()) {
            int level = INFO;
            // Try running through all parent classes to find an existing level
            string tmpname = name;
            bool triedAllParents = false;

            while (!triedAllParents) {
                // Is there a default level?
                if (defaultLevels().find(tmpname) != defaultLevels().end()) {
                    level = defaultLevels().find(tmpname)->second;
                    break;
                }

                // Is there already such a logger? (NB. tmpname != name)
                if (existingLogs().find(tmpname) != existingLogs().end()) {
                    level = existingLogs().find(tmpname)->second->getLevel();
                    break;
                }

                // Crop the string back to the next parent level
                size_t lastDot = tmpname.find_last_of('.');

                if (lastDot != string::npos) {
                    tmpname = tmpname.substr(0, lastDot);
                } else {
                    triedAllParents = true;
                }
            }
            int maxCount = 10;
            // now update data
            lock_guard<mutex> guard(lock());
            if (defaultMaxWarnings().find(tmpname) == defaultMaxWarnings().end()) {
                defaultMaxWarnings()[tmpname] = maxCount;
            } else {
                maxCount = defaultMaxWarnings()[tmpname];
            }

            // for (LevelMap::const_iterator l = defaultLevels().begin(); l != defaultLevels().end(); ++l) {
            //
            // }
            existingLogs().insert(make_pair(name, unique_ptr<Log>(new Log(name, level, maxCount))));
        }

        return *existingLogs()[name];
    }


    string Log::getLevelName(int level) {
        /// @todo Do the map::upper_limit thing to find nearest level...
        switch (level) {
        case TRACE:
            return "TRACE";

        case DEBUG:
            return "DEBUG";

        case INFO:
            return "INFO";

        case WARN:
            return "WARN";

        case ERROR:
            return "ERROR";

        default:
            return "";
        }

        // throw Error("Enum value was not a valid log level. How did that happen?");
    }


    string Log::getColorCode(int level) {
        if (!useShellColors) {
            return "";
        }

        // If the codes haven't been initialized, do so now.
        if (colorCodes().empty()) {
            lock_guard<mutex> guard(lock());
            // If stdout is a valid tty, try to use the appropriate codes.
            if (isatty(1) != 0) {
                /// @todo Test for VT100 compliance?
                colorCodes()[TRACE] = "\033[0;36m";
                colorCodes()[DEBUG] = "\033[0;34m";
                colorCodes()[INFO] = "\033[0;32m";
                colorCodes()[WARN] = "\033[0;33m";
                colorCodes()[ERROR] = "\033[0;31m";
                endColorCode() = "\033[0m";
            } else {
                colorCodes()[TRACE] = "";
                colorCodes()[DEBUG] = "";
                colorCodes()[INFO] = "";
                colorCodes()[WARN] = "";
                colorCodes()[ERROR] = "";
            }
        }

        // Return the appropriate code from the colour map.
        /// @todo Do the map::upper_limit thing to find nearest level...
        return colorCodes()[level];
    }


    Log::Level Log::getLevelFromName(const string& level) {
        if (level == "TRACE") {
            return TRACE;
        }

        if (level == "DEBUG") {
            return DEBUG;
        }

        if (level == "INFO") {
            return INFO;
        }

        if (level == "WARN") {
            return WARN;
        }

        if (level == "ERROR") {
            return ERROR;
        }

        throw Error("Couldn't create a log level from string '" + level + "'");
    }


    string Log::formatMessage(int level, const string& message) {
        string out;

        if (useShellColors) {
            out += getColorCode(level);
        }

        if (showLoggerName) {
            out += getName();
            out += ": ";
        }

        if (showLogLevel) {
            out += getLevelName(level);
            out += " ";
        }

        if (showTimestamp) {
            time_t rawtime;
            time(&rawtime);
            char* timestr = ctime(&rawtime);
            timestr[24] = ' ';
            out += timestr;
            out += " ";
        }

        if (useShellColors) {
            out += endColorCode();
        }

        out += " ";
        out += message;
        return out;
    }


    void Log::log(int level, const string& message) {
        if (isActive(level)) {
            int max_warning = getMaxWarning();
            if (level == WARN) {
                if (_warnCounter <= max_warning + 2) {
                    lock_guard<mutex> guard(lock());
                    ++_warnCounter;
                }
                if (_warnCounter > max_warning) {
                    return;
                }
            }
            std::string tmp = formatMessage(level, message);
            {
                lock_guard<mutex> guard(displayLock());
                cout << tmp << '\n';
                if (level == WARN && _warnCounter == max_warning) {
                    cout << formatMessage(level,
                                          "Maximum number of warnings reached. Further warning output will be "
                                          "suppressed for this class.")
                         << '\n';
                }
            }
        }
    }


    ostream& operator<<(Log& log, int level) {
        if (log.isActive(level)) {
            int max_warning = log.getMaxWarning();
            if (level == Log::WARN) {
                if (log._warnCounter <= max_warning + 2) {
                    lock_guard<mutex> guard(Hammer::Log::lock());
                    ++log._warnCounter;
                }
                if (log._warnCounter > max_warning) {
                    if (log._warnCounter == max_warning + 1) {
                        lock_guard<mutex> guard(Hammer::Log::displayLock());
                        cout << log.formatMessage(level,
                                                  "Maximum number of warnings reached. Further warning output will be "
                                                  "suppressed for this class.")
                             << '\n';
                    }
                    return *(log._nostream);
                }
            }
            {
                lock_guard<mutex> guard(Hammer::Log::displayLock());
                cout << log.formatMessage(level, "");
            }
            return cout;
        }
        return *(log._nostream);
    }

} // namespace Hammer
