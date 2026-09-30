///
/// @file  Logging.hh
/// @brief Message logging routines
/// @brief based on Logging from Rivet v.2.1.1
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_LOGGING_HH
#define HAMMER_LOGGING_HH

#include <cstdint>
#include <string>
#include <map>
#include <iostream>
#include <mutex>
#include <memory>

namespace Hammer {


    /// @brief Logging class
    ///
    /// Defines basic interface for logging and displaying messages
    /// Manages a static dictionary of different named loggers
    ///
    /// @ingroup Tools
    class Log {

    public:

        friend std::ostream& operator<<(Log& log, int level);


        /// @name Type definitions
        //@{

        /// Log priority levels.
        enum Level : std::uint8_t {
            TRACE = 0,
            DEBUG = 10,
            INFO = 20,
            WARN = 30,
            WARNING = 30,
            ERROR = 40,
            CRITICAL = 50,
            ALWAYS = 50
        };

        /// Typedef for a collection of named logs.
        using LogMap = std::map<std::string, std::unique_ptr<Log>>;

        /// Typedef for a collection of named log levels.
        using LevelMap = std::map<std::string, int>;

        /// Typedef for a counting the number of warnings in order to turn off warnings above a certain threshold
        using WarningCountMap = std::map<std::string, int>;

        /// Typedef for a collection of shell color codes, accessed by log level.
        using ColorCodes = std::map<int, std::string>;

        //@}

    private:

        /// @name Static data members access functions
        //@{

        /// A static map of existing logs: we don't make more loggers than necessary.
        /// @return  the existing logs
        static LogMap& existingLogs();

        /// A static map of default log levels.
        /// @return  the default levels
        static LevelMap& defaultLevels();

        /// A static map for counting how many warnings have been issued for each logger
        /// @return  the warning counts
        static WarningCountMap& defaultMaxWarnings();

        /// A static map of shell color codes for the log levels.
        /// @return  the color codes
        static ColorCodes& colorCodes();

        /// Shell color code for the end of the log levels.
        /// @return  the string
        static std::string& endColorCode();

        /// Show timestamp?
        static bool showTimestamp;

        /// Show log level?
        static bool showLogLevel;

        /// Show logger name?
        static bool showLoggerName;

        /// Use shell colour escape codes?
        static bool useShellColors;

    protected:

        /// Mutex to access screen (protected because external streaming operator uses it)
        /// @return  the mutex
        static std::mutex& displayLock();

        /// Mutex to modify global elements (protected because external streaming operator uses it)
        /// @return  the mutex
        static std::mutex& lock();

        //@}

    public:

        /// @name Public static interface
        //@{

        /// Get a logger with the given name. The level will be taken from the
        /// "requestedLevels" static map or will be INFO by default.
        /// @param[in]  name  the name of the logger
        /// @return  a reference to the logger
        static Log& getLog(const std::string& name);

        /// Set the log levels
        /// @param[in]  name   the name of the logger
        /// @param[in]  level  the verbosity level
        static void setLevel(const std::string& name, int level);

        /// Set the log levels all at once
        /// @param[in]  logLevels   the verbosity level map
        static void setLevels(const LevelMap& logLevels);

        /// toggle whether to include the timestamp in a logging message
        /// @param[in]  showTime   the flag value
        static void setShowTimestamp(bool showTime = true) {
            std::lock_guard<std::mutex> guard(lock());
            showTimestamp = showTime;
        }

        /// toggle whether to display the verbosity level in a logging message
        /// @param[in]  showLevel   the flag value
        static void setShowLevel(bool showLevel = true) {
            std::lock_guard<std::mutex> guard(lock());
            showLogLevel = showLevel;
        }

        /// toggle whether to display the logger name in a logging message
        /// @param[in]  showName   the flag value
        static void setShowLoggerName(bool showName = true) {
            std::lock_guard<std::mutex> guard(lock());
            showLoggerName = showName;
        }

        /// toggle whether to display colorful logging messages
        /// @param[in]  useColors   the flag value
        static void setUseColors(bool useColors = true) {
            std::lock_guard<std::mutex> guard(lock());
            useShellColors = useColors;
        }

        /// set the maximum number of warnings for a given logger
        /// @param[in]  name       the logger name
        /// @param[in]  maxCount   the maximum number of warnings
        static void setWarningMaxCount(const std::string& name, int maxCount);

        /// reset the warning counters for all loggers
        static void resetWarningCounters();

        //@}

    protected:

        /// @name Constructors
        //@{

        /// Constructor by name
        /// @param[in]  name   logger name
        explicit Log(std::string name);

        /// Constructor by name and verbosity level
        /// @param[in]  name   logger name
        /// @param[in]  level  verbosity level
        Log(std::string name, int level);

        /// Constructor by name, verbosity level and warning count
        /// @param[in]  name      logger name
        /// @param[in]  level     verbosity level
        /// @param[in]  maxCount  upper limit on the number of warnings
        Log(std::string name, int level, int maxCount);

        //@}


        /// @name Internal methods
        //@{

        /// updates the verbosity levels of all the loggers
        static void updateLevels();

        /// updates the max warning numbers of all the loggers
        static void updateCounters();

        /// provide the escape string for the color associated to a give log level
        /// @param[in]  level  the log level
        /// @return the escape string
        static std::string getColorCode(int level);

        /// Get the maximum number of warnings for this logger.
        /// @return the warning cutoff
        [[nodiscard]] int getMaxWarning() const {
            return _maxWarning;
        }

        /// Set the priority level of this logger.
        /// @param[in]  level the verbosity level
        /// @return  a reference to itself
        Log& setLevel(int level) {
            std::lock_guard<std::mutex> guard(lock());
            _level = level;
            return *this;
        }

        /// Set the maximum number of warnings of this logger.
        /// @param[in]  numCount the number of warnings
        /// @return  a reference to itself
        Log& setWarnCount(int numCount) {
            std::lock_guard<std::mutex> guard(lock());
            _maxWarning = numCount;
            return *this;
        }

        /// Get a log level enum from a string.
        /// @param[in]  level   the level name
        /// @return the level enum
        static Level getLevelFromName(const std::string& level);

        /// Get the std::string representation of a log level.
        /// @param[in]  level   the level enum
        /// @return the level name
        static std::string getLevelName(int level);

        /// Get the name of this logger.
        /// @return the logger name
        [[nodiscard]] const std::string& getName() const {
            return _name;
        }

        /// Set the name of this logger.
        /// @param[in]  name   the logger name
        /// @return  a reference to itself
        Log& setName(const std::string& name) {
            std::lock_guard<std::mutex> guard(lock());
            _name = name;
            return *this;
        }

        /// Write a message at a particular level.
        /// @param[in]  level  the log level
        /// @param[in]  message   the log message
        void log(int level, const std::string& message);

        /// Turn a message string into the current log format.
        /// @param[in]  level  the log level
        /// @param[in]  message   the log message
        /// @return the formatted message
        std::string formatMessage(int level, const std::string& message);

        //@}

    public:

        /// @name Non-static public interface
        //@{

        /// Will this log level produce output on this logger at the moment?
        [[nodiscard]] bool isActive(int level) const {
            return (level >= _level);
        }

        /// log a trace message
        /// @param[in] message   the message text
        void trace(const std::string& message) {
            log(TRACE, message);
        }

        /// log a debug message
        /// @param[in] message   the message text
        void debug(const std::string& message) {
            log(DEBUG, message);
        }

        /// log an information message
        /// @param[in] message   the message text
        void info(const std::string& message) {
            log(INFO, message);
        }

        /// log a warning message
        /// @param[in] message   the message text
        void warn(const std::string& message) {
            log(WARN, message);
        }

        /// log an error message
        /// @param[in] message   the message text
        void error(const std::string& message) {
            log(ERROR, message);
        }

        /// Get the priority level of this logger.
        /// @return the verbosity level
        [[nodiscard]] int getLevel() const {
            return _level;
        }

        //@}

    protected:

        /// A null output stream, used for piping discarded output to nowhere.
        std::unique_ptr<std::ostream> const _nostream;

    private:

        /// @name Data members
        //@{

        /// Threshold level for this logger.
        int _level;

        /// number of warnings issued by this logger.
        int _warnCounter;

        /// maximum number of warnings for this logger.
        int _maxWarning;

        /// This logger's name
        std::string _name;

        //@}
    };

    /// Streaming output to a logger must have a Level/int as its first argument.
    /// @param[in]  log   the logger
    /// @param[in]  level   the log level
    /// @return the output stream
    /// @ingroup Tools
    std::ostream& operator<<(Log& log, int level);

} // namespace Hammer


// Neat CPU-conserving logging macros. Use by preference!
// NB. Only usable in classes where a getLog() method is provided
#define MSG_LVL(lvl, x)                     \
    [&]() {                                 \
        if (getLog().isActive(lvl)) {       \
            getLog() << (lvl) << x << "\n"; \
        }                                   \
    }()

#define MSG_TRACE(x) MSG_LVL(Log::TRACE, x)
#define MSG_DEBUG(x) MSG_LVL(Log::DEBUG, x)
#define MSG_INFO(x) MSG_LVL(Log::INFO, x)
#define MSG_WARNING(x) MSG_LVL(Log::WARNING, x)
#define MSG_ERROR(x) MSG_LVL(Log::ERROR, x)


#endif
