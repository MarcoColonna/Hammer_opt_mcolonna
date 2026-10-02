#**** This file is a part of the HAMMER library
#**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
#**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
#**** Please note the MCnet academic guidelines; see GUIDELINES for details

## \file  cppdefs.pxd
#  \brief Cython class to import classes from hammer library
#

from libcpp.map cimport map
from libcpp.set cimport set as cset
from libcpp.string cimport string
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from libcpp.unordered_map cimport unordered_map
from libcpp cimport bool
from libcpp.complex cimport complex as cpp_complex

from libc.stdint cimport uintptr_t, uint32_t, uint8_t, uint16_t

cdef extern from "<array>" namespace "std" nogil:
    cdef cppclass array4 "std::array<double, 4>":
        array4() except +
        double& operator[](size_t)

    cdef cppclass array3 "std::array<double, 3>":
        array3() except +
        double& operator[](size_t)

cdef extern from "<iostream>" namespace "std":
    cdef cppclass ostream:
        ostream& write(const char*, int) except +

    cdef cppclass istream:
        istream& write(const char*, int) except +

# obviously std::ios_base isn't a namespace, but this lets
# Cython generate the correct C++ code
cdef extern from "<iostream>" namespace "std::ios_base":
    cdef cppclass open_mode:
        pass
    cdef open_mode binary
    # you can define other constants as needed

cdef extern from "<fstream>" namespace "std":
    cdef cppclass ofstream(ostream):
        # constructors
        ofstream(const char*) except +
        ofstream(const char*, open_mode) except+

    cdef cppclass ifstream(istream):
        # constructors
        ifstream(const char*) except +
        ifstream(const char*, open_mode) except+

cdef extern from "Hammer/Tools/Utils.hh" namespace "Hammer":

    cdef cppclass UMap[T, U]:
        ctypedef T key_type
        ctypedef U mapped_type
        ctypedef pair[const T, U] value_type
        cppclass iterator:
            pair[T, U]& operator*()
            iterator operator++()
            iterator operator--()
            bint operator==(iterator)
            bint operator!=(iterator)
        cppclass reverse_iterator:
            pair[T, U]& operator*()
            iterator operator++()
            iterator operator--()
            bint operator==(reverse_iterator)
            bint operator!=(reverse_iterator)
        cppclass const_iterator(iterator):
            pass
        cppclass const_reverse_iterator(reverse_iterator):
            pass
        UMap() except +
        UMap(UMap&) except +
        #UMap(key_compare&)
        U& operator[](T&)
        #UMap& operator=(UMap&)
        bint operator==(UMap&, UMap&)
        bint operator!=(UMap&, UMap&)
        bint operator<(UMap&, UMap&)
        bint operator>(UMap&, UMap&)
        bint operator<=(UMap&, UMap&)
        bint operator>=(UMap&, UMap&)
        U& at(const T&)
        const U& const_at "at"(const T&)
        iterator begin()
        const_iterator const_begin "begin"()
        void clear()
        size_t count(T&)
        bint empty()
        iterator end()
        const_iterator const_end "end"()
        pair[iterator, iterator] equal_range(T&)
        pair[const_iterator, const_iterator] const_equal_range "equal_range"(const T&)
        iterator erase(iterator)
        iterator erase(iterator, iterator)
        size_t erase(T&)
        iterator find(T&)
        const_iterator const_find "find"(T&)
        pair[iterator, bint] insert(pair[T, U]) # XXX pair[T,U]&
        iterator insert(iterator, pair[T, U]) # XXX pair[T,U]&
        iterator insert(iterator, iterator)
        #key_compare key_comp()
        iterator lower_bound(T&)
        const_iterator const_lower_bound "lower_bound"(T&)
        size_t max_size()
        reverse_iterator rbegin()
        const_reverse_iterator const_rbegin "rbegin"()
        reverse_iterator rend()
        const_reverse_iterator const_rend "rend"()
        size_t size()
        void swap(UMap&)
        iterator upper_bound(T&)
        const_iterator const_upper_bound "upper_bound"(T&)
        #value_compare value_comp()
        void max_load_factor(float)
        float max_load_factor()
        void rehash(size_t)
        void reserve(size_t)
        size_t bucket_count()
        size_t max_bucket_count()
        size_t bucket_size(size_t)
        size_t bucket(const T&)

cdef extern from "Hammer/Math/FourMomentum.hh" namespace "Hammer":

    ## wrapper of `Hammer::FourMomentum` class
    cdef cppclass FourMomentum:
        FourMomentum() except +
        FourMomentum(array4) except +
        FourMomentum(double, double, double, double) except +
        @staticmethod
        FourMomentum fromPM(double, double, double, double) except +
        @staticmethod
        FourMomentum fromPtEtaPhiM(double, double, double, double) except +
        @staticmethod
        FourMomentum fromEtaPhiME(double, double, double, double) except +
        void setE(double)
        void setPx(double)
        void setPy(double)
        void setPz(double)
        double E()
        double px()
        double py()
        double pz()
        double mass() except +
        double mass2()
        double p()
        double p2()
        double pt()
        double rapidity()
        double phi()
        double eta()
        double theta()
        array3 pVec()
        double dot(const FourMomentum&)
        double gamma() except +
        double beta() except +
        array3 boostVector() except +


cdef extern from "Hammer/Particle.hh" namespace "Hammer":

    ## wrapper of `Hammer::Particle` class
    cdef cppclass Particle:
        Particle() except +
        Particle(FourMomentum, int) except +
        Particle& setMomentum(FourMomentum)
        Particle& setPdgId(int)
        int pdgId()
        FourMomentum& momentum()
        FourMomentum& p()

cdef extern from "Hammer/Process.hh" namespace "Hammer":

    ## wrapper of `Hammer::Process` class
    cdef cppclass Process:
        Process() except +
        size_t addParticle(Particle)
        void addVertex(size_t, vector[size_t])
        void removeVertex(size_t, bool)
        vector[size_t]& getDaughtersIds(size_t) except +
        size_t getParentId(size_t)
        vector[size_t]& getSiblingsIds(size_t)
        Particle& getParticle(size_t) except +
        vector[Particle] getDaughters(size_t, bool) except +
        vector[Particle] getSiblings(size_t, bool)
        Particle& getParent(size_t)
        bool isParent(size_t)
        size_t getFirstVertex()
        size_t getId()
        cset[size_t] fullId()
        pair[Particle, vector[Particle]] getParticlesByVertex(int, vector[int])
        pair[Particle, vector[Particle]] getParticlesByVertex(string)
        size_t getVertexId(string)
        size_t numParticles(bool)

cdef extern from "Hammer/Tools/Logging.hh" namespace "Hammer":

    ## wrapper of `Hammer::Log` class
    cdef cppclass Log:
        @staticmethod
        void setLevel(string, int)
        @staticmethod
        void setLevels(map[string, int])
        @staticmethod
        void setShowTimestamp(bool)
        @staticmethod
        void setShowLevel(bool)
        @staticmethod
        void setShowLoggerName(bool)
        @staticmethod
        void setUseColors(bool)
        @staticmethod
        void setWarningMaxCount(string, int)
        @staticmethod
        void resetWarningCounters()


cdef extern from "Hammer/Tools/Utils.hh" namespace "Hammer":
    string version()

cdef extern from "Hammer/Tools/IOTypes.hh" namespace "Hammer":

    ## wrapper of `Hammer::IOBuffer` struct
    cdef cppclass IOBuffer:
        char kind
        uint32_t length
        uint8_t* start
        void init(size_t)
        void save(ostream)
        void load(istream)

    ## wrapper of `Hammer::BinContents` struct
    cdef cppclass BinContents:
        double sumWi
        double sumWi2
        size_t n

    ctypedef enum RecordType:
        UNDEFINED = 117 # b'u'
        HEADER = 98 # b'b'
        EVENT = 101 # b'e'
        HISTOGRAM = 104 # b'h'
        RATE = 114 # b'r'
        HISTOGRAM_DEFINITION = 100 # b'd'

    ## wrapper of `Hammer::IOBuffers` class
    cdef cppclass IOBuffers:
        IOBuffers() except +
        IOBuffers(IOBuffers&&)
        cppclass iterator:
            IOBuffer operator*()
            iterator operator++()
            bint operator==(iterator)
            bint operator!=(iterator)
        iterator begin()
        iterator end()
        IOBuffer& at(size_t)
        IOBuffer& front()
        IOBuffer& back()
        size_t size()
        bool empty()
        void clear()

    ## wrapper of `Hammer::HistoInfo` struct
    cdef cppclass HistoInfo:
        string name
        string scheme
        string specialization
        cset[cset[size_t]] eventGroupId


cdef extern from "Hammer/IndexTypes.hh" namespace "Hammer":

    ## wrapper of `Hammer::WTerm` enum class
    ctypedef enum WTerm:
        COMMON "Hammer::WTerm::COMMON"
        NUMERATOR "Hammer::WTerm::NUMERATOR"
        DENOMINATOR "Hammer::WTerm::DENOMINATOR"



from libcpp.memory cimport allocator

cdef extern from "boost/container/small_vector.hpp" namespace "boost::container":
    cdef cppclass SmallVector16 "boost::container::small_vector<uint16_t, 20, std::allocator<uint16_t>, void>":
        void push_back(uint16_t)
        size_t size()
        uint16_t& operator[](size_t)

ctypedef SmallVector16 IndexList



cdef extern from "Hammer/Hammer.hh" namespace "Hammer":

    ## wrapper of `Hammer::WTerm` enum class
    ctypedef enum PAction:
        ALL "Hammer::PAction::ALL"
        WEIGHTS "Hammer::PAction::WEIGHTS"
        HISTOGRAMS "Hammer::PAction::HISTOGRAMS"

    ## wrapper of `Hammer::Hammer` class
    cdef cppclass Hammer:
        Hammer() except +
        void initRun()
        void initEvent(double)
        size_t addProcess(Process)
        void removeProcess(size_t)
        void setEventHistogramBin(string, IndexList)
        void fillEventHistogram(string, vector[double])
        void setEventBaseWeight(double)
        void processEvent(PAction)
        bool loadEventWeights(IOBuffer, bool)
        IOBuffer saveEventWeights()
        bool loadRunHeader(IOBuffer, bool)
        IOBuffer saveRunHeader()
        string loadHistogramDefinition(IOBuffer, bool)
        HistoInfo loadHistogram(IOBuffer, bool)
        IOBuffers saveHistogram(string)
        IOBuffers saveHistogram(string, string)
        IOBuffers saveHistogram(string, string, string)
        IOBuffers saveHistogram(string, cset[cset[size_t]])
        IOBuffers saveHistogram(string, string, cset[cset[size_t]])
        IOBuffers saveHistogram(string, string, string, cset[cset[size_t]])
        IOBuffers saveHistogram(HistoInfo)
        bool loadRates(IOBuffer, bool)
        IOBuffer saveRates()
        void readCards(string, string)
        void saveOptionCard(string, bool)
        void saveHeaderCard(string)
        void saveReferences(string)
        void setOptions(string)
        void setHeader(string)
        void addTotalSumOfWeights(bool, bool)
        void addHistogram(string, IndexList, bool, vector[pair[double, double]])
        void addHistogram(string, vector[vector[double]], bool)
        void collapseProcessesInHistogram(string)
        void keepErrorsInHistogram(string, bool)
        void reconcileSpecializations()
        void saveEventGeneralWeights(bool)
        void reconcileGeneralWeightSetting()
        void reconcileGeneralHistogramSetting(string)
        void createWCSpecialization(string, string, vector[string])
        void setWCSpecializationOrigin(string, string, map[string, double complex])
        void setWCSpecializationOrigin(string, string, vector[double complex])
        void setWCSpecializationCoord(string, string, string, map[string, double complex])
        void setWCSpecializationCoord(string, string, string, vector[double complex])
        void setWCSpecializationBasis(string, string,  vector[ map[string,  cpp_complex[double]] ])
        void removeWCSpecialization(string, bool)
        void removeAllWCSpecializations(bool)
        cset[string] availableWCSpecializationIds()
        void applyWCSpecializationInWeights(string)
        void removeWCSpecializationInWeights(string)
        cset[string] appliedWCSpecializationsInWeights()
        void applyWCSpecializationInHistogram(string, string)
        void applyWCSpecializationInAllHistograms(string)
        void removeWCSpecializationInHistogram(string, string)
        cset[string] appliedWCSpecializationsInHistogram(string)
        ## void specializeWCInWeights(string, vector[double complex])
        ## void specializeWCInWeights(string, map[string, double complex])
        ## void specializeWCInHistogram(string, string, vector[double complex])
        ## void specializeWCInHistogram(string, string, map[string, double complex])
        void specializeFFInHistogram(string, string, string, vector[double])
        void specializeFFInHistogram(string, string, string, map[string, double])
        ## void resetSpecializationInHistogram(string)
        void removeFFSpecializationInHistogram(string)
        void createProjectedHistogram(string, string, cset[uint16_t])
        void removeHistogram(string)
        void addFFScheme(string, map[string, string])
        void setFFInputScheme(map[string, string])
        void removeFFScheme(string)
        vector[string] getFFSchemeNames()
        void showAvailableFFParams()
        void showAvailableFFParams(string)
        void includeDecay(vector[string])
        void includeDecay(string)
        void forbidDecay(vector[string])
        void forbidDecay(string)
        void addPurePSVertices(cset[string], WTerm)
        void clearPurePSVertices(WTerm)
        void setUnits(string)
        void renameFFEigenvectors(string, string, vector[string])
        void setWilsonCoefficients(string, vector[double complex], WTerm)
        void setWilsonCoefficients(string, map[string, double complex], WTerm)
        void setWilsonCoefficientsLocal(string, vector[double complex])
        void setWilsonCoefficientsLocal(string, map[string, double complex])
        void resetWilsonCoefficients(string, WTerm)
        map[string, double complex] retrieveWilsonCoefficients(string, WTerm)
        void setFFEigenvectors(string, string, vector[double])
        void setFFEigenvectors(string, string, map[string, double])
        void setFFEigenvectorsLocal(string, string, vector[double])
        void setFFEigenvectorsLocal(string, string, map[string, double])
        void resetFFEigenvectors(string, string)
        map[string, double] retrieveFFEigenvectors(string, string)
        double getWeight(string, vector[size_t], string)
        double getWeight(string, string)
        double getWeight(string, vector[vector[string]], string)
        map[size_t, double] getWeights(string, string)
        double getRate(size_t, string)
        double getRate(size_t, string, string)
        double getRate(int, vector[int], string)
        double getRate(int, vector[int], string, string)
        double getRate(string, string)
        double getRate(string, string, string)
        double getDenominatorRate(size_t)
        double getDenominatorRate(int, vector[int])
        double getDenominatorRate(string)
        vector[BinContents] getHistogram(string, string, string)
        UMap[cset[cset[size_t]], vector[BinContents]] getHistograms(string, string, string)
        cset[cset[size_t]] getHistogramEventIds(string, string, string)
        vector[vector[double]] getHistogramBinEdges(string)
        IndexList getHistogramShape(string)
        bool histogramHasUnderOverFlows(string)

