///
/// @file  Histogram.hh
/// @brief Hammer histogram class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_MATH_HISTOGRAM
#define HAMMER_MATH_HISTOGRAM

#include <cstdint>
#include <vector>
#include <string>
#include <memory>

#include "Hammer/Math/MultiDim/LabeledIndexing.hh"
#include "Hammer/Math/MultiDim/AlignedIndexing.hh"
#include "Hammer/Math/Tensor.hh"
#include "Hammer/Math/HistogramDefinition.hh"
#include "Hammer/Tools/HammerSerial.fhh"
#include "Hammer/Tools/IOTypes.hh"

namespace Hammer {

    class Log;

    /// @brief Multidimensional histogram class with Tensor as cell bins
    ///
    /// Contains the tensor version of an histogram defined by the user
    ///
    /// @ingroup Math
    class Histogram {

    protected:

        /// @brief Bin class with Tensor contents
        ///
        /// Contains an histogram bin, with value, error, etc.
        ///
        /// @ingroup Math
        class Bin {

        public:

            /// @brief constructor based on a tensor (for determining shape and labels)
            /// @param[in] t  the model tensor
            /// @param[in] withErrors whether errors should also be initialized (squared tensor weights)
            explicit Bin(const Tensor& t, bool withErrors = false);

            /// @brief serialization constructor
            /// @param[in] w the reader for the sum of tensor weights
            /// @param[in] w2 the reader for the sum of the squared tensor weights (if available)
            /// @param[in] num the number of entries
            explicit Bin(const Serial::FBTensor* w, const Serial::FBTensor* w2, size_t num);

            Bin(const Bin& other) = default;
            Bin& operator=(const Bin& other) = default;
            Bin(Bin&& other) = default;
            Bin& operator=(Bin&& other) = default;
            ~Bin() = default;

            /// @brief clears the contents of the bin (including the shapes of the tensors)
            void clear();

            /// @brief  add the values of a tensor to the current bin. Throws if the shapes do not match
            /// @param[in] t the tensor
            void addWeight(const Tensor& t);

            /// @brief  the sum of tensor weights contained in the bin
            /// @return  the weights
            const Tensor& getWeight() const;

            /// @brief the sum of (outer) squared tensor weights contained in the bin,
            ///        or an empty tensor if the errors are not tracked
            /// @return the squared weights
            const Tensor& getSquaredWeight() const;

            /// @brief the number of entries in the bin (number of times `addWeight` has been called)
            /// @return the number of events
            size_t getNumEvents() const;

            /// @brief
            /// @param[in] other
            /// @return
            Bin& operator+=(const Bin& other);

        private:

            Tensor _weight;       ///< tensor sum of weights
            Tensor _weight2;      ///< tensor sum of squared weights (in the outer square sense)
            size_t _numEvents{0}; ///< number of entries
            bool _withErrors{false};
        };

    public:

        /// @name Constructors
        //@{

        /// @brief constructor
        /// @param[in] withErrors whether to track squared weights
        Histogram();


        Histogram(const Histogram& other) = default;
        Histogram& operator=(const Histogram& other) = default;
        Histogram(Histogram&& other) = default;
        Histogram& operator=(Histogram&& other) = default;

        /// destructor
        virtual ~Histogram();

        //@}

        using DataType = std::vector<Bin>;
        using ElementType = Bin;
        using reference = ElementType&;
        using const_reference = const ElementType&;

        using iterator = DataType::iterator;
        using const_iterator = DataType::const_iterator;


        iterator begin();
        [[nodiscard]] const_iterator begin() const;

        iterator end();
        [[nodiscard]] const_iterator end() const;

        reference operator[](PositionType pos);
        const_reference operator[](PositionType pos) const;

        /// @brief  dimensionality of the histograms
        /// @return the number of components
        [[nodiscard]] size_t rank() const;

        /// @brief  dimension of a specific coordinate
        /// @param[in] index   the coordinate index
        /// @return  the dimension
        [[nodiscard]] IndexType dim(IndexType index) const;

        /// @brief  get the dimensions of all the coordinates at once
        /// @return  the list of dimensions
        [[nodiscard]] const IndexList& dims() const;

        /// @brief  get the histogram name
        /// @return  the name
        [[nodiscard]] const std::string& name() const;

        /// @brief set the value of a histogram bin based on the bin indices
        /// @param[in] value   value of the bin
        /// @param[in] rest    variadic argument with remaining indices
        template <typename... Args>
        void fillBin(const Tensor& value, Args... rest);

        /// @brief set the value of a histogram bin based on the bin indices
        /// @param[in] value   value of the bin
        /// @param[in] position list of coordinates
        void fillBin(const Tensor& value, const IndexList& position);

        /// @brief reset the histogram
        void clear();

        /// @brief reset the histogram contents but not the shape of the histogram
        void clearData();

        /// @brief get the value of a histogram bin based on an index of the bin
        /// @param[in] rest    variadic argument with remaining indices
        /// @return   value of the bin
        template <typename... Args>
        const Tensor& getWeight(Args... rest) const;

        /// @brief get the value of a histogram bin based on an index of the bin
        /// @param[in] rest    variadic argument with remaining indices
        /// @return   value of the bin
        template <typename... Args>
        const Tensor& getSquaredWeight(Args... rest) const;

        /// @brief get the value of a histogram bin based on an index of the bin
        /// @param[in] rest    variadic argument with remaining indices
        /// @return   value of the bin
        template <typename... Args>
        size_t getNumEvents(Args... rest) const;

        /// @brief
        /// @param[in] other
        /// @return
        Histogram& operator+=(const Histogram& other);

        /// @brief convert the indices into the position in the flattened
        ///        contents into vector organized as row-major
        /// @param[in] indices  index of coordinate in multidimensional histogram
        /// @return  the absolute bin position
        // PositionType indicesToPos(const IndexList& indices) const;

    protected:

        /// @brief access an element given its coordinates
        /// @param[in] indices  the list of coordinates, 0-based
        /// @return a reference to the bin
        Bin& element(const IndexList& indices);

        /// @brief access an element given its coordinates (const version)
        /// @param[in] indices  the list of coordinates, 0-based
        /// @return a reference to the bin
        [[nodiscard]] const Bin& element(const IndexList& indices) const;

        /// @brief  the number of elements (product of all the dimensions)
        [[nodiscard]] size_t numValues() const;

    public:

        /// @brief write the contents of the histogram for serialization
        /// @param[in] msgwriter  the writer object
        std::unique_ptr<Serial::FBHistogramBuilder> write(flatbuffers::FlatBufferBuilder* msgwriter) const;

        /// @brief read the contents of the histogram for serialization
        /// @param[in] msgreader the reader object
        /// @param[in] def  the location where to write
        void read(const Serial::FBHistogram* msgreader, const HistogramDefinition& def);

        /// @brief
        /// @return
        [[nodiscard]] const LabelsList& getLabelVector() const;

        /// @brief
        /// @param[in] externalData
        /// @return
        [[nodiscard]] IOHistogram evalToList(const Tensor& externalData) const;


        [[nodiscard]] std::unique_ptr<Histogram>
        collapseIntoNew(const std::string& newName, const HistogramDefinition& newDef,
                        const std::set<uint16_t>& collapsedIndexPositions) const;

    protected:

        /// @brief logging facility
        /// @return   stream to be used for logging
        static Log& getLog();

    private:

        friend std::unique_ptr<Histogram> makeHistogram(const std::string& /*name*/, const HistogramDefinition& /*def*/,
                                                        const Tensor& /*defaultValue*/);
        friend std::unique_ptr<Histogram> makeHistogram(const std::string& /*name*/, const HistogramDefinition& /*def*/,
                                                        const IndexList& /*defaultTensorDims*/,
                                                        const LabelsList& /*defaultTensorLabels*/);

        /// @brief set the dimensions of the various coordinates
        /// @param[in] def the dimensions list (its size specifies the histogram dimensionality)
        void setDefinition(const HistogramDefinition& def);

        /// @brief   give the histogram a name
        /// @param[in] name   the name
        void setName(const std::string& name);

        /// @brief set the value to be used for out-of-bounds coordinates if extrapolation is off
        /// @param[in] value the default value
        void setDefaultValue(const Tensor& value);
        void setDefaultValue(const IndexList& defaultDims, const LabelsList& defaultLabels);

        /// @brief allocate and initialize the histogram
        void initHistogram();

        const HistogramDefinition* _definition{nullptr};
        DataType _data;

        std::string _name;

        MultiDimensional::LabeledIndexing<MultiDimensional::AlignedIndexing> _defaultTensorIndexing;
        std::string _defaultTensorName;

        /// whether the grid is initialized
        bool _initialized{false};
    };

    /// @brief
    /// \overload makeHistogram
    /// @param[in] name
    /// @param[in] def
    /// @param[in] defaultValue
    /// @return
    std::unique_ptr<Histogram> makeHistogram(const std::string& name, const HistogramDefinition& def,
                                             const Tensor& defaultValue);

    /// @brief
    /// \overload makeHistogram
    /// @param[in] name
    /// @param[in] def
    /// @param[in] defaultTensorDims
    /// @param[in] defaultTensorLabels
    /// @return
    std::unique_ptr<Histogram> makeHistogram(const std::string& name, const HistogramDefinition& def,
                                             const IndexList& defaultTensorDims, const LabelsList& defaultTensorLabels);

    /// @brief
    /// @param[in] a
    /// @param[in] b
    /// @return
    Histogram operator+(const Histogram& first, const Histogram& second);

} // namespace Hammer

#include "Hammer/Math/HistogramDefs.hh"

#endif
