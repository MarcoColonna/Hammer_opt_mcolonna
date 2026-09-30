///
/// @file  Tensor.hh
/// @brief Hammer tensor class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_MATH_TENSOR
#define HAMMER_MATH_TENSOR

#include <complex>
#include <set>
#include <string>
#include <vector>

#include "Hammer/Math/MultiDimensional.fhh"
#include "Hammer/Tools/HammerSerial.fhh"

namespace Hammer {

    class Log;

    /// @brief Multidimensional tensor class with complex numbers as elements
    ///
    /// Contains a complex tensor flattened to a vector as raw-major
    ///
    /// @ingroup Math
    class Tensor {

    public:

        Tensor();

        explicit Tensor(std::string name);

        Tensor(std::string name, MultiDimensional::TensorData container);

        Tensor(std::string name, std::vector<std::pair<MultiDimensional::SharedTensorData, bool>>&& data);

        Tensor(const Tensor& other);
        Tensor& operator=(const Tensor& other);
        Tensor(Tensor&& other) noexcept;
        Tensor& operator=(Tensor&& other) noexcept;
        ~Tensor();

        /// @brief set the value of a specific tensor element
        /// @param[in] value  the value
        /// @param[in] rest   the rest of the indices (0-based)
        template <typename... Args>
        void setValue(std::complex<double> value, Args... rest);

        /// @brief returns the value of a specific tensor element
        /// @param[in] rest   the rest of the indices  (0-based)
        /// @return the element value
        template <typename... Args>
        std::complex<double> value(Args... rest) const;

        /// @brief access an element given its indices
        /// @param[in] indices  the list of indices, 0-based
        /// @return a reference to the element
        std::complex<double>& element(const IndexList& indices = {});

        /// @brief access an element given its indices (const version)
        /// @param[in] indices  the list of indices, 0-based
        /// @return a reference to the element
        std::complex<double> element(const IndexList& indices = {}) const;

        /// @brief  rank of the tensor
        /// @return the number of components
        size_t rank() const;

        /// @brief  get the dimensions of all the indices at once
        /// @return  the list of dimensions
        IndexList dims() const;

        /// @brief   get the labels of all the indices at once
        /// @return  the list of labels
        LabelsList labels() const;

        /// @brief   checks if Tensor has indices in the WC range
        /// @return  bool
        bool hasWCLabels() const;

        /// @brief   checks if Tensor has indices in the FF range
        /// @return  bool
        bool hasFFLabels() const;

        /// @brief   checks if Tensor has indices in the FF Var range
        /// @return  bool
        bool hasFFVarLabels() const;

        /// @brief  get the tensor name
        /// @return  the name
        const std::string& name() const;

        // protected:

        //     /// @brief fills specific tensor elements at once
        //     /// @param[in] values a list of (indices, value) for the elements to set
        //     /// @return a reference to itself
        //     void fillTensorElements(const std::vector<std::pair<std::vector<size_t>, std::complex<double>>>& values);

        bool isEqualTo(const Tensor& other) const;

        bool isSameLabelShape(const Tensor& other) const;

        /// @brief contract this tensor with another and stores the result in this tensor
        /// @param[in] other the othe tensor
        /// @param[in] indices restrict the indices to be contracted, identified by labels.
        ///                    If no labels are passed all the contractable indices are contracted.
        ///                    If a set of labels is passed, the indices to be contracted will be the
        ///                    intersection between the provided set and the set of contractable labels
        /// @return a reference to itself
        Tensor& dot(const Tensor& other, const UniqueLabelsList& indices = {});

        /// @brief  trace some of the indices of this tensor
        /// @return a reference to itself
        Tensor& spinSum();

        /// @brief trace over the traceable spin indices and divide by the product of the
        ///        dimensions of the traced indices (equal to \f$ 2s_i + 1 \f$)
        /// @return a reference to itself
        Tensor& spinAverage();

        /// @brief  creates a tensor with twice the rank
        ///         by multiplying the tensor with its hermitean conjugate
        /// @return a reference to itself
        Tensor& outerSquare();

        /// @brief  forces conversion of a tensor to vector type
        /// @return a reference to itself
        Tensor& toVector();


        /// @brief  forces conversion of a tensor to sparse type
        /// @return a reference to itself
        Tensor& toSparse();

        /// @brief multiply all the elements of the tensor by a real number
        /// @param[in] val the value
        /// @return a reference to itself
        Tensor& operator*=(double val);

        /// @brief multiply all the elements of the tensor by a complex number
        /// @param[in] val  the value
        /// @return a reference to itself
        Tensor& operator*=(std::complex<double> val);

        /// @brief sums another tensor to itself
        /// @param[in] other the other tensor
        /// @return a reference to itself. Throws if tensors have different shapes
        Tensor& operator+=(const Tensor& other);

        /// @brief multiply two tensors element by element and stores the result in this tensor
        /// @param[in] other the other tensor
        /// @return a reference to itself. Throws if tensors have different shapes
        Tensor& elementMultiplyBy(const Tensor& other);

        /// @brief divide two tensors element by element and stores the result in this tensor
        /// @param[in] other the other tensor
        /// @return a reference to itself. Throws if tensors have different shapes or if
        ///         a division by zero is encountered
        Tensor& elementDivideBy(const Tensor& other);

        /// @brief add a tensor of rank N-1 to a specific position in a specific coordinate
        ///        the dimension of the tensor should match those of this tensor in all the other dimensions
        /// @param[in] t  the other tensor
        /// @param[in] coord  the cooridnate, identified by its label. If multiple dimensions have the
        ///                   same label, the first one is considered.
        /// @param[in] position  the position along the dimension
        /// @return  a reference to itself. Throws if the dimensions of the two tensors are incompatible
        Tensor& addAt(const Tensor& t, IndexLabel coord, IndexType position);

        /// @brief check whether any tensor element is a NaN
        bool hasNaNs() const;

        /// @brief sets all the elements to 0
        void clearData();

        /// @brief write the contents of the tensor for serialization
        /// @param[in] msgwriter  the writer object
        /// @param[in] msg  the location where to write
        void write(flatbuffers::FlatBufferBuilder* msgwriter, flatbuffers::Offset<Serial::FBTensor>* msg) const;

        /// @brief read the contents of the tensor for serialization
        /// @param[in] msgreader the reader object
        void read(const Serial::FBTensor* msgreader);

    protected:

        /// @brief logging facility
        /// @return   stream to be used for logging
        Log& getLog() const;

    private:

        std::string _name; ///< the tensor name
        mutable MultiDimensional::TensorData _data;
    };

    /// @brief  tensor dot product
    /// \overload dot
    /// @param[in] first   the first tensor
    /// @param[in] second  the second tensor
    /// @param[in] indices restrict the indices to be contracted, identified by labels.
    ///                    If no labels are passed all the contractable indices are contracted.
    ///                    If a set of labels is passed, the indices to be contracted will be the
    ///                    intersection between the provided set and the set of contractable labels
    /// @return  the result of the dot product
    Tensor dot(const Tensor& first, const Tensor& second, const std::set<IndexLabel>& indices = {});

    /// @brief  trace a tensor
    /// @param[in] first    the tensor to be traced
    /// @return  the result of the trace
    Tensor spinSum(const Tensor& first);

    /// @brief trace a tensor over the traceable spin indices and divide by the product of the
    ///        dimensions of the traced indices (equal to \f$ 2s_i + 1 \f$)
    /// @param[in] first   the tensor to be spin averaged
    /// @return the result of the spin average
    Tensor spinAverage(const Tensor& first);

    /// @brief left multiplies a tensor by a real constant
    /// @param[in] first the tensor
    /// @param[in] val   the real number
    /// @return the result of the multiplication
    Tensor operator*(const Tensor& first, double val);

    /// @brief right multiplies a tensor by a real constant
    /// @param[in] val   the real number
    /// @param[in] first the tensor
    /// @return the result of the multiplication
    Tensor operator*(double val, const Tensor& first);

    /// @brief left multiplies a tensor by a complex constant
    /// @param[in] first the tensor
    /// @param[in] val   the complex number
    /// @return the result of the multiplication
    Tensor operator*(const Tensor& first, std::complex<double> val);

    /// @brief right multiplies a tensor by a complex constant
    /// @param[in] val   the complex number
    /// @param[in] first the tensor
    /// @return the result of the multiplication
    Tensor operator*(std::complex<double> val, const Tensor& first);

    /// @brief adds two tensors of the same rank and same dimensions
    /// @param[in] first  the first tensor
    /// @param[in] second  the second tensor
    /// @return the result of the sum.
    ///         Throws an exception if the tensors don't have the same dimensions
    Tensor operator+(const Tensor& first, const Tensor& second);

    /// @brief  creates a tensor with twice the rank
    ///         by multiplying the tensor with it's hermitean conjugate
    /// @param[in] first  the tensor
    /// @return  the result
    Tensor outerSquare(const Tensor& first);

    /// @brief multiplies two tensors of the same rank and same dimensions element by element
    /// @param[in] first  the first tensor
    /// @param[in] second  the second tensor
    /// @return the result of the multiplications
    ///         Throws an exception if the tensors don't have the same dimensions
    Tensor elementMultiply(const Tensor& first, const Tensor& second);

    /// @brief divides two tensors of the same rank and same dimensions element by element
    /// @param[in] first  the numerator tensor
    /// @param[in] second  the denominator tensor
    /// @return the result of the divisions
    ///         Throws an exception if the tensors don't have the same dimensions or if a
    ///         division by zero is encountered
    Tensor elementDivide(const Tensor& first, const Tensor& second);

} // namespace Hammer

#include "Hammer/Math/TensorDefs.hh"

#endif
