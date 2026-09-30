///
/// @file  FourMomentum.hh
/// @brief Hammer four momentum class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_MATH_FOURMOMENTUM
#define HAMMER_MATH_FOURMOMENTUM

#include <array>

#include <boost/operators.hpp>

#include "Hammer/Tools/HammerRoot.fhh"

namespace Hammer {

    /// @brief 4-momentum class
    ///
    /// Defines a 4-momentum and its basic properties and operations
    ///
    /// @ingroup Math
    class FourMomentum
        : public boost::addable<
              FourMomentum,                                                                         // fm + fm
              boost::subtractable<FourMomentum,                                                     // fm - fm
                                  boost::dividable2<FourMomentum, double,                           // fm / double
                                                    boost::multipliable2<FourMomentum, double>>>> { // fm * double,
                                                                                                    // double * fm

    public:

        /// @brief default constructor
        FourMomentum();

        /// @brief constructor from an array of numbers
        /// @param[in] other the array, defined as \f$ \left(E, p_x, p_y, p_z \right) \f$
        explicit FourMomentum(const std::array<double, 4>& other);

        /// @brief constructor from specific components
        /// @param[in] E  the energy
        /// @param[in] px  the x momentum component
        /// @param[in] py  the y momentum component
        /// @param[in] pz  the z momentum component
        FourMomentum(double E, double px, double py, double pz) noexcept;

        /// @brief  builder method based on 3-momentum and invariant mass
        /// @param[in] px  the x momentum component
        /// @param[in] py  the y momentum component
        /// @param[in] pz  the z momentum component
        /// @param[in] mass the invariant mass
        /// @return  the built 4-momentum
        static FourMomentum fromPM(double px, double py, double pz, double mass);

        /// @brief  builder method based on collider-friendly variables
        /// @param[in] pt  the transverse momentum
        /// @param[in] eta  the pseudorapidity
        /// @param[in] phi  the azimuthal angle
        /// @param[in] mass  the invariant mass
        /// @return  the built 4-momentum
        static FourMomentum fromPtEtaPhiM(double pt, double eta, double phi, double mass);

        /// @brief  builder method based on collider-friendly variables
        /// @param[in] eta  the pseudorapidity
        /// @param[in] phi  the azimuthal angle
        /// @param[in] mass  the invariant mass
        /// @param[in] E the energy
        /// @return  the built 4-momentum
        static FourMomentum fromEtaPhiME(double eta, double phi, double mass, double E);

#ifdef HAVE_ROOT

        /// @brief  builder method based on collider-friendly variables
        /// @param[in] v  the 4-vector
        /// @return  the built 4-momentum
        static FourMomentum fromRoot(const TLorentzVector& v);

        /// @brief  builder method based on collider-friendly variables
        /// @param[in] v  the 4-vector
        /// @return  the built 4-momentum
        template <class T>
        static FourMomentum fromRoot(const ROOT::Math::LorentzVector<T>& v);

#endif

        /// @brief set the 4-momentum energy
        /// @param[in] E the new energy
        /// @return  a reference to itself
        FourMomentum& setE(double E);

        /// @brief set the 4-momentum momentum x component
        /// @param[in] px the new momentum x component
        /// @return  a reference to itself
        FourMomentum& setPx(double px);

        /// @brief set the 4-momentum momentum y component
        /// @param[in] py the new momentum y component
        /// @return  a reference to itself
        FourMomentum& setPy(double py);

        /// @brief set the 4-momentum momentum z component
        /// @param[in] pz the new momentum z component
        /// @return  a reference to itself
        FourMomentum& setPz(double pz);

        /// @brief returns the energy
        [[nodiscard]] double E() const;

        /// @brief returns the momentum x component
        [[nodiscard]] double px() const;

        /// @brief returns the momentum y component
        [[nodiscard]] double py() const;

        /// @brief returns the momentum z component
        [[nodiscard]] double pz() const;

        /// @brief returns the invariant mass
        ///        if the invariant mass squared is negative
        ///        returns \f$ - \sqrt{|m^2|} \f$
        [[nodiscard]] double mass() const;

        /// @brief returns the squared invariant mass
        [[nodiscard]] double mass2() const;

        /// @brief returns the absolute value of the 3-momentum
        [[nodiscard]] double p() const;

        /// @brief returns the squared 3-momentum
        [[nodiscard]] double p2() const;

        /// @brief returns the transverse momentum
        [[nodiscard]] double pt() const;

        /// @brief returns the rapidity (along the z-axis)
        [[nodiscard]] double rapidity() const;

        /// @brief returns the azimuthal angle
        [[nodiscard]] double phi() const;

        /// @brief returns the pseudorapidity (along the z-axis)
        [[nodiscard]] double eta() const;

        /// @brief returns the polar angle (along the z-axis)
        [[nodiscard]] double theta() const;

        /// @brief returns the spatial 3-vector as an array
        [[nodiscard]] std::array<double, 3> pVec() const;

        /// @brief contracts the 4-momentum with another
        /// @param[in] v the other 4-momentum
        /// @return the dot product
        [[nodiscard]] double dot(const FourMomentum& v) const;

        /// @brief multiplies the components by a constant
        /// @param[in] a the real constant
        /// @return a reference to itself
        FourMomentum& operator*=(double a);

        /// @brief divides the components by a constant
        /// @param[in] a the real constant
        /// @return a reference to itself. Throws if division-by-zero.
        FourMomentum& operator/=(double a);

        /// @brief adds the another 4-momentum to itself
        /// @param[in] v the other 4-momentum
        /// @return a reference to itself
        FourMomentum& operator+=(const FourMomentum& v);

        /// @brief subtract the 4-momentum by another
        /// @param[in] v the other 4-momentum
        /// @return a reference to itself
        FourMomentum& operator-=(const FourMomentum& v);

        /// @brief construct a copy of itself with all the signs of the components flipped
        /// @return a new 4-vector with the components flipped
        [[nodiscard]] FourMomentum operator-() const;

        /// @brief construct a copy of itself with all the signs of the spatial components flipped
        /// @return a new 4-vector with the spatial components flipped
        [[nodiscard]] FourMomentum PFlip() const;

        /// @brief returns the \f$ \gamma = 1/\sqrt{1-\beta^2} \f$ for the Lorentz transformation to the rest frame
        /// @return the value. Throws for massless particles.
        [[nodiscard]] double gamma() const;

        /// @brief returns the value of the boost to the rest frame
        /// @return the value. Throws if zero energy.
        [[nodiscard]] double beta() const;

        /// @brief returns the boost 3-vector
        /// @return the vector. Throws if zero energy.
        [[nodiscard]] std::array<double, 3> boostVector() const;

        /// @brief returns a boosted 4-vector
        /// @param[in] v
        /// @return the boosted 4-vector.
        FourMomentum& boostToRestFrameOf(const FourMomentum& v);

        /// @brief
        /// @param[in] v
        /// @return
        FourMomentum& boostFromRestFrameOf(const FourMomentum& v);

        /// @brief
        /// @param[in] v
        /// @return
        FourMomentum& boostToRestFrameOf(const std::array<double, 3>& v);

    private:

        /// @brief
        /// @param[in] v
        /// @return
        void boostBy(const std::array<double, 3>& v);

        std::array<double, 4> _vec; ///< the contents of the 4-momentum in the notation \f$ (E, p_x, p_y, p_z) \f$.
    };

    /// @brief contracts two 4-momenta
    /// \overload dot
    /// @param[in] v the first 4-momentum
    /// @param[in] w the second 4-momentum
    double dot(const FourMomentum& v, const FourMomentum& w);

    /// @brief contracts two 4-momenta
    /// @param[in] v the first 4-momentum
    /// @param[in] w the second 4-momentum
    double operator*(const FourMomentum& v, const FourMomentum& w);

    /// @brief computes the angle between the spatial components of two 4-momenta
    /// @param[in] v the first 4-momentum
    /// @param[in] w the second 4-momentum
    [[nodiscard]] double angle(const FourMomentum& v, const FourMomentum& w);

    /// @brief computes \f$ \Delta R = \sqrt{\Delta \eta^2 + \Delta \phi^2} \f$
    ///        between two 4-momenta
    /// @param[in] v the first 4-momentum
    /// @param[in] w the second 4-momentum
    [[nodiscard]] double deltaR(const FourMomentum& v, const FourMomentum& w);

    /// @brief computes the difference between the azimuthal angles of two 4-momenta
    /// @param[in] v the first 4-momentum
    /// @param[in] w the second 4-momentum
    [[nodiscard]] double deltaPhi(const FourMomentum& v, const FourMomentum& w);

    /// @brief computes the difference between the pseudorapidities of two 4-momenta
    /// @param[in] v the first 4-momentum
    /// @param[in] w the second 4-momentum
    [[nodiscard]] double deltaEta(const FourMomentum& v, const FourMomentum& w);

    /// @brief contracts four 4-momenta with an 4D epsilon tensor.
    ///        The convetion used is \f$ \epsilon^{0123} = -1 \f$
    /// @param[in] a the first 4-momentum
    /// @param[in] b the second 4-momentum
    /// @param[in] c the third 4-momentum
    /// @param[in] d the fourth 4-momentum
    /// @return the contraction
    [[nodiscard]] double epsilon(const FourMomentum& a, const FourMomentum& b, const FourMomentum& c,
                                 const FourMomentum& d);

    [[nodiscard]] FourMomentum boostToRestFrameOf(const FourMomentum& mom, double vx, double vy, double vz);

    [[nodiscard]] FourMomentum boostToRestFrameOf(const FourMomentum& mom, const std::array<double, 3>& v);

    [[nodiscard]] FourMomentum boostToRestFrameOf(const FourMomentum& mom, const FourMomentum& v);

    /// @brief
    /// \overload dot
    /// @param[in] a the first 4-momentum
    /// @param[in] b the second 4-momentum
    [[nodiscard]] double dot(const std::array<double, 3>& a, const std::array<double, 3>& b);

    [[nodiscard]] double costheta(const std::array<double, 3>& a, const std::array<double, 3>& b);

} // namespace Hammer

#include "Hammer/Math/FourMomentumDefs.hh"

#endif
