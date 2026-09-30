///
/// @file  Pdg.hh
/// @brief Hammer particle data class
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_TOOLS_PDG
#define HAMMER_TOOLS_PDG

#ifndef PID_FRIENDS
#define PID_FRIENDS typedef bool DummyPIDFriend
#endif

#include <cstdint>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Hammer/Tools/Pdg.fhh"
#include "Hammer/IndexTypes.hh"

namespace Hammer {

    /// @brief Hammer class for dealing with particle data
    ///
    /// Provides PDG code information retrieval, particle masses and properties.
    /// Organized as a singleton accessed by the `instance()` method
    ///
    /// @ingroup Tools
    class PID {

        PID_FRIENDS;

    private:

        PID() = default;
        ~PID();

    public:

        PID(const PID&) = delete;
        PID& operator=(const PID&) = delete;
        PID(PID&&) = delete;
        PID& operator=(PID&&) = delete;


        /// @brief convert a particle name to its PDG code.
        ///        if the name corresponds to a set of particles (e.g. "D*"),
        ///        the first is returned
        /// @param[in] name the particle name, see the code for a list of names
        [[nodiscard]] PdgId toPdgCode(const std::string& name) const;

        /// @brief return a list of PDG codes corresponding to a particle name
        ///        E.g. "D*" returns the PDG codes for \f$ D^{*0}, D^{*+} \f$ and their conjugates
        ///        while "Pi0" returns a list with only one PDG code
        /// @param[in] name the particle name, see the code for a list of names
        [[nodiscard]] std::vector<PdgId> toPdgList(const std::string& name) const;

        /// @brief particle mass from a PDG code
        /// @param[in] id the PDG code
        [[nodiscard]] double getMass(PdgId id) const;

        /// @brief particle mass from a particle name
        /// @param[in] name the particle name, see the code for a list of names
        [[nodiscard]] double getMass(const std::string& name) const;

        /// @brief set the particle mass by PDG code
        /// @param[in] id the PDG code
        /// @param[in] value the mass
        void setMass(PdgId id, double value);

        /// @brief set the particle mass by particle name
        /// @param[in] name the particle name, see the code for a list of names
        /// @param[in] value the mass
        void setMass(const std::string& name, double value);

        /// @brief particle width from a PDG code
        /// @param[in] id the PDG code
        [[nodiscard]] double getWidth(PdgId id) const;

        /// @brief particle width from a particle name
        /// @param[in] name the particle name, see the code for a list of names
        [[nodiscard]] double getWidth(const std::string& name) const;

        /// @brief set the particle width by PDG code
        /// @param[in] id the PDG code
        /// @param[in] value the width
        void setWidth(PdgId id, double value);

        /// @brief set the particle width by particle name
        /// @param[in] name the particle name, see the code for a list of names
        /// @param[in] value the width
        void setWidth(const std::string& name, double value);

        /// @brief return \f$ 2s+1 \f$ where s is the particle spin from a PDG code
        /// @param[in] id the PDG code
        [[nodiscard]] static size_t getSpinMultiplicity(PdgId id);

        /// @brief return \f$ 2s+1 \f$ where s is the particle spin from a particle name
        /// @param[in] name the particle name, see the code for a list of names
        [[nodiscard]] size_t getSpinMultiplicity(const std::string& name) const;

        /// @brief return \f$ \prod_i (2s_i+1) \f$ where \f$ s_i  \f$ is the spin of particle i
        ///       from a list of particle PDG codes
        /// @param[in] ids the PDG codes
        [[nodiscard]] static size_t getSpinMultiplicities(const std::vector<PdgId>& ids);

        /// @brief return \f$ \prod_i (2s_i+1) \f$ where \f$ s_i  \f$ is the spin of particle i
        ///       from a list of particle names
        /// @param[in] names the names
        [[nodiscard]] size_t getSpinMultiplicities(const std::vector<std::string>& names) const;

        /// @brief return \f$ 3*q \f$ where q is the particle electric charge from a PDG code
        /// @param[in] id the PDG code
        static int getThreeCharge(PdgId id);

        /// @brief return 3 times the total electric charge from a list of particle PDG codes
        /// @param[in] ids the PDG codes
        static int getThreeCharge(const std::vector<PdgId>& ids);

        /// @brief the particle lepton number from a PDG code
        /// @param[in] id the PDG code
        static int getLeptonNumber(PdgId id);

        /// @brief return the total lepton number from a list of particle PDG codes
        /// @param[in] ids the PDG codes
        static int getLeptonNumber(const std::vector<PdgId>& ids);

        /// @brief the particle baryon number from a PDG code
        /// @param[in] id the PDG code
        static int getBaryonNumber(PdgId id);

        /// @brief return the total baryon number from a list of particle PDG codes
        /// @param[in] ids the PDG codes
        static int getBaryonNumber(const std::vector<PdgId>& ids);

        /// @brief the particle lepton numbers for each flavor from a PDG code
        /// @param[in] id the PDG code
        /// @return the triplet \f$ (L_e, L_\mu, L_\tau) \f$
        static std::tuple<int, int, int> getLeptonFlavorNumber(PdgId id);

        /// @brief return the total lepton numbers for each flabor from a list of particle PDG codes
        /// @param[in] ids the PDG codes
        /// @return the triplet \f$ (L_e, L_\mu, L_\tau) \f$
        static std::tuple<int, int, int> getLeptonFlavorNumber(const std::vector<PdgId>& ids);

        /// @brief
        /// @param[in] name
        /// @param[in] hadOnly
        /// @return
        [[nodiscard]] std::vector<HashId> expandToValidVertexUIDs(const std::string& name,
                                                                  const bool& hadOnly = false, const bool& flipped = false) const;

        /// @brief
        /// @param[in] name
        /// @return
        [[nodiscard]] std::vector<std::pair<PdgId, std::vector<PdgId>>>
        expandToValidVertices(const std::string& name) const;

        /// @brief get partial Widths
        [[nodiscard]] std::map<HashId, double> getPartialWidths();

        /// @brief
        /// @return
        static PID& instance();

    protected:

        /// @brief
        /// @return
        static PID* getPIDInstance();

        /// @brief initializes the data tables
        void init();

        /// @brief add BRs
        void addBR(PdgId parent, const std::vector<PdgId>& daughters, double br);


        enum location : std::uint8_t { nj = 1, nq3, nq2, nq1, nl, nr, n, n8, n9, n10 };

    public:

        /// @brief check whether a PDG code corresponds to a hadron
        /// @param[in] id the PDG code
        static bool isHadron(PdgId id);

        /// @brief check whether a PDG code corresponds to a meson
        /// @param[in] id the PDG code
        static bool isMeson(PdgId id);

        /// @brief check whether a PDG code corresponds to a baryon
        /// @param[in] id the PDG code
        static bool isBaryon(PdgId id);

        /// @brief check whether a PDG code corresponds to a charged lepton
        /// @param[in] id the PDG code
        static bool isLepton(PdgId id);

        /// @brief check whether a PDG code corresponds to a neutrino
        /// @param[in] id the PDG code
        static bool isNeutrino(PdgId id);

    protected:

        /// @brief check whether a PDG code corresponds to a diquark
        /// @param[in] id the PDG code
        static bool isDiQuark(PdgId id);

        /// @brief check whether a PDG code corresponds to a pentaquark
        /// @param[in] id the PDG code
        static bool isPentaquark(PdgId id);

        /// @brief
        /// @param[in] id the PDG code
        /// @return
        static PdgId fundamentalID(PdgId id);

        /// @brief
        /// @param[in] id the PDG code
        /// @return
        static PdgId extraBits(PdgId id);

        /// @brief takes the absolute value of a PDG code
        /// @param[in] id the PDG code
        /// @return the code stripped by the sign
        static PdgId abspid(PdgId id);

        /// @brief extracts a digit from a PDG code
        /// @param[in] loc digit position
        /// @param[in] pid PDG code
        /// @return the digit value
        static unsigned short digit(location loc, PdgId pid);


    private:

        static PID* _thePID;

        std::map<PdgId, double> _masses;

        std::map<PdgId, double> _widths;

        std::map<HashId, std::pair<double, PdgId>> _brs;

        std::map<std::string, std::vector<PdgId>> _names;

    public:

        /// @name Charged leptons
        //@{
        static constexpr PdgId ELECTRON = 11;
        static constexpr PdgId POSITRON = -ELECTRON;
        static constexpr PdgId EMINUS = ELECTRON;
        static constexpr PdgId EPLUS = POSITRON;
        static constexpr PdgId MUON = 13;
        static constexpr PdgId ANTIMUON = -MUON;
        static constexpr PdgId TAU = 15;
        static constexpr PdgId ANTITAU = -TAU;
        //@}

        /// @name Neutrinos
        //@{
        static constexpr PdgId NU_E = 12;
        static constexpr PdgId NU_EBAR = -NU_E;
        static constexpr PdgId NU_MU = 14;
        static constexpr PdgId NU_MUBAR = -NU_MU;
        static constexpr PdgId NU_TAU = 16;
        static constexpr PdgId NU_TAUBAR = -NU_TAU;
        //@}

        /// @name Bosons
        //@{
        static constexpr PdgId PHOTON = 22;
        static constexpr PdgId GAMMA = PHOTON;
        static constexpr PdgId WBOSON = 24;
        //@}

        /// @name Nucleons
        //@{
        static constexpr PdgId PROTON = 2212;
        static constexpr PdgId ANTIPROTON = -PROTON;
        static constexpr PdgId PBAR = ANTIPROTON;
        static constexpr PdgId NEUTRON = 2112;
        static constexpr PdgId ANTINEUTRON = -NEUTRON;
        //@}

        /// @name Light mesons
        //@{
        static constexpr PdgId PI0 = 111;
        static constexpr PdgId PIPLUS = 211;
        static constexpr PdgId PIMINUS = -PIPLUS;
        static constexpr PdgId RHO0 = 113;
        static constexpr PdgId RHOPLUS = 213;
        static constexpr PdgId RHOMINUS = -RHOPLUS;
        static constexpr PdgId K0L = 130;
        static constexpr PdgId K0S = 310;
        static constexpr PdgId K0 = 311;
        static constexpr PdgId KPLUS = 321;
        static constexpr PdgId KMINUS = -KPLUS;
        static constexpr PdgId ETA = 221;
        static constexpr PdgId ETAPRIME = 331;
        static constexpr PdgId PHI = 333;
        static constexpr PdgId OMEGA = 223;
        //@}

        /// @name Charmonia
        //@{
        static constexpr PdgId ETAC = 441;
        static constexpr PdgId JPSI = 443;
        static constexpr PdgId PSI2S = 100443;
        //@}

        /// @name Charm mesons
        //@{
        static constexpr PdgId D0 = 421;
        static constexpr PdgId DPLUS = 411;
        static constexpr PdgId DMINUS = -DPLUS;
        static constexpr PdgId DSTAR = 423;
        static constexpr PdgId DSTARPLUS = 413;
        static constexpr PdgId DSTARMINUS = -DSTARPLUS;
        static constexpr PdgId DSSD0STAR = 10421;
        static constexpr PdgId DSSD0STARPLUS = 10411;
        static constexpr PdgId DSSD0STARMINUS = -DSSD0STARPLUS;
        static constexpr PdgId DSSD1STAR = 20423;
        static constexpr PdgId DSSD1STARPLUS = 20413;
        static constexpr PdgId DSSD1STARMINUS = -DSSD1STARPLUS;
        static constexpr PdgId DSSD1 = 10423;
        static constexpr PdgId DSSD1PLUS = 10413;
        static constexpr PdgId DSSD1MINUS = -DSSD1PLUS;
        static constexpr PdgId DSSD2STAR = 425;
        static constexpr PdgId DSSD2STARPLUS = 415;
        static constexpr PdgId DSSD2STARMINUS = -DSSD2STARPLUS;
        //@}
        //// @name Charm strange mesons
        //@{
        static constexpr PdgId DSPLUS = 431;
        static constexpr PdgId DSMINUS = -DSPLUS;
        static constexpr PdgId DSSTARPLUS = 433;
        static constexpr PdgId DSSTARMINUS = -DSSTARPLUS;
        static constexpr PdgId DSSDS0STARPLUS = 10431;
        static constexpr PdgId DSSDS0STARMINUS = -DSSDS0STARPLUS;
        static constexpr PdgId DSSDS1STARPLUS = 20433;
        static constexpr PdgId DSSDS1STARMINUS = -DSSDS1STARPLUS;
        static constexpr PdgId DSSDS1PLUS = 10433;
        static constexpr PdgId DSSDS1MINUS = -DSSDS1PLUS;
        static constexpr PdgId DSSDS2STARPLUS = 435;
        static constexpr PdgId DSSDS2STARMINUS = -DSSDS2STARPLUS;
        //@}

        /// @name Bottomonia
        //@{
        static constexpr PdgId ETAB = 551;
        static constexpr PdgId UPSILON1S = 553;
        static constexpr PdgId UPSILON2S = 100553;
        static constexpr PdgId UPSILON3S = 200553;
        static constexpr PdgId UPSILON4S = 300553;
        //@}

        /// @name b mesons
        //@{
        static constexpr PdgId BZERO = 511;
        static constexpr PdgId BPLUS = 521;
        static constexpr PdgId BMINUS = -BPLUS;
        static constexpr PdgId BS = 531;
        static constexpr PdgId BCPLUS = 541;
        static constexpr PdgId BCMINUS = -BCPLUS;
        //@}

        /// @name Baryons
        //@{
        static constexpr PdgId LAMBDA = 3122;
        static constexpr PdgId SIGMA0 = 3212;
        static constexpr PdgId SIGMAPLUS = 3222;
        static constexpr PdgId SIGMAMINUS = 3112;
        static constexpr PdgId LAMBDACPLUS = 4122;
        static constexpr PdgId LAMBDACMINUS = -LAMBDACPLUS;
        static constexpr PdgId LAMBDACSTAR12PLUS = 14122;
        static constexpr PdgId LAMBDACSTAR12MINUS = -LAMBDACSTAR12PLUS;
        static constexpr PdgId LAMBDACSTAR32PLUS = 4124;
        static constexpr PdgId LAMBDACSTAR32MINUS = -LAMBDACSTAR32PLUS;
        static constexpr PdgId LAMBDAB = 5122;
        static constexpr PdgId XI0 = 3322;
        static constexpr PdgId XIMINUS = 3312;
        static constexpr PdgId XIPLUS = -XIMINUS;
        static constexpr PdgId OMEGAMINUS = 3334;
        static constexpr PdgId OMEGAPLUS = -OMEGAMINUS;
        //@}
    };

} // namespace Hammer

#endif
