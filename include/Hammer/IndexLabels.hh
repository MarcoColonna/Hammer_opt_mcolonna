///
/// @file  IndexLabels.hh
/// @brief Tensor indices label definitions
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

// -*- C++ -*-
#ifndef HAMMER_INDEXLABELS
#define HAMMER_INDEXLABELS

#include <cstdint>
#include <utility>
#include <string>

namespace Hammer {

    ///@brief label identifiers of tensor indices
    ///       they are used to determine which indices
    ///       can be contracted together or traced over
    ///       and which indices are quantum numbers, Wilson
    ///       coefficient indices, form factor indices, etc.
    ///       Same labels in different tensors can be contracted
    ///       together, labels with opposite sign in the same tensor
    ///       can be traced over when squaring matrix elements.
    ///       The opposite sign convention determine the hermitean
    ///       conjugate index of a given type.
    // NOLINTBEGIN(performance-enum-size)
    enum IndexLabel : int64_t {
        NONE = 0,

        SPIN_INDEX_START = 0,

        SPIN_TAUP = 1,
        SPIN_TAUP_HC = -SPIN_TAUP,
        SPIN_TAUM = 2,
        SPIN_TAUM_HC = -SPIN_TAUM,
        SPIN_MUP = 3,
        SPIN_MUP_HC = -SPIN_MUP,
        SPIN_MUM = 4,
        SPIN_MUM_HC = -SPIN_MUM,
        SPIN_EP = 5,
        SPIN_EP_HC = -SPIN_EP,
        SPIN_EM = 6,
        SPIN_EM_HC = -SPIN_EM,

        SPIN_NUTAU = 7,
        SPIN_NUTAU_HC = -SPIN_NUTAU,
        SPIN_NUTAU_BAR = 8,
        SPIN_NUTAU_BAR_HC = -SPIN_NUTAU_BAR,
        SPIN_NUMU = 9,
        SPIN_NUMU_HC = -SPIN_NUMU,
        SPIN_NUMU_BAR = 10,
        SPIN_NUMU_BAR_HC = -SPIN_NUMU_BAR,
        SPIN_NUE = 11,
        SPIN_NUE_HC = -SPIN_NUE,
        SPIN_NUE_BAR = 12,
        SPIN_NUE_BAR_HC = -SPIN_NUE_BAR,

        SPIN_GAMMA = 13,
        SPIN_GAMMA_HC = -SPIN_GAMMA,
        SPIN_DSTAR = 14,
        SPIN_DSTAR_HC = -SPIN_DSTAR,
        SPIN_DSSD1STAR = 15,
        SPIN_DSSD1STAR_HC = -SPIN_DSSD1STAR,
        SPIN_DSSD1 = 16,
        SPIN_DSSD1_HC = -SPIN_DSSD1,
        SPIN_DSSD2STAR = 17,
        SPIN_DSSD2STAR_HC = -SPIN_DSSD2STAR,

        SPIN_DSSTAR = 24,
        SPIN_DSSTAR_HC = -SPIN_DSSTAR,
        SPIN_DSSDS1STAR = 25,
        SPIN_DSSDS1STAR_HC = -SPIN_DSSDS1STAR,
        SPIN_DSSDS1 = 26,
        SPIN_DSSDS1_HC = -SPIN_DSSDS1,
        SPIN_DSSDS2STAR = 27,
        SPIN_DSSDS2STAR_HC = -SPIN_DSSDS2STAR,

        SPIN_JPSI = 51,
        SPIN_JPSI_HC = -SPIN_JPSI,

        SPIN_TAUP_REF = 101,
        SPIN_TAUP_HC_REF = -SPIN_TAUP_REF,
        SPIN_TAUM_REF = 102,
        SPIN_TAUM_HC_REF = -SPIN_TAUM_REF,
        SPIN_MUP_REF = 103,
        SPIN_MUP_HC_REF = -SPIN_MUP_REF,
        SPIN_MUM_REF = 104,
        SPIN_MUM_HC_REF = -SPIN_MUM_REF,
        SPIN_EP_REF = 105,
        SPIN_EP_HC_REF = -SPIN_EP_REF,
        SPIN_EM_REF = 106,
        SPIN_EM_HC_REF = -SPIN_EM_REF,
        SPIN_NUTAU_REF = 107,
        SPIN_NUTAU_HC_REF = -SPIN_NUTAU_REF,
        SPIN_NUTAU_BAR_REF = 108,
        SPIN_NUTAU_BAR_HC_REF = -SPIN_NUTAU_BAR_REF,
        SPIN_NUMU_REF = 109,
        SPIN_NUMU_HC_REF = -SPIN_NUMU_REF,
        SPIN_NUMU_BAR_REF = 110,
        SPIN_NUMU_BAR_HC_REF = -SPIN_NUMU_BAR_REF,
        SPIN_NUE_REF = 111,
        SPIN_NUE_HC_REF = -SPIN_NUE_REF,
        SPIN_NUE_BAR_REF = 112,
        SPIN_NUE_BAR_HC_REF = -SPIN_NUE_BAR_REF,

        SPIN_RHO = 202,
        SPIN_RHO_HC = -SPIN_RHO,
        SPIN_OMEGA = 203,
        SPIN_OMEGA_HC = -SPIN_OMEGA,

        SPIN_LB = 301,
        SPIN_LB_HC = -SPIN_LB,
        SPIN_LC = 351,
        SPIN_LC_HC = -SPIN_LC,
        SPIN_LCSTAR12 = 352,
        SPIN_LCSTAR12_HC = -SPIN_LCSTAR12,
        SPIN_LCSTAR32 = 353,
        SPIN_LCSTAR32_HC = -SPIN_LCSTAR32,

        SPIN_INDEX_END = 999,


        WC_INDEX_START = 1000,

        WILSON_BCTAUNU = 1001,
        WILSON_BCTAUNU_HC = -WILSON_BCTAUNU,

        WILSON_BCMUNU = 1002,
        WILSON_BCMUNU_HC = -WILSON_BCMUNU,

        WILSON_BCENU = 1003,
        WILSON_BCENU_HC = -WILSON_BCENU,

        WILSON_BUTAUNU = 1004,
        WILSON_BUTAUNU_HC = -WILSON_BUTAUNU,

        WILSON_BUMUNU = 1005,
        WILSON_BUMUNU_HC = -WILSON_BUMUNU,

        WILSON_BUENU = 1006,
        WILSON_BUENU_HC = -WILSON_BUENU,

        WC_INDEX_END = 1999,

        FF_INDEX_START = 2000,

        FF_BD = 2001,
        FF_BD_HC = -FF_BD,

        FF_BDSTAR = 2002,
        FF_BDSTAR_HC = -FF_BDSTAR,

        FF_BDSSD0STAR = 2003,
        FF_BDSSD0STAR_HC = -FF_BDSSD0STAR,

        FF_BDSSD1STAR = 2004,
        FF_BDSSD1STAR_HC = -FF_BDSSD1STAR,

        FF_BDSSD1 = 2005,
        FF_BDSSD1_HC = -FF_BDSSD1,

        FF_BDSSD2STAR = 2006,
        FF_BDSSD2STAR_HC = -FF_BDSSD2STAR,

        FF_BSDS = 2011,
        FF_BSDS_HC = -FF_BSDS,

        FF_BSDSSTAR = 2012,
        FF_BSDSSTAR_HC = -FF_BSDSSTAR,

        FF_BSDSSDS0STAR = 2013,
        FF_BSDSSDS0STAR_HC = -FF_BSDSSDS0STAR,

        FF_BSDSSDS1STAR = 2014,
        FF_BSDSSDS1STAR_HC = -FF_BSDSSDS1STAR,

        FF_BSDSSDS1 = 2015,
        FF_BSDSSDS1_HC = -FF_BSDSSDS1,

        FF_BSDSSDS2STAR = 2016,
        FF_BSDSSDS2STAR_HC = -FF_BSDSSDS2STAR,

        FF_LBLC = 2051,
        FF_LBLC_HC = -FF_LBLC,

        FF_LBLCSTAR12 = 2052,
        FF_LBLCSTAR12_HC = -FF_LBLCSTAR12,

        FF_LBLCSTAR32 = 2053,
        FF_LBLCSTAR32_HC = -FF_LBLCSTAR32,

        FF_BCJPSI = 2101,
        FF_BCJPSI_HC = -FF_BCJPSI,

        FF_BPI = 2201,
        FF_BPI_HC = -FF_BPI,

        FF_BRHO = 2202,
        FF_BRHO_HC = -FF_BRHO,

        FF_BOMEGA = 2203,
        FF_BOMEGA_HC = -FF_BOMEGA,

        FF_BSK = 2211,
        FF_BSK_HC = -FF_BSK,

        FF_DSSD1STAR = 2404,
        FF_DSSD1STAR_HC = -FF_DSSD1STAR,

        FF_DSSD1 = 2405,
        FF_DSSD1_HC = -FF_DSSD1,

        FF_DSSD2STAR = 2406,
        FF_DSSD2STAR_HC = -FF_DSSD2STAR,

        FF_DSSDS1STAR = 2414,
        FF_DSSDS1STAR_HC = -FF_DSSDS1STAR,

        FF_DSSDS1 = 2415,
        FF_DSSDS1_HC = -FF_DSSDS1,

        FF_DSSDS2STAR = 2416,
        FF_DSSDS2STAR_HC = -FF_DSSDS2STAR,

        FF_TAU3PI = 2501,
        FF_TAU3PI_HC = -FF_TAU3PI,

        FF_INDEX_END = 2999,

        FF_VAR_INDEX_START = 3000,

        FF_BD_VAR = 3001,
        FF_BD_VAR_HC = -FF_BD_VAR,

        FF_BDSTAR_VAR = 3002,
        FF_BDSTAR_VAR_HC = -FF_BDSTAR_VAR,

        FF_BDSSD0STAR_VAR = 3003,
        FF_BDSSD0STAR_VAR_HC = -FF_BDSSD0STAR_VAR,

        FF_BDSSD1STAR_VAR = 3004,
        FF_BDSSD1STAR_VAR_HC = -FF_BDSSD1STAR_VAR,

        FF_BDSSD1_VAR = 3005,
        FF_BDSSD1_VAR_HC = -FF_BDSSD1_VAR,

        FF_BDSSD2STAR_VAR = 3006,
        FF_BDSSD2STAR_VAR_HC = -FF_BDSSD2STAR_VAR,

        FF_BSDS_VAR = 3011,
        FF_BSDS_VAR_HC = -FF_BSDS_VAR,

        FF_BSDSSTAR_VAR = 3012,
        FF_BSDSSTAR_VAR_HC = -FF_BSDSSTAR_VAR,

        FF_BSDSSDS0STAR_VAR = 3013,
        FF_BSDSSDS0STAR_VAR_HC = -FF_BSDSSDS0STAR_VAR,

        FF_BSDSSDS1STAR_VAR = 3014,
        FF_BSDSSDS1STAR_VAR_HC = -FF_BSDSSDS1STAR_VAR,

        FF_BSDSSDS1_VAR = 3015,
        FF_BSDSSDS1_VAR_HC = -FF_BSDSSDS1_VAR,

        FF_BSDSSDS2STAR_VAR = 3016,
        FF_BSDSSDS2STAR_VAR_HC = -FF_BSDSSDS2STAR_VAR,

        FF_LBLC_VAR = 3051,
        FF_LBLC_VAR_HC = -FF_LBLC_VAR,

        FF_LBLCSTAR12_VAR = 3052,
        FF_LBLCSTAR12_VAR_HC = -FF_LBLCSTAR12_VAR,

        FF_LBLCSTAR32_VAR = 3053,
        FF_LBLCSTAR32_VAR_HC = -FF_LBLCSTAR32_VAR,

        FF_BCJPSI_VAR = 3101,
        FF_BCJPSI_VAR_HC = -FF_BCJPSI_VAR,

        FF_BPI_VAR = 3201,
        FF_BPI_VAR_HC = -FF_BPI_VAR,

        FF_BRHO_VAR = 3202,
        FF_BRHO_VAR_HC = -FF_BRHO_VAR,

        FF_BOMEGA_VAR = 3203,
        FF_BOMEGA_VAR_HC = -FF_BOMEGA_VAR,

        FF_BSK_VAR = 3211,
        FF_BSK_VAR_HC = -FF_BSK_VAR,

        FF_VAR_INDEX_END = 3999,

        INTEGRATION_INDEX = 4000,

        SPECIALIZATION_PAD = 32768

    };
    // NOLINTEND(performance-enum-size)

    /// @brief determine if a given label is corresponds to a spin index (or reference spin index)
    /// @param[in] val the label
    /// @return  true if label corresponds to a (reference) spin index
    bool isSpinIndex(IndexLabel val);

    /// @brief determine if a given label is corresponds to a form factor index
    /// @param[in] val the label
    /// @return  true if label corresponds to a form factor index
    bool isFFIndex(IndexLabel val);

    /// @brief determine if a given label is corresponds to a Wilson Coefficient index
    /// @param[in] val the label
    /// @return  true if label corresponds to a WC index
    bool isWCIndex(IndexLabel val);

    /// @brief determine if a given label is corresponds to a specialized index
    /// @param[in] val the label
    /// @return  true if label corresponds to a specialized index
    bool isSpecializedIndex(IndexLabel val);


    /// @brief break down a specialized label index into its components
    /// @param[in] val the label
    /// @return  a pair with the pad and the base index label
    std::pair<uint32_t, int> getSpecializedLabelComponents(IndexLabel val);

    IndexLabel specializeLabel(IndexLabel val, uint32_t pad = 1ul);

    uint32_t generateSpecializedPad(const std::string& id);

} // namespace Hammer

#endif
