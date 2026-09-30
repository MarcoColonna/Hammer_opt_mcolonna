message(STATUS "Preparing Hammer FormFactors library.")

add_library(FormFactors)
target_sources(FormFactors
        PRIVATE
        AsCoeffs.cc

        #
        ./BLPR/FFBLPRBase.cc
        ./BLPR/FFBLPRVarBase.cc
        ./BLPR/FFBtoDBLPR.cc
        ./BLPR/FFBtoDBLPRVar.cc
        ./BLPR/FFBtoDstarBLPR.cc
        ./BLPR/FFBtoDstarBLPRVar.cc

        #
        ./BLPRXP/FFBLPRXPBase.cc
        ./BLPRXP/FFBLPRXPVarBase.cc
        ./BLPRXP/FFBtoDBLPRXP.cc
        ./BLPRXP/FFBtoDBLPRXPVar.cc
        ./BLPRXP/FFBtoDstarBLPRXP.cc
        ./BLPRXP/FFBtoDstarBLPRXPVar.cc

        #
        ./BLRS/FFBLRSBase.cc
        ./BLRS/FFLbtoLcBLRS.cc
        ./BLRS/FFLbtoLcBLRSVar.cc

        #
        ./BLRSXP/FFBLRSXPBase.cc
        ./BLRSXP/FFLbtoLcBLRSXP.cc
        ./BLRSXP/FFLbtoLcBLRSXPVar.cc

        #
        ./BLR/FFBLRBase.cc
        ./BLR/FFBtoD0starBLR.cc
        ./BLR/FFBtoD0starBLRVar.cc
        ./BLR/FFBtoD1starBLR.cc
        ./BLR/FFBtoD1starBLRVar.cc
        ./BLR/FFBtoD1BLR.cc
        ./BLR/FFBtoD1BLRVar.cc
        ./BLR/FFBtoD2starBLR.cc
        ./BLR/FFBtoD2starBLRVar.cc

        #
        ./LLSW/FFLLSWBase.cc
        ./LLSW/FFBtoD0starLLSW.cc
        ./LLSW/FFBtoD1starLLSW.cc
        ./LLSW/FFBtoD1LLSW.cc
        ./LLSW/FFBtoD2starLLSW.cc

        #
        ./ISGW2/FFISGW2Base.cc
        ./ISGW2/FFBtoDISGW2.cc
        ./ISGW2/FFBtoDstarISGW2.cc
        ./ISGW2/FFBtoD0starISGW2.cc
        ./ISGW2/FFBtoD1starISGW2.cc
        ./ISGW2/FFBtoD1ISGW2.cc
        ./ISGW2/FFBtoD2starISGW2.cc
        ./ISGW2/FFBtoRhoOmegaISGW2.cc
        ./ISGW2/FFBtoPiISGW2.cc

        #
        ./CLN/FFCLNBase.cc
        ./CLN/FFBtoDstarCLN.cc
        ./CLN/FFBtoDCLN.cc
        ./CLN/FFBtoDstarCLNVar.cc
        ./CLN/FFBtoDCLNVar.cc

        #
        ./BGL/FFBGLBase.cc
        ./BGL/FFBtoDBGL.cc
        ./BGL/FFBtoDstarBGL.cc
        ./BGL/FFBtoDBGLVar.cc
        ./BGL/FFBtoDstarBGLVar.cc
        ./BGL/FFBctoJpsiBGL.cc
        ./BGL/FFBctoJpsiBGLVar.cc
	    ./BGL/FFBGLXBase.cc
        ./BGL/FFBtoDstarBGLX.cc
        ./BGL/FFBtoDstarBGLXVar.cc

        #
        ./PCR/FFPCRBase.cc
        ./PCR/FFLbtoLcPCR.cc
        ./PCR/FFLbtoLcstar12PCR.cc
        ./PCR/FFLbtoLcstar32PCR.cc

        #
        ./LSPR/FFLSPRBase.cc
        ./LSPR/FFLSPRVarBase.cc
        ./LSPR/FFLbtoLcstar12LSPR.cc
        ./LSPR/FFLbtoLcstar32LSPR.cc
        ./LSPR/FFLbtoLcstar12LSPRVar.cc
        ./LSPR/FFLbtoLcstar32LSPRVar.cc

        #
        ./EFG/FFEFGBase.cc
        ./EFG/FFBctoJpsiEFG.cc

        #
        ./Kiselev/FFKiselevBase.cc
        ./Kiselev/FFBctoJpsiKiselev.cc

        #
        ./BCL/FFBCLBase.cc
        ./BCL/FFBtoPiBCL.cc
        ./BCL/FFBtoPiBCLVar.cc
        ./BCL/FFBstoKBCL.cc
        ./BCL/FFBstoKBCLVar.cc
	./BCL/FFBtoDBCL.cc
	./BCL/FFBtoDBCLVar.cc
        ./BCL/FFBstoDsBCL.cc
        ./BCL/FFBstoDsBCLVar.cc

        #
        ./BSZ/FFBSZBase.cc
        ./BSZ/FFBtoRhoBSZ.cc
        ./BSZ/FFBtoOmegaBSZ.cc
        ./BSZ/FFBtoRhoBSZVar.cc
        ./BSZ/FFBtoOmegaBSZVar.cc

        #
        ./GKvD/FFBtoPiGKvD.cc

        #
        ./PW/FFPWBase.cc
        ./PW/FFD1DstarPiPW.cc
        ./PW/FFD1starDstarPiPW.cc
        ./PW/FFD2starDstarPiPW.cc
        ./PW/FFD2starDPiPW.cc

        #
        ./RCT/FFRCTBase.cc
        ./RCT/FFTauto3PiRCT.cc
)

set(FormFactors_INCLUDE_DIR "../../include/Hammer/FormFactors")

target_sources(FormFactors
        PRIVATE
        ${FormFactors_INCLUDE_DIR}/AsCoeffs.hh

        #
        ${FormFactors_INCLUDE_DIR}/BLPR/FFBLPRBase.hh
        ${FormFactors_INCLUDE_DIR}/BLPR/FFBLPRVarBase.hh
        ${FormFactors_INCLUDE_DIR}/BLPR/FFBtoDBLPR.hh
        ${FormFactors_INCLUDE_DIR}/BLPR/FFBtoDBLPRVar.hh
        ${FormFactors_INCLUDE_DIR}/BLPR/FFBtoDstarBLPR.hh
        ${FormFactors_INCLUDE_DIR}/BLPR/FFBtoDstarBLPRVar.hh

        #
        ${FormFactors_INCLUDE_DIR}/BLPRXP/FFBLPRXPBase.hh
        ${FormFactors_INCLUDE_DIR}/BLPRXP/FFBLPRXPVarBase.hh
        ${FormFactors_INCLUDE_DIR}/BLPRXP/FFBtoDBLPRXP.hh
        ${FormFactors_INCLUDE_DIR}/BLPRXP/FFBtoDBLPRXPVar.hh
        ${FormFactors_INCLUDE_DIR}/BLPRXP/FFBtoDstarBLPRXP.hh
        ${FormFactors_INCLUDE_DIR}/BLPRXP/FFBtoDstarBLPRXPVar.hh

        #
        ${FormFactors_INCLUDE_DIR}/BLRS/FFBLRSBase.hh
        ${FormFactors_INCLUDE_DIR}/BLRS/FFLbtoLcBLRS.hh
        ${FormFactors_INCLUDE_DIR}/BLRS/FFLbtoLcBLRSVar.hh

        #
        ${FormFactors_INCLUDE_DIR}/BLRSXP/FFBLRSXPBase.hh
        ${FormFactors_INCLUDE_DIR}/BLRSXP/FFLbtoLcBLRSXP.hh
        ${FormFactors_INCLUDE_DIR}/BLRSXP/FFLbtoLcBLRSXPVar.hh

        #
        ${FormFactors_INCLUDE_DIR}/BLR/FFBLRBase.hh
        ${FormFactors_INCLUDE_DIR}/BLR/FFBtoD0starBLR.hh
        ${FormFactors_INCLUDE_DIR}/BLR/FFBtoD0starBLRVar.hh
        ${FormFactors_INCLUDE_DIR}/BLR/FFBtoD1starBLR.hh
        ${FormFactors_INCLUDE_DIR}/BLR/FFBtoD1starBLRVar.hh
        ${FormFactors_INCLUDE_DIR}/BLR/FFBtoD1BLR.hh
        ${FormFactors_INCLUDE_DIR}/BLR/FFBtoD1BLRVar.hh
        ${FormFactors_INCLUDE_DIR}/BLR/FFBtoD2starBLR.hh
        ${FormFactors_INCLUDE_DIR}/BLR/FFBtoD2starBLRVar.hh

        #
        ${FormFactors_INCLUDE_DIR}/LLSW/FFLLSWBase.hh
        ${FormFactors_INCLUDE_DIR}/LLSW/FFBtoD0starLLSW.hh
        ${FormFactors_INCLUDE_DIR}/LLSW/FFBtoD1starLLSW.hh
        ${FormFactors_INCLUDE_DIR}/LLSW/FFBtoD1LLSW.hh
        ${FormFactors_INCLUDE_DIR}/LLSW/FFBtoD2starLLSW.hh

        #
        ${FormFactors_INCLUDE_DIR}/ISGW2/FFISGW2Base.hh
        ${FormFactors_INCLUDE_DIR}/ISGW2/FFBtoDISGW2.hh
        ${FormFactors_INCLUDE_DIR}/ISGW2/FFBtoDstarISGW2.hh
        ${FormFactors_INCLUDE_DIR}/ISGW2/FFBtoD0starISGW2.hh
        ${FormFactors_INCLUDE_DIR}/ISGW2/FFBtoD1starISGW2.hh
        ${FormFactors_INCLUDE_DIR}/ISGW2/FFBtoD1ISGW2.hh
        ${FormFactors_INCLUDE_DIR}/ISGW2/FFBtoD2starISGW2.hh
        ${FormFactors_INCLUDE_DIR}/ISGW2/FFBtoRhoOmegaISGW2.hh
        ${FormFactors_INCLUDE_DIR}/ISGW2/FFBtoPiISGW2.hh

        #
        ${FormFactors_INCLUDE_DIR}/CLN/FFCLNBase.hh
        ${FormFactors_INCLUDE_DIR}/CLN/FFBtoDstarCLN.hh
        ${FormFactors_INCLUDE_DIR}/CLN/FFBtoDCLN.hh
        ${FormFactors_INCLUDE_DIR}/CLN/FFBtoDstarCLNVar.hh
        ${FormFactors_INCLUDE_DIR}/CLN/FFBtoDCLNVar.hh

        #
        ${FormFactors_INCLUDE_DIR}/BGL/FFBGLBase.hh
        ${FormFactors_INCLUDE_DIR}/BGL/FFBtoDBGL.hh
        ${FormFactors_INCLUDE_DIR}/BGL/FFBtoDstarBGL.hh
        ${FormFactors_INCLUDE_DIR}/BGL/FFBtoDBGLVar.hh
        ${FormFactors_INCLUDE_DIR}/BGL/FFBtoDstarBGLVar.hh
        ${FormFactors_INCLUDE_DIR}/BGL/FFBctoJpsiBGL.hh
        ${FormFactors_INCLUDE_DIR}/BGL/FFBctoJpsiBGLVar.hh

        #
        ${FormFactors_INCLUDE_DIR}/PCR/FFPCRBase.hh
        ${FormFactors_INCLUDE_DIR}/PCR/FFLbtoLcPCR.hh
        ${FormFactors_INCLUDE_DIR}/PCR/FFLbtoLcstar12PCR.hh
        ${FormFactors_INCLUDE_DIR}/PCR/FFLbtoLcstar32PCR.hh

        #
        ${FormFactors_INCLUDE_DIR}/LSPR/FFLSPRBase.hh
        ${FormFactors_INCLUDE_DIR}/LSPR/FFLSPRVarBase.hh
        ${FormFactors_INCLUDE_DIR}/LSPR/FFLbtoLcstar12LSPR.hh
        ${FormFactors_INCLUDE_DIR}/LSPR/FFLbtoLcstar32LSPR.hh
        ${FormFactors_INCLUDE_DIR}/LSPR/FFLbtoLcstar12LSPRVar.hh
        ${FormFactors_INCLUDE_DIR}/LSPR/FFLbtoLcstar32LSPRVar.hh

        #
        ${FormFactors_INCLUDE_DIR}/EFG/FFEFGBase.hh
        ${FormFactors_INCLUDE_DIR}/EFG/FFBctoJpsiEFG.hh

        #
        ${FormFactors_INCLUDE_DIR}/Kiselev/FFKiselevBase.hh
        ${FormFactors_INCLUDE_DIR}/Kiselev/FFBctoJpsiKiselev.hh

        #
        ${FormFactors_INCLUDE_DIR}/BCL/FFBCLBase.hh
        ${FormFactors_INCLUDE_DIR}/BCL/FFBtoPiBCL.hh
        ${FormFactors_INCLUDE_DIR}/BCL/FFBtoPiBCLVar.hh
        ${FormFactors_INCLUDE_DIR}/BCL/FFBstoKBCL.hh
        ${FormFactors_INCLUDE_DIR}/BCL/FFBstoKBCLVar.hh
        ${FormFactors_INCLUDE_DIR}/BCL/FFBstoDsBCL.hh
        ${FormFactors_INCLUDE_DIR}/BCL/FFBstoDsBCLVar.hh

        #
        ${FormFactors_INCLUDE_DIR}/BSZ/FFBSZBase.hh
        ${FormFactors_INCLUDE_DIR}/BSZ/FFBtoRhoBSZ.hh
        ${FormFactors_INCLUDE_DIR}/BSZ/FFBtoOmegaBSZ.hh
        ${FormFactors_INCLUDE_DIR}/BSZ/FFBtoRhoBSZVar.hh
        ${FormFactors_INCLUDE_DIR}/BSZ/FFBtoOmegaBSZVar.hh

        #
        ${FormFactors_INCLUDE_DIR}/GKvD/FFBtoPiGKvD.hh

        #
        ${FormFactors_INCLUDE_DIR}/PW/FFPWBase.hh
        ${FormFactors_INCLUDE_DIR}/PW/FFD1DstarPiPW.hh
        ${FormFactors_INCLUDE_DIR}/PW/FFD1starDstarPiPW.hh
        ${FormFactors_INCLUDE_DIR}/PW/FFD2starDstarPiPW.hh
        ${FormFactors_INCLUDE_DIR}/PW/FFD2starDPiPW.hh

        #
        ${FormFactors_INCLUDE_DIR}/RCT/FFRCTBase.hh
        ${FormFactors_INCLUDE_DIR}/RCT/FFTauto3PiRCT.hh
)

target_include_directories(FormFactors PRIVATE ${PROJECT_SOURCE_DIR}/include)
set_target_properties(FormFactors PROPERTIES
        POSITION_INDEPENDENT_CODE ${Hammer_PIC}
        VERSION ${Hammer_VERSION}
        SOVERSION ${Hammer_VERSION_MAJOR})
target_compile_features(FormFactors PUBLIC cxx_std_17)

if(VERBOSE_DEBUG)
        target_compile_options(FormFactors PRIVATE ${Hammer_VerboseOptions})
endif()

target_link_libraries(FormFactors PUBLIC HammerTools HammerBase HammerTensor HammerMath)

if(SANITIZE)
        target_link_libraries(FormFactors PRIVATE HammerSanitizers)
endif()

install(TARGETS FormFactors
        EXPORT HammerTargets
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        COMPONENT Hammer_Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/Hammer
        COMPONENT Hammer_Runtime
        NAMELINK_COMPONENT Hammer_Development
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/Hammer
        COMPONENT Hammer_Development
)