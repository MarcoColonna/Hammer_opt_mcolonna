///
/// @file  TestProcRequirements.cc
/// @brief Tests for ProcRequirements
///

//**** This file is a part of the HAMMER library
//**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
//**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
//**** Please note the MCnet academic guidelines; see GUIDELINES for details

#include "Hammer/Process.hh"
#include "Hammer/ProcGraph.hh"
#include "Hammer/ProcRequirements.hh"
#include "Hammer/ProvidersRepo.hh"
#include "Hammer/ExternalData.hh"
#include "Hammer/DictionaryManager.hh"
#include "Hammer/PurePhaseSpaceDefs.hh"
#include "Hammer/SchemeDefinitions.hh"
#include "Hammer/AmplitudeBase.hh"
#include "Hammer/Amplitudes/AmplBToQLepNuBase.hh"
#include "Hammer/IndexLabels.hh"
#include "Hammer/Tools/ParticleUtils.hh"
#include "Hammer/Math/Utils.hh"
#include "Hammer/Math/PhaseSpace.hh"
#include "Hammer/Math/MultiDim/SparseContainer.hh"
#include "Hammer/Tools/SettingsConsumer.hh"
#include "Hammer/Tools/Pdg.hh"
#include "gtest/gtest.h"
#include <string>

using namespace std;


namespace Hammer {

    namespace MD = MultiDimensional;

    class Ampl1ToN : public AmplBToQLepNuBase {

    public:

        Ampl1ToN(const PdgId parent, const vector<PdgId>& daughters, const LabelsList& labels, const size_t mult) {
            // For simplicity, everything has same dimensions: rank 2 tensor
            IndexList dims = {2, 2};
            string fs;
            for (int daughter : daughters) {
                fs.append(to_string(daughter) + ":");
            }
            string name{to_string(parent) + "To" + fs};
            addProcessSignature(parent, daughters);
            addTensor(Tensor{name, MD::makeEmptySparse(dims, labels)});

            setSignatureIndex();
            _multiplicity = mult;
        }

        void eval(const Particle& /*parent*/, const ParticleList& /*daughters*/,
                  const ParticleList& /*references*/) override {
            Tensor& t = getTensor();
            t.clearData();

            t.element({0, 0}) = 0.5 - 0.5i;
            t.element({0, 1}) = 1.5;
            t.element({1, 0}) = -1.5 + 1.i;
            t.element({1, 1}) = -1. + 0.5i;
        }


        void defineSettings() override {
        }

        void updateWilsonCoeffLabelPrefix() override {
        }

    private:
    };

    class DummyExternalData : public ExternalData {
    public:

        DummyExternalData() : ExternalData{nullptr, nullptr} {
        }

        MD::SharedTensorData getWilsonCoefficients(PdgId /*parent*/, const std::vector<PdgId>& /*daughters*/,
                                                   const std::vector<PdgId>& /*granddaughters*/,
                                                   WTerm /*what*/) const override {
            return {};
        }

        MD::SharedTensorData getWilsonCoefficients(const AmplitudeBase* /*amp*/, WTerm /*what*/) const override {
            return {};
        }
    };

    class DummyPurePhaseSpaceDefs : public PurePhaseSpaceDefs {
    public:

        DummyPurePhaseSpaceDefs(const bool PSon2, const bool PSon1) : _mPSon2{PSon2}, _mPSon1{PSon1} {
            // Parent     P
            // Tests  c-c=o=o=o
            //             |
            //             c
        }

        ~DummyPurePhaseSpaceDefs() override = default;

        [[nodiscard]] NumDenPair<bool> isPurePhaseSpace(PdgId parent,
                                                        const std::vector<PdgId>& daughters) const override {
            if (_mPSon2) {
                vector<PdgId> vd = {5, 6};
                if (parent == 2 && (daughters == vd)) {
                    return {true, false};
                }
                vd = {7, 8};
                if (parent == 6 && (daughters == vd)) {
                    return {true, false};
                }
            }
            if (_mPSon1) {
                vector<PdgId> vd2 = {2, 3, 4};
                if (parent == 1 && (daughters == vd2)) {
                    return {true, false};
                }
                vd2 = {11, 12};
                if (parent == 3 && (daughters == vd2)) {
                    return {true, false};
                }
            }
            return {false, false};
        }

    private:

        bool _mPSon2;
        bool _mPSon1;
    };

    class DummySchemeDefinitions : public SchemeDefinitions {
    public:

        DummySchemeDefinitions() = default;
        ~DummySchemeDefinitions() override = default;

        [[nodiscard]] std::set<FFIndex> getFormFactorIndices(HashId /*processId*/) const override {
            return {};
        }

        [[nodiscard]] FFIndex getDenominatorFormFactor(HashId /*processId*/) const override {
            return 0;
        }
    };

    class DummyProvidersRepo : public ProvidersRepo {
    public:

        DummyProvidersRepo()
            : ProvidersRepo{nullptr},
              // Parent     P
              // Tests  c-c=o=o=o
              //             |
              //             c
              // Edge 1 -> (2 -> 5,6),3,4: Edge o=o (Give parent a 'fake' mult of 2)
              _amp1{new Ampl1ToN(1, {6, 5, 4, 3}, {SPIN_DSTAR, SPIN_TAUM}, 2)},
              // Edge 2 -> 5(6 -> 7,8): Edge o=o
              _amp2{new Ampl1ToN(2, {8, 7, 5}, {SPIN_DSSD2STAR, SPIN_GAMMA}, 5)},
              // Vertex 5 -> 9,10: Vertex c
              _amp3{new Ampl1ToN(5, {10, 9}, {SPIN_MUP, SPIN_MUM}, 2)},
              // Edge 1 -> 2(3 -> 11,12)4: Edge o=c
              _amp4{new Ampl1ToN(1, {12, 11, 4, 2}, {SPIN_DSSD2STAR, SPIN_NUTAU}, 2)},
              // Vertex 3 -> 11,12: Vertex c
              _amp5{new Ampl1ToN(3, {12, 11}, {SPIN_TAUM, SPIN_NUTAU}, 2)},
              // Vertex 12 -> 13,14: Vertex c
              _amp6{new Ampl1ToN(12, {14, 13}, {SPIN_NUTAU, SPIN_NUE}, 3)}, _pw{&pw} {
        }

        [[nodiscard]] AmplitudeBase* getAmplitude(PdgId parent, const std::vector<PdgId>& daughters,
                                                  const std::vector<PdgId>& granddaughters = {}) const override {
            vector<PdgId> vd = {2, 3, 4};
            vector<PdgId> vgd = {5, 6};
            if (parent == 1 && (daughters == vd) && (granddaughters == vgd)) {
                return _amp1.get();
            }
            vd = {5, 6};
            vgd = {7, 8};
            if (parent == 2 && (daughters == vd) && (granddaughters == vgd)) {
                return _amp2.get();
            }
            vd = {9, 10};
            vgd = {};
            if (parent == 5 && (daughters == vd) && (granddaughters == vgd)) {
                return _amp3.get();
            }
            vd = {2, 3, 4};
            vgd = {11, 12};
            if (parent == 1 && (daughters == vd) && (granddaughters == vgd)) {
                return _amp4.get();
            }
            vd = {11, 12};
            vgd = {};
            if (parent == 3 && (daughters == vd) && (granddaughters == vgd)) {
                return _amp5.get();
            }
            vd = {13, 14};
            vgd = {};
            if (parent == 12 && (daughters == vd) && (granddaughters == vgd)) {
                return _amp6.get();
            }
            return nullptr;
        }

        [[nodiscard]] RateBase* getRate(PdgId /*parent*/, const std::vector<PdgId>& /*daughters*/,
                                        const std::vector<PdgId>& /*granddaughters*/) const override {
            return nullptr;
        }

        [[nodiscard]] const double* getPartialWidth(PdgId /*parent*/, const std::vector<PdgId>& /*daughters*/,
                                                    const std::vector<PdgId>& /*granddaughters*/ = {}) const override {
            return _pw;
        }

        [[nodiscard]] std::vector<FormFactorBase*> getFormFactor(HashId /*processId*/) const override {
            return {};
        }

    private:

        // 5 amplitudes
        unique_ptr<AmplitudeBase> _amp1;
        unique_ptr<AmplitudeBase> _amp2;
        unique_ptr<AmplitudeBase> _amp3;
        unique_ptr<AmplitudeBase> _amp4;
        unique_ptr<AmplitudeBase> _amp5;
        unique_ptr<AmplitudeBase> _amp6;
        double pw{1.};
        double* _pw;
    };

    class DummyDictionaryManager : public DictionaryManager {

    public:

        DummyDictionaryManager(const bool PSon2, const bool PSon5) : _purePS{PSon2, PSon5} {
        }

        const ProvidersRepo& providers() const override {
            return _prov;
        }

        PurePhaseSpaceDefs& purePSDefs() override {
            return _purePS;
        }

        const PurePhaseSpaceDefs& purePSDefs() const override {
            return _purePS;
        }

        SchemeDefinitions& schemeDefs() override {
            return _schemes;
        }

        const SchemeDefinitions& schemeDefs() const override {
            return _schemes;
        }

        const ExternalData& externalData() const override {
            return _ext;
        }

        ExternalData& externalData() override {
            return _ext;
        }


    private:

        DummyProvidersRepo _prov;
        DummyPurePhaseSpaceDefs _purePS;
        DummySchemeDefinitions _schemes;
        DummyExternalData _ext;
    };

    TEST(ProcRequirementsTest, Construction) {
        SettingsHandler sh;
        sh.addSetting<string>("Hammer", "Units", "GeV");

        Process proc;
        proc.setSettingsHandler(sh);
        // masses set to parent/(num daughters*2)
        Particle p1{{1., 0., 0., 0.}, 1};
        Particle p2{{1. / 6., 0., 0., 0.}, 2};
        Particle p3{{1. / 6., 0., 0., 0.}, 3};
        Particle p4{{1. / 6., 0., 0., 0.}, 4};
        Particle p5{{1. / 6. / 4., 0., 0., 0.}, 5};
        Particle p6{{1. / 6. / 4., 0., 0., 0.}, 6};
        Particle p7{{1. / 6. / 4. / 4., 0., 0., 0.}, 7};
        Particle p8{{1. / 6. / 4. / 4., 0., 0., 0.}, 8};
        Particle p9{{1. / 6. / 4. / 4., 0., 0., 0.}, 9};
        Particle p10{{1. / 6. / 4. / 4., 0., 0., 0.}, 10};
        Particle p11{{1. / 6. / 4., 0., 0., 0.}, 11};
        Particle p12{{1. / 6. / 4., 0., 0., 0.}, 12};
        Particle p13{{1. / 6. / 4. / 4., 0., 0., 0.}, 13};
        Particle p14{{1. / 6. / 4. / 4., 0., 0., 0.}, 14};
        Particle p15{{1. / 6. / 4. / 4. / 4., 0., 0., 0.}, 15};
        Particle p16{{1. / 6. / 4. / 4. / 4., 0., 0., 0.}, 16};

        auto idx1 = proc.addParticle(p1);
        auto idx2 = proc.addParticle(p2);
        auto idx3 = proc.addParticle(p3);
        auto idx4 = proc.addParticle(p4);
        auto idx5 = proc.addParticle(p5);
        auto idx6 = proc.addParticle(p6);
        auto idx7 = proc.addParticle(p7);
        auto idx8 = proc.addParticle(p8);
        auto idx9 = proc.addParticle(p9);
        auto idx10 = proc.addParticle(p10);
        auto idx11 = proc.addParticle(p11);
        auto idx12 = proc.addParticle(p12);
        auto idx13 = proc.addParticle(p13);
        auto idx14 = proc.addParticle(p14);
        auto idx15 = proc.addParticle(p15);
        auto idx16 = proc.addParticle(p16);

        PID& pdg = PID::instance();
        pdg.setMass(p1.pdgId(), p1.p().mass());
        pdg.setMass(p2.pdgId(), p2.p().mass());
        pdg.setMass(p3.pdgId(), p3.p().mass());
        pdg.setMass(p4.pdgId(), p4.p().mass());
        pdg.setMass(p5.pdgId(), p5.p().mass());
        pdg.setMass(p6.pdgId(), p6.p().mass());
        pdg.setMass(p7.pdgId(), p7.p().mass());
        pdg.setMass(p8.pdgId(), p8.p().mass());
        pdg.setMass(p9.pdgId(), p9.p().mass());
        pdg.setMass(p10.pdgId(), p10.p().mass());
        pdg.setMass(p11.pdgId(), p11.p().mass());
        pdg.setMass(p12.pdgId(), p12.p().mass());
        pdg.setMass(p13.pdgId(), p13.p().mass());
        pdg.setMass(p14.pdgId(), p14.p().mass());
        pdg.setMass(p15.pdgId(), p15.p().mass());
        pdg.setMass(p16.pdgId(), p16.p().mass());

        proc.addVertex(idx1, {idx2, idx3, idx4});
        proc.addVertex(idx2, {idx5, idx6});
        proc.addVertex(idx6, {idx7, idx8});
        proc.addVertex(idx5, {idx9, idx10});
        proc.addVertex(idx3, {idx11, idx12});
        proc.addVertex(idx12, {idx13, idx14});
        proc.addVertex(idx13, {idx15, idx16}); // extra guy not implemented in hammer!

        DummyDictionaryManager dummyHammer(false, false);
        dummyHammer.setSettingsHandler(sh);
        proc.initialize();
        ProcGraph pg;
        pg.initialize(&dummyHammer, &proc);
        ProcRequirements pr;
        pr.initialize(&dummyHammer, &proc, &pg);

        // vector of tuples: {parent index, edge daughter index (or 0),  pair<num amplitude ptr, type>, pair<den
        // amplitude prt, type>, num corrmult factor, den corrmult factor, num corrmass factor, den corrmass factor}
        auto procReqAmpl = pr.generatedAmplsMultsPS();
        auto sizeRA = procReqAmpl.size();

        EXPECT_EQ(sizeRA, 4);


        // should contain edges: parent idx1 daughter idx3; parent idx2 daughter idx6; and: vertex idx5; vertex idx12
        // corrresponds to amplitudes with ordering in _requiredAmplitudes: 1 -> 12 11 4 2; 2 -> 8 7 5; 5 -> 10 9; 12 ->
        // 13 14
        vector<tuple<size_t, size_t>> verReqAmpl{tuple<size_t, size_t>{idx1, idx3}, tuple<size_t, size_t>{idx2, idx6},
                                                 tuple<size_t, size_t>{idx5, 0}, tuple<size_t, size_t>{idx12, 0}};
        vector<string> nameReqAmpl = {"1To12:11:4:2:", "2To8:7:5:", "5To10:9:", "12To14:13:"};
        for (size_t ix = 0; ix < sizeRA; ++ix) {
            EXPECT_EQ(get<0>(procReqAmpl[ix]), get<0>(verReqAmpl[ix]));
            EXPECT_EQ(get<1>(procReqAmpl[ix]), get<1>(verReqAmpl[ix]));
            EXPECT_EQ(nameReqAmpl[ix].compare(((get<2>(procReqAmpl[ix])).numerator.ptr)->getTensor().name()), 0);
        }

        // check amplitudes by labels
        vector<LabelsList> labReqAmpl = {
            {SPIN_DSSD2STAR, SPIN_NUTAU}, {SPIN_DSSD2STAR, SPIN_GAMMA}, {SPIN_MUP, SPIN_MUM}, {SPIN_NUTAU, SPIN_NUE}};
        for (size_t ix = 0; ix < sizeRA; ++ix) {
            auto procAmplElemLabs = ((get<2>(procReqAmpl[ix])).numerator.ptr)->getTensor().labels();
            EXPECT_EQ(procAmplElemLabs[0], labReqAmpl[ix][0]);
            EXPECT_EQ(procAmplElemLabs[1], labReqAmpl[ix][1]);
        }

        // Now check PS. Set edge 2 -> 8 7 5 to PS by setting both vertices
        DummyDictionaryManager dummyHammerPS(true, false);
        dummyHammerPS.setSettingsHandler(sh);
        pg.initialize(&dummyHammerPS, &proc);
        pr.initialize(&dummyHammerPS, &proc, &pg);

        // vector of tuples: {parent index, edge daughter index (or 0),  pair<num amplitude ptr, type>, pair<den
        // amplitude prt, type>, num corrmult factor, den corrmult factor, num corrmass factor, den corrmass factor}
        auto procReqAmplPS = pr.generatedAmplsMultsPS();
        auto sizeRAPS = procReqAmplPS.size();

        EXPECT_EQ(sizeRAPS, 4);

        // check nullptr in numerator of  2 -> 8 7 5
        EXPECT_EQ((get<2>(procReqAmplPS[1])).numerator.ptr == nullptr, true);
        EXPECT_EQ((get<2>(procReqAmplPS[1])).denominator.ptr == nullptr, false);

        // check correction factors.
        // Only 5 -> 10 9 (mult 2) should be mult corrected. Initial parent 1/2 factors should be present.
        // Only edge 2 -> 5(6 -> 7,8) should have mass corrections
        vector<vector<double>> reqCorrs = {{1. / 2., 1. / 2.}, {1., 1.}, {1. / 2., 1.}, {1., 1.}};
        double m2m6 = pow(p2.p().mass(), 6. - 2. * 2.) * pow(p6.p().mass(), 6. - 2. * 2.);
        vector<vector<double>> reqMasses = {{1., 1.}, {m2m6, 1.}, {1., 1.}, {1., 1.}};
        for (size_t ix = 0; ix < sizeRAPS; ++ix) {
            EXPECT_NEAR(get<3>(procReqAmplPS[ix]).numerator, reqCorrs[ix][0], 1e-5);
            EXPECT_NEAR(get<3>(procReqAmplPS[ix]).denominator, reqCorrs[ix][1], 1e-5);
            EXPECT_NEAR(get<4>(procReqAmplPS[ix]).numerator, reqMasses[ix][0], 1e-5);
            EXPECT_NEAR(get<4>(procReqAmplPS[ix]).denominator, reqMasses[ix][1], 1e-5);
        }

        // Check PSrates 2 -> 5(6->7,8) edge
        // double p2rateNum = dummyHammerPS.getVertexRate(processID(2,{6,5}), "ISGW2");
        // double p2rateDen = dummyHammerPS.getVertexRate(processID(2,{6,5}), "Denominator");
        // double p2rateNumExp = 0.5*phaseSpaceNBody(p2.p().mass(), {p5.p().mass(),
        // p6.p().mass()})*pow(p2.p().mass(),5.-2.*2.); double p2rateDenExp = 1.; double p6rateNum =
        // dummyHammerPS.getVertexRate(processID(6,{8,7}), "ISGW2"); double p6rateDen =
        // dummyHammerPS.getVertexRate(processID(6,{8,7}), "Denominator"); double p6rateNumExp =
        // 0.5*phaseSpaceNBody(p6.p().mass(), {p8.p().mass(), p7.p().mass()})*pow(p6.p().mass(),5.-2.*2.); double
        // p6rateDenExp = 1.; EXPECT_NEAR(p2rateNum/p2rateNumExp, 1., 1e-4); EXPECT_NEAR(p2rateDen/p2rateDenExp, 1.,
        // 1e-4); EXPECT_NEAR(p6rateNum/p6rateNumExp, 1., 1e-4); EXPECT_NEAR(p6rateDen/p6rateDenExp, 1., 1e-4);

        // Now check PS on two adjacent sites, with non-trivial ancestor on the 12 -> 13 14 vertex
        DummyDictionaryManager dummyHammerPS2(true, true);
        dummyHammerPS2.setSettingsHandler(sh);
        pg.initialize(&dummyHammerPS2, &proc);
        pr.initialize(&dummyHammerPS2, &proc, &pg);

        // vector of tuples: {parent index, edge daughter index (or 0), pair<num amplitude ptr, type>, pair<den
        // amplitude prt, type>, num corr factor, den corr factor}
        auto procReqAmplPS2 = pr.generatedAmplsMultsPS();
        auto sizeRAPS2 = procReqAmplPS2.size();

        EXPECT_EQ(sizeRAPS2, 4);

        // check nullptr in numerator of  2 -> 8 7 5 and 1 -> 12 11 4 2
        EXPECT_EQ((get<2>(procReqAmplPS2[1])).numerator.ptr == nullptr, true);
        EXPECT_EQ((get<2>(procReqAmplPS2[1])).denominator.ptr == nullptr, false);
        EXPECT_EQ((get<2>(procReqAmplPS2[0])).numerator.ptr == nullptr, true);
        EXPECT_EQ((get<2>(procReqAmplPS2[0])).denominator.ptr == nullptr, false);

        // check correction factors.
        // 5 -> 10 9 (mult 2) should be corrected as should 12 -> 13 14 (mult 3).
        // Initial parent 1/2 should be absent in numerator
        // Edge 1 -> 2(3 -> 11,12)4 and edge 2 -> 5(6 -> 7,8) should have mass corrections
        vector<vector<double>> reqCorrs2 = {{1., 1. / 2.}, {1., 1.}, {1. / 2., 1.}, {1. / 3., 1.}};
        double m1m3 = pow(p1.p().mass(), 6. - 2. * 3.) * pow(p3.p().mass(), 6. - 2. * 2.);
        vector<vector<double>> reqMasses2 = {{m1m3, 1.}, {m2m6, 1.}, {1., 1.}, {1., 1.}};
        for (size_t ix = 0; ix < sizeRAPS2; ++ix) {
            EXPECT_NEAR(get<3>(procReqAmplPS2[ix]).numerator, reqCorrs2[ix][0], 1e-5);
            EXPECT_NEAR(get<3>(procReqAmplPS2[ix]).denominator, reqCorrs2[ix][1], 1e-5);
            EXPECT_NEAR(get<4>(procReqAmplPS2[ix]).numerator, reqMasses2[ix][0], 1e-5);
            EXPECT_NEAR(get<4>(procReqAmplPS2[ix]).denominator, reqMasses2[ix][1], 1e-5);
        }

        // Check PSrates 1 -> 2,4(3->11,12) edge
        //  double p1rateNum = dummyHammerPS2.getVertexRate(processID(1,{4,3,2}), "ISGW2");
        //  double p1rateDen = dummyHammerPS2.getVertexRate(processID(1,{4,3,2}), "Denominator");
        //  double p1rateNumExp = 0.5*phaseSpaceNBody(p1.p().mass(), {p4.p().mass(), p3.p().mass(),
        //  p2.p().mass()})*pow(p1.p().mass(),5.-2.*3.); double p1rateDenExp = 1.; double p3rateNum =
        //  dummyHammerPS2.getVertexRate(processID(3,{12,11}), "ISGW2"); double p3rateDen =
        //  dummyHammerPS2.getVertexRate(processID(3,{12,11}), "Denominator"); double p3rateNumExp =
        //  0.5*phaseSpaceNBody(p3.p().mass(), {p11.p().mass(), p12.p().mass()})*pow(p3.p().mass(),5.-2.*2.); double
        //  p3rateDenExp = 1.; EXPECT_NEAR(p1rateNum/p1rateNumExp, 1., 1e-4); EXPECT_NEAR(p1rateDen/p1rateDenExp, 1.,
        //  1e-4); EXPECT_NEAR(p3rateNum/p3rateNumExp, 1., 1e-4); EXPECT_NEAR(p3rateDen/p3rateDenExp, 1., 1e-4);
    }

    //////////////////////////////////////////////
    ////////  Second Test for Phase Space  ///////
    //////////////////////////////////////////////

    class DummyPurePhaseSpaceDefs2 : public PurePhaseSpaceDefs {
    public:

        DummyPurePhaseSpaceDefs2(const bool PSon2, const bool PSon5) : _mPSon2{PSon2}, _mPSon5{PSon5} {
            // Parent   P
            // Tests    c
            //         |
            //       c-o=c=c
            //         ^ ^
            //     PSon2 PSon5
        }

        ~DummyPurePhaseSpaceDefs2() override = default;

        [[nodiscard]] NumDenPair<bool> isPurePhaseSpace(PdgId parent,
                                                        const std::vector<PdgId>& daughters) const override {
            if (_mPSon2) {
                vector<PdgId> vd = {4, 5};
                if (parent == 2 && (daughters == vd)) {
                    return {true, false};
                }
            }
            if (_mPSon5) {
                vector<PdgId> vd2 = {8, 9};
                if (parent == 5 && (daughters == vd2)) {
                    return {true, false};
                }
            }
            return {false, false};
        }

    private:

        bool _mPSon2;
        bool _mPSon5;
    };


    class DummyProvidersRepo2 : public ProvidersRepo {
    public:

        DummyProvidersRepo2()
            : ProvidersRepo{nullptr},
              // Parent   P
              // Tests    c
              //         |
              //       c-o=c=c
              //         ^ ^
              //     PSon2 PSon5
              // Vertex 1 -> 2,3: parent c
              _amp1{new Ampl1ToN(1, {3, 2}, {SPIN_TAUM, SPIN_TAUP}, 1)},
              // Edge 2 -> 4(5 -> 8,9): Edge o=c
              _amp2{new Ampl1ToN(2, {9, 8, 4}, {SPIN_TAUP, SPIN_DSTAR}, 2)},
              // Vertex 4 -> 6,7: Vertex c
              _amp3{new Ampl1ToN(4, {7, 6}, {SPIN_MUP, SPIN_MUM}, 2)},
              // Edge 5 -> (8 -> 10,11)9: Edge c=c
              _amp4{new Ampl1ToN(5, {11, 10, 9, 8}, {SPIN_DSSD2STAR, SPIN_GAMMA}, 5)},
              // Vertex 5 -> 8,9: Vertex c
              _amp5{new Ampl1ToN(5, {9, 8}, {SPIN_DSSD2STAR, SPIN_DSTAR}, 5)},
              // Vertex 8 -> 10,11: Vertex c
              _amp6{new Ampl1ToN(8, {11, 10}, {SPIN_DSTAR, SPIN_GAMMA}, 3)}, _pw{&pw} {
        }

        [[nodiscard]] AmplitudeBase* getAmplitude(PdgId parent, const std::vector<PdgId>& daughters,
                                                  const std::vector<PdgId>& granddaughters = {}) const override {
            vector<PdgId> vd = {2, 3};
            vector<PdgId> vgd = {};
            if (parent == 1 && (daughters == vd) && (granddaughters == vgd)) {
                return _amp1.get();
            }
            vd = {4, 5};
            vgd = {8, 9};
            if (parent == 2 && (daughters == vd) && (granddaughters == vgd)) {
                return _amp2.get();
            }
            vd = {6, 7};
            vgd = {};
            if (parent == 4 && (daughters == vd) && (granddaughters == vgd)) {
                return _amp3.get();
            }
            vd = {8, 9};
            vgd = {10, 11};
            if (parent == 5 && (daughters == vd) && (granddaughters == vgd)) {
                return _amp4.get();
            }
            vd = {8, 9};
            vgd = {};
            if (parent == 5 && (daughters == vd) && (granddaughters == vgd)) {
                return _amp5.get();
            }
            vd = {10, 11};
            vgd = {};
            if (parent == 8 && (daughters == vd) && (granddaughters == vgd)) {
                return _amp6.get();
            }
            return nullptr;
        }

        [[nodiscard]] RateBase* getRate(PdgId /*parent*/, const std::vector<PdgId>& /*daughters*/,
                                        const std::vector<PdgId>& /*granddaughters*/) const override {
            return nullptr;
        }

        [[nodiscard]] const double* getPartialWidth(PdgId /*parent*/, const std::vector<PdgId>& /*daughters*/,
                                                    const std::vector<PdgId>& /*granddaughters*/ = {}) const override {
            return _pw;
        }

        [[nodiscard]] std::vector<FormFactorBase*> getFormFactor(HashId /*processId*/) const override {
            return {};
        }


    private:

        // 5 amplitudes
        unique_ptr<AmplitudeBase> _amp1;
        unique_ptr<AmplitudeBase> _amp2;
        unique_ptr<AmplitudeBase> _amp3;
        unique_ptr<AmplitudeBase> _amp4;
        unique_ptr<AmplitudeBase> _amp5;
        unique_ptr<AmplitudeBase> _amp6;
        double pw{1.};
        double* _pw;
    };

    class DummyDictionaryManager2 : public DictionaryManager {

    public:

        DummyDictionaryManager2(const bool PSon2, const bool PSon5) : _purePS{PSon2, PSon5} {
        }

        const ProvidersRepo& providers() const override {
            return _prov;
        }

        PurePhaseSpaceDefs& purePSDefs() override {
            return _purePS;
        }

        const PurePhaseSpaceDefs& purePSDefs() const override {
            return _purePS;
        }

        SchemeDefinitions& schemeDefs() override {
            return _schemes;
        }

        const SchemeDefinitions& schemeDefs() const override {
            return _schemes;
        }

        const ExternalData& externalData() const override {
            return _ext;
        }

        ExternalData& externalData() override {
            return _ext;
        }

    private:

        DummyProvidersRepo2 _prov;
        DummyPurePhaseSpaceDefs2 _purePS;
        DummySchemeDefinitions _schemes;
        DummyExternalData _ext;
    };


    TEST(ProcRequirementsTest, PhaseSpace) {
        SettingsHandler sh;
        sh.addSetting<string>("Hammer", "Units", "GeV");

        Process proc;
        proc.setSettingsHandler(sh);

        // masses set to parent/(num daughters*2)
        Particle p1{{1., 0., 0., 0.}, 1};
        Particle p2{{1. / 4., 0., 0., 0.}, 2};
        Particle p3{{1. / 4., 0., 0., 0.}, 3};
        Particle p4{{1. / 4. / 4., 0., 0., 0.}, 4};
        Particle p5{{1. / 4. / 4., 0., 0., 0.}, 5};
        Particle p6{{1. / 4. / 4. / 4., 0., 0., 0.}, 6};
        Particle p7{{1. / 4. / 4. / 4., 0., 0., 0.}, 7};
        Particle p8{{1. / 4. / 4. / 4., 0., 0., 0.}, 8};
        Particle p9{{1. / 4. / 4. / 4., 0., 0., 0.}, 9};
        Particle p10{{1. / 4. / 4. / 4. / 4., 0., 0., 0.}, 10};
        Particle p11{{1. / 4. / 4. / 4. / 4., 0., 0., 0.}, 11};

        PID& pdg = PID::instance();
        pdg.setMass(p1.pdgId(), p1.p().mass());
        pdg.setMass(p2.pdgId(), p2.p().mass());
        pdg.setMass(p3.pdgId(), p3.p().mass());
        pdg.setMass(p4.pdgId(), p4.p().mass());
        pdg.setMass(p5.pdgId(), p5.p().mass());
        pdg.setMass(p6.pdgId(), p6.p().mass());
        pdg.setMass(p7.pdgId(), p7.p().mass());
        pdg.setMass(p8.pdgId(), p8.p().mass());
        pdg.setMass(p9.pdgId(), p9.p().mass());
        pdg.setMass(p10.pdgId(), p10.p().mass());
        pdg.setMass(p11.pdgId(), p11.p().mass());

        auto idx1 = proc.addParticle(p1);
        auto idx2 = proc.addParticle(p2);
        auto idx3 = proc.addParticle(p3);
        auto idx4 = proc.addParticle(p4);
        auto idx5 = proc.addParticle(p5);
        auto idx6 = proc.addParticle(p6);
        auto idx7 = proc.addParticle(p7);
        auto idx8 = proc.addParticle(p8);
        auto idx9 = proc.addParticle(p9);
        auto idx10 = proc.addParticle(p10);
        auto idx11 = proc.addParticle(p11);

        proc.addVertex(idx1, {idx2, idx3});
        proc.addVertex(idx2, {idx4, idx5});
        proc.addVertex(idx4, {idx6, idx7});
        proc.addVertex(idx5, {idx8, idx9});
        proc.addVertex(idx8, {idx10, idx11});

        // set just the o to PS
        proc.initialize();
        DummyDictionaryManager2 dummyHammer(true, false);
        dummyHammer.setSettingsHandler(sh);
        ProcGraph pg;
        pg.initialize(&dummyHammer, &proc);
        ProcRequirements pr;
        pr.initialize(&dummyHammer, &proc, &pg);

        // vector of tuples: {parent index, edge daughter index (or 0), pair<num amplitude ptr, type>, pair<den
        // amplitude prt, type>, num corr factor, den corr factor}
        auto procReqAmpl = pr.generatedAmplsMultsPS();
        auto sizeRA = procReqAmpl.size();
        EXPECT_EQ(sizeRA, 4);

        // Should have types V, D, V, V in num and V, E, V, V in den
        // corrresponds to amplitudes with ordering in _requiredAmplitudes: 1 -> 23; 2 -> 4 8 9; 4 -> 6 7; 8 -> 10 11
        vector<AmplType> typeRAnum{AmplType::VERTEX, AmplType::DAUGHTEREDGE, AmplType::VERTEX, AmplType::VERTEX};
        vector<AmplType> typeRAden{AmplType::VERTEX, AmplType::FULLEDGE, AmplType::VERTEX, AmplType::VERTEX};
        for (size_t idx = 0; idx < sizeRA; ++idx) {
            EXPECT_EQ((get<2>(procReqAmpl[idx])).numerator.type, typeRAnum[idx]);
        }
        for (size_t idx = 0; idx < sizeRA; ++idx) {
            EXPECT_EQ((get<2>(procReqAmpl[idx])).denominator.type, typeRAden[idx]);
        }

        // Correction factors in num should be on 4 -> 6,7 (mult 2) and 5 -> 8,9 (mult 5) inside daughter edge 2 -> 4 8
        // 9 Mass corrections on 2 -> 4,5 vertex inside daughter edge 2 -> 4 8 9
        vector<vector<double>> reqCorrs = {{1., 1.}, {1. / 5., 1.}, {1. / 2., 1.}, {1., 1.}};
        double m2 = pow(p2.p().mass(), 6. - 2. * 2.);
        vector<vector<double>> reqMasses = {{1., 1.}, {m2, 1.}, {1., 1.}, {1., 1.}};
        for (size_t ix = 0; ix < sizeRA; ++ix) {
            EXPECT_NEAR(get<3>(procReqAmpl[ix]).numerator, reqCorrs[ix][0], 1e-5);
            EXPECT_NEAR(get<3>(procReqAmpl[ix]).denominator, reqCorrs[ix][1], 1e-5);
            EXPECT_NEAR(get<4>(procReqAmpl[ix]).numerator, reqMasses[ix][0], 1e-5);
            EXPECT_NEAR(get<4>(procReqAmpl[ix]).denominator, reqMasses[ix][1], 1e-5);
        }

        // Check PSrates on 2 -> 4,5
        //  double p2rateNum = dummyHammer.getVertexRate(processID(2,{5,4}), "ISGW2");
        //  double p2rateDen = dummyHammer.getVertexRate(processID(2,{5,4}), "Denominator");
        //  double p2rateNumExp = 0.5*phaseSpaceNBody(p2.p().mass(), {p4.p().mass(),
        //  p5.p().mass()})*pow(p2.p().mass(),5.-2.*2.); double p2rateDenExp = 1.;
        //  EXPECT_NEAR(p2rateNum/p2rateNumExp, 1., 1e-4);
        //  EXPECT_NEAR(p2rateDen/p2rateDenExp, 1., 1e-4);

        // Now test second with PS vertex on too, setting o=c to PS.
        DummyDictionaryManager2 dummyHammer2(true, true);
        dummyHammer2.setSettingsHandler(sh);
        pg.initialize(&dummyHammer2, &proc);
        pr.initialize(&dummyHammer2, &proc, &pg);

        auto procReqAmpl2 = pr.generatedAmplsMultsPS();
        auto sizeRA2 = procReqAmpl2.size();
        EXPECT_EQ(sizeRA, 4);

        // Should have types V, E==null, V, V in num and V, E, V, V in den
        // corrresponds to amplitudes with ordering in _requiredAmplitudes: 1 -> 23; 2 -> 4 8 9; 4 -> 6 7; 8 -> 10 11
        vector<AmplType> typeRAnum2{AmplType::VERTEX, AmplType::FULLEDGE, AmplType::VERTEX, AmplType::VERTEX};
        vector<AmplType> typeRAden2{AmplType::VERTEX, AmplType::FULLEDGE, AmplType::VERTEX, AmplType::VERTEX};
        for (size_t idx = 0; idx < sizeRA2; ++idx) {
            EXPECT_EQ((get<2>(procReqAmpl2[idx])).numerator.type, typeRAnum2[idx]);
        }
        for (size_t idx = 0; idx < sizeRA; ++idx) {
            EXPECT_EQ((get<2>(procReqAmpl2[idx])).denominator.type, typeRAden2[idx]);
        }
        EXPECT_EQ((get<2>(procReqAmpl2[1])).numerator.ptr == nullptr, true);

        // Correction factors in num should be on 4 -> 6,7 (mult 2) and 8 -> 11,10 (mult 3)
        // Mass corrections on 2 -> 4(5-> 8,9) edge
        vector<vector<double>> reqCorrs2 = {{1., 1.}, {1., 1.}, {1. / 2., 1.}, {1. / 3., 1.}};
        double m2m5 = pow(p2.p().mass(), 6. - 2. * 2.) * pow(p5.p().mass(), 6. - 2. * 2.);
        vector<vector<double>> reqMasses2 = {{1., 1.}, {m2m5, 1.}, {1., 1.}, {1., 1.}};
        for (size_t ix = 0; ix < sizeRA; ++ix) {
            EXPECT_NEAR(get<3>(procReqAmpl2[ix]).numerator, reqCorrs2[ix][0], 1e-5);
            EXPECT_NEAR(get<3>(procReqAmpl2[ix]).denominator, reqCorrs2[ix][1], 1e-5);
            EXPECT_NEAR(get<4>(procReqAmpl2[ix]).numerator, reqMasses2[ix][0], 1e-5);
            EXPECT_NEAR(get<4>(procReqAmpl2[ix]).denominator, reqMasses2[ix][1], 1e-5);
        }

        // Check PSrates on 2 -> 4(5-> 8,9)
        //  p2rateNum = dummyHammer2.getVertexRate(processID(2,{5,4}), "ISGW2");
        //  p2rateDen = dummyHammer2.getVertexRate(processID(2,{5,4}), "Denominator");
        //  double p5rateNum = dummyHammer2.getVertexRate(processID(5,{9,8}), "ISGW2");
        //  double p5rateDen = dummyHammer2.getVertexRate(processID(5,{9,8}), "Denominator");
        //  double p5rateNumExp = 0.5*phaseSpaceNBody(p5.p().mass(), {p9.p().mass(),
        //  p8.p().mass()})*pow(p5.p().mass(),5.-2.*2.); double p5rateDenExp = 1.;
        //  EXPECT_NEAR(p2rateNum/p2rateNumExp, 1., 1e-4);
        //  EXPECT_NEAR(p2rateDen/p2rateDenExp, 1., 1e-4);
        //  EXPECT_NEAR(p5rateNum/p5rateNumExp, 1., 1e-4);
        //  EXPECT_NEAR(p5rateDen/p5rateDenExp, 1., 1e-4);
    }

} // namespace Hammer
