set(Hammer_Math_TEST_LIBS HammerMath yaml-cpp::yaml-cpp)
set(Hammer_INT_TEST_LIBS HammerCore Rates)

new_test(FourMomentum "${Hammer_Math_TEST_LIBS}")
new_test(Tensor "${Hammer_Math_TEST_LIBS}")
new_test(Histogram "${Hammer_Math_TEST_LIBS}")
new_test(Integrator "${Hammer_INT_TEST_LIBS}")
new_test(PhaseSpace "${Hammer_INT_TEST_LIBS}")
new_test(Differentiator "${Hammer_INT_TEST_LIBS}")
new_test(MathUtils "${Hammer_Math_TEST_LIBS}")
new_test(Units "${Hammer_Math_TEST_LIBS}")

add_subdirectory(MultiDim)
