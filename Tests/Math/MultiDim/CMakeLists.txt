set(Hammer_Tensor_TEST_LIBS HammerTensor)

new_test(AlignedIndexing "${Hammer_Tensor_TEST_LIBS}")
new_test(SequentialIndexing "${Hammer_Tensor_TEST_LIBS}")
new_test(LabeledIndexing "${Hammer_Tensor_TEST_LIBS}")
new_test(BlockIndexing "${Hammer_Tensor_TEST_LIBS}")
new_test(BruteForceIterator "${Hammer_Tensor_TEST_LIBS}")

new_test(VectorContainer "${Hammer_Tensor_TEST_LIBS}")
new_test(SparseContainer "${Hammer_Tensor_TEST_LIBS}")
new_test(ScalarContainer "${Hammer_Tensor_TEST_LIBS}")
new_test(OuterContainer "${Hammer_Tensor_TEST_LIBS}")

add_subdirectory(Ops)
