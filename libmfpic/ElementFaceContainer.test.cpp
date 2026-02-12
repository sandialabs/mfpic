#include <libmfpic/ElementFaceContainer.hpp>

#include <gtest/gtest.h>
#include <stdexcept>

namespace {

using namespace mfpic;

using TestType = int;
using TestContainer = ElementFaceContainer<TestType>;
using InsertResult = std::pair<TestContainer::UnderlyingContainerType::iterator, bool>;

TEST(ElementFaceContainer, CanRetrieveOneInsertedObject) {
  TestContainer container;

  constexpr TestType thing_to_insert = 3;
  InsertResult insert_result = container.insert(0, 1, thing_to_insert);

  EXPECT_TRUE(insert_result.second);
  EXPECT_EQ(container.at(0, 1), thing_to_insert);
}

TEST(ElementFaceContainer, DoubleInsertionFails) {
  TestContainer container;

  constexpr TestType thing_to_insert = 3;
  InsertResult first_insert_result = container.insert(0, 1, thing_to_insert);
  InsertResult second_insert_result = container.insert(0, 1, thing_to_insert);

  EXPECT_TRUE(first_insert_result.second);
  EXPECT_FALSE(second_insert_result.second);
}

TEST(ElementFaceContainer, AtThrowsIfNoSuchElementExists) {
  TestContainer container;

  EXPECT_THROW(container.at(0, 0), std::out_of_range);
}

} // namespace
