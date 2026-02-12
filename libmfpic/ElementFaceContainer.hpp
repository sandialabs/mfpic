#pragma once

#include <unordered_map>

namespace mfpic {

/// Stores an object on element-local faces.
template <typename T>
class ElementFaceContainer {
public:
  /**
   * @brief Insert a thing on a face on an element.
   *
   * @param[in] element         Element to which desired face belongs.
   * @param[in] element_face    Element-local face.
   * @param[in] thing_to_insert Thing to insert.
   *
   * @returns Same thing as std::unordered_map::insert.
   */
  auto insert(int element, int element_face, const T& thing_to_insert);

  /**
   * @brief Request a thing on a face on an element, throwing an exception if no such thing exists.
   *
   * @param[in] element         Element to which desired face belongs.
   * @param[in] element_face    Element-local face.
   *
   * @returns Same thing as std::unordered_map::at.
   */
  const auto& at(int element, int element_face) const;

  /**
   * @brief Request a thing on a face on an element, throwing an exception if no such thing exists.
   *
   * @param[in] element         Element to which desired face belongs.
   * @param[in] element_face    Element-local face.
   *
   * @returns Same thing as std::unordered_map::at.
   */
  auto& at(int element, int element_face);

  /**
   * @brief Checks whether something exists on the given element and face.
   *
   * @param[in] element         Element to which desired face belongs.
   * @param[in] element_face    Element-local face.
   *
   * @returns Whether something exists on the given element and face.
   */
  bool contains(int element, int element_face) const;

  /// Hashes a pair of ints.
  struct IntPairHasher {
    /**
     * @brief Hash a pair of ints.
     *
     * @param[in] pair Pair of ints.
     *
     * @returns Hash of pair of ints.
     */
    std::size_t operator() (const std::pair<int, int>& pair) const;
  };

  /// Container type to which this class is an interface.
  using UnderlyingContainerType = std::unordered_map<std::pair<int, int>, T, IntPairHasher>;

private:
  /// Maps element-face pairs to objects.
  UnderlyingContainerType data_;
};

// template definitions
template <typename T>
auto ElementFaceContainer<T>::insert(int element, int element_face, const T& thing_to_insert) {
  return data_.insert(std::make_pair(std::make_pair(element, element_face), thing_to_insert));
}

template <typename T>
const auto& ElementFaceContainer<T>::at(int element, int element_face) const {
  return data_.at(std::make_pair(element, element_face));
}

template <typename T>
auto& ElementFaceContainer<T>::at(int element, int element_face) {
  return data_.at(std::make_pair(element, element_face));
}

template <typename T>
bool ElementFaceContainer<T>::contains(int element, int element_face) const {
  return data_.contains(std::make_pair(element, element_face));
}

template <typename T>
std::size_t ElementFaceContainer<T>::IntPairHasher::operator()(const std::pair<int, int> &pair) const {
  std::hash<int> int_hasher;
  return int_hasher(pair.first) ^ int_hasher(pair.second);
}

} // namespace mfpic
