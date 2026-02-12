#pragma once

namespace mfpic {

#define EXPECT_NEAR_RELATIVE( v1, v2, relative_tolerance ) \
  EXPECT_LT(std::abs(v1 - v2), relative_tolerance * v1)

}