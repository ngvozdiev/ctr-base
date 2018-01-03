#include "geo.h"

#include <gtest/gtest.h>
#include <limits>
#include <tuple>

#include "ncode/common.h"

namespace nc {
namespace geo {
namespace {

class LocalizerTest : public ::testing::Test {
 protected:
  Localizer localizer_;
};

TEST_F(LocalizerTest, Init) { ASSERT_LT(0ul, localizer_.AllCities().size()); }

TEST_F(LocalizerTest, EmptyLookup) {
  FindCityRequest request;

  // Since we have not specified the city name and by default the filters are .*
  // this should return the city with the shortest name that has the most
  // population.
  const CityData* city = localizer_.FindCityOrNull(request);
  CHECK_NE(nullptr, city);

  const std::vector<CityData>& all_cities = localizer_.AllCities();
  const CityData* shortest_named_city = nullptr;
  size_t name_len = std::numeric_limits<size_t>::max();
  size_t population = 0;
  for (const CityData& city_data : all_cities) {
    if (city_data.name.size() < name_len) {
      shortest_named_city = &city_data;
      name_len = city_data.name.size();
      population = city_data.population;
    } else if (city_data.name.size() == name_len) {
      if (city_data.population > population) {
        shortest_named_city = &city_data;
        name_len = city_data.name.size();
        population = city_data.population;
      }
    }
  }

  ASSERT_EQ(city, shortest_named_city);
}

TEST_F(LocalizerTest, SimpleLookup) {
  FindCityRequest request;
  request.ascii_name = "London";

  const CityData* city = localizer_.FindCityOrNull(request);
  CHECK_NE(nullptr, city);
  ASSERT_EQ("London", city->name);
  ASSERT_EQ("GB", city->country_code);

  // There is also a London in the US, but with a smaller population.
  request.regex_country_code_filter = "US";
  city = localizer_.FindCityOrNull(request);
  CHECK_NE(nullptr, city);
  ASSERT_EQ("London", city->name);
  ASSERT_EQ("US", city->country_code);
}

TEST_F(LocalizerTest, CaseInsensitive) {
  FindCityRequest request;

  request.ascii_name = "NewYork";
  const CityData* c1 = localizer_.FindCityOrNull(request);

  request.ascii_name = "New York";
  const CityData* c2 = localizer_.FindCityOrNull(request);

  request.ascii_name = "New york";
  const CityData* c3 = localizer_.FindCityOrNull(request);

  request.ascii_name = "newYork";
  const CityData* c4 = localizer_.FindCityOrNull(request);

  request.ascii_name = "newyork";
  const CityData* c5 = localizer_.FindCityOrNull(request);

  ASSERT_TRUE(AllEqual(c1, c2, c3, c4, c5));
}

TEST_F(LocalizerTest, XYProjection) {
  WebMercator projection;

  FindCityRequest request;
  request.ascii_name = "London";
  const CityData* city = localizer_.FindCityOrNull(request);

  double x, y;
  std::tie(x, y) = projection.Project(city->longitude, city->latitude);

  ASSERT_LT(0, x);
  ASSERT_LT(0, y);
}
}
}
}
