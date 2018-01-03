#ifndef NCODE_GEO_H
#define NCODE_GEO_H

#include <stddef.h>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <mutex>

#include "ncode/substitute.h"
#include "ncode/logging.h"

namespace nc {
namespace geo {

struct CityData {
  CityData(const std::string& name, const std::string& country_code,
           size_t population, double longitude, double latitude)
      : name(name),
        country_code(country_code),
        population(population),
        longitude(longitude),
        latitude(latitude) {}

  std::string ToString() const {
    return StrCat("name: ", name, " country: ", country_code, " pop: ",
                  population);
  }

  // Returns the distance between this city and another one.
  double DistanceKm(const CityData& other_city) const;

  // The ASCII name of the city.
  const std::string name;

  // Country code.
  const std::string country_code;

  // Population.
  const uint64_t population;

  // Coordinates of the city in decimal degrees.
  const double longitude;
  const double latitude;
};

std::ostream& operator<<(std::ostream& out, const CityData& city);

// A generic interface for projections.
class ProjectionInterface {
 public:
  virtual ~ProjectionInterface() {}

  // Returns the x,y coordinates of a lon,lat point. Both should be in decimal
  // degrees.
  virtual std::pair<double, double> Project(double longitude,
                                            double latitude) const = 0;
};

class WebMercator : public ProjectionInterface {
 public:
  WebMercator(double zoom_level = 1.0) : zoom_level_(zoom_level) {}

  std::pair<double, double> Project(double longitude,
                                    double latitude) const override {
    double lon_rad = (M_PI / 180.0) * longitude;
    double lat_rad = (M_PI / 180.0) * latitude;

    double c = (128.0 / M_PI) * std::pow(2, zoom_level_);
    double x = c * (lon_rad + M_PI);
    double y = c * (M_PI - std::log(std::tan(M_PI / 4.0 + lat_rad / 2.0)));

    DCHECK(x == x);
    DCHECK(y == y);
    return std::make_pair(x, y);
  }

 private:
  double zoom_level_;
};

// Given a city populates its x,y coordinates according to some projection.
bool CityXY(const CityData& city, const ProjectionInterface& projection,
            double* x, double* y);

struct FindCityRequest {
  // The name of the city. Out of all non-filtered cities the one that whose
  // name has the closest distance (Livenshtein) to this string will be
  // returned.
  std::string ascii_name;

  // A regex filter to be applied to names before shortest-distance matching.
  std::string regex_ascii_name_filter = ".*";

  // Only cities with country codes that pass this filter will be considered.
  std::string regex_country_code_filter = ".*";

  bool operator<(const FindCityRequest& other) const {
    return std::tie(ascii_name, regex_ascii_name_filter,
                    regex_country_code_filter) <
           std::tie(other.ascii_name, other.regex_ascii_name_filter,
                    other.regex_country_code_filter);
  }
};

class Localizer {
 public:
  static constexpr char kDefaultCitiesFileLocation[] =
      "geonames/cities5000.txt";

  Localizer(
      const std::string& world_cities_db_file = kDefaultCitiesFileLocation);

  // Finds a city by name.
  const CityData* FindCityOrNull(const FindCityRequest& request);

  // Returns all cities in no particular order.
  const std::vector<CityData>& AllCities() const { return cities_; }

 private:
  static constexpr size_t kCacheSize = 10000;

  // Only cities with non-0 population are retained.
  std::vector<CityData> cities_;

  // Cache of requests.
  std::map<FindCityRequest, const CityData*> request_to_city_;

  // Protects the cache.
  std::mutex mu_;
};

}  // namespace geo
}  // namespace nc
#endif
