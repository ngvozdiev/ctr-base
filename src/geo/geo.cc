#include "geo.h"

#include <chrono>
#include <cstdint>
#include <functional>
#include <limits>
#include <regex>
#include <iostream>
#include <fstream>

#include "ncode_common/src/map_util.h"
#include "ncode_common/src/strutil.h"

namespace nc {
namespace geo {

static constexpr size_t kEntriesPerRow = 19;
static constexpr size_t kNameIndex = 2;
static constexpr size_t kPopulationIndex = 14;
static constexpr size_t kLatIndex = 4;
static constexpr size_t kLonIndex = 5;
static constexpr size_t kCCIndex = 8;

constexpr char Localizer::kDefaultCitiesFileLocation[];
static constexpr double kEarthRadiusKm = 6371.0;

static void ReadLinesWithCallback(
    const std::string& name,
    std::function<void(const std::string& line)> callback) {
  std::ifstream infile(name);
  std::string line;
  while (std::getline(infile, line)) {
    callback(line);
  }
}

static bool ParseRow(const std::string& row, std::vector<CityData>* out) {
  std::vector<std::string> pieces;
  SplitStringAllowEmpty(row, "\t", &pieces);
  if (pieces.size() != kEntriesPerRow) {
    return false;
  }

  uint32_t population;
  double latitude;
  double longitude;
  if (!safe_strtou32(pieces[kPopulationIndex], &population) ||
      !safe_strtod(pieces[kLatIndex], &latitude) ||
      !safe_strtod(pieces[kLonIndex], &longitude)) {
    return false;
  }

  const std::string& ascii_name = pieces[kNameIndex];
  const std::string& country_code = pieces[kCCIndex];
  if (ascii_name.empty() || country_code.empty()) {
    return false;
  }

  out->emplace_back(ascii_name, country_code, population, longitude, latitude);
  return true;
}

Localizer::Localizer(const std::string& world_cities_db_file) {
  using namespace std::chrono;
  size_t num_errors = 0;
  size_t num_rows = 0;
  auto start_time_point = high_resolution_clock::now();
  ReadLinesWithCallback(world_cities_db_file, [this, &num_errors, &num_rows](
                                                  const std::string& line) {
    if (!ParseRow(line, &cities_)) {
      ++num_errors;
    }
    ++num_rows;
  });
  auto end_time_point = high_resolution_clock::now();
  auto delta = duration_cast<milliseconds>(end_time_point - start_time_point);
  LOG(INFO) << Substitute("Processed $0 rows in $1ms ($2 errors)", num_rows,
                          delta.count(), num_errors);
}

const CityData* Localizer::FindCityOrNull(const FindCityRequest& request) {
  const CityData* in_cache = FindPtrOrNull(request_to_city_, request);
  if (in_cache != nullptr) {
    return in_cache;
  }

  std::regex cc_regex(request.regex_country_code_filter,
                      std::regex_constants::icase);
  std::regex name_regex(request.regex_ascii_name_filter,
                        std::regex_constants::icase);

  std::vector<const CityData*> cities_to_consider;
  for (const auto& city : cities_) {
    if (std::regex_search(city.country_code, cc_regex) &&
        std::regex_search(city.name, name_regex)) {
      cities_to_consider.emplace_back(&city);
    }
  }

  size_t best_population = 0;
  int best_score = std::numeric_limits<int>::max();
  const CityData* best_city = nullptr;
  for (const CityData* city : cities_to_consider) {
    int score = StrDistanceCaseInsensitive(city->name, request.ascii_name);
    if (score < best_score) {
      best_city = city;
      best_score = score;
      best_population = city->population;
    } else if (score == best_score) {
      if (city->population > best_population) {
        best_city = city;
        best_score = score;
        best_population = city->population;
      }
    }
  }

  if (request_to_city_.size() == kCacheSize) {
    request_to_city_.erase(request_to_city_.begin());
  }
  request_to_city_.emplace(request, best_city);

  return best_city;
}

std::ostream& operator<<(std::ostream& out, const CityData& city) {
  out << Substitute("name: $0, cc: $1, population: $2, long: $3, lat: $4",
                    city.name, city.country_code, city.population,
                    city.longitude, city.latitude);
  return out;
}

// This function converts decimal degrees to radians
static double Deg2rad(double deg) { return (deg * M_PI / 180); }

// Stolen from
// https://stackoverflow.com/questions/10198985/calculating-the-distance-between-2-latitudes-and-longitudes-that-are-saved-in-a
static double DistanceEarth(double lat1d, double lon1d, double lat2d,
                            double lon2d) {
  double lat1r, lon1r, lat2r, lon2r, u, v;
  lat1r = Deg2rad(lat1d);
  lon1r = Deg2rad(lon1d);
  lat2r = Deg2rad(lat2d);
  lon2r = Deg2rad(lon2d);
  u = sin((lat2r - lat1r) / 2);
  v = sin((lon2r - lon1r) / 2);
  return 2.0 * kEarthRadiusKm *
         asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));
}

double CityData::DistanceKm(const CityData& other_city) const {
  return DistanceEarth(latitude, longitude, other_city.latitude,
                       other_city.longitude);
}

}  // namespace geo
}  // namespace nc
