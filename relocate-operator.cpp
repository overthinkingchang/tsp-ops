#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

using City = std::pair<unsigned int, unsigned int>;
using Cities = std::vector<City>;
using Tour = std::vector<unsigned int>;
using Distances = std::vector<double>;

double calculate_distance(const City& city1, const City& city2) {
  return std::sqrt((city2.first - city1.first) * (city2.first - city1.first) +
                   (city2.second - city1.second) *
                       (city2.second - city1.second));
}

double calculate_total_distance(const Tour& tour, const Distances& distance) {
  double total_distance = 0;
  std::size_t num_cities = tour.size();

  for (std::size_t i = 0; i < num_cities; ++i) {
    total_distance +=
        distance[tour[i] * num_cities + tour[(i + 1) % num_cities]];
  }

  return total_distance;
}

std::pair<double, Tour> relocate_operator(const Tour& tour,
                                          const Distances& distance) {
  double best_distance = calculate_total_distance(tour, distance);
  Tour best_tour = tour;
  std::size_t num_cities = tour.size();

  for (std::size_t i = 0; i < num_cities - 1; ++i) {
    for (std::size_t j = i + 2; j < num_cities; ++j) {
      // Remove city i + 1 and insert it after city j
      Tour new_tour(num_cities);

      // copy from 0 to i
      auto new_tour_iterator =
          std::copy_n(tour.begin(), i + 1, new_tour.begin());

      // copy from i + 2 to j
      // j + 1 - (i + 2) = j - i - 1
      new_tour_iterator =
          std::copy_n(tour.begin() + i + 2, j - i - 1, new_tour_iterator);

      // assign i + 1 to j + 1
      *new_tour_iterator = tour[i + 1];
      ++new_tour_iterator;

      // copy from j + 1 to end
      std::copy(tour.begin() + j + 1, tour.end(), new_tour_iterator);

      double new_distance = calculate_total_distance(new_tour, distance);

      if (new_distance < best_distance) {
        best_distance = new_distance;
        best_tour = new_tour;
      }
    }
  }

  return std::make_pair(best_distance, best_tour);
}

Cities import_tsp_data(const std::string& filename) {
  Cities cities;

  std::ifstream is(filename);
  for (std::string line; std::getline(is, line);) {
    if (std::isdigit(line[0])) {
      std::size_t first_whitespace = line.find(' ');
      std::size_t second_whitespace = line.find(' ', first_whitespace + 1);

      unsigned int x = std::stoul(line.substr(
          first_whitespace + 1, second_whitespace - first_whitespace - 1));

      unsigned int y = std::stoul(line.substr(second_whitespace + 1));

      cities.emplace_back(x, y);
    }
  }

  return cities;
}

int main() {
  Cities cities = import_tsp_data("pbn423.tsp");
  std::size_t num_cities = cities.size();

  std::vector<double> distances(num_cities * num_cities);

  for (std::size_t i = 0; i < num_cities; ++i) {
    for (std::size_t j = i + 1; j < num_cities; ++j) {
      distances[i * num_cities + j] = calculate_distance(cities[i], cities[j]);
      distances[j * num_cities + i] = distances[i * num_cities + j];
    }
  }

  Tour best_tour(num_cities);
  std::iota(best_tour.begin(), best_tour.end(), 0);
  double best_distance = calculate_total_distance(best_tour, distances);

  bool improved = true;
  while (improved) {
    auto result = relocate_operator(best_tour, distances);
    double new_distance = result.first;
    Tour new_tour = result.second;

    if (new_distance < best_distance) {
      best_tour = new_tour;
      best_distance = new_distance;
      std::cout << "Improved distance: " << best_distance << "\n";
    } else {
      improved = false;
    }
  }

  std::cout << "Optimized tour: [";
  for (std::size_t i = 0; i < best_tour.size(); ++i) {
    std::cout << best_tour[i];
    if (i < best_tour.size() - 1) {
      std::cout << ", ";
    }
  }
  std::cout << "]\n";

  std::cout << "Total distance: " << best_distance << "\n";

  return 0;
}
