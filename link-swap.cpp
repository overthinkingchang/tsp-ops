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

std::pair<double, Tour> link_swap(const Tour& tour, double distance,
                                  const Distances& distances) {
  std::size_t num_cities = tour.size();
  double best_distance = distance;

  unsigned int reverse_method = 0;
  unsigned int reverse_index = 0;

  for (std::size_t i = 1; i < num_cities - 1; ++i) {
    // i-th edge: connect i-th node and (i + 1)-th node
    // num_edges = num_cities - 1
    // for example 0 -> 1 -> 2 -> 3 -> 4 -> 5 -> 6 -> 7 -> 8 -> 0
    // remove 0-th edge => the same distance
    // 0 -> 8 -> 7 -> 6 -> 5 -> 4 -> 3 -> 2 -> 1 -> 0
    // remove 1-th edge => only one neighbor
    // 0 -> 1 /-> 8 -> 7 -> 6 -> 5 -> 4 -> 3 -> 2 /-> 0 (reverse second)
    // from 2-th edge => 3 neighbors
    // 0 /-> 1 -> 2 /-> 8 -> 7 -> 6 -> 5 -> 4 -> 3 /-> 0 (reverse second)
    // 0 /-> 2 -> 1 /-> 3 -> 4 -> 5 -> 6 -> 7 -> 8 /-> 0 (reverse first)
    // 0 /-> 2 -> 1 /-> 8 -> 7 -> 6 -> 5 -> 4 -> 3 /-> 0 (reverse both)
    // remove 3-th edge => example
    // 0 /-> 1 -> 2 -> 3 /-> 8 -> 7 -> 6 -> 5 -> 4 /-> 0 (reverse second)
    // 0 /-> 3 -> 2 -> 1 /-> 4 -> 5 -> 6 -> 7 -> 8 /-> 0 (reverse first)
    // 0 /-> 3 -> 2 -> 1 /-> 8 -> 7 -> 6 -> 5 -> 4 /-> 0 (reverse both)
    // reverse_method
    // 0: no reverse
    // 1: reverse second
    // 2: reverse first
    // 3: reverse both

    // tour[0] always = 0

    // d1 = d0 - d_{i, i + 1} - d_{last, 0} + d_{i, last} + d_{i + 1, 0}
    double d1 = distance - distances[tour[i] * num_cities + tour[i + 1]] -
                distances[tour[num_cities - 1]] +
                distances[tour[i] * num_cities + tour[num_cities - 1]] +
                distances[tour[i + 1]];

    // d2 = d0 - d_{i, i + 1} - d_{0, 1} + d_{i, 0} + d_{1, i + 1}
    double d2 = distance - distances[tour[i] * num_cities + tour[i + 1]] -
                distances[tour[1]] + distances[tour[i]] +
                distances[tour[1] * num_cities + tour[i + 1]];

    // d3 = d0 - d_{i, i + 1} - d_{0, 1} - d_{last, 0}
    // + d_{i, 0} + d_{1, last} + d_{i + 1, 0}
    double d3 = distance - distances[tour[i] * num_cities + tour[i + 1]] -
                distances[tour[1]] - distances[tour[num_cities - 1]] +
                distances[tour[i]] +
                distances[tour[1] * num_cities + tour[num_cities - 1]] +
                distances[tour[i + 1]];

    if (best_distance - d1 > 1e-10) {
      best_distance = d1;
      reverse_method = 1;
      reverse_index = i;
    }
    if (best_distance - d2 > 1e-10) {
      best_distance = d2;
      reverse_method = 2;
      reverse_index = i;
    }
    if (best_distance - d3 > 1e-10) {
      best_distance = d3;
      reverse_method = 3;
      reverse_index = i;
    }
  }

  Tour best_tour = tour;
  if (reverse_method == 1) {
    std::reverse(best_tour.begin() + reverse_index + 1, best_tour.end());
  } else if (reverse_method == 2) {
    std::reverse(best_tour.begin() + 1, best_tour.begin() + reverse_index + 1);
  } else if (reverse_method == 3) {
    std::reverse(best_tour.begin() + 1, best_tour.begin() + reverse_index + 1);
    std::reverse(best_tour.begin() + reverse_index + 1, best_tour.end());
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
    auto result = link_swap(best_tour, best_distance, distances);
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
