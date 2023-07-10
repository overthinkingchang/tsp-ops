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

double calculate_total_distance(const Tour& tour, const Distances& distances) {
  double total_distance = 0;
  std::size_t num_cities = tour.size();

  for (std::size_t i = 0; i < num_cities; ++i) {
    total_distance +=
        distances[tour[i] * num_cities + tour[(i + 1) % num_cities]];
  }

  return total_distance;
}

double two_opt(Tour& tour, const Distances& distances) {
  std::size_t num_cities = tour.size();
  bool improvement = true;
  double best_distance = calculate_total_distance(tour, distances);

  while (improvement) {
    improvement = false;

    for (std::size_t i = 1; i < num_cities - 1; ++i) {
      for (std::size_t j = i + 1; j < num_cities; ++j) {
        // copy current tour to new tour
        // keep element from 0 to i - 1
        // and element from j + 1 to num_cities - 1
        // reverse the order of elements from i to j
        Tour new_tour = tour;
        std::reverse(new_tour.begin() + i, new_tour.begin() + j + 1);

        // D1: a -> b -> b1 -> b2 -> b3 -> c -> d
        // D2: a -> c -> b3 -> b2 -> b1 -> b -> d
        // D2 =  D1 - D_{ab} - D_{cd} + D_{ac} + D_{bd}

        double new_distance =
            best_distance - distances[tour[i - 1] * num_cities + tour[i]] -
            distances[tour[j] * num_cities + tour[(j + 1) % num_cities]] +
            distances[tour[i - 1] * num_cities + tour[j]] +
            distances[tour[i] * num_cities + tour[(j + 1) % num_cities]];

        if (best_distance - new_distance > 1e-10) {
          tour = new_tour;
          best_distance = new_distance;
          improvement = true;
        }
      }
    }

    if (improvement) {
      std::cout << "Improved distance: " << best_distance << "\n";
    }
  }

  return best_distance;
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

  std::vector<double> distances(num_cities * num_cities, 0);

  for (std::size_t i = 0; i < num_cities; ++i) {
    for (std::size_t j = i + 1; j < num_cities; ++j) {
      distances[i * num_cities + j] = calculate_distance(cities[i], cities[j]);
      distances[j * num_cities + i] = distances[i * num_cities + j];
    }
  }

  Tour tour(num_cities);
  std::iota(tour.begin(), tour.end(), 0);

  double total_distance = two_opt(tour, distances);

  std::cout << "Optimized tour: [";
  for (std::size_t i = 0; i < tour.size(); ++i) {
    std::cout << tour[i];
    if (i < tour.size() - 1) {
      std::cout << ", ";
    }
  }
  std::cout << "]\n";

  std::cout << "Total distance: " << total_distance << "\n";

  return 0;
}
