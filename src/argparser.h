#ifndef __ARG_PARSER_H__
#define __ARG_PARSER_H__

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>

// ====================================================================
// ====================================================================

inline void separatePathAndFile(const std::string &input, std::string &path,
                                std::string &file) {
  // we need to separate the filename from the path
  // (we assume the vertex & fragment shaders are in the same directory)
  // first, locate the last '/' in the filename
  size_t last = std::string::npos;
  while (1) {
    int next = input.find('/', last + 1);
    if (next != (int)std::string::npos) {
      last = next;
      continue;
    }
    next = input.find('\\', last + 1);
    if (next != (int)std::string::npos) {
      last = next;
      continue;
    }
    break;
  }
  if (last == std::string::npos) {
    // if there is no directory in the filename
    file = input;
    path = ".";
  } else {
    // separate filename & path
    file = input.substr(last + 1, input.size() - last - 1);
    path = input.substr(0, last);
  }
}

// ====================================================================
// ====================================================================

class ArgParser {

public:
  ArgParser() { DefaultValues(); }

  ArgParser(int argc, char *argv[]) {
    DefaultValues();
    // parse the command line arguments
    for (int i = 1; i < argc; i++) {
      if (argv[i] == std::string("-fluid")) {
        i++;
        assert(i < argc);
        separatePathAndFile(argv[i], path, fluid_file);
      } else if (argv[i] == std::string("-size")) {
        i++;
        assert(i < argc);
        width = height = atoi(argv[i]);
      } else if (argv[i] == std::string("-timestep")) {
        i++;
        assert(i < argc);
        timestep = atof(argv[i]);
        assert(timestep > 0);
      } else if (argv[i] == std::string("-fluid_type")) {
        i++;
        assert(i < argc);
        fluid_type = argv[i];
      } else if (argv[i] == std::string("-enable_jets")) {
        enable_jets = true;
      } else if (argv[i] == std::string("-record")) {
        record = true;
      } else {
        std::cout << "ERROR: unknown command line argument " << i << ": '"
                  << argv[i] << "'" << std::endl;
        exit(1);
      }
    }
  }

  double rand() {
    static std::random_device rd;
    static std::mt19937 engine(rd());
    static std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(engine);
  }

  void DefaultValues() {
    width = 2000;
    height = 1000;

    timestep = 0.01;
    animate = false;

    particles = true;
    velocity = true;

    face_velocity = 0;
    dense_velocity = 0;

    surface = false;
    isosurface = 0.7;

    bounding_box = true;
    cubes = false;
    pressure = false;
    fluid_type = "other";
    enable_jets = false;
    record = false;

    gravity = glm::vec3(0, -9.8, 0);

    perspective = true;
  }

  // ==============
  // REPRESENTATION
  // all public! (no accessors)

  std::string fluid_file;
  std::string path;
  int width;
  int height;

  // animation control
  double timestep;
  bool animate;
  glm::vec3 gravity;

  // display option toggles
  // (used by both)
  bool particles;
  bool velocity;
  bool surface;
  bool bounding_box;
  std::string fluid_type;
  bool enable_jets;
  bool record;

  // used by fluid
  int face_velocity;
  int dense_velocity;
  double isosurface;
  bool cubes;
  bool pressure;

  bool perspective;
};

#endif
