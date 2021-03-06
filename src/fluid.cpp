#include "glCanvas.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "argparser.h"
#include "boundingbox.h"
#include "fluid.h"
#include "marching_cubes.h"
#include "utils.h"

#define BETA_0 1.7
#define EPSILON 0.0001

// ==============================================================
// ==============================================================
// CONSTRUCTOR
// ==============================================================
// ==============================================================

Fluid::Fluid(ArgParser *_args) {
  steps = 0;
  args = _args;
  Load();
  marchingCubes = new MarchingCubes(nx + 1, ny + 1, nz + 1, dx, dy, dz, args->fluid_type);
  SetEmptySurfaceFull();
}

Fluid::~Fluid() {
  delete[] cells;
  delete marchingCubes;
  cleanupVBOs();
}

// ==============================================================

void Fluid::Load() {

  // open the file
  std::ifstream istr(std::string(args->path + '/' + args->fluid_file).c_str());
  assert(istr.good());
  std::string token, token2, token3;

  // load in the grid size & dimensions
  istr >> token >> nx >> ny >> nz;
  assert(token == "grid");
  assert(nx > 0 && ny > 0 && nz > 0);
  istr >> token >> dx >> dy >> dz;
  assert(token == "cell_dimensions");
  cells = new Cell[(nx + 2) * (ny + 2) * (nz + 2)];

  // simulation parameters
  istr >> token >> token2;
  assert(token == "flow");
  if (token2 == "compressible")
    compressible = true;
  else {
    assert(token2 == "incompressible");
    compressible = false;
  }
  istr >> token >> token2;
  assert(token == "xy_boundary");
  if (token2 == "free_slip")
    xy_free_slip = true;
  else {
    assert(token2 == "no_slip");
    xy_free_slip = false;
  }
  istr >> token >> token2;
  assert(token == "yz_boundary");
  if (token2 == "free_slip")
    yz_free_slip = true;
  else {
    assert(token2 == "no_slip");
    yz_free_slip = false;
  }
  istr >> token >> token2;
  assert(token == "zx_boundary");
  if (token2 == "free_slip")
    zx_free_slip = true;
  else {
    assert(token2 == "no_slip");
    zx_free_slip = false;
  }
  istr >> token >> viscosity;
  assert(token == "viscosity");
  double gravity;
  istr >> token >> gravity;
  assert(token == "gravity");
  args->gravity = glm::vec3(0, -9.8, 0) * float(gravity);

  // initialize marker particles
  istr >> token >> token2 >> token3;
  assert(token == "initial_particles");
  istr >> token >> density;
  assert(token == "density");
  GenerateParticles(token2, token3);

  // initialize velocities
  istr >> token >> token2;
  assert(token == "initial_velocity");
  if (token2 == "zero") {
    // default is zero
  } else {
    assert(token2 == "random");
    int i, j, k;
    double max_dim = std::max(dx, std::max(dy, dz));
    for (i = -1; i <= nx; i++) {
      for (j = -1; j <= ny; j++) {
        for (k = -1; k <= nz; k++) {
          getCell(i, j, k)->set_u_plus((2 * args->rand() - 1) * max_dim);
          getCell(i, j, k)->set_v_plus((2 * args->rand() - 1) * max_dim);
          getCell(i, j, k)->set_w_plus((2 * args->rand() - 1) * max_dim);
        }
      }
    }
  }
  // read in custom velocities
  while (istr >> token) {
    if (token == "barrier_cells" || token == "sources")
      break;

    int i, j, k;
    double velocity;
    assert(token == "u" || token == "v" || token == "w");
    istr >> i >> j >> k >> velocity;
    assert(i >= 0 && i < nx);
    assert(j >= 0 && j < ny);
    assert(k >= 0 && k < nz);
    jets.emplace_back(token, i, j, k, velocity);
    if (token == "u")
      getCell(i, j, k)->set_u_plus(velocity);
    else if (token == "v")
      getCell(i, j, k)->set_v_plus(velocity);
    else if (token == "w")
      getCell(i, j, k)->set_w_plus(velocity);
    else
      assert(0);
  }

  if (token == "barrier_cells") {
    while (istr >> token) {
      if (token == "sources")
        break;

      int i, j, k;

      istr >> i >> j >> k;
      assert(i >= 0 && i < nx);
      assert(j >= 0 && j < ny);
      assert(k >= 0 && k < nz);

      barrierCells.emplace(i, j, k);
    }
  }

  if (token == "sources") {
    while (istr >> token) {
      std::string shape;
      std::string placement;

      istr >> shape >> placement;

      sources.emplace_back(shape, placement);
    }
  }

  std::cout << barrierCells.size() << '\n';

  SetBoundaryVelocities();
}

// ==============================================================

bool Fluid::inShape(glm::vec3 &pos, const std::string &shape) {
  // return true if this point is inside the "shape"
  // defined procedurally (using an implicit surface)
  if (shape == "everywhere") {
    return true;
  } else if (shape == "left") {
    // a blob of particles on the lower left (for the dam)
    return (pos.x < 0.2 * nx * dx && pos.y < 0.5 * ny * dy);
  } else if (shape == "top_left") {
    return (pos.x < 0.1 * nx * dx && pos.y > 0.85 * ny * dy);
  } else if (shape == "drop") {
    // a shallow pool of particles on the bottom
    double h = ny * dy / 6.0;
    if (pos.y < 2 * h)
      return true;
    // and a sphere of particles above
    glm::vec3 center = glm::vec3(nx * dx * 0.5, 5 * h, nz * dz * 0.5);
    double length = glm::length(center - pos);
    if (length < 0.8 * h)
      return true;
    return false;
  } else if (shape == "middle") {
    double h = ny * dy / 6.0;

    glm::vec3 center = glm::vec3(nx * dx * 0.5, 5 * h, nz * dz * 0.5);
    double length = glm::length(center - pos);
    if (length < 0.8 * h)
      return true;
    return false;
  } else if (shape == "bottom") {
    double h = ny * dy / 6.0;

    glm::vec3 center = glm::vec3(nx * dx * 0.5, h, nz * dz * 0.5);
    double length = glm::length(center - pos);
    if (length < 0.8 * h)
      return true;
    return false;
  } else {
    std::cout << "unknown shape: " << shape << std::endl;
    exit(0);
  }
}

// ==============================================================

void Fluid::GenerateParticles(const std::string &shape,
                              const std::string &placement) {
  // create a set of points according to the "placement" token,
  // then check whether they are inside of the "shape"
  if (placement == "uniform") {
    int dens = (int)pow(density, 0.334);
    assert(dens * dens * dens == density);
    // the uniform grid spacing
    double spacing = 1 / double(dens);
    for (double x = 0.5 * spacing * dx; x < nx * dx; x += spacing * dx) {
      for (double y = 0.5 * spacing * dy; y < ny * dy; y += spacing * dy) {
        for (double z = 0.5 * spacing * dz; z < nz * dz; z += spacing * dz) {
          glm::vec3 pos = glm::vec3(x, y, z);
          if (inShape(pos, shape)) {
            Cell *cell = getCell(int(x / dx), int(y / dy), int(z / dz));
            FluidParticle *p = new FluidParticle();
            p->setPosition(pos);
            cell->addParticle(p);
          }
        }
      }
    }
  } else {
    assert(placement == "random");
    // note: we don't necessarily have the same number of particles in each cell
    for (int n = 0; n < nx * ny * nz * density; n++) {
      glm::vec3 pos = glm::vec3(args->rand() * nx * dx, args->rand() * ny * dy,
                                args->rand() * nz * dz);
      if (inShape(pos, shape)) {
        Cell *cell = getCell(int(pos.x / dx), int(pos.y / dy), int(pos.z / dz));
        FluidParticle *p = new FluidParticle();
        p->setPosition(pos);
        cell->addParticle(p);
      }
    }
  }
}

// ==============================================================
// ==============================================================
// ANIMATION
// ==============================================================
// ==============================================================

void Fluid::Animate() {

  // the animation manager:  this is what gets done each timestep!
  int new_particles = 3000;
  int new_barriers = 2000;

  if (steps % new_particles == 0) {
    for (const auto& s : sources) {
      GenerateParticles(std::get<0>(s), std::get<1>(s));
    }
  }

  ComputeNewVelocities();
  SetBoundaryVelocities();
  SetBarrierVelocities();
  SetJetVelocities();

  // compressible / incompressible flow
  if (compressible == false) {
    for (int iters = 0; iters < 20; iters++) {
      double max_divergence = AdjustForIncompressibility();
      SetBoundaryVelocities();
      SetBarrierVelocities();
      SetJetVelocities();
      if (max_divergence < EPSILON)
        break;
    }
  }

  UpdatePressures();

  // if (steps % new_barriers == 0 && steps > 0)
  //   CreateNewBarriers();

  CopyVelocities();

  // advanced the particles through the fluid
  MoveParticles();
  ReassignParticles();
  SetEmptySurfaceFull();

  setupVBOs();

  steps++;
}

// ==============================================================

void Fluid::ComputeNewVelocities() {
  double dt = args->timestep;
  int i, j, k;

  // using the formulas from Foster & Metaxas

  for (i = 0; i < nx - 1; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        Cell *cell = getCell(i, j, k);
        double new_u_plus =
            get_u_plus(i, j, k) +
            dt * ((1 / dx) * (square(get_u_avg(i, j, k)) -
                              square(get_u_avg(i + 1, j, k))) +
                  (1 / dy) * (get_uv_plus(i, j - 1, k) - get_uv_plus(i, j, k)) +
                  (1 / dz) * (get_uw_plus(i, j, k - 1) - get_uw_plus(i, j, k)) +
                  args->gravity.x +
                  (1 / dx) * (getPressure(i, j, k) - getPressure(i + 1, j, k)) +
                  (viscosity / square(dx)) *
                      (get_u_plus(i + 1, j, k) - 2 * get_u_plus(i, j, k) +
                       get_u_plus(i - 1, j, k)) +
                  (viscosity / square(dy)) *
                      (get_u_plus(i, j + 1, k) - 2 * get_u_plus(i, j, k) +
                       get_u_plus(i, j - 1, k)) +
                  (viscosity / square(dz)) *
                      (get_u_plus(i, j, k + 1) - 2 * get_u_plus(i, j, k) +
                       get_u_plus(i, j, k - 1)));
        cell->set_new_u_plus(new_u_plus);
      }
    }
  }

  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny - 1; j++) {
      for (k = 0; k < nz; k++) {
        Cell *cell = getCell(i, j, k);
        double new_v_plus =
            get_v_plus(i, j, k) +
            dt * ((1 / dx) * (get_uv_plus(i - 1, j, k) - get_uv_plus(i, j, k)) +
                  (1 / dy) * (square(get_v_avg(i, j, k)) -
                              square(get_v_avg(i, j + 1, k))) +
                  (1 / dz) * (get_vw_plus(i, j, k - 1) - get_vw_plus(i, j, k)) +
                  args->gravity.y +
                  (1 / dy) * (getPressure(i, j, k) - getPressure(i, j + 1, k)) +
                  (viscosity / square(dx)) *
                      (get_v_plus(i + 1, j, k) - 2 * get_v_plus(i, j, k) +
                       get_v_plus(i - 1, j, k)) +
                  (viscosity / square(dy)) *
                      (get_v_plus(i, j + 1, k) - 2 * get_v_plus(i, j, k) +
                       get_v_plus(i, j - 1, k)) +
                  (viscosity / square(dz)) *
                      (get_v_plus(i, j, k + 1) - 2 * get_v_plus(i, j, k) +
                       get_v_plus(i, j, k - 1)));
        cell->set_new_v_plus(new_v_plus);
      }
    }
  }

  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz - 1; k++) {
        Cell *cell = getCell(i, j, k);
        double new_w_plus =
            get_w_plus(i, j, k) +
            dt * ((1 / dx) * (get_uw_plus(i - 1, j, k) - get_uw_plus(i, j, k)) +
                  (1 / dy) * (get_vw_plus(i, j - 1, k) - get_vw_plus(i, j, k)) +
                  (1 / dz) * (square(get_w_avg(i, j, k)) -
                              square(get_w_avg(i, j, k + 1))) +
                  args->gravity.z +
                  (1 / dz) * (getPressure(i, j, k) - getPressure(i, j, k + 1)) +
                  (viscosity / square(dx)) *
                      (get_w_plus(i + 1, j, k) - 2 * get_w_plus(i, j, k) +
                       get_w_plus(i - 1, j, k)) +
                  (viscosity / square(dy)) *
                      (get_w_plus(i, j + 1, k) - 2 * get_w_plus(i, j, k) +
                       get_w_plus(i, j - 1, k)) +
                  (viscosity / square(dz)) *
                      (get_w_plus(i, j, k + 1) - 2 * get_w_plus(i, j, k) +
                       get_w_plus(i, j, k - 1)));
        cell->set_new_w_plus(new_w_plus);
      }
    }
  }
}

// ==============================================================
void Fluid::SetJetVelocities() {
  if (args->enable_jets) {
    for (const auto& jet: jets) {
      std::string direction;
      int i, j, k;
      double velocity;
      std::tie(direction, i, j, k, velocity) = jet;

      if (direction == "u")
        set_new_u_plus(i, j, k, velocity);
      else if (direction == "v")
        set_new_v_plus(i, j, k, velocity);
      else if (direction == "w")
        set_new_w_plus(i, j, k, velocity);
    }
  }
}

// ==============================================================

void Fluid::SetBoundaryVelocities() {

  // zero out flow perpendicular to the boundaries (no sources or sinks)
  for (int j = -1; j <= ny; j++) {
    for (int k = -1; k <= nz; k++) {
      getCell(-1, j, k)->set_u_plus(0);
      getCell(nx - 1, j, k)->set_u_plus(0);
      getCell(nx, j, k)->set_u_plus(0);
    }
  }
  for (int i = -1; i <= nx; i++) {
    for (int k = -1; k <= nz; k++) {
      getCell(i, -1, k)->set_v_plus(0);
      getCell(i, ny - 1, k)->set_v_plus(0);
      getCell(i, ny, k)->set_v_plus(0);
    }
  }
  for (int i = -1; i <= nx; i++) {
    for (int j = -1; j <= ny; j++) {
      getCell(i, j, -1)->set_w_plus(0);
      getCell(i, j, nz - 1)->set_w_plus(0);
      getCell(i, j, nz)->set_w_plus(0);
    }
  }

  // free slip or no slip boundaries (friction with boundary)
  double xy_sign = (xy_free_slip) ? 1 : -1;
  double yz_sign = (yz_free_slip) ? 1 : -1;
  double zx_sign = (zx_free_slip) ? 1 : -1;
  for (int i = 0; i < nx; i++) {
    for (int j = -1; j <= ny; j++) {
      getCell(i, j, -1)->set_u_plus(xy_sign * getCell(i, j, 0)->get_u_plus());
      getCell(i, j, nz)->set_u_plus(xy_sign *
                                    getCell(i, j, nz - 1)->get_u_plus());
    }
    for (int k = -1; k <= nz; k++) {
      getCell(i, -1, k)->set_u_plus(zx_sign * getCell(i, 0, k)->get_u_plus());
      getCell(i, ny, k)->set_u_plus(zx_sign *
                                    getCell(i, ny - 1, k)->get_u_plus());
    }
  }
  for (int j = 0; j < ny; j++) {
    for (int i = -1; i <= nx; i++) {
      getCell(i, j, -1)->set_v_plus(xy_sign * getCell(i, j, 0)->get_v_plus());
      getCell(i, j, nz)->set_v_plus(xy_sign *
                                    getCell(i, j, nz - 1)->get_v_plus());
    }
    for (int k = -1; k <= nz; k++) {
      getCell(-1, j, k)->set_v_plus(yz_sign * getCell(0, j, k)->get_v_plus());
      getCell(nx, j, k)->set_v_plus(yz_sign *
                                    getCell(nx - 1, j, k)->get_v_plus());
    }
  }
  for (int k = 0; k < nz; k++) {
    for (int i = -1; i <= nx; i++) {
      getCell(i, -1, k)->set_w_plus(zx_sign * getCell(i, 0, k)->get_w_plus());
      getCell(i, ny, k)->set_w_plus(zx_sign *
                                    getCell(i, ny - 1, k)->get_w_plus());
    }
    for (int j = -1; j <= ny; j++) {
      getCell(-1, j, k)->set_w_plus(yz_sign * getCell(0, j, k)->get_w_plus());
      getCell(nx, j, k)->set_w_plus(yz_sign *
                                    getCell(nx - 1, j, k)->get_w_plus());
    }
  }
}

// ==============================================================

void Fluid::SetBarrierVelocities() {
  for (const auto& barrier : barrierCells) {
    int i, j, k;
    std::tie(i, j, k) = barrier;

    getCell(i, j, k)->set_new_u_plus(0);
    getCell(i, j, k)->set_new_v_plus(0);
    getCell(i, j, k)->set_new_w_plus(0);

    getCell(i - 1, j, k)->set_new_u_plus(0);
    getCell(i, j - 1, k)->set_new_v_plus(0);
    getCell(i, j, k - 1)->set_new_w_plus(0);
  }
}

// ==============================================================

void Fluid::EmptyVelocities(int i, int j, int k) {
  Cell *c = getCell(i, j, k);
  if (c->getStatus() != CELL_EMPTY)
    return;
  Cell *ciplus = getCell(i + 1, j, k);
  Cell *cjplus = getCell(i, j + 1, k);
  Cell *ckplus = getCell(i, j, k + 1);
  if (ciplus->getStatus() == CELL_EMPTY)
    c->set_new_u_plus(0);
  if (cjplus->getStatus() == CELL_EMPTY)
    c->set_new_v_plus(0);
  if (ckplus->getStatus() == CELL_EMPTY)
    c->set_new_w_plus(0);
}

// move to new timestep
void Fluid::CopyVelocities() {
  double dt = args->timestep;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        Cell *c = getCell(i, j, k);

        EmptyVelocities(i, j, k);

        c->copyVelocity();
        if (fabs(c->get_u_plus()) > 0.5 * dx / dt ||
            fabs(c->get_v_plus()) > 0.5 * dy / dt ||
            fabs(c->get_w_plus()) > 0.5 * dz / dt) {
          // velocity has exceeded reasonable threshhold
          std::cout << "velocity has exceeded reasonable threshhold, stopping "
                       "animation"
                    << std::endl;
          // args->animate = false;
        }
      }
    }
  }
}

// ==============================================================

double Fluid::ComputeDivergence(int i, int j, int k) const {
  double divergence =
      -((1 / dx) * (get_new_u_plus(i, j, k) - get_new_u_plus(i - 1, j, k)) +
        (1 / dy) * (get_new_v_plus(i, j, k) - get_new_v_plus(i, j - 1, k)) +
        (1 / dz) * (get_new_w_plus(i, j, k) - get_new_w_plus(i, j, k - 1)));

  return divergence;
}

// ==============================================================

double Fluid::ComputeDP(double divergence) const {
  double dt = args->timestep;
  double beta =
      BETA_0 / ((2 * dt) * (1 / square(dx) + 1 / square(dy) + 1 / square(dz)));
  double dp = beta * divergence;

  return dp;
}

// ==============================================================

double Fluid::AdjustForIncompressibility() {
  double max_divergence = 0;

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        Cell *c = getCell(i, j, k);
        if (c->getStatus() != CELL_EMPTY) {
          // compute divergence and increment/decrement pressure
          double divergence = ComputeDivergence(i, j, k);
          max_divergence = std::max(divergence, max_divergence);

          double dt = args->timestep;
          double dp = ComputeDP(divergence);

          c->adjust_new_u_plus(dt * dp / dx);
          c->adjust_new_v_plus(dt * dp / dy);
          c->adjust_new_w_plus(dt * dp / dz);

          adjust_new_u_plus(i - 1, j, k, -dt * dp / dx);
          adjust_new_v_plus(i, j - 1, k, -dt * dp / dy);
          adjust_new_w_plus(i, j, k - 1, -dt * dp / dz);
        }
      }
    }
  }

  // return the maximum divergence
  // (will iterate for specified # of iterations or until divergence is near
  // zero)
  return max_divergence;
}

// ==============================================================

void Fluid::UpdatePressures() {
  for (int i = -1; i <= nx; i++) {
    for (int j = -1; j <= ny; j++) {
      for (int k = -1; k <= nz; k++) {
        Cell *c = getCell(i, j, k);
        if (i >= 0 && i < nx && j >= 0 && j < ny && k >= 0 && k < nz) {
          // compute divergence and increment/decrement pressure
          double pressure = c->getPressure();
          double divergence = ComputeDivergence(i, j, k);

          double dp = ComputeDP(divergence);
          c->setPressure(pressure + dp);
        } else {
          // zero out boundary cells (just in case)
          c->setPressure(0);
        }

        // =======================================
        // HACK? From Foster 2001 paper?
        // zero out empty cells
        if (c->getStatus() == CELL_EMPTY) {
          c->setPressure(0);
        }
        // ========================================
      }
    }
  }
}

// ==============================================================

void Fluid::MoveParticles() {
  double dt = args->timestep;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        Cell *cell = getCell(i, j, k);
        std::vector<FluidParticle *> &particles = cell->getParticles();
        for (unsigned int iter = 0; iter < particles.size(); iter++) {
          FluidParticle *p = particles[iter];
          glm::vec3 pos = p->getPosition();
          glm::vec3 vel = getInterpolatedVelocity(pos);
          glm::vec3 pos2 = pos + float(dt) * vel;
#if 0
          // euler integration
          p->setPosition(pos2);
#else
          // trapezoid integration
          glm::vec3 vel2 = getInterpolatedVelocity(pos2);
          glm::vec3 pos3 = pos + float(0.5 * dt) * (vel + vel2);
          p->setPosition(pos3);
#endif
        }
      }
    }
  }
}

// ==============================================================

void Fluid::ReassignParticles() {
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        Cell *cell = getCell(i, j, k);
        std::vector<FluidParticle *> &particles = cell->getParticles();
        for (unsigned int iter = 0; iter < particles.size(); iter++) {
          FluidParticle *p = particles[iter];
          glm::vec3 pos = p->getPosition();
          int i2 =
              (int)std::min(double(nx - 1), std::max(0.0, floor(pos.x / dx)));
          int j2 =
              (int)std::min(double(ny - 1), std::max(0.0, floor(pos.y / dy)));
          int k2 =
              (int)std::min(double(nz - 1), std::max(0.0, floor(pos.z / dz)));
          // if the particle has crossed one of the cell faces
          // assign it to the new cell
          if (i != i2 || j != j2 || k != k2) {
            cell->removeParticle(p);
            getCell(i2, j2, k2)->addParticle(p);
          }
        }
      }
    }
  }
}

// ==============================================================

void Fluid::SetEmptySurfaceFull() {
  int i, j, k;
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        Cell *cell = getCell(i, j, k);
        // if (cell->numParticles() == 0 && barrierCells.count({i, j, k}) == 0)
        if (cell->numParticles() == 0)
          cell->setStatus(CELL_EMPTY);
        else
          cell->setStatus(CELL_FULL);
      }
    }
  }

  // pick out the boundary cells
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        Cell *cell = getCell(i, j, k);
        if (cell->getStatus() == CELL_FULL &&
            (getCell(i - 1, j, k)->getStatus() == CELL_EMPTY ||
             getCell(i + 1, j, k)->getStatus() == CELL_EMPTY ||
             getCell(i, j - 1, k)->getStatus() == CELL_EMPTY ||
             getCell(i, j + 1, k)->getStatus() == CELL_EMPTY ||
             getCell(i, j, k - 1)->getStatus() == CELL_EMPTY ||
             getCell(i, j, k + 1)->getStatus() == CELL_EMPTY)) {
          cell->setStatus(CELL_SURFACE);
        }
      }
    }
  }
}

// ==============================================================

glm::vec3 Fluid::getInterpolatedVelocity(const glm::vec3 &pos) const {
  int i = int(floor(pos.x / dx));
  int j = int(floor(pos.y / dy));
  int k = int(floor(pos.z / dz));

  if (i < 0 || i >= nx || j < 0 || j >= ny || k < 0 || k >= nz) {
    i = std::min(std::max(0, i), nx - 1);
    j = std::min(std::max(0, j), ny - 1);
    k = std::min(std::max(0, k), nz - 1);

    return glm::vec3(get_u_avg(i, j, k), get_v_avg(i, j, k),
                     get_w_avg(i, j, k));
  }

  double u = 0;
  double v = 0;
  double w = 0;

  for (int di = -1; di <= 1; di++) {
    for (int dj = -1; dj <= 1; dj++) {
      for (int dk = -1; dk <= 1; dk++) {
        Cell *cell = getCell(i + di, j + dj, k + dk);
        double cx = dx * (i + di) + (dx / 2);
        double cy = dy * (j + dj) + (dy / 2);
        double cz = dz * (k + dk) + (dz / 2);

        double width, height, depth;

        // Calculating u
        width = dx - abs(pos.x - cx - (dx / 2));
        height = dy - abs(pos.y - cy);
        depth = dz - abs(pos.z - cz);

        if (width >= 0 && height >= 0 && depth >= 0)
          u += cell->get_u_plus() * width * height * depth;

        // Calculating v
        width = dx - abs(pos.x - cx);
        height = dy - abs(pos.y - cy - (dy / 2));
        depth = dz - abs(pos.z - cz);

        if (width >= 0 && height >= 0 && depth >= 0)
          v += cell->get_v_plus() * width * height * depth;

        // Calculating w
        width = dx - abs(pos.x - cx);
        height = dy - abs(pos.y - cy);
        depth = dz - abs(pos.z - cz - (dz / 2));

        if (width >= 0 && height >= 0 && depth >= 0)
          w += cell->get_w_plus() * width * height * depth;
      }
    }
  }

  glm::vec3 interpolated(u, v, w);
  interpolated /= (dx * dy * dz);

  return interpolated;
}

// ==============================================================

void Fluid::CreateNewBarriers() {
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        Cell *c = getCell(i, j, k);

        if (c->getStatus() != CELL_EMPTY) {
          float norm1 = glm::length(glm::vec3(c->get_new_u_plus(), c->get_new_v_plus(), c->get_new_w_plus()));
          float norm2 = glm::length(glm::vec3(get_new_u_plus(i - 1, j, k), get_new_v_plus(i, j - 1, k), get_new_w_plus(i, j, k - 1)));
          if (norm1 < EPSILON * 100 && norm2 < EPSILON * 100) {
            c->setStatus(CELL_EMPTY);
            c->deleteAllParticles();
            barrierCells.emplace(i, j, k);
          }
        }

        // if (c->getStatus() == CELL_FULL) {
        //   c->setStatus(CELL_SURFACE);
        //   c->deleteAllParticles();
        //   barrierCells.emplace(i, j, k);
        // } else if (c->getStatus() == CELL_SURFACE) {
        //   c->setStatus(CELL_EMPTY);
        //   c->deleteAllParticles();
        // }
      }
    }
  }
}