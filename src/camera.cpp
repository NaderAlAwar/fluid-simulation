#include "camera.h"
#include "glCanvas.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>

// ====================================================================
// CONSTRUCTORS
// ====================================================================

Camera::Camera(const glm::vec3 &c, const glm::vec3 &poi, const glm::vec3 &u) {
  camera_position = c;
  point_of_interest = poi;
  up = glm::normalize(u);
}

OrthographicCamera::OrthographicCamera(const glm::vec3 &c, const glm::vec3 &poi,
                                       const glm::vec3 &u, float s)
    : Camera(c, poi, u) {
  size = s;
}

PerspectiveCamera::PerspectiveCamera(const glm::vec3 &c, const glm::vec3 &poi,
                                     const glm::vec3 &u, float a)
    : Camera(c, poi, u) {
  angle = a;
}

// ====================================================================
// Construct the ViewMatrix & ProjectionMatrix for GL Rendering
// ====================================================================

void OrthographicCamera::glPlaceCamera() {
  glfwGetWindowSize(GLCanvas::window, &width, &height);
  float aspect = width / (float)height;
  float w;
  float h;
  // handle non square windows
  if (aspect < 1.0) {
    w = size / 2.0;
    h = w / aspect;
  } else {
    h = size / 2.0;
    w = h * aspect;
  }
  ProjectionMatrix = glm::ortho<float>(-w, w, -h, h, 0.1f, 100.0f);
  ViewMatrix = glm::lookAt(camera_position, point_of_interest, getScreenUp());
}

void PerspectiveCamera::glPlaceCamera() {
  glfwGetWindowSize(GLCanvas::window, &width, &height);
  float aspect = width / (float)height;
  // must convert angle from degrees to radians
  ProjectionMatrix =
      glm::perspective<float>(glm::radians(angle), aspect, 0.1f, 1000.0f);
  ViewMatrix = glm::lookAt(camera_position, point_of_interest, getScreenUp());
}

// ====================================================================
// dollyCamera: Move camera along the direction vector
// ====================================================================

void Camera::dollyCamera(float dist) {
  glm::vec3 diff = camera_position - point_of_interest;
  float d = glm::length(diff);
  glm::vec3 translate = float(0.004 * d * dist) * getDirection();
  camera_position += translate;
}

// ====================================================================
// zoomCamera: Change the field of view/angle
// ====================================================================

void OrthographicCamera::zoomCamera(float factor) {
  size *= pow(1.003, factor);
}

void PerspectiveCamera::zoomCamera(float dist) {
  angle *= pow(1.002, dist);
  // put some reasonable limits on the camera angle (in degrees)
  if (angle < 5)
    angle = 5;
  if (angle > 175)
    angle = 175;
}

// ====================================================================
// truckCamera: Translate camera perpendicular to the direction vector
// ====================================================================

void Camera::truckCamera(float dx, float dy) {
  glm::vec3 diff = camera_position - point_of_interest;
  float d = glm::length(diff);
  glm::vec3 translate =
      (d * 0.0007f) * (getHorizontal() * float(dx) + getScreenUp() * float(dy));
  camera_position += translate;
  point_of_interest += translate;
}

// ====================================================================
// rotateCamera: Rotate around the up and horizontal vectors
// ====================================================================

// adjust this number if desired
#define ROTATE_SPEED 0.05

void Camera::rotateCamera(float rx, float ry) {

  // this version of rotate doesn't let the model flip "upside-down"

  glm::vec4 tmp(camera_position, 1);
  glm::vec4 tmpUp(camera_position + up, 1);

  glm::vec3 axis = -ry * getHorizontal() + rx * up;
  tmpUp = glm::rotate<GLfloat>(ROTATE_SPEED, axis) * tmpUp;
  tmp = glm::rotate<GLfloat>(ROTATE_SPEED, axis) * tmp;
  tmpUp -= tmp;

  for (int i = 0; i < 3; i++) {
    up[i] = tmpUp[i];
    camera_position[i] = tmp[i];
  }
}

// ====================================================================
// ====================================================================

inline std::ostream &operator<<(std::ostream &ostr, const glm::vec3 &v) {
  ostr << "<" << v.x << "," << v.y << "," << v.z << ">";
  return ostr;
}

inline std::istream &operator>>(std::istream &istr, glm::vec3 &v) {
  char c;
  istr >> c;
  assert(c == '<');
  istr >> v.x;
  istr >> c;
  assert(c == ',');
  istr >> v.y;
  istr >> c;
  assert(c == ',');
  istr >> v.z;
  istr >> c;
  assert(c == '>');
  return istr;
}

std::ostream &operator<<(std::ostream &ostr, const Camera &c) {
  const Camera *cp = &c;
  if (dynamic_cast<const OrthographicCamera *>(cp)) {
    const OrthographicCamera *ocp = (const OrthographicCamera *)cp;
    ostr << *ocp << std::endl;
  } else if (dynamic_cast<const PerspectiveCamera *>(cp)) {
    const PerspectiveCamera *pcp = (const PerspectiveCamera *)cp;
    ostr << *pcp << std::endl;
  }
  return ostr;
}

std::ostream &operator<<(std::ostream &ostr, const OrthographicCamera &c) {
  ostr << "OrthographicCamera {" << std::endl;
  ostr << "    camera_position   " << c.camera_position << std::endl;
  ostr << "    point_of_interest " << c.point_of_interest << std::endl;
  ostr << "    up                " << c.up << std::endl;
  ostr << "    size              " << c.size << std::endl;
  ostr << "}" << std::endl;
  return ostr;
}

std::ostream &operator<<(std::ostream &ostr, const PerspectiveCamera &c) {
  ostr << "PerspectiveCamera {" << std::endl;
  ostr << "  camera_position    " << c.camera_position << std::endl;
  ostr << "  point_of_interest  " << c.point_of_interest << std::endl;
  ostr << "  up                 " << c.up << std::endl;
  ostr << "  angle              " << c.angle << std::endl;
  ostr << "}" << std::endl;
  return ostr;
}

std::istream &operator>>(std::istream &istr, OrthographicCamera &c) {
  std::string token;
  istr >> token;
  assert(token == "{");
  istr >> token;
  assert(token == "camera_position");
  istr >> c.camera_position;
  istr >> token;
  assert(token == "point_of_interest");
  istr >> c.point_of_interest;
  istr >> token;
  assert(token == "up");
  istr >> c.up;
  istr >> token;
  assert(token == "size");
  istr >> c.size;
  istr >> token;
  assert(token == "}");
  return istr;
}

std::istream &operator>>(std::istream &istr, PerspectiveCamera &c) {
  std::string token;
  istr >> token;
  assert(token == "{");
  istr >> token;
  assert(token == "camera_position");
  istr >> c.camera_position;
  istr >> token;
  assert(token == "point_of_interest");
  istr >> c.point_of_interest;
  istr >> token;
  assert(token == "up");
  istr >> c.up;
  istr >> token;
  assert(token == "angle");
  istr >> c.angle;
  istr >> token;
  assert(token == "}");
  return istr;
}
