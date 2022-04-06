#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <stdio.h>
#include <vector>

const double EPS = 1e-7;

class Point {
public:
  double x, y;
  
  Point() = default;
  Point(double x1, double y1) {
    x = x1;
    y = y1;
  };
  bool operator == (const Point& another) const{
    if (abs(x - another.x) < EPS && abs(y - another.y) < EPS) return 1;
    return 0;
    
  }
  
  bool operator != (const Point& another) {
    return !(*this == another);
  }
  
  Point rotate(Point center, double angle) {
    angle = round(angle * 10000000) / 10000000;
    //angle = angle * M_PI / 180;
    double new_point_x = center.x + (x - center.x) * cos(angle) - (y - center.y) * sin(angle);
    double new_point_y = center.y + (x - center.x) * sin(angle) + (y - center.y) * cos(angle);
    return Point(round(new_point_x * 10000000) / 10000000, round(new_point_y * 10000000) / 10000000);
  }
  
  void reflex(Point center) {
    x = 2 * center.x - x;
    y = 2 * center.y - y;
  }
  
  void scale(Point center, double coefficient) {
    x = center.x + coefficient * (x - center.x);
    y = center.y + coefficient * (y - center.y);
    x = round(x * 10000000) / 10000000;
  }
};

double dist_points(const Point& a, const Point& b) {
  return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

class Vector {
public:
  double x, y;
  double len;
  
  Vector() = default;
  Vector(double x_new, double y_new): x(x_new), y(y_new) {
    len = dist_points(Point(x, y), Point(0, 0));
  }
  explicit Vector(const Point& a): x(a.x), y(a.y) {
    len = dist_points(a, Point(0, 0));
  }
  Vector(const Point& a, const Point& b) {
    x = b.x - a.x;
    y = b.y - a.y;
    len = dist_points(a, b);
  }
  Vector rotate(double degree) {
    double t_x = Point(x, y).rotate(Point(0, 0), degree).x;
    double t_y = Point(x, y).rotate(Point(0, 0), degree).y;
    return Vector(Point(t_x, t_y));
  }
};

double scalar_product(const Vector& a, const Vector& b) {
  return (a.x * b.x + a.y * b.y) / (a.len * b.len);
}

double degree_vectors(const Vector& a, const Vector& b) {
  return (a.x * b.x + a.y * b.y) / (a.len * b.len);
}

class Line {
public:
  double A, B, C;
  
  Line() = default;
  Line(double A_new, double B_new, double C_new): A(A_new), B(B_new), C(C_new)
  {
  }
  Line(const Point& point_a, const Point& point_b) {
    A = point_b.y - point_a.y;
    B = point_a.x - point_b.x;
    C = point_a.y * point_b.x - point_a.x * point_b.y;
    if (A < 0 || (A == 0 && B < 0)) {
      A *= -1;
      B *= -1;
      C *= -1;
    }
  }
  Line(double k, double b) {
    A = -k;
    B = 1;
    C = -b;
    if (A < 0 || (A == 0 && B < 0)) {
      A *= -1;
      B *= -1;
      C *= -1;
    }
  }
  Line(const Point& point, double k) {
    // b = y - kx
    A = -k;
    B = 1;
    C = -(point.y - k * point.x);
    if (A < 0 || (A == 0 && B < 0)) {
      A *= -1;
      B *= -1;
      C *= -1;
    }
  }
  
  bool operator == (const Line& another) {
    double x1 = 0;
    double y11 = (-C - A * x1) / B, y12 = (-another.C - another.A * x1) / another.B;
    double x2 = 1;
    double y21 = (-C - A * x2) / B, y22 = (-another.C - another.A * x2) / another.B;
    if (abs(y11 - y12) < EPS && abs(y21 - y22) < EPS) {
      return 1;
    }
    return 0;
  }
  
  bool operator != (const Line& another) {
    return !(*this == another);
  }
  
};

Point reflex_point(Point a, Line axis) {
  double new_point_x = ((axis.B * axis.B - axis.A * axis.A) * a.x - 2 * axis.A * axis.B * a.y - 2 * axis.A * axis.C) / (axis.A * axis.A + axis.B * axis.B);
  double new_point_y = ((axis.A * axis.A - axis.B * axis.B) * a.y - 2 * axis.A * axis.B * a.x - 2 * axis.B * axis.C) / (axis.A * axis.A + axis.B * axis.B);
  return Point(new_point_x, new_point_y);
}

Point line_intersection(const Line& a, const Line& b) {
  return Point((a.B * b.C - b.B * a.C) / (a.A * b.B - a.B * b.A), (a.C * b.A - b.C * a.A) / (a.A * b.B - a.B * b.A));
}

double dist_point_line(const Point& a, const Line& b) {
  return abs((b.A * a.x + b.B * a.y + b.C) / sqrt(b.A * b.A + b.B * b.B));
}

class Shape {
public:
  virtual double perimeter() const = 0;
  
  virtual double area() const = 0;
  
  virtual bool operator==(const Shape& another) const;
  
  virtual bool operator!=(const Shape& another) const{
    return !(*this == another);
  }
  
  virtual bool isCongruentTo(const Shape& another);
  
  virtual bool isSimilarTo(const Shape& another);
  
  virtual bool containsPoint(Point point) = 0;
  
  virtual void rotate(Point center, double angle) = 0;
  
  virtual void reflex(Point center) = 0;
  
  virtual void reflex(Line axis) = 0;
  
  virtual void scale(Point center, double coefficient) = 0;
  
  virtual ~Shape() = 0;
};

Shape::~Shape() {}

class Polygon: public Shape {
protected:
  std::vector<Point> points_;
public:
  Polygon() = default;
  explicit Polygon(const std::vector<Point>& points) {

    points_ = points;
  }
  
  explicit Polygon(const std::initializer_list<Point>& points) {

    points_ = points;
  }
  
  size_t verticesCount() {
    return this->points_.size();
  }
  
  std::vector<Point> getVertices() {
    return this->points_;
  }
  
  bool isConvex() {
    int sign = 0;
    std::vector<Point> points_first_last;
    points_first_last.push_back(points_[points_.size() - 1]);
    for (auto point: points_)
      points_first_last.push_back(point);
    points_first_last.push_back(points_[0]);
    for (size_t i = 1; i < points_first_last.size() - 1; ++i) {
      Point vec1(points_[i].x - points_[i - 1].x,
        points_[i].y - points_[i - 1].y);
      Point vec2(points_[i + 1].x - points_[i].x,
        points_[i + 1].y - points_[i].y);
      double product = vec1.x * vec2.y - vec1.y * vec2.x;
      if (sign == 0) {
        if (product < 0) sign = -1;
        else sign = 1;
      }
      else {
        if (sign * product < 0) {
          return 0;
        }
      }
    }
    return 1;
  }
  
  virtual double area() const {
    double ans = 0;
    for (size_t i = 0; i < points_.size() - 1; ++i) {
      ans += points_[i].x * points_[i + 1].y;
    }
    ans += points_[points_.size() - 1].x * points_[0].y;
    for (size_t i = 0; i < points_.size() - 1; ++i) {
      ans -= points_[i + 1].x * points_[i].y;
    }
    ans -= points_[points_.size() - 1].y * points_[0].x;
    return abs(ans) * 0.5;
  }
  
  virtual double perimeter() const{
    double ans = 0;
    for (size_t i = 0; i < points_.size() - 1; ++i) {
      ans += dist_points(points_[i], points_[i + 1]);
    }
    ans += dist_points(points_[0], points_[points_.size() - 1]);
    return ans;
  }
  
  void rotate(Point center, double angle) {
    for (size_t i = 0; i < points_.size(); ++i) {
       points_[i] = points_[i].rotate(center, angle * M_PI / 180);
    }
  }
  
  void reflex(Line axis) {
    for (size_t i = 0; i < points_.size(); ++i) {
      points_[i] = reflex_point(points_[i], axis);
    }
  }
  
  void reflex(Point center) {
    for (size_t i = 0; i < points_.size(); i++) {
      points_[i].reflex(center);
    }
  }
  
  void scale(Point center, double coefficient) {
    for (size_t i = 0; i < points_.size(); i++) {
      points_[i].scale(center, coefficient);
    }
  }
  
  virtual bool operator==(const Polygon& another) const {
    if (points_.size() != another.points_.size()) return 0;
    int idx_first_point = -1, idx_second_point = -1;
    for (size_t i = 0; i < points_.size(); ++i) {
      for (size_t j = 0; j < another.points_.size(); ++j) {
        if (points_[i] == another.points_[j]) {
          idx_first_point = int(i);
          idx_second_point = int(j);
          break;
        }
      }
    }
    if (idx_first_point < 0) return 0;
    std::vector<Point> tmp1_points_, tmp2_points_, tmp3_points_;
    for (size_t i = idx_first_point; i < points_.size(); ++i) {
      tmp1_points_.push_back(points_[i]);
    }
    for (size_t i = 0; int(i) < idx_first_point; ++i) {
      tmp1_points_.push_back(points_[i]);
    }
    for (size_t i = idx_second_point; i < another.points_.size(); ++i) {
      tmp2_points_.push_back(another.points_[i]);
    }
    for (size_t i = 0; int(i) < idx_second_point; ++i) {
      tmp2_points_.push_back(another.points_[i]);
    }
    for (size_t i = idx_second_point + 1; i > 0; --i) {
      tmp3_points_.push_back(another.points_[i - 1]);
    }
    for (size_t i = another.points_.size(); int(i) > idx_second_point + 1; --i) {
      tmp3_points_.push_back(another.points_[i - 1]);
    }
    if (tmp1_points_ == tmp2_points_ || tmp1_points_ == tmp3_points_) {
      return 1;
    }
    return 0;
  }
  
  bool containsPoint(Point point);
  
  
  virtual bool isCongruentTo(const Polygon& another) const{
    if (another.points_.size() != points_.size()) {
      return 0;
    }
    std::vector<double> len1, len2;
    size_t n = points_.size();
    for (size_t i = 0; i < n; ++i) {
      len1.push_back(Vector(points_[i], points_[(i + 1) % n]).len);
      len2.push_back(Vector(another.points_[i], another.points_[(i + 1) % n]).len);
    }
    for (size_t i = 0; i < len1.size(); ++i) {
      for (size_t j = 0; j < len2.size(); ++j) {
        if (abs(len1[i] - len2[j]) < 1e-4) {
          bool flag_is_congruent = 1;
          for (size_t k = 0; k < n; ++k) {
            if (abs(len1[(i + k) % n] - len2[(j + k) % n]) > 1e-4) {
              flag_is_congruent = 0;
              break;
            }
          }
          if(flag_is_congruent) return 1;
          flag_is_congruent = 1;
          for (size_t k = 0; k < n; ++k) {
            if (abs(len1[(i + k) % n] - len2[(j - k + n) % n]) > 1e-4) {
              flag_is_congruent = 0;
              break;
            }
          }
          if(flag_is_congruent) return 1;
        }
      }
    }
    return 0;
  }
  
  virtual bool isSimilarTo(const Polygon& another) const{
    if (another.points_.size() != points_.size()) {
      return 0;
    }
    std::vector<Vector> edges1, edges2;
    Vector edge1(points_[0], points_[1]);
    for (size_t i = 0; i < points_.size(); ++i) {
      edges1.push_back(Vector(points_[i], points_[(i + 1) % points_.size()]));
      edges2.push_back(Vector(another.points_[i], another.points_[(i + 1) % points_.size()]));
    }
    int n = int(edges1.size());
    for (size_t i = 0; i < edges2.size(); ++i) {
      int idx_first_point = 0;
      int idx_second_point = int(i);
      double coeff_similarity = edges1[idx_first_point].len / edges2[idx_second_point].len;
      int kol_correct_edges = 1;
      while (kol_correct_edges < n) {
        idx_first_point++;
        idx_second_point++;
        idx_first_point %= n;
        idx_second_point %= n;
        if (abs(coeff_similarity - edges1[idx_first_point].len / edges2[idx_second_point].len) < EPS) {
          kol_correct_edges++;
        }
        else {
          break;
        }
      }
      if (kol_correct_edges == n) {
        int kol_correct_angle = 0;
        while(kol_correct_angle < n) {
          double angle1 = scalar_product(edges1[idx_first_point], edges1[(idx_first_point + 1) % n]);
          double angle2 = scalar_product(edges1[idx_second_point], edges1[(idx_second_point + 1) % n]);
          if (abs(angle2 - angle1) < EPS) {
            kol_correct_angle++;
          }
          else {
            break;
          }
          idx_first_point++;
          idx_second_point++;
          idx_first_point %= n;
          idx_second_point %= n;
        }
        if (kol_correct_angle == n) {
          return 1;
        }
      }
      //--------------------------------
      idx_first_point = 0;
      idx_second_point = int(i);
      coeff_similarity = edges1[idx_first_point].len / edges2[idx_second_point].len;
      kol_correct_edges = 1;
      while (kol_correct_edges < n) {
        idx_first_point++;
        idx_second_point--;
        if (idx_second_point < 0) idx_second_point += n;
        idx_first_point %= n;
        idx_second_point %= n;
        if (abs(coeff_similarity - edges1[idx_first_point].len / edges2[idx_second_point].len) < EPS) {
          kol_correct_edges++;
        }
        else {
          break;
        }
      }
      if (kol_correct_edges == n) {
        int kol_correct_angle = 0;
        while(kol_correct_angle < n) {
          double angle1 = scalar_product(edges1[idx_first_point], edges1[(idx_first_point + 1) % n]);
          double angle2 = scalar_product(edges1[idx_second_point], edges1[(idx_second_point + 1) % n]);
          if (abs(angle2 - angle1) < EPS) {
            kol_correct_angle++;
          }
          else {
            break;
          }
          idx_first_point++;
          idx_second_point--;
          if (idx_second_point < 0) idx_second_point += n;
          idx_first_point %= n;
          idx_second_point %= n;
        }
        if (kol_correct_angle == n) {
          return 1;
        }
      }
    }
    return 0;
  }
  
  virtual ~Polygon(){
    points_.clear();
  };
};

class Ellipse: public Shape {
protected:
  double a, b;
  
  Point c1 = Point(0, 0);
  Point c2 = Point(0, 0);
  Point center_ = Point(0, 0);
public:
  Ellipse() = default;
  Ellipse(const Point& focus1, const Point& focus2, double dist) {

    Point cnt = Point((focus1.x + focus2.x) / 2, (focus1.y + focus2.y) / 2);
    Point f1_new = Point(focus1.x - cnt.x, focus1.y - cnt.y);
    Point f2_new = Point(focus2.x - cnt.x, focus2.y - cnt.y);
    c1 = focus1;
    c2 = focus2;
    center_.x = cnt.x;
    center_.y = cnt.y;
    Line focuses(f1_new, f2_new);
    a = dist / 2;
    b = sqrt(a * a - ((f1_new.x * f1_new.x) + (f1_new.y * f1_new.y)));
  }
  
  std::pair<Point,Point> focuses() {
    return std::make_pair(Point(center_.x - sqrt(a * a - b * b), center_.y), Point(center_.x + sqrt(a * a - b * b), center_.y));
  }
  
  std::pair<Line, Line> directrices() {
    return std::make_pair(Line(1, 0, center_.x -(a * a) / (sqrt(a * a - b * b))), Line(1, 0, center_.x + (a * a) / (sqrt(a * a - b * b))));
  }
  
  double eccentricity() const{
    return sqrt(1 - (b / a) * (b / a));
  }
  
  Point center() {
    return center_;
  }
  
  virtual double area() const{
    return M_PI * a * b;
  }
  
  virtual double perimeter() const{
    return M_PI * (3 * (a + b) - sqrt((3 * a + b) * (3 * b + a)));
  }
  
  virtual bool containsPoint(Point point) {
    if (dist_points(point, c1) + dist_points(point, c2) <= EPS + 2 * (a > b ? a: b)) {
      return 1;
    }
    return 0;
  }
  
  void rotate(Point center, double angle) {
    center_ = center_.rotate(center, angle * M_PI / 180);
    c1 = c1.rotate(center, angle * M_PI / 180);
    c2 = c2.rotate(center, angle * M_PI / 180);
  }
  
  void reflex(Point center) {
    center.reflex(center);
    c1.reflex(center);
    c2.reflex(center);
  }
  
  void reflex(Line axis) {
    center_ = reflex_point(center_, axis);
    c1 = reflex_point(c1, axis);
    c2 = reflex_point(c2, axis);
  }
  
  void scale(Point center, double coefficient) {
    center_.scale(center, coefficient);
    c1.scale(center, coefficient);
    c2.scale(center, coefficient);
    a *= coefficient;
    b *= coefficient;
  }
  
  virtual bool operator==(const Ellipse& another) const{
    if (center_ == another.center_ && a * a == another.a * another.a && b * b == another.b * another.b) {
      return 1;
    }
    return 0;
  }
  
  virtual bool isCongruentTo(const Ellipse& another) const{
    if (abs(a - another.a) < EPS && abs(b - another.b) < EPS) return 1;
    return 0;
  }
  
  virtual bool isSimilarTo(const Ellipse& another) const{
    if (abs(a / b - another.a / another.b) < EPS) {
      return 1;
    }
    return 0;
  }
};


class Circle: public Ellipse {
protected:
  double r;
  
public:
  Circle() = default;
  Circle(const Point& center_new, double radius): Ellipse(center_new, center_new, 2 * radius) {
    r = radius;
  }
  
  double radius() {
    return r;
  }
  
  Point get_center() {
    return center_;
  }
  
  
  virtual double perimeter() {
    return 2 * M_PI * r;
  }
  
  bool operator==(const Circle& another) const{
    if (center_ == another.center_ && r == another.r) {
      return 1;
    }
    return 0;
  }
  
  virtual bool isSimilarTo(const Circle& another) const{
    return 1 + 0 * another.r;
  }
  
  virtual bool isCongruentTo(const Circle& another) const{
    if (abs(r - another.r) < EPS) return 1;
    return 0;
  }
  
};


class Rectangle: public Polygon {
public:
  Rectangle() = default;
  Rectangle(const Point& a, const Point& b, double relation_near_side) {
    points_.resize(4);
    Point center_ = Point((a.x + b.x) / 2, (a.y + b.y) / 2);
    double tga = relation_near_side;
    if (relation_near_side < 1) relation_near_side = 1 / relation_near_side;
    double cosa = sqrt(1 / (1 + tga * tga));
    double cos2a = 2 * cosa * cosa - 1;
    points_[0] = a;
    points_[2] = b;
    Vector rotat(center_, b);
    
    points_[1] = Point(center_.x + rotat.rotate(acos(cos2a)).x, center_.y + rotat.rotate(acos(cos2a)).y);
    points_[3] = Point(center_.x - rotat.rotate(acos(cos2a)).x, center_.y - rotat.rotate(acos(cos2a)).y);
  }
  Point center() {
    return Point((points_[0].x + points_[2].x) / 2, (points_[0].y + points_[2].y) / 2);
  }
  std::pair<Line, Line> diagonals() {
    return std::make_pair(Line(points_[0], points_[2]), Line(points_[1], points_[3]));
  }
  
};

class Square: public Rectangle {
public:
  Square() = default;
  Square(const Point& a, const Point& b): Rectangle(a, b, 1){
    
  }
  Circle circumscribedCircle() {
    Point vertex1 = points_[0];
    Point vertex3 = points_[2];
    Point center_ = Point((vertex1.x + vertex3.x) / 2, (vertex1.y + vertex3.y) / 2);
    return Circle(center_, dist_points(vertex1, vertex3) / 2);
  }
  
  Circle inscribedCircle() {
    Point vertex1 = points_[0];
    Point vertex2 = points_[1];
    
    Point vertex3 = points_[2];
    Point center_ = Point((vertex1.x + vertex3.x) / 2, (vertex1.y + vertex3.y) / 2);
    return Circle(center_, dist_points(vertex1, vertex2) / 2);
  }
  
};

class Triangle: public Polygon {
public:
  Triangle() = default;
  explicit Triangle(Point a, Point b, Point c){
    points_ = {a, b, c};
  };
  Circle circumscribedCircle() {
    Point vertex1 = points_[0];
    Point vertex2 = points_[1];
    Point vertex3 = points_[2];
    
    Vector tmp(vertex2, vertex3);
    Vector serper = tmp.rotate(M_PI / 2);
    Point one_on_sp = Point((vertex2.x + vertex3.x) / 2, (vertex2.y + vertex3.y) / 2);
    Line serper1(one_on_sp, Point(one_on_sp.x + serper.x, one_on_sp.y + serper.y));
    tmp = Vector(vertex1, vertex2);
    serper = tmp.rotate(M_PI / 2);
    one_on_sp = Point((vertex1.x + vertex2.x) / 2, (vertex1.y + vertex2.y) / 2);
    Line serper2(one_on_sp, Point(one_on_sp.x + serper.x, one_on_sp.y + serper.y));
    Point center = line_intersection(serper1, serper2);

    return Circle(center, dist_points(center, vertex1));
  }
  
  Circle inscribedCircle() {
    Point vertex1 = points_[0];
    Point vertex2 = points_[1];
    Point vertex3 = points_[2];
    Vector st12(vertex1, vertex2);
    Vector st23(vertex2, vertex3);
    Vector st31(vertex3, vertex1);
    double radius = 2 * area() / (st12.len + st23.len + st31.len);
    Point incenter;
    incenter.x = (st12.len * vertex3.x + st23.len * vertex1.x + st31.len * vertex2.x) / (st12.len + st23.len + st31.len);
    incenter.y = (st12.len * vertex3.y + st23.len * vertex1.y + st31.len * vertex2.y) / (st12.len + st23.len + st31.len);
    
    return Circle(incenter, radius);
  }
  
  Point centroid() {
    Point vertex1 = points_[0];
    Point vertex2 = points_[1];
    Point vertex3 = points_[2];
    return line_intersection(Line(vertex1, Point((vertex2.x + vertex3.x) / 2, (vertex2.y + vertex3.y) / 2)), Line(vertex2, Point((vertex1.x + vertex3.x) / 2, (vertex1.y + vertex3.y) / 2)));
  }
  
  Point orthocenter() {
    Point vertex1 = points_[0];
    Point vertex2 = points_[1];
    Point vertex3 = points_[2];
    Line ab(vertex1, vertex2);
    Line bc(vertex2, vertex3);
    Line ca(vertex3, vertex1);
    Vector ab_v(vertex1, vertex2);
    Vector bc_v(vertex2, vertex3);
    Vector ca_v(vertex3, vertex1);
    
    Line h_a(bc_v.x, bc_v.y, -bc_v.x * vertex1.x - bc_v.y * vertex1.y);
    Line h_b(ca_v.x, ca_v.y, -ca_v.x * vertex2.x - ca_v.y * vertex2.y);
    return line_intersection(h_a, h_b);
  }
  
  Line EulerLine() {
    Circle a = circumscribedCircle();
    return Line(a.get_center(), centroid());
  }
  
  Circle ninePointsCircle() {
    Point center((circumscribedCircle().get_center().x + orthocenter().x) / 2, (circumscribedCircle().get_center().y + orthocenter().y) / 2);
    return Circle(center, circumscribedCircle().radius() / 2);
  }
  
};

bool Polygon::containsPoint(Point point) {
  double our_area = area();
  double triangulation_area = Triangle(point, points_[0], points_[points_.size() - 1]).area();
  for (size_t i = 0; i < points_.size() - 1; ++i) {
    triangulation_area += Triangle(point, points_[i], points_[i + 1]).area();
  }
  if (abs(triangulation_area - our_area) < EPS) return 1;
  return 0;
}

bool Shape::operator==(const Shape& another) const{
  if (dynamic_cast<const Polygon*>(&another) && dynamic_cast<const Polygon*>(this)) {
    const Polygon& first_polygon = *dynamic_cast<const Polygon*>(&another);
    const Polygon& second_polygon = *dynamic_cast<const Polygon*>(this);
    if (first_polygon == second_polygon) {
      return 1;
    }
  }
  if (dynamic_cast<const Circle*>(&another) && dynamic_cast<const Circle*>(this)) {
    const Circle& first_circle = *dynamic_cast<const Circle*>(&another);
    const Circle& second_circle = *dynamic_cast<const Circle*>(this);
    if (first_circle == second_circle) {
      return 1;
    }
  }
  if (dynamic_cast<const Ellipse*>(&another) && dynamic_cast<const Ellipse*>(this)) {
    const Ellipse& first_ellipse = *dynamic_cast<const Ellipse*>(&another);
    const Ellipse& second_ellipse = *dynamic_cast<const Ellipse*>(this);
    if (first_ellipse == second_ellipse) {
      return 1;
    }
  }
  return 0;
}

bool Shape::isSimilarTo(const Shape& another) {
  if (dynamic_cast<const Polygon*>(&another) && dynamic_cast<const Polygon*>(this)) {
    const Polygon& first_polygon = *dynamic_cast<const Polygon*>(&another);
    const Polygon& second_polygon = *dynamic_cast<const Polygon*>(this);
    if (first_polygon.isSimilarTo(second_polygon)) {
      return 1;
    }
  }
  if (dynamic_cast<const Circle*>(&another) && dynamic_cast<const Circle*>(this)) {
    const Circle& first_circle = *dynamic_cast<const Circle*>(&another);
    const Circle& second_circle = *dynamic_cast<const Circle*>(this);
    if (first_circle.isSimilarTo(second_circle)) {
      return 1;
    }
  }
  if (dynamic_cast<const Ellipse*>(&another) && dynamic_cast<const Ellipse*>(this)) {
    const Ellipse& first_ellipse = *dynamic_cast<const Ellipse*>(&another);
    const Ellipse& second_ellipse = *dynamic_cast<const Ellipse*>(this);
    if (first_ellipse.isSimilarTo(second_ellipse)) {
      return 1;
    }
  }
  return 0;
}

bool Shape::isCongruentTo(const Shape& another) {
  if (dynamic_cast<const Polygon*>(&another) && dynamic_cast<const Polygon*>(this)) {
    const Polygon& first_polygon = *dynamic_cast<const Polygon*>(&another);
    const Polygon& second_polygon = *dynamic_cast<const Polygon*>(this);
    if (first_polygon.isCongruentTo(second_polygon)) {
      return 1;
    }
  }
  if (dynamic_cast<const Circle*>(&another) && dynamic_cast<const Circle*>(this)) {
    const Circle& first_circle = *dynamic_cast<const Circle*>(&another);
    const Circle& second_circle = *dynamic_cast<const Circle*>(this);
    if (first_circle.isCongruentTo(second_circle)) {
      return 1;
    }
  }
  if (dynamic_cast<const Ellipse*>(&another) && dynamic_cast<const Ellipse*>(this)) {
    const Ellipse& first_ellipse = *dynamic_cast<const Ellipse*>(&another);
    const Ellipse& second_ellipse = *dynamic_cast<const Ellipse*>(this);
    if (first_ellipse.isCongruentTo(second_ellipse)) {
      return 1;
    }
  }
  return 0;
}

