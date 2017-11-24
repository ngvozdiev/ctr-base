#include "sfml.h"

#include <algorithm>
#include <cmath>
#include <tuple>

#include "ncode_common/src/common.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/strutil.h"

namespace ctr {
namespace sfml {

template <typename T>
T DotProduct(sf::Vector2<T> lhs, sf::Vector2<T> rhs) {
  return lhs.x * rhs.x + lhs.y * rhs.y;
}

template <typename T>
T SquaredLength(sf::Vector2<T> vector) {
  return DotProduct(vector, vector);
}

template <typename T>
T Length(sf::Vector2<T> vector) {
  return std::sqrt(SquaredLength(vector));
}

template <typename T>
sf::Vector2<T> UnitVector(sf::Vector2<T> vector) {
  return vector / Length(vector);
}

template <typename T>
sf::Vector2<T> PerpendicularVector(sf::Vector2<T> vector) {
  return sf::Vector2<T>(-vector.y, vector.x);
}

template <typename T>
T LinearDistance(sf::Vector2<T> from, sf::Vector2<T> to) {
  return std::sqrt((from.x - to.x) * (from.x - to.x) +
                   (from.y - to.y) * (from.y - to.y));
}

static constexpr float PI = 3.141592653;
static constexpr float RadToDeg(float rad) { return 180 / PI * rad; }

sf::ConvexShape Line(const sf::Vector2f& point, const LineStyle& style) {
  sf::Vector2f perpendicular =
      UnitVector(PerpendicularVector(point)) * 0.5f * style.thickness;

  sf::ConvexShape line;
  line.setFillColor(style.color);
  line.setPointCount(4);
  line.setPoint(0, -perpendicular);
  line.setPoint(1, perpendicular);
  line.setPoint(2, point + perpendicular);
  line.setPoint(3, point - perpendicular);
  return line;
}

Arrow::Arrow(sf::Vector2f direction, float length, const ArrowStyle& style)
    : length_(length), direction_(direction) {
  Update(style);
}

void Arrow::draw(sf::RenderTarget& target, sf::RenderStates states) const {
  states.transform *= getTransform();
  target.draw(line_, states);
  target.draw(triangle_, states);
}

void Arrow::Update(const ArrowStyle& style) {
  float arrowhead_height = style.arrowhead_height;
  float thickness = style.line_style.thickness;

  // Update the line.
  sf::Vector2f arrow_direction =
      (length_ - arrowhead_height) * UnitVector(direction_);
  line_ = Line(arrow_direction, style.line_style);

  // Update the triangle
  const float end = std::max(length_, arrowhead_height);
  const float begin = end - arrowhead_height;

  triangle_.setFillColor(style.line_style.color);
  triangle_.setRotation(RadToDeg(std::atan2(direction_.y, direction_.x)));

  triangle_.setPointCount(3);
  triangle_.setPoint(0, sf::Vector2f(end, 0.f));
  triangle_.setPoint(1, sf::Vector2f(begin, 1.5f * thickness));
  triangle_.setPoint(2, sf::Vector2f(begin, -1.5f * thickness));
}

void AxisBase::draw(sf::RenderTarget& target, sf::RenderStates states) const {
  states.transform *= getTransform();

  target.draw(arrow_, states);
  target.draw(label_, states);
  for (const auto& tick : ticks_) {
    target.draw(tick, states);
  }

  for (const auto& tick_label : tick_labels_) {
    target.draw(tick_label, states);
  }
}

HAxis::HAxis(float length, const std::string& label,
             const std::vector<std::pair<double, std::string>>& ticks,
             const AxisStyle& style)
    : AxisBase(Arrow(sf::Vector2f(length, 0), length, style.arrow_style), label,
               style) {
  for (const auto& fraction_and_label : ticks) {
    const std::string& tick_label = fraction_and_label.second;
    double fraction_of_range = fraction_and_label.first;
    float location = length * fraction_of_range;

    DrawTick(location, tick_label, style);
  }

  float w_offset = label_.getLocalBounds().width / 2;
  float h_offset = label_.getLocalBounds().height / 2;
  label_.setPosition(length / 2 - w_offset, style.label_offset + h_offset);
  label_.setColor(style.tick_lines_style.color);
}

void HAxis::DrawTick(float location, const std::string& label,
                     const AxisStyle& style) {
  sf::ConvexShape tick_line =
      Line(sf::Vector2f(0, style.tick_length), style.tick_lines_style);
  tick_line.setPosition(location, 0);
  ticks_.emplace_back(tick_line);

  sf::Text tick_text(label, style.tick_font, style.tick_font_size);
  float w_offset = tick_text.getLocalBounds().width / 2;
  float h_offset = tick_text.getLocalBounds().height / 2;
  tick_text.setPosition(location - w_offset, style.tick_length + h_offset);
  tick_text.setColor(style.tick_lines_style.color);
  tick_labels_.emplace_back(tick_text);
}

VAxis::VAxis(float length, const std::string& label,
             const std::vector<std::pair<double, std::string>>& ticks,
             const AxisStyle& style)
    : AxisBase(Arrow(sf::Vector2f(0, -length), length, style.arrow_style),
               label, style) {
  float max_tick_label_width = 0;
  for (const auto& fraction_and_label : ticks) {
    const std::string& tick_label = fraction_and_label.second;
    double fraction_of_range = fraction_and_label.first;
    float location = length * fraction_of_range;

    float tick_label_width = DrawTick(location, tick_label, style);
    max_tick_label_width = std::max(max_tick_label_width, tick_label_width);
  }

  float w_offset = label_.getLocalBounds().width / 2;
  float h_offset = label_.getLocalBounds().height / 2;
  label_.rotate(90);
  label_.setPosition(-style.label_offset - h_offset - max_tick_label_width,
                     -length / 2 - w_offset);
  label_.setColor(style.tick_lines_style.color);
}

float VAxis::DrawTick(float location, const std::string& label,
                      const AxisStyle& style) {
  float tick_length = style.tick_length;
  sf::ConvexShape tick_line =
      Line(sf::Vector2f(tick_length, 0), style.tick_lines_style);
  tick_line.setPosition(-tick_length, -location);
  ticks_.emplace_back(tick_line);

  sf::Text tick_text(label, style.tick_font, style.tick_font_size);
  float w_offset = tick_text.getLocalBounds().width;
  float h_offset = tick_text.getLocalBounds().height;
  tick_text.setPosition(-tick_length - w_offset - 10, -location - h_offset);
  tick_text.setColor(style.tick_lines_style.color);
  tick_labels_.emplace_back(tick_text);
  return tick_text.getLocalBounds().width;
}

static std::string FormatNumber(double num) {
  if (num >= 100000000) {
    return nc::StrCat(nc::ToStringMaxDecimals(num / 1000000.0, 1), "M");
  }

  if (num >= 1000000) {
    return nc::StrCat(nc::ToStringMaxDecimals(num / 1000000.0, 2), "M");
  }

  if (num >= 100000) {
    return nc::StrCat(nc::ToStringMaxDecimals(num / 1000.0, 1), "K");
  }

  if (num >= 1000) {
    return nc::StrCat(nc::ToStringMaxDecimals(num / 1000.0, 2), "K");
  }

  return nc::ToStringMaxDecimals(num, 2);
}

// Given a range and number of ticks will produce ticks along the range.
static std::vector<std::pair<double, std::string>> GetTicks(const Range& range,
                                                            double step = 0) {
  double range_min, range_max;
  std::tie(range_min, range_max) = range;
  CHECK(range_max > range_min);

  if (step == 0) {
    // Will try to auto-pick reasonable step.
    double range_len = range_max - range_min;
    step = range_len / 5;
  }

  CHECK(step > 0);
  std::vector<std::pair<double, std::string>> out;
  for (double value = range_min; value <= range_max; value += step) {
    double fraction_of_range = (value - range_min) / (range_max - range_min);
    out.emplace_back(fraction_of_range, FormatNumber(value));
  }

  return out;
}

TimePlot::TimePlot(float width, float height, const PlotParams& params,
                   const PlotStyle& style)
    : width_(width),
      height_(height),
      params_(params),
      x_axis_(width, "time (ms)", GetTicks({0, params.time_span.count()}),
              style.x_axis_style),
      y_axis_(height, params.y_label,
              GetTicks(params.range_y, params.y_tick_step),
              style.y_axis_style) {
  double time_span = params_.time_span.count();
  double time_step = params_.time_step.count();
  max_values_count_ = time_span / time_step;
}

void TimePlot::AddData(double value) {
  double time_span = params_.time_span.count();
  double time_step = params_.time_step.count();
  float step_size = width_ * time_step / time_span;

  CHECK(value >= 0);
  double max_y = params_.range_y.second;
  value = std::min(value, max_y);

  values_.push_front(TransformY(value));
  if (values_.size() > max_values_count_) {
    values_.pop_back();
  }

  std::vector<sf::Vector2f> points;
  for (size_t i = 0; i < values_.size(); ++i) {
    float x = i * step_size;
    float y = values_[i];
    points.emplace_back(x, y);
  }

  plot_ = GenerateTrianglesStrip(points, params_.plot_line_style);
}

sf::VertexArray TimePlot::GenerateTrianglesStrip(
    const std::vector<sf::Vector2f>& points, const LineStyle& line_style) {
  using namespace sf;

  float thickness = line_style.thickness;
  const sf::Color& color = line_style.color;

  float ht = thickness / 2.0f;
  auto array = VertexArray(PrimitiveType::TrianglesStrip);
  for (int i = 0; i < points.size(); ++i) {
    Vector2f v0 =
        (i == 0) ? points[i] + points[i] - points[i + 1] : points[i - 1];
    Vector2f v1 = points[i];
    Vector2f v2 = (i == points.size() - 1)
                      ? points[i] + points[i] - points[i - 1]
                      : points[i + 1];

    Vector2f v01 = UnitVector(v1 - v0);
    Vector2f v12 = UnitVector(v2 - v1);
    Vector2f d = PerpendicularVector(v01 + v12);
    Vector2f n0 = PerpendicularVector(v01) * ht;
    Vector2f n1 = PerpendicularVector(v12) * ht;
    float fdot = DotProduct(d, n0);

    // special treatment for angles smaller 2*(90°-arccos(0.87)) = 120°
    if (fdot < 0.87f * Length(d)) {
      array.append(Vertex(v1 - n0, color));
      array.append(Vertex(v1 + n0, color));
      if (DotProduct(n0, v12) > 0) {
        // turn right
        array.append(Vertex(v1 - UnitVector(d) * ht, color));
        array.append(Vertex(v1 + n1, color));
      } else {
        // turn left
        array.append(Vertex(v1 - n1, color));
        array.append(Vertex(v1 + UnitVector(d) * ht, color));
      }
      array.append(Vertex(v1 - n1, color));
      array.append(Vertex(v1 + n1, color));
    } else {
      d *= ht / fdot;
      array.append(Vertex(v1 - d, color));
      array.append(Vertex(v1 + d, color));
    }
  }
  return array;
}

float TimePlot::TransformY(double y) {
  double max_y, min_y;
  std::tie(min_y, max_y) = params_.range_y;
  double fraction = (y - min_y) / (max_y - min_y);
  return -height_ * fraction;
}

void TimePlot::draw(sf::RenderTarget& target, sf::RenderStates states) const {
  states.transform *= getTransform();

  target.draw(x_axis_, states);
  target.draw(y_axis_, states);
  target.draw(plot_, states);
}

DashedLine::DashedLine(float total_len, const DashedLineStyle& style) {
  size_t count = total_len / (style.dash_len + style.skip_len);
  float offset = 0;
  for (size_t i = 0; i < count; ++i) {
    sf::ConvexShape line =
        Line(sf::Vector2f(style.dash_len, 0), style.line_style);
    line.setPosition(sf::Vector2f(offset, 0));
    lines_.emplace_back(line);

    offset += style.dash_len + style.skip_len;
  }
}

void DashedLine::draw(sf::RenderTarget& target, sf::RenderStates states) const {
  states.transform *= getTransform();
  for (const auto& line : lines_) {
    target.draw(line, states);
  }
}

VisualLink::VisualLink(sf::Vector2f from, sf::Vector2f to, float height,
                       bool flip_v, std::chrono::milliseconds time_step,
                       const nc::net::GraphLinkBase& link,
                       const VisualLinkStyle& style)
    : width_(LinearDistance(from, to)),
      height_(height),
      link_(link),
      style_(style),
      time_step_(time_step),
      base_line_(Line(sf::Vector2f(width_, 0), style.base_line_style)),
      upper_line_(width_, style.top_line_style),
      flip_v_(flip_v) {
  using namespace std::chrono;

  // Will create as many deques as there are tags.
  for (const auto& tag_and_color : style.plot_colors) {
    values_[tag_and_color.first];
  }

  if (flip_v) {
    upper_line_.setPosition(0, height_);
    rect_.setSize(sf::Vector2f(width_, height_));
  } else {
    upper_line_.setPosition(0, -height_);
    rect_.setSize(sf::Vector2f(width_, -height_));
  }

  max_values_count_ =
      duration_cast<milliseconds>(link.delay()).count() / time_step.count();
  rect_.setFillColor(style_.background_color);

  basic_transform_.rotate(RadToDeg(atan2(to.y - from.y, to.x - from.x)), from);
  basic_transform_.translate(from);
}

void VisualLink::AddValues(
    const std::map<nc::htsim::PacketTag, double>& new_values) {
  double max_y = link_.bandwidth().Mbps();

  double value_cumulative = 0;
  for (auto& tag_and_deque : values_) {
    std::deque<float>& values = tag_and_deque.second;
    nc::htsim::PacketTag tag = tag_and_deque.first;

    double new_value = 0;
    if (nc::ContainsKey(new_values, tag)) {
      new_value = nc::FindOrDie(new_values, tag);
    }
    CHECK(new_value >= 0);

    value_cumulative += new_value;
    value_cumulative = std::min(value_cumulative, max_y);

    values.push_front(TransformY(value_cumulative));
    if (values.size() > max_values_count_) {
      values.pop_back();
    }
  }

  plot_.clear();
  const std::deque<float>* bottom = nullptr;
  for (const auto& tag_and_deque : values_) {
    const std::deque<float>& values = tag_and_deque.second;
    sf::Color color =
        nc::FindOrDieNoPrint(style_.plot_colors, tag_and_deque.first);
    sf::VertexArray strip = GenerateTrianglesStrip(color, values, bottom);
    plot_.emplace_back(strip);
    bottom = &values;
  }
}

sf::VertexArray VisualLink::GenerateTrianglesStrip(
    const sf::Color& color, const std::deque<float>& values,
    const std::deque<float>* bottom) {
  using namespace sf;
  if (bottom != nullptr) {
    CHECK(values.size() == bottom->size());
  }

  std::chrono::milliseconds time_span =
      std::chrono::duration_cast<std::chrono::milliseconds>(link_.delay());
  time_span = std::max(time_span, std::chrono::milliseconds(1));

  float step_size = width_ * time_step_.count() / time_span.count();
  VertexArray array(PrimitiveType::TrianglesStrip);
  if (values.size() < 2) {
    return array;
  }

  // The first add the initial vertex.
  float first_bottom = bottom == nullptr ? 0 : (*bottom)[0];
  array.append(Vertex(Vector2f(0, first_bottom), color));

  float offset = 0;
  for (size_t i = 0; i < values.size(); ++i) {
    float value = values[i];
    array.append(Vertex(Vector2f(offset, value), color));

    float point_y;
    if (i != values.size() - 1) {
      offset += step_size;
      point_y = bottom == nullptr ? 0 : (*bottom)[i + 1];
    } else {
      point_y = bottom == nullptr ? 0 : (*bottom)[i];
    }

    array.append(Vertex(Vector2f(offset, point_y), color));
  }

  return array;
}

float VisualLink::TransformY(double y) {
  double min_y = 0;
  double max_y = link_.bandwidth().Mbps();
  double fraction = (y - min_y) / (max_y - min_y);

  if (flip_v_) {
    return height_ * fraction;
  }
  return -height_ * fraction;
}

void VisualLink::draw(sf::RenderTarget& target, sf::RenderStates states) const {
  states.transform *= getTransform();
  states.transform *= basic_transform_;
  target.draw(rect_, states);
  for (const auto& strip : plot_) {
    target.draw(strip, states);
  }
  target.draw(upper_line_, states);
  target.draw(base_line_, states);
}

Node::Node(float radius, sf::Vector2f location, const NodeStyle& style)
    : circle_(radius) {
  circle_.setFillColor(style.fill);
  circle_.setOutlineColor(style.outline);
  circle_.setOutlineThickness(style.outline_thickness);

  circle_.setPosition(location.x - radius, location.y - radius);
}

}  // namespace sfml
}  // namespace ctr
