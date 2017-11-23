#include <SFML/Graphics.hpp>
#include <algorithm>
#include <cmath>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <queue>

#include "ncode_common/src/common.h"
#include "ncode_common/src/logging.h"
#include "ncode_common/src/strutil.h"

template <typename T>
T DotProduct(const sf::Vector2<T>& lhs, const sf::Vector2<T>& rhs) {
  return lhs.x * rhs.x + lhs.y * rhs.y;
}

template <typename T>
T SquaredLength(const sf::Vector2<T>& vector) {
  return DotProduct(vector, vector);
}

template <typename T>
T Length(const sf::Vector2<T>& vector) {
  return std::sqrt(SquaredLength(vector));
}

template <typename T>
sf::Vector2<T> UnitVector(const sf::Vector2<T>& vector) {
  return vector / Length(vector);
}

template <typename T>
sf::Vector2<T> PerpendicularVector(const sf::Vector2<T>& vector) {
  return sf::Vector2<T>(-vector.y, vector.x);
}

static constexpr float PI = 3.141592653;
static constexpr float RadToDeg(float rad) { return 180 / PI * rad; }
static constexpr float DegToRad(float deg) { return PI / 180 * deg; }

class LineStyle {
 public:
  void set_color(sf::Color color) { color_ = color; }

  void set_thickness(float thickness) { thickness_ = thickness; }

  const sf::Color& color() const { return color_; }

  float thickness() const { return thickness_; }

 private:
  sf::Color color_ = sf::Color::White;
  float thickness_ = 5.0;
};

// Line from the origin to a point.
sf::ConvexShape Line(const sf::Vector2f& point, const LineStyle& style) {
  sf::Vector2f perpendicular =
      UnitVector(PerpendicularVector(point)) * 0.5f * style.thickness();

  sf::ConvexShape line;
  line.setFillColor(style.color());
  line.setPointCount(4);
  line.setPoint(0, -perpendicular);
  line.setPoint(1, perpendicular);
  line.setPoint(2, point + perpendicular);
  line.setPoint(3, point - perpendicular);
  return line;
}

class ArrowStyle {
 public:
  float arrowhead_height() const { return arrowhead_height_; }

  void set_arrowhead_height(float arrowhead_height) {
    arrowhead_height_ = arrowhead_height;
  }

  const LineStyle& line_style() const { return line_style_; }

  void set_line_style(LineStyle line_style) { line_style_ = line_style; }

 private:
  LineStyle line_style_;
  float arrowhead_height_ = 15.0;
};

// A 2D arrow.
class Arrow : public sf::Drawable, public sf::Transformable {
 public:
  Arrow(sf::Vector2f direction, float length, const ArrowStyle& style)
      : length_(length), direction_(direction) {
    Update(style);
  }

 private:
  void draw(sf::RenderTarget& target, sf::RenderStates states) const {
    states.transform *= getTransform();
    target.draw(line_, states);
    target.draw(triangle_, states);
  }

  void Update(const ArrowStyle& style) {
    float arrowhead_height = style.arrowhead_height();
    float thickness = style.line_style().thickness();

    // Update the line.
    sf::Vector2f arrow_direction =
        (length_ - arrowhead_height) * UnitVector(direction_);
    line_ = Line(arrow_direction, style.line_style());

    // Update the triangle
    const float end = std::max(length_, arrowhead_height);
    const float begin = end - arrowhead_height;

    triangle_.setFillColor(style.line_style().color());
    triangle_.setRotation(RadToDeg(std::atan2(direction_.y, direction_.x)));

    triangle_.setPointCount(3);
    triangle_.setPoint(0, sf::Vector2f(end, 0.f));
    triangle_.setPoint(1, sf::Vector2f(begin, 1.5f * thickness));
    triangle_.setPoint(2, sf::Vector2f(begin, -1.5f * thickness));
  }

  float length_;
  sf::Vector2f direction_;

  sf::ConvexShape line_;
  sf::ConvexShape triangle_;
};

// Styles an axis.
struct AxisStyle {
  const ArrowStyle& arrow_style() const { return arrow_style_; }

  void set_arrow_style(const ArrowStyle& arrow_style) {
    arrow_style_ = arrow_style;
  }

  const sf::Font& label_font() const { return label_font_; }

  void set_label_font(const sf::Font& label_font) { label_font_ = label_font; }

  float label_font_size() const { return label_font_size_; }

  void set_label_font_size(float label_font_size_Size) {
    label_font_size_ = label_font_size_Size;
  }

  float label_offset() const { return label_offset_; }

  void set_label_offset(float label_offset) { label_offset_ = label_offset; }

  const sf::Font& tick_font() const { return tick_font_; }

  void set_tick_font(const sf::Font& tick_font) { tick_font_ = tick_font; }

  float tick_font_size() const { return tick_font_size_; }

  void set_tick_font_size(float tick_font_size) {
    tick_font_size_ = tick_font_size;
  }

  float tick_length() const { return tick_length_; }

  void set_tick_length(float tick_length) { tick_length_ = tick_length; }

  const LineStyle& tick_lines_style() const { return tick_lines_style_; }

  void set_tick_lines_style(const LineStyle& tick_lines_style) {
    tick_lines_style_ = tick_lines_style;
  }

 private:
  // Style of axis line.
  ArrowStyle arrow_style_;

  // Style of ticks.
  sf::Font tick_font_;
  LineStyle tick_lines_style_;
  float tick_font_size_ = 10;
  float tick_length_ = 10;

  // Style of label.
  sf::Font label_font_;
  float label_font_size_ = 15;
  float label_offset_ = 15;
};

// Axes inherit from this class.
class AxisBase : public sf::Drawable, public sf::Transformable {
 protected:
  AxisBase(const Arrow& arrow, const std::string& label, const AxisStyle& style)
      : arrow_(arrow),
        label_(label, style.label_font(), style.label_font_size()) {}

  // The arrow.
  Arrow arrow_;

  // The label.
  sf::Text label_;

  // Ticks.
  std::vector<sf::ConvexShape> ticks_;
  std::vector<sf::Text> tick_labels_;

 private:
  void draw(sf::RenderTarget& target, sf::RenderStates states) const {
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
};

// A horizontal 1D axis.
class HAxis : public AxisBase {
 public:
  HAxis(float length, const std::string& label,
        const std::vector<std::pair<double, std::string>>& ticks,
        const AxisStyle& style)
      : AxisBase(Arrow(sf::Vector2f(length, 0), length, style.arrow_style()),
                 label, style) {
    for (const auto& fraction_and_label : ticks) {
      const std::string& tick_label = fraction_and_label.second;
      double fraction_of_range = fraction_and_label.first;
      float location = length * fraction_of_range;

      DrawTick(location, tick_label, style);
    }

    float w_offset = label_.getLocalBounds().width / 2;
    float h_offset = label_.getLocalBounds().height / 2;
    label_.setPosition(length / 2 - w_offset, style.label_offset() + h_offset);
    label_.setColor(style.tick_lines_style().color());
  }

  void DrawTick(float location, const std::string& label,
                const AxisStyle& style) {
    sf::ConvexShape tick_line =
        Line(sf::Vector2f(0, style.tick_length()), style.tick_lines_style());
    tick_line.setPosition(location, 0);
    ticks_.emplace_back(tick_line);

    sf::Text tick_text(label, style.tick_font(), style.tick_font_size());
    float w_offset = tick_text.getLocalBounds().width / 2;
    float h_offset = tick_text.getLocalBounds().height / 2;
    tick_text.setPosition(location - w_offset, style.tick_length() + h_offset);
    tick_text.setColor(style.tick_lines_style().color());
    tick_labels_.emplace_back(tick_text);
  }
};

// A vertical 1D axis.
class VAxis : public AxisBase {
 public:
  VAxis(float length, const std::string& label,
        const std::vector<std::pair<double, std::string>>& ticks,
        const AxisStyle& style)
      : AxisBase(Arrow(sf::Vector2f(0, -length), length, style.arrow_style()),
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
    label_.setPosition(-style.label_offset() - h_offset - max_tick_label_width,
                       -length / 2 - w_offset);
    label_.setColor(style.tick_lines_style().color());
  }

  float DrawTick(float location, const std::string& label,
                 const AxisStyle& style) {
    float tick_length = style.tick_length();
    sf::ConvexShape tick_line =
        Line(sf::Vector2f(tick_length, 0), style.tick_lines_style());
    tick_line.setPosition(-tick_length, -location);
    ticks_.emplace_back(tick_line);

    sf::Text tick_text(label, style.tick_font(), style.tick_font_size());
    float w_offset = tick_text.getLocalBounds().width;
    float h_offset = tick_text.getLocalBounds().height;
    tick_text.setPosition(-tick_length - w_offset - 10, -location - h_offset);
    tick_text.setColor(style.tick_lines_style().color());
    tick_labels_.emplace_back(tick_text);
    return tick_text.getLocalBounds().width;
  }
};

// Just a range oof numbers.
using Range = std::pair<double, double>;

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

class PlotStyle {
 public:
  const AxisStyle& x_axis_style() const { return x_axis_style_; }
  const AxisStyle& y_axis_style() const { return y_axis_style_; }

  void set_x_axis_style(const AxisStyle& axis_style) {
    x_axis_style_ = axis_style;
  }

  void set_y_axis_style(const AxisStyle& axis_style) {
    y_axis_style_ = axis_style;
  }

  AxisStyle* mutable_x_axis_style() { return &x_axis_style_; }
  AxisStyle* mutable_y_axis_style() { return &y_axis_style_; }

  // Helper function to set all fonts of all axes to some value.
  void SetFonts(const sf::Font& font) {
    x_axis_style_.set_label_font(font);
    x_axis_style_.set_tick_font(font);
    y_axis_style_.set_label_font(font);
    y_axis_style_.set_tick_font(font);
  }

 private:
  AxisStyle x_axis_style_;
  AxisStyle y_axis_style_;
};

struct PlotParams {
 public:
  const Range& range_y() const { return range_y_; }

  void set_range_y(const Range& range_y) { range_y_ = range_y; }

  double y_tick_step() const { return y_tick_step_; }

  void set_y_tick_step(double tick_step) { y_tick_step_ = tick_step; }

  const std::string& y_label() const { return y_label_; }

  void set_y_label(const std::string& label) { y_label_ = label; }

  const std::chrono::milliseconds& time_span() const { return time_span_; }

  void set_time_span(const std::chrono::milliseconds& time_span) {
    time_span_ = time_span;
  }

  const std::chrono::milliseconds& time_step() const { return time_step_; }

  void set_time_step(const std::chrono::milliseconds& time_step) {
    time_step_ = time_step;
  }

 private:
  std::chrono::milliseconds time_step_;
  std::chrono::milliseconds time_span_;

  Range range_y_;
  double y_tick_step_ = 0;
  std::string y_label_;
};

class TimePlot : public sf::Drawable, public sf::Transformable {
 public:
  TimePlot(float width, float height, const PlotParams& params,
           const PlotStyle& style)
      : width_(width),
        height_(height),
        params_(params),
        x_axis_(width, "time (ms)", GetTicks({0, params.time_span().count()}),
                style.x_axis_style()),
        y_axis_(height, params.y_label(),
                GetTicks(params.range_y(), params.y_tick_step()),
                style.y_axis_style()) {
    double time_span = params_.time_span().count();
    double time_step = params_.time_step().count();
    max_values_count_ = time_span / time_step;
  }

  void AddData(double value) {
    double time_span = params_.time_span().count();
    double time_step = params_.time_step().count();
    float step_size = width_ * time_step / time_span;

    CHECK(value >= 0);
    double max_y = params_.range_y().second;
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

    LineStyle line_style;
    sf::Color c = sf::Color::White;
    line_style.set_color(c);
    line_style.set_thickness(2);

    plot_ = GenerateTrianglesStrip(points, line_style);
  }

 private:
  // adapted from
  // https://github.com/SFML/SFML/wiki/Source%3A-line-with-thickness
  static sf::VertexArray GenerateTrianglesStrip(
      const std::vector<sf::Vector2f>& points, const LineStyle& line_style) {
    using namespace sf;

    float thickness = line_style.thickness();
    const sf::Color& color = line_style.color();

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

  float TransformY(double y) {
    double max_y, min_y;
    std::tie(min_y, max_y) = params_.range_y();
    double fraction = (y - min_y) / (max_y - min_y);
    return -height_ * fraction;
  }

  void draw(sf::RenderTarget& target, sf::RenderStates states) const {
    states.transform *= getTransform();

    target.draw(x_axis_, states);
    target.draw(y_axis_, states);
    target.draw(plot_, states);
  }

  // Plot w and h
  const float width_;
  const float height_;

  // The plot parameters.
  const PlotParams params_;

  HAxis x_axis_;
  VAxis y_axis_;

  // The current values.
  std::deque<float> values_;

  // Max number of values.
  size_t max_values_count_;

  // The actual plot.
  sf::VertexArray plot_;
};

std::unique_ptr<ctr::manual::NetworkContainer> BootstrapNetwork(
    std::map<ctr::AggregateId, std::unique_ptr<ctr::BinSequence>>&& sequences,
    const nc::net::GraphStorage& graph, nc::EventQueue* event_queue) {
  auto network_container =
      nc::make_unique<ctr::manual::NetworkContainer>(&graph, event_queue);
  ctr::MockSimDeviceFactory device_factory(enter_port, event_queue);
  for (auto& id_and_bin_sequence : sequences) {
    const ctr::AggregateId& id = id_and_bin_sequence.first;
    std::unique_ptr<ctr::BinSequence>& bin_sequence =
        id_and_bin_sequence.second;

    nc::htsim::MatchRuleKey key_for_aggregate =
        network_container->AddAggregate(id);

    // Will add the bin sequence to the source device.
    const std::string& src_device_id = graph.GetNode(id.src())->id();
    device_factory.AddBinSequence(src_device_id, key_for_aggregate,
                                  std::move(bin_sequence));
  }
}

// File to use for all fonts.
static constexpr char kDefaultFontFile[] = "DejaVuSans.ttf";

int main() {
  sf::ContextSettings settings;
  settings.antialiasingLevel = 4;

  sf::RenderWindow window(sf::VideoMode(800, 600), "SFML works!",
                          sf::Style::Default, settings);
  window.setFramerateLimit(60);

  ctr::controller::NetworkContainer container;

  sf::Font font;
  CHECK(font.loadFromFile(kDefaultFontFile));
  PlotStyle plot_style;
  plot_style.SetFonts(font);

  PlotParams plot_params;
  plot_params.set_time_span(std::chrono::milliseconds(10000));
  plot_params.set_time_step(std::chrono::milliseconds(100));
  plot_params.set_range_y({0, 5000});
  plot_params.set_y_label("some y label");

  TimePlot plot(300, 150, plot_params, plot_style);
  plot.move(200, 200);

  std::mt19937 rnd(1);
  std::uniform_real_distribution<> y_dist(0, 5000);

  sf::Clock clock;
  while (window.isOpen()) {
    sf::Event event;
    while (window.pollEvent(event)) {
      if (event.type == sf::Event::Closed) window.close();
    }

    plot.AddData(y_dist(rnd));
    sf::Time time = clock.restart();
    LOG(INFO) << "FPS" << 1 / time.asSeconds();

    window.clear();
    window.draw(plot);
    window.display();
  }

  return 0;
}
