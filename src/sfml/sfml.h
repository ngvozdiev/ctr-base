#ifndef CTR_SFML_H
#define CTR_SFML_H

#include <stddef.h>
#include <SFML/Graphics.hpp>
#include <chrono>
#include <deque>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "ncode_common/src/htsim/packet.h"
#include "ncode_common/src/net/net_common.h"

namespace ctr {
namespace sfml {

// Styles a line.
struct LineStyle {
  sf::Color color = sf::Color::White;
  float thickness = 5.0;
};

// Draws a line from (0,0) to a given point.
sf::ConvexShape Line(const sf::Vector2f& point, const LineStyle& style);

// Styles a dashed line.
struct DashedLineStyle {
  LineStyle line_style;
  float dash_len = 20;
  float skip_len = 10;
};

// Dashed line from (0,0) to (len, 0).
class DashedLine : public sf::Drawable, public sf::Transformable {
 public:
  DashedLine(float total_len, const DashedLineStyle& style);

 private:
  void draw(sf::RenderTarget& target, sf::RenderStates states) const;

  std::vector<sf::ConvexShape> lines_;
};

// Styles an arrow.
struct ArrowStyle {
  LineStyle line_style;
  float arrowhead_height = 15.0;
};

// A combination of a line and an arrow at the end. Drawn from (0,0) to point at
// a given direction (point) and be of a given length.
class Arrow : public sf::Drawable, public sf::Transformable {
 public:
  Arrow(sf::Vector2f direction, float length, const ArrowStyle& style);

 private:
  void draw(sf::RenderTarget& target, sf::RenderStates states) const;

  // Updates the arrow. Called once upon construction.
  void Update(const ArrowStyle& style);

  float length_;
  sf::Vector2f direction_;

  sf::ConvexShape line_;
  sf::ConvexShape triangle_;
};

// Styles an axis.
struct AxisStyle {
  // Style of axis line.
  ArrowStyle arrow_style;

  // Style of ticks.
  sf::Font tick_font;
  LineStyle tick_lines_style;
  float tick_font_size = 10;
  float tick_length = 10;

  // Style of label.
  sf::Font label_font;
  float label_font_size = 15;
  float label_offset = 15;
};

// Axes inherit from this class.
class AxisBase : public sf::Drawable, public sf::Transformable {
 protected:
  AxisBase(const Arrow& arrow, const std::string& label, const AxisStyle& style)
      : arrow_(arrow), label_(label, style.label_font, style.label_font_size) {}

  // The arrow.
  Arrow arrow_;

  // The label.
  sf::Text label_;

  // Ticks.
  std::vector<sf::ConvexShape> ticks_;
  std::vector<sf::Text> tick_labels_;

 private:
  void draw(sf::RenderTarget& target, sf::RenderStates states) const;
};

// A horizontal 1D axis.
class HAxis : public AxisBase {
 public:
  HAxis(float length, const std::string& label,
        const std::vector<std::pair<double, std::string>>& ticks,
        const AxisStyle& style);

  void DrawTick(float location, const std::string& label,
                const AxisStyle& style);
};

// A vertical 1D axis.
class VAxis : public AxisBase {
 public:
  VAxis(float length, const std::string& label,
        const std::vector<std::pair<double, std::string>>& ticks,
        const AxisStyle& style);

  float DrawTick(float location, const std::string& label,
                 const AxisStyle& style);
};

// Styles a plot.
struct PlotStyle {
  // Helper function to set all fonts of all axes to some value.
  void SetFonts(const sf::Font& font) {
    x_axis_style.label_font = font;
    x_axis_style.tick_font = font;
    y_axis_style.label_font = font;
    y_axis_style.tick_font = font;
  }

  AxisStyle x_axis_style;
  AxisStyle y_axis_style;
};

// Just a range of numbers.
using Range = std::pair<double, double>;

struct PlotParams {
  std::chrono::milliseconds time_step;
  std::chrono::milliseconds time_span;
  LineStyle plot_line_style;

  Range range_y;
  double y_tick_step = 0;
  std::string y_label;
};

// A plot whose x axis is always time in milliseconds, and data is assumed to be
// added at regular intervals.
class TimePlot : public sf::Drawable, public sf::Transformable {
 public:
  TimePlot(float width, float height, const PlotParams& params,
           const PlotStyle& style);

  void AddData(double value);

  const PlotParams& params() const { return params_; }

 private:
  // adapted from
  // https://github.com/SFML/SFML/wiki/Source%3A-line-with-thickness
  static sf::VertexArray GenerateTrianglesStrip(
      const std::vector<sf::Vector2f>& points, const LineStyle& line_style);

  float TransformY(double y);

  void draw(sf::RenderTarget& target, sf::RenderStates states) const;

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

  // Plot of the traffic level on the link.
  sf::VertexArray plot_;
};

struct VisualLinkStyle {
  VisualLinkStyle() {
    background_color = sf::Color::White;
    background_color.a = 128;
  }

  sf::Color background_color;
  std::map<nc::htsim::PacketTag, sf::Color> plot_colors;
  LineStyle base_line_style;
  DashedLineStyle top_line_style;
};

class VisualLink : public sf::Drawable, public sf::Transformable {
 public:
  VisualLink(sf::Vector2f from, sf::Vector2f to, float height, bool flip_v,
             std::chrono::milliseconds time_step,
             const nc::net::GraphLinkBase& link, const VisualLinkStyle& style);

  void AddValues(const std::map<nc::htsim::PacketTag, double>& new_values);

  std::chrono::milliseconds time_step() const { return time_step_; }

 private:
  sf::VertexArray GenerateTrianglesStrip(const sf::Color& color,
                                         const std::deque<float>& values,
                                         const std::deque<float>* bottom);

  float TransformY(double y);

  void draw(sf::RenderTarget& target, sf::RenderStates states) const;

  // Width and heights of the visual representation of this link.
  const float width_;
  const float height_;

  // The link.
  const nc::net::GraphLinkBase link_;

  // The style.
  const VisualLinkStyle style_;

  // Time step.
  std::chrono::milliseconds time_step_;

  // The current values.
  std::map<nc::htsim::PacketTag, std::deque<float>> values_;

  // Max number of values.
  size_t max_values_count_;

  // Plot of the link's traffic level.
  std::vector<sf::VertexArray> plot_;

  // Outline in the background.
  sf::RectangleShape rect_;

  // Line at the base.
  sf::ConvexShape base_line_;

  // Line on top.
  DashedLine upper_line_;

  // Basic transform that is always applied.
  sf::Transform basic_transform_;

  // Flips the representation vertically.
  bool flip_v_;
};

struct NodeStyle {
  sf::Color outline = sf::Color::Yellow;
  sf::Color fill = sf::Color::White;
  float outline_thickness = 2;
};

// A node in the graph.
class Node : public sf::Drawable, public sf::Transformable {
 public:
  Node(float radius, sf::Vector2f location, const NodeStyle& style);

 private:
  void draw(sf::RenderTarget& target, sf::RenderStates states) const {
    states.transform *= getTransform();
    target.draw(circle_, states);
  }

  sf::CircleShape circle_;
};

struct CircleGaugeStyle {
  sf::Font font;
  size_t label_font_size = 30;
  size_t limits_font_size = 20;
  size_t annotation_font_size = 25;
};

// Interface for things that can change opacity.
class OpacityVariable {
 public:
  virtual ~OpacityVariable() {}

  OpacityVariable() : opacity_(1.0) {}

  // Changes the opacity. Value should be in range 0 (transparent) - 1 (opaque).
  virtual void SetOpacity(double opacity) {
    if (opacity == opacity_) {
      return;
    }

    opacity_ = opacity;
    ChangeOpacity(opacity);
  }

  double opacity() const { return opacity_; }

 protected:
  virtual void ChangeOpacity(double opacity) = 0;

  double opacity_;
};

class CircleGauge : public sf::Drawable,
                    public sf::Transformable,
                    public OpacityVariable {
 public:
  CircleGauge(const std::string& annotation, const std::string& units,
              size_t min, size_t max, const CircleGaugeStyle& style);

  void Update(size_t value);

 protected:
  void ChangeOpacity(double opacity) override;

 private:
  static constexpr size_t kNumPoints = 100;
  static constexpr char kTexture[] = "gauge.png";
  static constexpr float kGaugeBezel = 44;

  void RegenerateWithOpacity(double opacity);

  void draw(sf::RenderTarget& target, sf::RenderStates states) const override {
    states.transform *= getTransform();
    target.draw(triangle_fan_, states);
    target.draw(sprite_, states);
    target.draw(min_label_, states);
    target.draw(max_label_, states);
    target.draw(label_, states);
    target.draw(annotation_label_, states);
  }

  const size_t min_;
  const size_t max_;
  const CircleGaugeStyle style_;
  const std::string annotation_;
  const std::string units_;
  sf::Text min_label_;
  sf::Text max_label_;
  sf::Text label_;
  sf::Text annotation_label_;
  sf::Texture texture_;
  sf::VertexArray triangle_fan_;
  sf::Sprite sprite_;
};

}  // namespace sfml
}  // namespace ctr
#endif
