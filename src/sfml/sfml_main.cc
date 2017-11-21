#include "ncode_common/src/common.h"
#include <SFML/Graphics.hpp>

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

struct LineStyle {
  LineStyle(const sf::Color& color, float thickness)
      : color(color), thickness(thickness) {}

  sf::Color color;
  float thickness;
};

sf::ConvexShape Line(const sf::Vector2f& direction, const LineStyle& style) {
  sf::Vector2f perpendicular =
      UnitVector(PerpendicularVector(direction)) * 0.5f * style.thickness;

  sf::ConvexShape line;
  line.setFillColor(style.color);
  line.setPointCount(4);
  line.setPoint(0, -perpendicular);
  line.setPoint(1, perpendicular);
  line.setPoint(2, direction + perpendicular);
  line.setPoint(3, direction - perpendicular);
  return line;
}

struct ArrowStyle {
  ArrowStyle(const LineStyle& line_style, float arrowhead_height)
      : line_style(line_style), arrowhead_height(arrowhead_height) {}

  LineStyle line_style;
  float arrowhead_height;
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

  float length_;
  sf::Vector2f direction_;

  sf::ConvexShape line_;
  sf::ConvexShape triangle_;
};

struct AxisStyle {
  AxisStyle(const ArrowStyle& arrow_style, const std::string& tick_font_file,
            const LineStyle& tick_lines_style, float tick_font_size,
            float tick_length, const std::string& label_font_file,
            float label_font_size, float label_offset)
      : arrow_style(arrow_style),
        tick_lines_style(tick_lines_style),
        tick_font_size(tick_font_size),
        tick_length(tick_length),
        label_font_size(label_font_size),
        label_offset(label_offset) {
    CHECK(tick_font.loadFromFile(tick_font_file));
    CHECK(label_font.loadFromFile(label_font_file));
  }

  // Style of axis line.
  ArrowStyle arrow_style;

  // Style of ticks.
  sf::Font tick_font;
  LineStyle tick_lines_style;
  float tick_font_size;
  float tick_length;

  // Style of label.
  sf::Font label_font;
  float label_font_size;
  float label_offset;
};

// A horizontal 1D axis.
class HAxis : public sf::Drawable, public sf::Transformable {
 public:
  HAxis(float length, const std::string& label,
        const std::vector<std::pair<float, std::string>>& ticks,
        const AxisStyle& style)
      : arrow_(sf::Vector2f(length, 0), length, style.arrow_style),
        label_(label, style.label_font, style.label_font_size) {
    for (const auto& location_and_label : ticks) {
      DrawTick(location_and_label.first, location_and_label.second, style);
    }

    float w_offset = label_.getLocalBounds().width / 2;
    float h_offset = label_.getLocalBounds().height / 2;
    label_.setPosition(length / 2 - w_offset, style.label_offset + h_offset);
    label_.setColor(style.tick_lines_style.color);
  }

  void DrawTick(float location, const std::string& label,
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

  // The arrow.
  Arrow arrow_;

  // The label.
  sf::Text label_;

  // Ticks.
  std::vector<sf::ConvexShape> ticks_;
  std::vector<sf::Text> tick_labels_;
};

// A vertical 1D axis.
class VAxis : public sf::Drawable, public sf::Transformable {
 public:
  VAxis(float length, const std::string& label,
        const std::vector<std::pair<float, std::string>>& ticks,
        const AxisStyle& style)
      : arrow_(sf::Vector2f(0, -length), length, style.arrow_style),
        label_(label, style.label_font, style.label_font_size) {
    for (const auto& location_and_label : ticks) {
      DrawTick(location_and_label.first, location_and_label.second, style);
    }

    float w_offset = label_.getLocalBounds().width / 2;
    float h_offset = label_.getLocalBounds().height / 2;
    float offset_total = style.label_offset + w_offset;
    label_.setPosition(-offset_total - w_offset - 10, -length / 2 + h_offset);
    label_.setColor(style.tick_lines_style.color);
  }

  void DrawTick(float location, const std::string& label,
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
  }

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

  // The arrow.
  Arrow arrow_;

  // The label.
  sf::Text label_;

  // Ticks.
  std::vector<sf::ConvexShape> ticks_;
  std::vector<sf::Text> tick_labels_;
};

// File to use for all fonts.
static constexpr char kDefaultFontSize[] = "DejaVuSans.ttf";

int main() {
  sf::RenderWindow window(sf::VideoMode(800, 600), "SFML works!");

  LineStyle axis_line_style(sf::Color::White, 5);
  LineStyle tick_line_style(sf::Color::White, 2);
  ArrowStyle axis_arrow_style(axis_line_style, 15);
  AxisStyle axis_style(axis_arrow_style, kDefaultFontSize, tick_line_style, 10,
                       10, kDefaultFontSize, 15, 15);

  HAxis haxis(100.0f, "some label", {{10.0f, "10"}, {50.0f, "2"}}, axis_style);
  VAxis vaxis(100.0f, "some labez", {{10.0f, "1"}, {50.0f, "200"}}, axis_style);
  haxis.move(200, 200);
  vaxis.move(200, 200);

  while (window.isOpen()) {
    sf::Event event;
    while (window.pollEvent(event)) {
      if (event.type == sf::Event::Closed) window.close();
    }

    window.clear();
    window.draw(haxis);
    window.draw(vaxis);
    window.display();
  }

  return 0;
}
