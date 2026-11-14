"""Plot the fire-return-interval results from flytrap_ipm.py as an SVG."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def load_points(path: Path) -> list[tuple[int, float]]:
    with path.open(encoding="utf-8") as handle:
        report = json.load(handle)
    values = report.get("projected_lambda_by_fire_return_interval")
    if not isinstance(values, dict) or not values:
        raise ValueError("JSON has no projected_lambda_by_fire_return_interval values")
    points = sorted((int(interval), float(value)) for interval, value in values.items())
    if any(interval < 1 or not math.isfinite(value) or value < 0 for interval, value in points):
        raise ValueError("Intervals must be positive and multipliers finite and nonnegative")
    if len({interval for interval, _ in points}) != len(points):
        raise ValueError("Duplicate fire-return intervals")
    return points


def nice_step(span: float, target_ticks: int = 5) -> float:
    rough = span / target_ticks
    magnitude = 10 ** math.floor(math.log10(rough))
    for multiplier in (1, 2, 2.5, 5, 10):
        step = multiplier * magnitude
        if step >= rough:
            return step
    return 10 * magnitude


def svg_chart(points: list[tuple[int, float]]) -> str:
    width, height = 920, 570
    left, right, top, bottom = 105, 55, 105, 100
    plot_w, plot_h = width - left - right, height - top - bottom
    xmin, xmax = points[0][0], points[-1][0]
    if xmin == xmax:
        xmin -= 1
        xmax += 1
    ymin = min(1.0, *(value for _, value in points))
    ymax = max(1.0, *(value for _, value in points))
    padding = max(0.025, 0.10 * (ymax - ymin))
    ymin = max(0.0, ymin - padding)
    ymax += padding
    if ymax - ymin < 1e-9:
        ymax = ymin + 0.1

    def x_position(value: float) -> float:
        return left + (value - xmin) / (xmax - xmin) * plot_w

    def y_position(value: float) -> float:
        return top + (ymax - value) / (ymax - ymin) * plot_h

    best_interval, best_value = max(points, key=lambda item: item[1])
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" role="img" '
        f'aria-label="Flytrap population growth by fire-return interval">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="105" y="44" font-family="Arial,sans-serif" font-size="27" '
        'font-weight="700" fill="#173a30">Venus flytrap population projection</text>',
        '<text x="105" y="72" font-family="Arial,sans-serif" font-size="15" '
        'fill="#52645d">Fixed observed mean climate · deterministic fire intervals</text>',
    ]
    step = nice_step(ymax - ymin)
    tick = math.ceil(ymin / step) * step
    while tick <= ymax + step * 1e-8:
        y = y_position(tick)
        parts.append(f'<line x1="{left}" y1="{y:.2f}" x2="{width-right}" y2="{y:.2f}" '
                     'stroke="#e2e9e5" stroke-width="1"/>')
        parts.append(f'<text x="{left-12}" y="{y+5:.2f}" text-anchor="end" '
                     'font-family="Arial,sans-serif" font-size="13" fill="#52645d">'
                     f'{tick:.2f}</text>')
        tick += step
    y_unity = y_position(1.0)
    parts.append(f'<line x1="{left}" y1="{y_unity:.2f}" x2="{width-right}" '
                 f'y2="{y_unity:.2f}" stroke="#a35728" stroke-width="2" '
                 'stroke-dasharray="7 6"/>')
    parts.append(f'<text x="{width-right-5}" y="{y_unity-8:.2f}" text-anchor="end" '
                 'font-family="Arial,sans-serif" font-size="12" fill="#a35728">'
                 'stable population (λ = 1)</text>')
    for interval, _ in points:
        x = x_position(interval)
        parts.append(f'<line x1="{x:.2f}" y1="{height-bottom}" x2="{x:.2f}" '
                     f'y2="{height-bottom+6}" stroke="#697b72"/>')
        parts.append(f'<text x="{x:.2f}" y="{height-bottom+26}" text-anchor="middle" '
                     f'font-family="Arial,sans-serif" font-size="13" fill="#3f5048">'
                     f'{interval}</text>')
    parts.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}" '
                 'stroke="#697b72" stroke-width="1.5"/>')
    parts.append(f'<line x1="{left}" y1="{height-bottom}" x2="{width-right}" '
                 f'y2="{height-bottom}" stroke="#697b72" stroke-width="1.5"/>')
    path = " ".join(f'{"M" if i == 0 else "L"} {x_position(interval):.2f} '
                    f'{y_position(value):.2f}' for i, (interval, value) in enumerate(points))
    parts.append(f'<path d="{path}" fill="none" stroke="#157b62" stroke-width="3" '
                 'stroke-linecap="round" stroke-linejoin="round"/>')
    for interval, value in points:
        x, y = x_position(interval), y_position(value)
        parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="5" fill="#157b62" '
                     'stroke="white" stroke-width="2"/>')
    bx, by = x_position(best_interval), y_position(best_value)
    label_x = min(bx + 15, width - right - 160)
    label_y = max(top + 18, by - 18)
    parts.append(f'<circle cx="{bx:.2f}" cy="{by:.2f}" r="9" fill="none" '
                 'stroke="#a35728" stroke-width="2"/>')
    parts.append(f'<text x="{label_x:.2f}" y="{label_y:.2f}" '
                 'font-family="Arial,sans-serif" font-size="14" font-weight="700" '
                 f'fill="#173a30">Highest tested: {best_interval} years '
                 f'(λ = {best_value:.3f})</text>')
    parts.append(f'<text x="{left+plot_w/2:.2f}" y="{height-28}" text-anchor="middle" '
                 'font-family="Arial,sans-serif" font-size="16" fill="#273d33">'
                 'Fire-return interval (years)</text>')
    parts.append(f'<text transform="translate(30,{top+plot_h/2:.2f}) rotate(-90)" '
                 'text-anchor="middle" font-family="Arial,sans-serif" font-size="16" '
                 'fill="#273d33">Annual population multiplier (λ)</text>')
    parts.append('</svg>')
    return "\n".join(parts) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="JSON report from flytrap_ipm.py")
    parser.add_argument("--output", type=Path, help="SVG filename (default: input stem + .svg)")
    args = parser.parse_args()
    output = args.output or args.input.with_suffix(".svg")
    points = load_points(args.input)
    output.write_text(svg_chart(points), encoding="utf-8")
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
