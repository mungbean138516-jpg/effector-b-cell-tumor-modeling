from __future__ import annotations

import math
import textwrap
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont, JpegImagePlugin


ROOT = Path("/Users/moyingli/USC/SP26/MACLEAN_LAB")
OUTDIR = ROOT / "presentation" / "long_term_ode_b6"
SLIDES_DIR = OUTDIR / "slides"
PDF_PATH = OUTDIR / "long_term_ode_b6_summary.pdf"

IMG_SUMMARY = ROOT / "4:15 after meeting modifications" / "step1_d2" / "b6_dense_step1_outputs_pretty" / "b6_dense_ode_summary_pretty.png"
IMG_ZOOM = ROOT / "4:15 after meeting modifications" / "step1_d2" / "b6_dense_step1_outputs_pretty" / "b6_dense_threshold_zoom_pretty.png"
IMG_TRAJ_THRESHOLD = ROOT / "4:15 after meeting modifications" / "step1_d2" / "b6_dense_step1_outputs_pretty" / "b6_threshold_trajectories_0_10000_pretty.png"


SLIDE_W = 1920
SLIDE_H = 1080
MARGIN_X = 84
TITLE_Y = 62
CONTENT_Y = 178
IMAGE_W = 1080
IMAGE_H = 740
RIGHT_X = 1210
RIGHT_W = 620
BG = "#F7F8FB"
INK = "#172033"
MUTED = "#56637A"
ACCENT = "#0E7490"
ACCENT_2 = "#D97706"
PANEL = "#FFFFFF"
PANEL_BORDER = "#D7DDE8"


def load_font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates = []
    if bold:
        candidates.extend(
            [
                "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
                "/System/Library/Fonts/Supplemental/Helvetica Bold.ttf",
                "/System/Library/Fonts/Supplemental/Trebuchet MS Bold.ttf",
            ]
        )
    else:
        candidates.extend(
            [
                "/System/Library/Fonts/Supplemental/Arial.ttf",
                "/System/Library/Fonts/Supplemental/Helvetica.ttf",
                "/System/Library/Fonts/Supplemental/Trebuchet MS.ttf",
            ]
        )
    for path in candidates:
        if Path(path).exists():
            return ImageFont.truetype(path, size=size)
    return ImageFont.load_default()


FONT_TITLE = load_font(40, bold=True)
FONT_SUBTITLE = load_font(23, bold=False)
FONT_BODY = load_font(28, bold=False)
FONT_BODY_BOLD = load_font(28, bold=True)
FONT_CAPTION = load_font(18, bold=False)
FONT_SMALL = load_font(20, bold=False)
FONT_COMPACT = load_font(22, bold=False)


def ensure_dirs() -> None:
    SLIDES_DIR.mkdir(parents=True, exist_ok=True)


def bullet_lines(text: str, width: int = 34) -> list[str]:
    wrapped = textwrap.wrap(text, width=width)
    if not wrapped:
        return []
    lines = [f"• {wrapped[0]}"]
    lines.extend(f"  {line}" for line in wrapped[1:])
    return lines


def draw_bullets(draw: ImageDraw.ImageDraw, bullets: list[str], x: int, y: int, width: int) -> None:
    line_gap = 16
    block_gap = 18
    cursor_y = y
    for bullet in bullets:
        lines = bullet_lines(bullet, width=max(24, int(width / 16)))
        for line in lines:
            draw.text((x, cursor_y), line, font=FONT_BODY, fill=INK)
            bbox = draw.textbbox((x, cursor_y), line, font=FONT_BODY)
            cursor_y = bbox[3] + line_gap
        cursor_y += block_gap


def draw_compact_bullets(draw: ImageDraw.ImageDraw, bullets: list[str], x: int, y: int, width: int) -> None:
    line_gap = 8
    block_gap = 8
    cursor_y = y
    for bullet in bullets:
        lines = bullet_lines(bullet, width=max(36, int(width / 17)))
        for line in lines:
            draw.text((x, cursor_y), line, font=FONT_COMPACT, fill=INK)
            bbox = draw.textbbox((x, cursor_y), line, font=FONT_COMPACT)
            cursor_y = bbox[3] + line_gap
        cursor_y += block_gap


def fit_image_obj(image: Image.Image, box_w: int, box_h: int) -> Image.Image:
    scale = min(box_w / image.width, box_h / image.height)
    new_size = (int(image.width * scale), int(image.height * scale))
    return image.resize(new_size, Image.Resampling.LANCZOS)


def fit_image(path: Path, box_w: int, box_h: int) -> Image.Image:
    image = Image.open(path).convert("RGB")
    return fit_image_obj(image, box_w, box_h)


def paste_image_obj(canvas: Image.Image, image: Image.Image, x: int, y: int, w: int, h: int) -> None:
    draw = ImageDraw.Draw(canvas)
    draw.rounded_rectangle(
        (x - 2, y - 2, x + w + 2, y + h + 2),
        radius=24,
        fill=PANEL,
        outline=PANEL_BORDER,
        width=2,
    )
    fitted = fit_image_obj(image, w - 22, h - 22)
    paste_x = x + (w - fitted.width) // 2
    paste_y = y + (h - fitted.height) // 2
    canvas.paste(fitted, (paste_x, paste_y))


def paste_image(canvas: Image.Image, img_path: Path, x: int, y: int, w: int, h: int) -> None:
    paste_image_obj(canvas, Image.open(img_path).convert("RGB"), x, y, w, h)


def crop_image(img_path: Path, crop_box: tuple[int, int, int, int]) -> Image.Image:
    return Image.open(img_path).convert("RGB").crop(crop_box)


def title_block(draw: ImageDraw.ImageDraw, title: str, subtitle: str | None = None) -> None:
    draw.text((MARGIN_X, TITLE_Y), title, font=FONT_TITLE, fill=INK)
    if subtitle:
        draw.text((MARGIN_X, TITLE_Y + 54), subtitle, font=FONT_SUBTITLE, fill=MUTED)
    draw.rounded_rectangle((MARGIN_X - 34, TITLE_Y - 20, MARGIN_X - 20, TITLE_Y + 30), radius=7, fill=ACCENT)


def footer(draw: ImageDraw.ImageDraw, text: str) -> None:
    draw.line((MARGIN_X, SLIDE_H - 62, SLIDE_W - MARGIN_X, SLIDE_H - 62), fill=PANEL_BORDER, width=2)
    draw.text((MARGIN_X, SLIDE_H - 46), text, font=FONT_CAPTION, fill=MUTED)


def new_canvas() -> tuple[Image.Image, ImageDraw.ImageDraw]:
    canvas = Image.new("RGB", (SLIDE_W, SLIDE_H), BG)
    draw = ImageDraw.Draw(canvas)
    return canvas, draw


def make_title_slide(path: Path) -> None:
    canvas, draw = new_canvas()
    draw.rounded_rectangle((MARGIN_X, 130, 900, 930), radius=28, fill=PANEL, outline=PANEL_BORDER, width=2)
    draw.rounded_rectangle((MARGIN_X, 154, MARGIN_X + 210, 194), radius=20, fill="#DFF3F7")
    draw.text((MARGIN_X + 20, 162), "MacLean Lab update", font=FONT_SMALL, fill=ACCENT)
    draw.text((MARGIN_X, 240), "Long-Term ODE Dense b6 Sweep", font=load_font(56, bold=True), fill=INK)
    draw.text(
        (MARGIN_X, 320),
        "Model fixed after the meeting; only b6 is varied and each run is extended to 10,000 days.",
        font=FONT_SUBTITLE,
        fill=MUTED,
    )
    bullets = [
        "Script in focus: b6_dense_step0_step1_pretty.jl",
        "Baseline setup: direct CTL, b5 = 1e-4, baseline initial condition",
        "Main question: when does stronger B-cell killing create durable control instead of only delayed escape?",
    ]
    draw_bullets(draw, bullets, MARGIN_X, 420, 740)
    paste_image(canvas, IMG_SUMMARY, 980, 150, 860, 760)
    footer(draw, "Long-term ODE summary from 1e-7 to 3e-6 in b6")
    canvas.save(path)


def make_result_slide(
    path: Path,
    title: str,
    subtitle: str,
    image_path: Path,
    bullets: list[str],
    footer_text: str,
) -> None:
    canvas, draw = new_canvas()
    title_block(draw, title, subtitle)
    paste_image(canvas, image_path, MARGIN_X, CONTENT_Y, IMAGE_W, IMAGE_H)
    draw.rounded_rectangle(
        (RIGHT_X, CONTENT_Y, RIGHT_X + RIGHT_W, CONTENT_Y + IMAGE_H),
        radius=28,
        fill=PANEL,
        outline=PANEL_BORDER,
        width=2,
    )
    draw.rounded_rectangle((RIGHT_X + 28, CONTENT_Y + 28, RIGHT_X + 270, CONTENT_Y + 66), radius=18, fill="#EFF6FF")
    draw.text((RIGHT_X + 46, CONTENT_Y + 34), "Key takeaways", font=FONT_SMALL, fill=ACCENT)
    draw_bullets(draw, bullets, RIGHT_X + 34, CONTENT_Y + 110, RIGHT_W - 70)
    footer(draw, footer_text)
    canvas.save(path)


def make_big_visual_slide(
    path: Path,
    title: str,
    subtitle: str,
    image_path: Path,
    bullets: list[str],
    footer_text: str,
) -> None:
    canvas, draw = new_canvas()
    title_block(draw, title, subtitle)
    paste_image(canvas, image_path, 68, 154, 1784, 688)

    draw.rounded_rectangle((68, 864, 1852, 1000), radius=24, fill=PANEL, outline=PANEL_BORDER, width=2)
    draw.rounded_rectangle((94, 889, 248, 925), radius=18, fill="#EFF6FF")
    draw.text((118, 897), "Main points", font=FONT_SMALL, fill=ACCENT)
    draw_compact_bullets(draw, bullets, 286, 886, 1526)

    footer(draw, footer_text)
    canvas.save(path)


def make_dual_crop_slide(
    path: Path,
    title: str,
    subtitle: str,
    left_path: Path,
    left_crop: tuple[int, int, int, int],
    right_path: Path,
    right_crop: tuple[int, int, int, int],
    bullets: list[str],
    footer_text: str,
) -> None:
    canvas, draw = new_canvas()
    title_block(draw, title, subtitle)
    paste_image_obj(canvas, crop_image(left_path, left_crop), 68, 156, 872, 666)
    paste_image_obj(canvas, crop_image(right_path, right_crop), 980, 156, 872, 666)

    draw.rounded_rectangle((68, 848, 1852, 1000), radius=24, fill=PANEL, outline=PANEL_BORDER, width=2)
    draw.rounded_rectangle((96, 874, 276, 912), radius=18, fill="#FEF3C7")
    draw.text((120, 882), "Takeaway", font=FONT_SMALL, fill=ACCENT_2)
    draw_compact_bullets(draw, bullets, 314, 875, 1490)

    footer(draw, footer_text)
    canvas.save(path)


def make_discussion_slide(
    path: Path,
    title: str,
    subtitle: str,
    bullets: list[str],
    side_title: str,
    side_lines: list[str],
    image_path: Path,
    footer_text: str,
) -> None:
    canvas, draw = new_canvas()
    title_block(draw, title, subtitle)

    draw.rounded_rectangle((MARGIN_X, CONTENT_Y, 1040, 920), radius=28, fill=PANEL, outline=PANEL_BORDER, width=2)
    draw.rounded_rectangle((MARGIN_X + 30, CONTENT_Y + 28, MARGIN_X + 285, CONTENT_Y + 68), radius=18, fill="#FEF3C7")
    draw.text((MARGIN_X + 50, CONTENT_Y + 35), side_title, font=FONT_SMALL, fill=ACCENT_2)
    draw_bullets(draw, bullets, MARGIN_X + 34, CONTENT_Y + 110, 880)

    draw.rounded_rectangle((1120, CONTENT_Y, 1836, 920), radius=28, fill=PANEL, outline=PANEL_BORDER, width=2)
    paste_image(canvas, image_path, 1148, CONTENT_Y + 34, 660, 430)
    box_y = CONTENT_Y + 500
    for idx, line in enumerate(side_lines):
        yy = box_y + idx * 56
        draw.rounded_rectangle((1150, yy, 1800, yy + 42), radius=16, fill="#F8FAFC")
        draw.text((1172, yy + 8), line, font=FONT_SMALL, fill=INK)

    footer(draw, footer_text)
    canvas.save(path)


def build_slides() -> list[Path]:
    slides = []

    slide1 = SLIDES_DIR / "slide_01_title.png"
    make_title_slide(slide1)
    slides.append(slide1)

    slide2 = SLIDES_DIR / "slide_02_early_burden.png"
    make_big_visual_slide(
        slide2,
        "Finding 1: larger b6 strongly suppresses early tumor growth",
        "Full 10,000-day sweep summary",
        IMG_SUMMARY,
        [
            "T100 drops from 1.28e4 at b6 = 1e-7 to 47 at b6 = 3e-6.",
            "Time to 1e5 shifts from 136 days to not reached by 10,000 days.",
            "AUC(0-365) falls from 1.95e9 to 1.63e4, a drop of more than five orders of magnitude.",
        ],
        "Early burden and threshold-crossing times both improve monotonically as b6 increases.",
    )
    slides.append(slide2)

    slide3 = SLIDES_DIR / "slide_03_transition_window.png"
    make_big_visual_slide(
        slide3,
        "Finding 2: long-time behavior changes in a narrow b6 window",
        "Threshold zoom around the switch region",
        IMG_ZOOM,
        [
            "Tumor-dominant solutions persist up to about 6.56e-7.",
            "An intermediate or dormancy-like regime appears from about 6.69e-7 to 1.05e-6.",
            "Clear tumor-controlled solutions begin around 1.15e-6 and continue through the top of the sweep.",
        ],
        "The control boundary is not broad; it is concentrated near b6 ~ 10^-6.",
    )
    slides.append(slide3)

    slide4 = SLIDES_DIR / "slide_04_day365_vs_longtime.png"
    make_dual_crop_slide(
        slide4,
        "Finding 3: day 365 can look good even when long-time escape still occurs",
        "Large-panel view of the threshold-region curves",
        IMG_ZOOM,
        (55, 55, 650, 455),
        IMG_ZOOM,
        (650, 55, 1250, 455),
        [
            "At b6 ≈ 6.05e-7, T365 = 5.76e3 but T10000 rebounds to 9.74e6.",
            "T365 < 1e5 starts near 5.54e-7; T10000 < 1e5 starts later, near 6.69e-7.",
            "Long runs separate delay from durable control.",
        ],
        "This is the main reason the long-term ODE view is more informative than the 365-day snapshot alone.",
    )
    slides.append(slide4)

    slide5 = SLIDES_DIR / "slide_05_immune_state.png"
    make_dual_crop_slide(
        slide5,
        "Finding 4: crossing the threshold reshapes the immune steady state",
        "Early-burden collapse and long-time immune-state expansion",
        IMG_SUMMARY,
        (55, 455, 650, 900),
        IMG_SUMMARY,
        (650, 455, 1255, 900),
        [
            "B at day 10,000 rises from about 205 at b6 = 1e-7 to 1.42e5 at b6 = 3e-6.",
            "NK at day 10,000 rises from about 716 to 2.28e5 across the same range.",
            "CTL peaks near b6 = 6.69e-7 at about 960, then falls once tumor becomes very small.",
        ],
        "The long-time control state is associated with a high-B and high-NK background, not simply with maximal CTL.",
    )
    slides.append(slide5)

    slide6 = SLIDES_DIR / "slide_06_questions.png"
    make_discussion_slide(
        slide6,
        "Questions to ask the leading PhD",
        "Suggested discussion prompts for the meeting",
        [
            "Which endpoint should define control for this project: T365, T1000, T10000, or a threshold-crossing rule?",
            "Should the intermediate regime be interpreted as biologically meaningful dormancy, or do we treat it as a model artifact until validated?",
            "Is the estimated control window around b6 ~ 10^-6 biologically plausible for the intended B-cell mechanism?",
            "Do we want the next long-time scan to focus on b5 × b6, since B and NK change sharply across the threshold?",
        ],
        "Quick numerical anchors",
        [
            "First T10000 < 1e5: b6 ≈ 6.69e-7",
            "First T10000 < 1e3: b6 ≈ 1.15e-6",
            "Transition window: ~6.7e-7 to ~1.15e-6",
            "CTL peak: ~960 near b6 = 6.69e-7",
        ],
        IMG_ZOOM,
        "These questions help turn the figure set into concrete next decisions.",
    )
    slides.append(slide6)

    slide7 = SLIDES_DIR / "slide_07_next_steps.png"
    make_discussion_slide(
        slide7,
        "Future directions",
        "Reasonable next steps after the long-term ODE result",
        [
            "Run a long-time b5 × b6 sweep to see how the control boundary moves when MDSC suppression of B cells changes.",
            "Add focused SDE only around 6.5e-7 to 1.2e-6 to test whether the ODE threshold survives stochasticity.",
            "Check sensitivity to initial tumor seed and to the definition of control versus delay.",
            "If the threshold remains sharp, use continuation or bifurcation analysis to explain the switch mechanistically.",
        ],
        "Short recommendation",
        [
            "Best immediate next step: long-time b5 × b6",
            "Best robustness check: threshold-window SDE",
            "Best interpretation check: redefine control criteria",
            "Best mechanism study: bifurcation near b6 ~ 10^-6",
        ],
        IMG_SUMMARY,
        "This creates a clean bridge from descriptive plots to the next modelling task.",
    )
    slides.append(slide7)

    return slides


def build_pdf(slides: list[Path]) -> None:
    images = [Image.open(path).convert("RGB") for path in slides]
    first, rest = images[0], images[1:]
    first.save(PDF_PATH, save_all=True, append_images=rest, resolution=150.0)


def main() -> None:
    ensure_dirs()
    slides = build_slides()
    build_pdf(slides)
    print(f"Generated {len(slides)} slide images in {SLIDES_DIR}")
    print(f"Generated PDF at {PDF_PATH}")


if __name__ == "__main__":
    main()
