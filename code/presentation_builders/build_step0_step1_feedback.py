from __future__ import annotations

import textwrap
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont, JpegImagePlugin


ROOT = Path("/Users/moyingli/USC/SP26/MACLEAN_LAB")
OUTDIR = ROOT / "presentation" / "step0_step1_feedback"
SLIDES_DIR = OUTDIR / "slides"
PDF_PATH = OUTDIR / "step0_step1_feedback_summary.pdf"

IMG_CTL_SUMMARY = ROOT / "B6_CTL_compare" / "ctl_form_compare_near_b6_summary.png"
IMG_TRAJ_365 = ROOT / "4:15 after meeting modifications" / "step1_d2" / "b6_dense_step1_outputs_pretty" / "b6_dense_trajectories_0_365_pretty.png"
IMG_TRAJ_10000 = ROOT / "4:15 after meeting modifications" / "step1_d2" / "b6_dense_step1_outputs_pretty" / "b6_dense_trajectories_0_10000_pretty.png"
IMG_TRAJ_THRESH = ROOT / "4:15 after meeting modifications" / "step1_d2" / "b6_dense_step1_outputs_pretty" / "b6_threshold_trajectories_0_10000_pretty.png"
IMG_SUMMARY = ROOT / "4:15 after meeting modifications" / "step1_d2" / "b6_dense_step1_outputs_pretty" / "b6_dense_ode_summary_pretty.png"
IMG_ZOOM = ROOT / "4:15 after meeting modifications" / "step1_d2" / "b6_dense_step1_outputs_pretty" / "b6_dense_threshold_zoom_pretty.png"

SLIDE_W = 1920
SLIDE_H = 1080
MARGIN_X = 58
TITLE_Y = 54
CONTENT_Y = 148
BG = "#F7F8FB"
INK = "#1A2437"
MUTED = "#5B6B82"
ACCENT = "#0F766E"
ACCENT_2 = "#CA8A04"
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
FONT_SUBTITLE = load_font(24, bold=False)
FONT_BODY = load_font(24, bold=False)
FONT_SMALL = load_font(20, bold=False)
FONT_CAPTION = load_font(18, bold=False)


def ensure_dirs() -> None:
    SLIDES_DIR.mkdir(parents=True, exist_ok=True)


def new_canvas() -> tuple[Image.Image, ImageDraw.ImageDraw]:
    canvas = Image.new("RGB", (SLIDE_W, SLIDE_H), BG)
    draw = ImageDraw.Draw(canvas)
    return canvas, draw


def title_block(draw: ImageDraw.ImageDraw, title: str, subtitle: str | None = None) -> None:
    draw.rounded_rectangle((MARGIN_X - 28, TITLE_Y - 18, MARGIN_X - 14, TITLE_Y + 32), radius=7, fill=ACCENT)
    draw.text((MARGIN_X, TITLE_Y), title, font=FONT_TITLE, fill=INK)
    if subtitle:
        draw.text((MARGIN_X, TITLE_Y + 54), subtitle, font=FONT_SUBTITLE, fill=MUTED)


def footer(draw: ImageDraw.ImageDraw, text: str) -> None:
    draw.line((MARGIN_X, SLIDE_H - 62, SLIDE_W - MARGIN_X, SLIDE_H - 62), fill=PANEL_BORDER, width=2)
    draw.text((MARGIN_X, SLIDE_H - 46), text, font=FONT_CAPTION, fill=MUTED)


def bullet_lines(text: str, width: int) -> list[str]:
    wrapped = textwrap.wrap(text, width=width)
    if not wrapped:
        return []
    lines = [f"• {wrapped[0]}"]
    lines.extend(f"  {line}" for line in wrapped[1:])
    return lines


def draw_bullets(draw: ImageDraw.ImageDraw, bullets: list[str], x: int, y: int, width_px: int) -> None:
    line_gap = 8
    block_gap = 8
    cursor = y
    wrap_width = max(34, int(width_px / 15))
    for bullet in bullets:
        for line in bullet_lines(bullet, wrap_width):
            draw.text((x, cursor), line, font=FONT_BODY, fill=INK)
            bbox = draw.textbbox((x, cursor), line, font=FONT_BODY)
            cursor = bbox[3] + line_gap
        cursor += block_gap


def paste_image_obj(canvas: Image.Image, image: Image.Image, x: int, y: int, w: int, h: int) -> None:
    draw = ImageDraw.Draw(canvas)
    draw.rounded_rectangle((x - 2, y - 2, x + w + 2, y + h + 2), radius=24, fill=PANEL, outline=PANEL_BORDER, width=2)
    scale = min((w - 18) / image.width, (h - 18) / image.height)
    new_size = (int(image.width * scale), int(image.height * scale))
    fitted = image.resize(new_size, Image.Resampling.LANCZOS)
    px = x + (w - fitted.width) // 2
    py = y + (h - fitted.height) // 2
    canvas.paste(fitted, (px, py))


def crop(path: Path, box: tuple[int, int, int, int]) -> Image.Image:
    return Image.open(path).convert("RGB").crop(box)


def paste_image(canvas: Image.Image, path: Path, x: int, y: int, w: int, h: int) -> None:
    paste_image_obj(canvas, Image.open(path).convert("RGB"), x, y, w, h)


def bottom_strip(draw: ImageDraw.ImageDraw, label: str, bullets: list[str]) -> None:
    y0, y1 = 874, 998
    draw.rounded_rectangle((MARGIN_X, y0, SLIDE_W - MARGIN_X, y1), radius=24, fill=PANEL, outline=PANEL_BORDER, width=2)
    draw.rounded_rectangle((MARGIN_X + 18, y0 + 24, MARGIN_X + 170, y0 + 60), radius=18, fill="#FEF3C7")
    draw.text((MARGIN_X + 38, y0 + 31), label, font=FONT_SMALL, fill=ACCENT_2)
    draw_bullets(draw, bullets, MARGIN_X + 210, y0 + 18, SLIDE_W - MARGIN_X - 240)


def make_title_slide(path: Path) -> None:
    canvas, draw = new_canvas()
    draw.rounded_rectangle((MARGIN_X, 146, 896, 926), radius=28, fill=PANEL, outline=PANEL_BORDER, width=2)
    draw.rounded_rectangle((MARGIN_X + 22, 170, MARGIN_X + 270, 210), radius=20, fill="#DFF3F7")
    draw.text((MARGIN_X + 40, 178), "Step 0 + Step 1 only", font=FONT_SMALL, fill=ACCENT)
    draw.text((MARGIN_X + 24, 258), "Post-meeting feedback deck", font=load_font(54, bold=True), fill=INK)
    intro = (
        "This version follows Adam's request more literally:\n"
        "freeze the structure first, then make the long-term dense b6 ODE\n"
        "result clear through trajectories and summary plots."
    )
    draw.multiline_text((MARGIN_X + 24, 340), intro, font=FONT_SUBTITLE, fill=MUTED, spacing=8)
    draw_bullets(
        draw,
        [
            "Step 0: fix the simpler direct CTL formulation and stop changing model structure.",
            "Step 1: long-time dense 1D b6 ODE sweep with b5 fixed at 1e-4 and baseline initial condition.",
            "Presentation emphasis: trajectory path first, endpoint plots second.",
        ],
        MARGIN_X + 24,
        460,
        760,
    )
    paste_image(canvas, IMG_ZOOM, 960, 150, 860, 760)
    footer(draw, "Trajectory-first presentation for Step 0 and Step 1")
    canvas.save(path)


def make_full_image_slide(path: Path, title: str, subtitle: str, image_path: Path, bullets: list[str], footer_text: str) -> None:
    canvas, draw = new_canvas()
    title_block(draw, title, subtitle)
    paste_image(canvas, image_path, MARGIN_X, CONTENT_Y, 1804, 694)
    bottom_strip(draw, "Takeaway", bullets)
    footer(draw, footer_text)
    canvas.save(path)


def make_full_crop_slide(path: Path, title: str, subtitle: str, image: Image.Image, bullets: list[str], footer_text: str) -> None:
    canvas, draw = new_canvas()
    title_block(draw, title, subtitle)
    paste_image_obj(canvas, image, MARGIN_X, CONTENT_Y, 1804, 694)
    bottom_strip(draw, "Takeaway", bullets)
    footer(draw, footer_text)
    canvas.save(path)


def make_dual_crop_slide(
    path: Path,
    title: str,
    subtitle: str,
    left_image: Image.Image,
    right_image: Image.Image,
    bullets: list[str],
    footer_text: str,
    left_caption: str | None = None,
    right_caption: str | None = None,
) -> None:
    canvas, draw = new_canvas()
    title_block(draw, title, subtitle)
    paste_image_obj(canvas, left_image, MARGIN_X, CONTENT_Y, 876, 694)
    paste_image_obj(canvas, right_image, 986, CONTENT_Y, 876, 694)
    if left_caption:
        draw.rounded_rectangle((MARGIN_X + 28, CONTENT_Y + 18, MARGIN_X + 248, CONTENT_Y + 54), radius=18, fill="#EFF6FF")
        draw.text((MARGIN_X + 48, CONTENT_Y + 25), left_caption, font=FONT_SMALL, fill=ACCENT)
    if right_caption:
        draw.rounded_rectangle((1014, CONTENT_Y + 18, 1234, CONTENT_Y + 54), radius=18, fill="#EFF6FF")
        draw.text((1034, CONTENT_Y + 25), right_caption, font=FONT_SMALL, fill=ACCENT)
    bottom_strip(draw, "Takeaway", bullets)
    footer(draw, footer_text)
    canvas.save(path)


def make_text_slide(path: Path, title: str, subtitle: str, cards: list[tuple[str, list[str]]], footer_text: str) -> None:
    canvas, draw = new_canvas()
    title_block(draw, title, subtitle)

    card_w = 540
    xs = [MARGIN_X, 690, 1274]
    for x, (card_title, card_bullets) in zip(xs, cards):
        draw.rounded_rectangle((x, 204, x + card_w, 818), radius=28, fill=PANEL, outline=PANEL_BORDER, width=2)
        draw.rounded_rectangle((x + 24, 232, x + 260, 270), radius=18, fill="#EFF6FF")
        draw.text((x + 42, 240), card_title, font=FONT_SMALL, fill=ACCENT)
        draw_bullets(draw, card_bullets, x + 28, 324, card_w - 56)

    footer(draw, footer_text)
    canvas.save(path)


def build_slides() -> list[Path]:
    slides: list[Path] = []

    slide1 = SLIDES_DIR / "slide_01_title.png"
    make_title_slide(slide1)
    slides.append(slide1)

    slide2 = SLIDES_DIR / "slide_02_step0_freeze.png"
    ctl_top_left = crop(IMG_CTL_SUMMARY, (12, 48, 532, 405))
    ctl_top_right = crop(IMG_CTL_SUMMARY, (544, 48, 1084, 405))
    make_dual_crop_slide(
        slide2,
        "Step 0: freeze the model structure before testing robustness",
        "Direct versus coupled CTL form near the b6 threshold",
        ctl_top_left,
        ctl_top_right,
        [
            "Tumor day 100 and tumor day 365 are almost identical for the two CTL formulations near the threshold.",
            "So Step 1 fixes the simpler direct CTL formulation instead of changing structure again.",
            "All Step 1 runs then keep baseline initial condition and b5 = 1e-4 while varying only b6.",
        ],
        "Step 0 is a modelling decision slide, not a new parameter search slide.",
        left_caption="Tumor day 365",
        right_caption="Tumor day 100",
    )
    slides.append(slide2)

    slide3 = SLIDES_DIR / "slide_03_traj_365.png"
    traj365_top = crop(IMG_TRAJ_365, (0, 0, 1100, 1090))
    traj365_bottom = crop(IMG_TRAJ_365, (0, 1090, 1100, 2200))
    make_dual_crop_slide(
        slide3,
        "Step 1A: trajectories over 0-365 days already separate by b6",
        "Lower b6 cases on the left, threshold-to-high b6 cases on the right",
        traj365_top,
        traj365_bottom,
        [
            "Below about 6e-7, the first-year trajectories remain clearly tumor-dominant.",
            "Around 7e-7 and above, tumor drops toward the 10^3 scale by day 365.",
            "This is why the trajectory path matters: it shows how the response develops, not just the endpoint value.",
        ],
        "Adam likes trajectories because they reveal the path, not only T365.",
    )
    slides.append(slide3)

    slide4 = SLIDES_DIR / "slide_04_traj_10000.png"
    traj10000_top = crop(IMG_TRAJ_10000, (0, 0, 1100, 1090))
    traj10000_bottom = crop(IMG_TRAJ_10000, (0, 1090, 1100, 2200))
    make_dual_crop_slide(
        slide4,
        "Step 1B: by 10,000 days, delayed growth is not always durable control",
        "Same representative b6 values, now followed to long time",
        traj10000_top,
        traj10000_bottom,
        [
            "The b6 ≈ 6.05e-7 case still escapes by long time even though it looks partially suppressed at day 365.",
            "The ~7e-7 region settles near an intermediate or dormancy-like state instead of rebounding to 10^7.",
            "The 2e-6 and 3e-6 cases stay controlled through the full 10,000-day window.",
        ],
        "Long-time trajectories are what turn a threshold hint into a robustness story.",
    )
    slides.append(slide4)

    slide5 = SLIDES_DIR / "slide_05_threshold_lower.png"
    thresh_lower_left = crop(IMG_TRAJ_THRESH, (0, 0, 1098, 660))
    thresh_lower_right = crop(IMG_TRAJ_THRESH, (0, 660, 1098, 1320))
    make_dual_crop_slide(
        slide5,
        "Step 1C: lower side of the threshold shows delayed escape",
        "Threshold-region trajectories, roughly 4e-7 through 6.56e-7",
        thresh_lower_left,
        thresh_lower_right,
        [
            "From 4e-7 to about 6.56e-7, all cases still move toward tumor-dominant behavior.",
            "What changes is the escape timing: growth becomes slower and slower as b6 increases.",
            "This is the pre-switch side of the result, not true long-time control yet.",
        ],
        "Threshold-region trajectory figure is included because the switch is the main Step 1 question.",
    )
    slides.append(slide5)

    slide6 = SLIDES_DIR / "slide_06_threshold_upper.png"
    thresh_upper_left = crop(IMG_TRAJ_THRESH, (0, 1080, 1098, 1760))
    thresh_upper_right = crop(IMG_TRAJ_THRESH, (0, 1760, 1098, 2420))
    make_dual_crop_slide(
        slide6,
        "Step 1D: upper side of the threshold stays low over long time",
        "Threshold-region trajectories, roughly 6.95e-7 through 9e-7",
        thresh_upper_left,
        thresh_upper_right,
        [
            "Across this narrow band, the long-time tumor level stays near the 10^3 scale rather than returning to 10^7.",
            "So the switch window is narrow, and it is visible directly in the trajectories.",
            "The real controlled regime begins later, around 1.15e-6, outside this threshold-only figure.",
        ],
        "This is the clearest trajectory evidence that the change is sharp, not broad.",
    )
    slides.append(slide6)

    slide7 = SLIDES_DIR / "slide_07_summary.png"
    summary_top = crop(IMG_SUMMARY, (0, 0, 1300, 470))
    summary_bottom = crop(IMG_SUMMARY, (0, 470, 1300, 950))
    make_dual_crop_slide(
        slide7,
        "Step 1E: ODE summary plots confirm the same story numerically",
        "Top row: tumor outputs and crossing times; bottom row: burden and immune state",
        summary_top,
        summary_bottom,
        [
            "T100, T365 and T10000 all move with b6, but the long-time readout is the strictest one.",
            "Time-to-1e5 and time-to-1e6 stretch sharply near the threshold window.",
            "The long-time control state comes with high B and high NK, while CTL peaks near the transition itself.",
        ],
        "Trajectories tell the path; summary plots anchor the thresholds and metrics.",
    )
    slides.append(slide7)

    slide8 = SLIDES_DIR / "slide_08_zoom.png"
    zoom_top = crop(IMG_ZOOM, (0, 0, 1300, 470))
    zoom_bottom = crop(IMG_ZOOM, (0, 470, 1300, 950))
    make_dual_crop_slide(
        slide8,
        "Step 1F: the threshold is narrow, and day 365 is not the whole story",
        "Zoomed view of the switch region",
        zoom_top,
        zoom_bottom,
        [
            "The first T10000 < 1e5 appears near 6.69e-7, while the first T10000 < 1e3 appears near 1.15e-6.",
            "T365 and T10000 do not define exactly the same threshold point.",
            "That is why the 10,000-day extension was necessary after the meeting feedback.",
        ],
        "The zoomed summary turns the trajectory story into explicit threshold numbers.",
    )
    slides.append(slide8)

    slide9 = SLIDES_DIR / "slide_09_feedback_response.png"
    make_text_slide(
        slide9,
        "How this Step 0 + Step 1 deck answers the meeting feedback",
        "Only the current work, not future extensions",
        [
            (
                "Step 0",
                [
                    "Freeze the structure instead of re-editing the CTL term again.",
                    "Use the simpler direct CTL form because tumor-level outputs are unchanged near threshold.",
                ],
            ),
            (
                "Step 1",
                [
                    "Use finer b6 grades with a dense threshold region.",
                    "Extend ODE runs to 10,000 days.",
                    "Keep trajectory plots as a main result, not a supplement.",
                ],
            ),
            (
                "Current discussion points",
                [
                    "Is the intermediate regime biologically meaningful?",
                    "Which endpoint should define control: day 365 or long-time behavior?",
                    "Is the narrow b6 switch strong enough to motivate the next analysis step?",
                ],
            ),
        ],
        "This version is built to show Adam the path, the threshold, and the modelling decision.",
    )
    slides.append(slide9)

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
