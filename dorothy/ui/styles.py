"""
Design system for Dorothy UI.

Centralizes all colors, typography, and widget styles.
Applied per-widget (no global stylesheet) to avoid cascade issues.
"""

# =============================================================================
# Colors
# =============================================================================

ACCENT = "#e05a2b"
ACCENT_HOVER = "#c94d24"
ACCENT_PRESSED = "#a84120"

TEXT = "#1a1a1a"
TEXT_SECONDARY = "#333333"
TEXT_MUTED = "#555555"

BG = "#ffffff"
BG_CARD = "#f0f0f0"       # card/container backgrounds
BG_HOVER = "#e4e4e4"

BORDER = "#aaaaaa"         # strong enough to see on white
BORDER_FOCUS = "#666666"

SUCCESS = "#2e7d32"
ERROR = "#c62828"

# =============================================================================
# Typography — minimum 11pt
# =============================================================================

FONT_FAMILY = "Atkinson Hyperlegible"

SIZE_TITLE = 18
SIZE_SUBTITLE = 14
SIZE_BODY = 12
SIZE_SMALL = 11


# =============================================================================
# Buttons — WHITE background so they work on both white screens AND gray cards
# =============================================================================

def btn_secondary() -> str:
    """White bg, dark text, visible border. Works on white AND gray backgrounds."""
    return f"""
        QPushButton {{
            background-color: white;
            color: {TEXT};
            border: 1px solid {BORDER};
            border-radius: 6px;
            padding: 6px 16px;
            font-size: {SIZE_BODY}pt;
            font-family: "{FONT_FAMILY}";
        }}
        QPushButton:hover {{
            background-color: {BG_CARD};
            border-color: {BORDER_FOCUS};
        }}
        QPushButton:pressed {{
            background-color: {BG_HOVER};
        }}
    """


def btn_primary() -> str:
    """Orange bg, white text."""
    return f"""
        QPushButton {{
            background-color: {ACCENT};
            color: white;
            border: none;
            border-radius: 8px;
            padding: 10px 24px;
            font-size: 13pt;
            font-weight: bold;
            font-family: "{FONT_FAMILY}";
        }}
        QPushButton:hover {{
            background-color: {ACCENT_HOVER};
        }}
        QPushButton:pressed {{
            background-color: {ACCENT_PRESSED};
        }}
    """


def btn_ghost() -> str:
    """Transparent bg, accent text."""
    return f"""
        QPushButton {{
            background-color: transparent;
            color: {ACCENT};
            border: none;
            border-radius: 6px;
            padding: 6px 16px;
            font-size: {SIZE_BODY}pt;
            font-weight: bold;
            font-family: "{FONT_FAMILY}";
        }}
        QPushButton:hover {{
            background-color: {BG_CARD};
        }}
    """


def btn_toggle() -> str:
    """Checkable pill toggle."""
    return f"""
        QPushButton {{
            background-color: transparent;
            color: {TEXT_SECONDARY};
            border: none;
            border-radius: 6px;
            padding: 6px 14px;
            font-size: {SIZE_BODY}pt;
            font-family: "{FONT_FAMILY}";
        }}
        QPushButton:checked {{
            background-color: {ACCENT};
            color: white;
            font-weight: bold;
        }}
        QPushButton:hover:!checked {{
            background-color: {BG_CARD};
        }}
    """


def btn_toolbar() -> str:
    """Compact uniform toolbar button. Works for regular and checkable buttons."""
    return f"""
        QPushButton {{
            background-color: white;
            color: {TEXT};
            border: 1px solid {BORDER};
            border-radius: 4px;
            padding: 4px 10px;
            font-size: {SIZE_SMALL}pt;
            font-family: "{FONT_FAMILY}";
        }}
        QPushButton:hover {{
            background-color: {BG_CARD};
            border-color: {BORDER_FOCUS};
        }}
        QPushButton:pressed {{
            background-color: {BG_HOVER};
        }}
        QPushButton:checked {{
            background-color: {ACCENT};
            color: white;
            border-color: {ACCENT};
        }}
        QPushButton:checked:hover {{
            background-color: {ACCENT_HOVER};
            border-color: {ACCENT_HOVER};
        }}
    """


def btn_icon() -> str:
    """Small square button for +/- zoom, arrows."""
    return f"""
        QPushButton {{
            background-color: white;
            color: {TEXT};
            border: 1px solid {BORDER};
            border-radius: 4px;
            padding: 2px;
            font-size: {SIZE_BODY}pt;
            font-family: "{FONT_FAMILY}";
            min-width: 28px;
            min-height: 28px;
            max-width: 28px;
            max-height: 28px;
        }}
        QPushButton:hover {{
            background-color: {BG_CARD};
            border-color: {BORDER_FOCUS};
        }}
        QPushButton:pressed {{
            background-color: {BG_HOVER};
        }}
    """


def btn_search() -> str:
    """Pill-shaped button matching the search input height/shape."""
    return f"""
        QPushButton {{
            background-color: white;
            color: {TEXT};
            border: 1px solid {BORDER};
            border-radius: 20px;
            padding: 6px 20px;
            font-size: {SIZE_BODY}pt;
            font-family: "{FONT_FAMILY}";
        }}
        QPushButton:hover {{
            background-color: {BG_CARD};
            border-color: {BORDER_FOCUS};
        }}
    """


# =============================================================================
# Input styles
# =============================================================================

def input_search() -> str:
    """Google-style rounded search field. No min-height — sized by font+padding."""
    return f"""
        QLineEdit {{
            background-color: white;
            color: {TEXT};
            border: 1px solid {BORDER};
            border-radius: 20px;
            padding: 8px 20px;
            font-size: {SIZE_SUBTITLE}pt;
            font-family: "{FONT_FAMILY}";
        }}
        QLineEdit:focus {{
            border-color: {BORDER_FOCUS};
        }}
    """


def combo_style() -> str:
    return f"""
        QComboBox {{
            background-color: white;
            color: {TEXT};
            border: 1px solid {BORDER};
            border-radius: 4px;
            padding: 4px 10px;
            font-size: {SIZE_BODY}pt;
            font-family: "{FONT_FAMILY}";
        }}
        QComboBox:hover {{
            border-color: {BORDER_FOCUS};
        }}
        QComboBox::drop-down {{
            border: none;
        }}
    """


def spinbox_style() -> str:
    return f"""
        QSpinBox {{
            background-color: white;
            color: {TEXT};
            border: 1px solid {BORDER};
            border-radius: 4px;
            padding: 4px 10px;
            font-size: {SIZE_BODY}pt;
            font-family: "{FONT_FAMILY}";
        }}
    """


def checkbox_style() -> str:
    return f"""
        QCheckBox {{
            color: {TEXT};
            font-size: {SIZE_BODY}pt;
            font-family: "{FONT_FAMILY}";
            spacing: 6px;
        }}
    """


# =============================================================================
# Containers
# =============================================================================

def card_style() -> str:
    """Gray card on white screen — buttons inside are white, so they stand out."""
    return f"""
        QFrame {{
            background: {BG_CARD};
            border: 1px solid {BORDER};
            border-radius: 8px;
            padding: 12px;
        }}
        QFrame:hover {{
            background: {BG_HOVER};
            border-color: {BORDER_FOCUS};
        }}
    """


def groupbox_style() -> str:
    return f"""
        QGroupBox {{
            font-family: "{FONT_FAMILY}";
            font-size: {SIZE_SMALL}pt;
            font-weight: bold;
            color: {TEXT_SECONDARY};
            border: 1px solid {BORDER};
            border-radius: 6px;
            margin-top: 12px;
            padding: 16px 8px 8px 8px;
        }}
        QGroupBox::title {{
            subcontrol-origin: margin;
            left: 10px;
            padding: 0 4px;
        }}
    """


def screen_bg() -> str:
    return f"background-color: {BG};"


def progress_bar() -> str:
    return f"""
        QProgressBar {{
            border: 1px solid {BORDER};
            border-radius: 4px;
            text-align: center;
            font-family: "{FONT_FAMILY}";
            color: {TEXT};
            background: {BG_CARD};
        }}
        QProgressBar::chunk {{
            background-color: {ACCENT};
            border-radius: 3px;
        }}
    """


# =============================================================================
# Labels
# =============================================================================

def label_title() -> str:
    return f"color: {TEXT};"

def label_secondary() -> str:
    return f"color: {TEXT_SECONDARY};"

def label_muted() -> str:
    return f"color: {TEXT_MUTED};"

def label_success() -> str:
    return f"color: {SUCCESS};"

def label_error() -> str:
    return f"color: {ERROR};"
