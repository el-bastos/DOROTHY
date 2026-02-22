"""
Dorothy - Main Application Entry Point
"""

import sys
from pathlib import Path
from PyQt6.QtWidgets import QApplication, QStyleFactory
from PyQt6.QtCore import QTranslator, QLocale, QLibraryInfo
from PyQt6.QtGui import QFontDatabase, QFont
from dorothy.ui.main_window import MainWindow
from dorothy.ui.styles import FONT_FAMILY, SIZE_BODY


def main():
    """Launch the Dorothy application."""
    app = QApplication(sys.argv)
    app.setStyle(QStyleFactory.create("Fusion"))
    app.setApplicationName("Dorothy")
    app.setOrganizationName("Dorothy")
    app.setApplicationVersion("0.1.0")

    # Load Atkinson Hyperlegible font
    font_dir = Path(__file__).parent.parent / "font"
    if font_dir.exists():
        for ttf in font_dir.glob("*.ttf"):
            QFontDatabase.addApplicationFont(str(ttf))

    # Set as application-wide default font
    app_font = QFont(FONT_FAMILY, SIZE_BODY)
    app.setFont(app_font)

    # Set up translations
    translator = QTranslator()
    locale = QLocale.system()

    # Try to load translation for system locale
    translations_path = Path(__file__).parent / "resources" / "translations"
    if translator.load(locale, "dorothy", "_", str(translations_path)):
        app.installTranslator(translator)

    # Load Qt's own translations for standard dialogs
    qt_translator = QTranslator()
    if qt_translator.load(locale, "qtbase", "_", QLibraryInfo.path(QLibraryInfo.LibraryPath.TranslationsPath)):
        app.installTranslator(qt_translator)

    # Create and show main window
    window = MainWindow()
    window.show()

    sys.exit(app.exec())


if __name__ == "__main__":
    main()
