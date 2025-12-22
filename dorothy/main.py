"""
Dorothy - Main Application Entry Point
"""

import sys
from pathlib import Path
from PyQt6.QtWidgets import QApplication
from PyQt6.QtCore import QTranslator, QLocale, QLibraryInfo
from dorothy.ui.main_window import MainWindow


def main():
    """Launch the Dorothy application."""
    app = QApplication(sys.argv)
    app.setApplicationName("Dorothy")
    app.setOrganizationName("Dorothy")
    app.setApplicationVersion("0.1.0")

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
