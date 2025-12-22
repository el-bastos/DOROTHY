# Translations

Dorothy uses Qt's translation system for internationalization.

## Current Languages

- English (default, built-in)

## Adding a New Language

1. Extract translatable strings:
   ```bash
   pylupdate6 dorothy/**/*.py -ts dorothy/resources/translations/dorothy_XX.ts
   ```
   Replace `XX` with the language code (e.g., `pt_BR` for Brazilian Portuguese, `es` for Spanish)

2. Open the `.ts` file in Qt Linguist to translate
   - Download Qt Linguist: https://doc.qt.io/qt-6/linguist-translators.html
   - Or edit the XML manually

3. Compile the translation:
   ```bash
   lrelease dorothy/resources/translations/dorothy_XX.ts
   ```
   This creates `dorothy_XX.qm` which the app will load automatically

## Language Codes

Common codes:
- `en` - English (default)
- `pt_BR` - Portuguese (Brazil)
- `es` - Spanish
- `fr` - French
- `de` - German
- `it` - Italian
- `ja` - Japanese
- `zh_CN` - Chinese (Simplified)

Dorothy will automatically detect the system language and load the appropriate translation if available.
