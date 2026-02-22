# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller spec for Dorothy."""

import sys
from pathlib import Path

block_cipher = None
base = Path(SPECPATH)

datas = [
    (str(base / 'font'), 'font'),
    (str(base / 'logo' / 'dorothy_logo.png'), 'logo'),
    (str(base / 'logo' / 'dorothy_logo_icon.png'), 'logo'),
    (str(base / 'examples'), 'examples'),
    (str(base / 'dorothy' / 'resources'), 'dorothy/resources'),
]

a = Analysis(
    [str(base / 'dorothy' / 'main.py')],
    pathex=[str(base)],
    binaries=[],
    datas=datas,
    hiddenimports=[
        'matplotlib.backends.backend_agg',
        'matplotlib.backends.backend_qtagg',
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=['tkinter'],
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='Dorothy',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    icon=str(base / 'logo' / 'dorothy_logo_icon.png'),
)

coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    name='Dorothy',
)

if sys.platform == 'darwin':
    app = BUNDLE(
        coll,
        name='Dorothy.app',
        icon=str(base / 'logo' / 'dorothy_logo_icon.png'),
        bundle_identifier='com.dorothy.app',
        info_plist={
            'CFBundleShortVersionString': '0.7.0',
            'CFBundlePackageType': 'APPL',
            'NSPrincipalClass': 'NSApplication',
            'NSHighResolutionCapable': True,
            'LSMinimumSystemVersion': '12.0',
        },
    )
