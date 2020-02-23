#!/bin/sh

python3 -m PyInstaller gagrank-cli.py -D -c --clean --noconfirm --add-data ../lib/GAGfragDB.db:GAGfragDB.db --exclude-module PyQt4 --exclude-module PyQt5 --exclude-module _tkinter