@echo off
REM MSPP Data Plotter - Desktop Launcher
REM Double-click this file to launch the app

setlocal enabledelayedexpansion

cd /d "c:\GitHub\BYU-MS-Core-Automative-Proteomics-Tools\programs\mspp_web"

echo ===============================================================
echo Starting MSPP Data Plotter...
echo ===============================================================
echo.

REM Run the Python launcher
"c:\GitHub\BYU-MS-Core-Automative-Proteomics-Tools\.venv\Scripts\python.exe" launch_app.py

if errorlevel 1 (
    echo.
    echo ===============================================================
    echo ERROR: Failed to start the app. Check the message above.
    echo ===============================================================
    pause
)
