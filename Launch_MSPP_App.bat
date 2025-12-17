@echo off
REM MSPP Data Plotter - Desktop Launcher
REM Double-click this file to launch the app

setlocal enabledelayedexpansion

cd /d "c:\Repositories\BYU-Core-MS-Lab\programs\mspp_web"

echo ===============================================================
echo Starting MSPP Data Plotter...
echo ===============================================================
echo.

REM Run the Python launcher
"c:\Repositories\BYU-Core-MS-Lab\.venv\Scripts\python.exe" launch_app.py

if errorlevel 1 (
    echo.
    echo ===============================================================
    echo ERROR: Failed to start the app. Check the message above.
    echo ===============================================================
    pause
)
