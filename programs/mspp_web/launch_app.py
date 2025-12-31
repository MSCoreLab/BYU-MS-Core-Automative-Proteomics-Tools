#!/usr/bin/env python3
"""
MSPP Data Plotter - Simple Desktop Launcher
Opens the app in your default browser - no extra dependencies needed!
"""

import sys
import threading
import time
import webbrowser
from pathlib import Path

# Add backend to path
backend_path = Path(__file__).parent / 'backend'
sys.path.insert(0, str(backend_path))

from backend.app import app


def start_flask():
    """Start Flask server in background thread."""
    app.run(port=5000, debug=False, use_reloader=False, threaded=True)


if __name__ == '__main__':
    print("=" * 60)
    print("🚀 MSPP Data Plotter - Desktop Application")
    print("=" * 60)

    # Start Flask in background
    print("\n⏳ Starting backend server...")
    flask_thread = threading.Thread(target=start_flask, daemon=True)
    flask_thread.start()

    # Wait for Flask to start
    time.sleep(2)

    print("✅ Backend ready on http://localhost:5000")
    print("🎨 Opening application in your default browser...\n")

    # Open in default browser only
    webbrowser.open('http://localhost:5000')

    print("=" * 60)
    print("✨ Application is running!")
    print("=" * 60)
    print("\n📌 The app should open in your browser automatically.")
    print("📌 If not, manually visit: http://localhost:5000")
    print("\n⚠️  Press CTRL+C to stop the server when done.\n")

    try:
        # Keep the server running
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        print("\n\n👋 Shutting down... Goodbye!")
        sys.exit(0)
