# Quick Reference Card

## 🚀 Common Commands

### Setup (First Time)
```powershell
# Automated setup
.\scripts\setup_dev.ps1

# Manual setup
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

### Daily Development
```powershell
# Activate environment
.\.venv\Scripts\Activate.ps1

# Run web app
python programs/mspp_web/launch_app.py

# Run desktop app
python programs/pyscripts/MSPP_data_plotter.py

# Deactivate when done
deactivate
```

### Frontend Development (Web App)
```powershell
cd programs/mspp_web/frontend

# Install dependencies
npm install

# Development mode (hot reload)
npm run dev

# Production build
npm run build
```

### Code Quality
```powershell
# Format code with Black
black .

# Lint with Ruff
ruff check .

# Run tests
pytest tests/
```

### Git Workflow
```powershell
# Create feature branch
git checkout -b feature/my-feature

# Stage changes
git add .

# Commit with conventional message
git commit -m "feat: Add new feature"

# Push to GitHub
git push origin feature/my-feature
```

## 📁 Key Files

| File | Purpose |
|------|---------|
| `pyproject.toml` | Python project metadata & dependencies |
| `requirements.txt` | Python package list |
| `.gitignore` | Files to exclude from git |
| `CONTRIBUTING.md` | Contribution guidelines |
| `CHANGELOG.md` | Version history |

## 🔧 Troubleshooting

### Python Environment Issues
```powershell
# Remove old environment
Remove-Item -Recurse -Force .venv

# Create fresh environment
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

### Frontend Build Issues
```powershell
cd programs/mspp_web/frontend

# Clear cache and reinstall
Remove-Item -Recurse -Force node_modules
Remove-Item package-lock.json
npm install
npm run build
```

### Git Issues
```powershell
# See what's being ignored
git status --ignored

# Force add if needed (careful!)
git add -f path/to/file
```

## 📊 Application URLs

| App | URL |
|-----|-----|
| Web App | http://localhost:5000 |
| Frontend Dev | http://localhost:5173 (if running separately) |

## 🎯 Quick Tips

- **Always activate venv** before running Python scripts
- **Commit `frontend/dist/`** for deployment
- **Don't commit** `node_modules/`, `.venv/`, `__pycache__/`
- **Use conventional commits**: `feat:`, `fix:`, `docs:`, `refactor:`
- **Test before pushing**: Run app to verify it works
- **Update CHANGELOG.md** when adding features

## 📞 Getting Help

1. Check `CONTRIBUTING.md` for guidelines
2. Check `README.md` for documentation
3. Check `programs/mspp_web/README.md` for web app details
4. Open an issue on GitHub

---

**Tip:** Bookmark this file for quick reference! 📑
