# Python System Package Cleanup Script
# This script removes unnecessary packages from the system Python installation
# Date: January 2, 2026

Write-Host "Python Package Cleanup Script" -ForegroundColor Cyan
Write-Host "==============================`n" -ForegroundColor Cyan

# Verify Python location
$pythonPath = python -c "import sys; print(sys.executable)"
Write-Host "Python Installation: $pythonPath`n" -ForegroundColor Yellow

# Packages to remove - organized by category
$packagesToRemove = @(
    # Azure AI Services
    "azure-ai-documentintelligence",
    "azure-core",
    "azure-identity",
    "msal",
    "msal-extensions",
    
    # Document Processing (keeping markitdown as requested)
    "mammoth",
    "python-pptx",
    "pdfminer.six",
    
    # Excel/Spreadsheet (unless actively used)
    "openpyxl",
    "xlrd",
    "xlsxwriter",
    
    # Audio/Media Processing
    "audioop-lts",
    "pydub",
    "SpeechRecognition",
    
    # YouTube API
    "youtube-transcript-api",
    
    # OpenAI API (unless actively using)
    "openai"
)

Write-Host "The following packages will be removed:`n" -ForegroundColor Yellow
$packagesToRemove | ForEach-Object { Write-Host "  - $_" }

Write-Host "`nTotal packages to remove: $($packagesToRemove.Count)`n" -ForegroundColor Yellow

# Prompt for confirmation
$confirmation = Read-Host "Do you want to proceed with uninstallation? (y/n)"

if ($confirmation -ne 'y') {
    Write-Host "`nOperation cancelled." -ForegroundColor Red
    exit
}

Write-Host "`nStarting package removal...`n" -ForegroundColor Green

$successCount = 0
$failCount = 0
$notFoundCount = 0

foreach ($package in $packagesToRemove) {
    Write-Host "Removing: $package..." -ForegroundColor Cyan
    
    try {
        $output = python -m pip uninstall -y $package 2>&1
        
        if ($LASTEXITCODE -eq 0) {
            Write-Host "  ✓ Successfully removed $package" -ForegroundColor Green
            $successCount++
        } else {
            if ($output -match "not installed") {
                Write-Host "  ⚠ Package not found: $package" -ForegroundColor Yellow
                $notFoundCount++
            } else {
                Write-Host "  ✗ Failed to remove $package" -ForegroundColor Red
                $failCount++
            }
        }
    } catch {
        Write-Host "  ✗ Error removing $package : $_" -ForegroundColor Red
        $failCount++
    }
    
    Start-Sleep -Milliseconds 100
}

Write-Host "`n==============================" -ForegroundColor Cyan
Write-Host "Cleanup Summary" -ForegroundColor Cyan
Write-Host "==============================" -ForegroundColor Cyan
Write-Host "Successfully removed: $successCount" -ForegroundColor Green
Write-Host "Not found: $notFoundCount" -ForegroundColor Yellow
Write-Host "Failed: $failCount" -ForegroundColor Red
Write-Host "`nDone!" -ForegroundColor Green
