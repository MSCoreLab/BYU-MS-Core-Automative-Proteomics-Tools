# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.1] - 2025-12-18

### Added
- MSPP Web App with React + TypeScript frontend and Flask backend
- Grouped fold change plot with regex pattern matching
- Instance-level data caching for improved performance
- Browser-based launcher with Chrome app mode support
- Dark mode UI for all visualizations
- Protein ID count labels directly on bar chart segments for easier reading
- CodeQL security analysis workflow for automated vulnerability scanning
- SECURITY.md policy document with supported versions and vulnerability reporting procedures
- `.bat` shortcut file for quick webapp launch

### Changed
- Refactored MSPP_data_plotter.py into three-layer architecture (DataProcessor, PlotGenerator, MSPPDataPlotter)
- Separated data processing logic from UI elements and plotting functions
- Migrated modular architecture to webapp backend for feature parity
- Implemented consensus protein calculation (only compares proteins present in both E25 and E100 samples)
- Updated Yeast protein labels to Y150/Y75 for accurate spike-in ratio representation
- Optimized backend API with data caching (5-10x faster subsequent requests)
- Extracted boxplot styling to helper methods (reduced code duplication)
- Rebuilt frontend to reflect updated backend API structure

### Fixed
- TypeScript build errors with CSS module declarations
- Organism identification pattern matching for E.coli variants
- Fold change "dilution" issue by implementing consensus protein filtering
- Mix identifier extraction to support 4-parameter pattern for scalability
- Flask debug mode security vulnerability (disabled debug in production)

### Security
- Disabled Flask debug mode in production environment to prevent information disclosure
- Added automated CodeQL analysis for continuous security monitoring

## [0.1.0] - 2025-12-08

### Added
- Initial repository structure
- MSPP Data Plotter desktop GUI (tkinter)
- FASTA filter GUI tool
- Protein ID bar chart visualization
- E.coli vs Yeast fold change analysis
- Organisms vs HeLa spike-in validation plots
- Dark mode matplotlib styling
- Vectorized organism identification for performance

### Documentation
- Project README
- FA Workflow Tutorial
- Reference literature collection

---

## Format Guide

### Added
New features

### Changed
Changes in existing functionality

### Deprecated
Soon-to-be removed features

### Removed
Removed features

### Fixed
Bug fixes

### Security
Vulnerability fixes
