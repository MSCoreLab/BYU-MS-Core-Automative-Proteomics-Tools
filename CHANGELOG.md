# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- MSPP Web App with React + TypeScript frontend and Flask backend
- Grouped fold change plot with regex pattern matching
- Instance-level data caching for improved performance
- Browser-based launcher with Chrome app mode support
- Dark mode UI for all visualizations

### Changed
- Refactored MSPP_data_plotter.py for optimal performance
- Optimized backend API with data caching (5-10x faster subsequent requests)
- Extracted boxplot styling to helper methods (reduced code duplication)

### Fixed
- TypeScript build errors with CSS module declarations
- Organism identification pattern matching for E.coli variants

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
