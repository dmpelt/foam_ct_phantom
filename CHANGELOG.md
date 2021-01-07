# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).


## [Unreleased]
### Added
- Multi-GPU support for cone-beam projections

### Changed
- Faster generation of phantoms

### Fixed
- Generating phantoms with the same seed now actually produces identical phantoms

### Removed

## [1.1.2]
### Added

### Changed
- Fix for Windows builds

### Fixed

### Removed

## [1.1.1]
### Added

### Changed
- Reduced phantom generation time by using a custom skiplist

### Fixed

### Removed

## [1.1.0]
### Added
- Support for nonuniform void distributions
- Support for creating 3D renderings of voids

### Changed

### Fixed
- Saving and loading projection datasets with many angles (8192 or more)

### Removed

## 1.0.0 - 2019-08-02

### Added
- Initial release

[Unreleased]: https://github.com/dmpelt/foam_ct_phantom/compare/v1.1.2...HEAD
[1.1.2]: https://github.com/dmpelt/foam_ct_phantom/compare/v1.1.1...v1.1.2
[1.1.1]: https://github.com/dmpelt/foam_ct_phantom/compare/v1.1.0...v1.1.1
[1.1.0]: https://github.com/dmpelt/foam_ct_phantom/compare/v1.0.0...v1.1.0
