# Implementation Plan: MAFFT to Mappy Migration

## Status: IN PROGRESS

**Last Updated:** January 28, 2026

This document tracks the implementation of the complete replacement of MAFFT with mappy/minimap2.

## Implementation Summary

The implementation completely replaces MAFFT with:
1. **Mappy** for ORF-by-ORF local alignment
2. **Mappy with context windows** for small regions (MSD, PSI, RRE) - \u00b1100bp around feature
3. **aligntools** library (CigarHit) for coordinate mapping

**Key Update (Jan 28):** Unified approach using only mappy. No BioPython PairwiseAligner needed.

## Files Created

✅ `DESIGN.md` - Complete design document (updated for unified mappy approach)  
✅ `src/cfeintact/orf_aligner.py` - Mappy-based alignment for ORFs and small regions
  - `align_orf()` for large ORFs
  - `align_region()` for small features with context window

## Files Removed

✅ `src/cfeintact/region_aligner.py` - No longer needed (mappy handles everything)

## Files To Modify/Remove

### High Priority
- [ ] `src/cfeintact/intact.py` - Main analysis logic
  - Remove AlignedSequence usage
  - Use OrfAligner.align_orf() for ORF detection
  - Use OrfAligner.align_region() for MSD/PSI/RRE (with context=100)
  
- [ ] `src/cfeintact/find_orf.py` - ORF detection logic
  - Replace AlignedSequence with OrfAlignment
  - Use CigarHit from aligntools
  
- [ ] `src/cfeintact/mapped_orf.py` - Add cigar_hit field
  - Include CigarHit for coordinate mapping
  
- [ ] `src/cfeintact/initialize_orf.py` - Update to not use AlignedSequence

### Files To Remove
- [ ] `src/cfeintact/aligned_sequence.py` - MAFFT-based alignment
- [ ] `src/cfeintact/wrappers.py` - Contains mafft() function

### Dependencies
- [ ] `pyproject.toml` - Add mappy, already has aligntools
- [ ] `Dockerfile` - Remove MAFFT installation
- [ ] `.github/workflows/main.yml` - Remove MAFFT from CI

### Documentation
- [ ] `docs/workflow.md` - Update to describe mappy approach with context windows
- [ ] `docs/installation.md` - Remove MAFFT instructions
- [ ] `README.md` - Update installation steps

### Tests
- [ ] `tests/test_orf_aligner.py` - Unit tests for ORF alignment
- [ ] `tests/test_region_aligner.py` - Unit tests for region alignment
- [ ] `tests/test_end_to_end.py` - Update integration tests
- [ ] Remove `use_mappy` parameter (no longer needed)

## Next Step

Need to understand the complete flow and refactor systematically. The code has circular dependencies that need careful handling.

Current blockers:
1. `initialize_orf` depends on `find_orf` which depends on `AlignedSequence`
2. `intact.py` uses `AlignedSequence` for coordinate mapping
3. Need to trace through complete usage to refactor properly

##Approach

1. First, create parallel implementations that don't remove old code
2. Add new ORF detection using mappy
3. Test thoroughly
4. Remove old MAFFT code once validated

