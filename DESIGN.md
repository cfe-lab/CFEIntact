# CFEIntact: Minimap2 (Mappy) Integration Design

**Date:** January 27, 2026  
**Version:** 2.0 (Complete Redesign)

---

## Executive Summary

This document describes the complete replacement of MAFFT with minimap2 (via mappy Python bindings) in CFEIntact. The key insight is to leverage minimap2's strength as a **local aligner** by aligning individual ORFs rather than performing whole-sequence alignment.

**Key Design Decisions:**
1. **ORF-by-ORF local alignment** using mappy for large regions (gag, pol, env, vif, vpr, tat, rev, vpu, nef)
2. **Context-window alignment** using mappy for small critical regions (MSD, PSI, RRE) - extract ±100bp around feature
3. **Use existing aligntools library** (CigarHit, Cigar, CoordinateMapping) instead of reimplementing
4. **Complete removal of MAFFT** - no hybrid approach, clean replacement
5. **Unified tool** - mappy for all alignment tasks

---

## Problem Analysis

### Current MAFFT Usage

MAFFT is currently used for:
1. **Whole sequence global alignment** between query and reference subtype
2. **Coordinate mapping** from reference (HXB2) positions to query positions
3. **Feature detection** (MSD, PSI, RRE) via aligned sequences

**Problems with MAFFT:**
- Slow (~5-10 seconds per 9kb HIV sequence)
- External dependency requiring system installation
- Overkill for the actual task (we only need specific region alignments)
- Global alignment when we need local alignment for ORFs

### Why Minimap2/Mappy?

**Advantages:**
- **10-50x faster** than MAFFT for local alignment
- **Native Python binding** (pip installable, no external binaries)
- **Designed for local alignment** - perfect for ORF detection
- **CIGAR string output** - compatible with aligntools library
- **Industry standard** - used widely in bioinformatics

**Minimap2's Strength:**
Minimap2 excels at finding **local alignments** between sequences. This is perfect for:
- Finding where gag, pol, env map between reference and query
- Handling insertions/deletions within ORFs
- Fast mapping of HIV-sized sequences

---

## Architecture Design

### Core Principle

**Align what you need, when you need it, at the appropriate granularity - using mappy for everything**

Instead of:
```
[Whole sequence alignment] → [Extract ORF regions] → [Analyze]
```

Do:
```
[Mappy: Local ORF alignment] → [Analyze ORF]
[Mappy: Small region with context] → [Check feature]
```

### Component Breakdown

**Unified Approach:** Use mappy for all alignments, adjusting the alignment window size based on the feature.

#### 1. ORF Alignment (Large Regions)

For large ORFs (gag, pol, env) and small ORFs (vif, vpr, tat, rev, vpu, nef):

```python
# Extract reference ORF sequence
ref_orf_seq = reference_sequence[orf_start:orf_end]

# Align query sequence to this ORF using mappy
aligner = mappy.Aligner(seq=ref_orf_seq, preset="asm20")
hits = list(aligner.map(str(query_sequence)))

# Get best hit
best_hit = max(hits, key=lambda h: h.mlen)

# Convert to aligntools CigarHit
cigar_hit = CigarHit(
    cigar=Cigar.coerce(best_hit.cigar_str),
    r_st=best_hit.r_st,
    r_ei=best_hit.r_en - 1,
    q_st=best_hit.q_st,
    q_ei=best_hit.q_en - 1
)

# Use coordinate mapping for analysis
coordinate_mapping = cigar_hit.coordinate_mapping
```

**Why this works:**
- Mappy finds the best local alignment of the ORF
- Returns CIGAR string describing the alignment
- aligntools CigarHit provides coordinate mapping
- Fast and accurate for ORF-sized regions (500-3000bp)

#### 2. Small Region Detection (MSD, PSI, RRE)

For small features (2-300bp), we know their position in the reference genome. Use mappy with a context window:

```python
# Example: MSD at position 743-744 in HXB2
MSD_POS = 743
CONTEXT = 100  # bp before and after

# Extract reference region with context
ref_start = MSD_POS - CONTEXT  # 643
ref_end = MSD_POS + 2 + CONTEXT  # 845 (2bp feature + context)
ref_region = reference_sequence[ref_start:ref_end]  # 202bp

# Align this region to query using mappy
aligner = mappy.Aligner(seq=ref_region, preset="asm20")
hits = list(aligner.map(str(query_sequence)))
best_hit = max(hits, key=lambda h: h.mlen)

# Convert to CigarHit
cigar_hit = CigarHit(
    cigar=Cigar.coerce(best_hit.cigar_str),
    r_st=best_hit.r_st,
    r_ei=best_hit.r_en - 1,
    q_st=best_hit.q_st,
    q_ei=best_hit.q_en - 1
)

# Map the feature position through the alignment
# MSD is at offset 100 in the extracted region (due to context)
feature_offset = CONTEXT
coord_map = cigar_hit.coordinate_mapping
query_msd_pos = coord_map.ref_to_query.right_min(feature_offset)

# Check if MSD is intact in query
msd_query = query_sequence[query_msd_pos:query_msd_pos + 2]
if msd_query.upper() != "GT":
    report_defect("MSD mutated")
```

**Why this works:**
- We know the exact position in the reference
- Context window (±100bp) provides anchor points for alignment
- Mappy finds where this region maps in the query
- Use coordinate mapping to locate the exact feature
- Consistent with ORF alignment approach

#### 3. Subtype Determination

Current approach uses BLAST for subtype determination. **Keep this unchanged** - BLAST works well for this purpose.

---

## Detailed Implementation Plan

### Phase 1: Core Infrastructure

#### File: `src/cfeintact/orf_aligner.py`

**This module handles both ORF alignment and small region alignment using mappy.**

```python
"""
Local alignment of ORFs using minimap2/mappy.
"""

import mappy  # type: ignore
from Bio.SeqRecord import SeqRecord
from aligntools.cigar_hit import CigarHit
from aligntools.cigar import Cigar
from typing import Optional, List
from dataclasses import dataclass

from cfeintact.original_orf import OriginalORF


@dataclass
class OrfAlignment:
    """Result of aligning an ORF."""
    cigar_hit: CigarHit
    reference_orf: OriginalORF
    query_start: int  # Position in full query sequence
    query_end: int
    mapq: int
    strand: int


class OrfAligner:
    """Aligns ORFs and small regions using minimap2."""
    
    def __init__(self, preset: str = "asm20"):
        self.preset = preset
    
    def align_orf(self, 
                  reference: SeqRecord,
                  query: SeqRecord,
                  orf: OriginalORF) -> Optional[OrfAlignment]:
        """
        Align a specific ORF region from reference to query.
        
        Args:
            reference: Reference sequence (e.g., HXB2 or subtype)
            query: Query sequence to analyze
            orf: ORF definition with start/end positions in reference
        
        Returns:
            OrfAlignment if found, None otherwise
        """
        # Extract ORF sequence from reference
        ref_orf_seq = str(reference.seq[orf.start:orf.end + 1])
        query_seq = str(query.seq)
        
        # Create aligner for this ORF
        aligner = mappy.Aligner(seq=ref_orf_seq, preset=self.preset)
        
        # Find alignments
        hits = list(aligner.map(query_seq))
        
        if not hits:
            return None
        
        # Get best hit (most matching bases)
        best_hit = max(hits, key=lambda h: h.mlen)
        
        # Convert to CigarHit
        cigar = Cigar.coerce(best_hit.cigar_str)
        cigar_hit = CigarHit(
            cigar=cigar,
            r_st=best_hit.r_st,
            r_ei=best_hit.r_en - 1,  # mappy uses exclusive end
            q_st=best_hit.q_st,
            q_ei=best_hit.q_en - 1
        )
        
        return OrfAlignment(
            cigar_hit=cigar_hit,
            reference_orf=orf,
            query_start=best_hit.q_st,
            query_end=best_hit.q_en - 1,
            mapq=best_hit.mapq,
            strand=best_hit.strand
        )
    
    def align_region(self,
                    reference: SeqRecord,
                    query: SeqRecord,
                    region_start: int,
                    region_end: int,
                    context: int = 100) -> Optional[OrfAlignment]:
        """
        Align a small region with context window.
        
        Args:
            reference: Reference sequence
            query: Query sequence
            region_start: Start position of feature in reference
            region_end: End position of feature in reference (inclusive)
            context: Number of bp to include before/after (default 100)
        
        Returns:
            OrfAlignment with coordinates in context window space
        """
        # Extract reference region with context
        ref_start = max(0, region_start - context)
        ref_end = min(len(reference.seq), region_end + 1 + context)
        ref_region = str(reference.seq[ref_start:ref_end])
        query_seq = str(query.seq)
        
        # Align to query
        aligner = mappy.Aligner(seq=ref_region, preset=self.preset)
        hits = list(aligner.map(query_seq))
        
        if not hits:
            return None
        
        best_hit = max(hits, key=lambda h: h.mlen)
        
        # Convert to CigarHit
        cigar = Cigar.coerce(best_hit.cigar_str)
        cigar_hit = CigarHit(
            cigar=cigar,
            r_st=best_hit.r_st,
            r_ei=best_hit.r_en - 1,
            q_st=best_hit.q_st,
            q_ei=best_hit.q_en - 1
        )
        
        # Create a pseudo-ORF for the region
        # Note: coordinates are in the context window space
        from cfeintact.original_orf import OriginalORF
        pseudo_orf = OriginalORF(
            name=f"region_{region_start}_{region_end}",
            start=ref_start,
            end=ref_end - 1,
            nucleotides=reference.seq[ref_start:ref_end],
            aminoacids="",
            protein="",
            max_deletions=0,
            max_insertions=0,
            max_distance=0.0,
            is_small=True
        )
        
        return OrfAlignment(
            cigar_hit=cigar_hit,
            reference_orf=pseudo_orf,
            query_start=best_hit.q_st,
            query_end=best_hit.q_en - 1,
            mapq=best_hit.mapq,
            strand=best_hit.strand
        )
```

### Phase 2: Refactor ORF Detection

Replace `find_orf.py` logic to use `OrfAligner`:

```python
# Old approach (via whole-sequence alignment):
aligned_sequence = AlignedSequence(this=query, reference=ref)
q_start = aligned_sequence.coordinate_mapping.ref_to_query.right_min(orf.start)
q_end = aligned_sequence.coordinate_mapping.ref_to_query.left_max(orf.end)

# New approach (direct ORF alignment):
orf_aligner = OrfAligner()
orf_alignment = orf_aligner.align_orf(reference, query, orf)
if orf_alignment:
    q_start = orf_alignment.query_start
    q_end = orf_alignment.query_end
    cigar_hit = orf_alignment.cigar_hit
    # Use cigar_hit.coordinate_mapping for detailed analysis
```

### Phase 3: Refactor Small Region Detection

For MSD, PSI, RRE - use OrfAligner with context:

```python
# Old approach (via whole-sequence alignment):
alignment = aligned_sequence.alignment
sd_begin = [m.start() for m in re.finditer(r"[^-]", str(alignment[0].seq))][pos]
sd = alignment[1].seq[sd_begin:sd_end]

# New approach (mappy with context window):
orf_aligner = OrfAligner()

# Example: MSD is at position 743 in HXB2
MSD_REF_POS = 743
CONTEXT = 100

# Align region with context
msd_alignment = orf_aligner.align_region(
    reference=reference,
    query=query,
    region_start=MSD_REF_POS,
    region_end=MSD_REF_POS + 1,  # 2bp feature (743-744)
    context=CONTEXT
)

if msd_alignment:
    # Map feature position through alignment
    # Feature is at offset CONTEXT in the extracted region
    coord_map = msd_alignment.cigar_hit.coordinate_mapping
    query_msd_start = coord_map.ref_to_query.right_min(CONTEXT)
    query_msd_end = coord_map.ref_to_query.left_max(CONTEXT + 1)
    
    # Check the actual sequence
    msd_seq = query.seq[msd_alignment.query_start + query_msd_start:
                        msd_alignment.query_start + query_msd_start + 2]
    
    if str(msd_seq).upper() != "GT":
        report_defect("MSD mutated")
else:
    report_defect("MSD not found")
```

### Phase 4: Remove MAFFT Completely

1. **Delete files:**
   - `src/cfeintact/wrappers.py` (contains mafft function)
   - `src/cfeintact/aligned_sequence.py` (uses mafft)

2. **Remove from:**
   - `Dockerfile` - remove mafft installation
   - `.github/workflows/main.yml` - remove mafft from CI
   - `docs/installation.md` - remove mafft instructions
   - `docs/workflow.md` - update to describe mappy approach

3. **Update imports:**
   - Remove `from cfeintact.aligned_sequence import AlignedSequence`
   - Replace with `from cfeintact.orf_aligner import OrfAligner`

---

## Integration with Existing Code

### Mapping to Current Structure

**Current OriginalORF:**
```python
@dataclass
class OriginalORF:
    name: str
    start: int  # Position in reference
    end: int
    nucleotides: Seq
    aminoacids: str
    protein: str
    max_deletions: int
    max_insertions: int
    max_distance: float
    is_small: bool
```

**New workflow:**
```python
# 1. Initialize ORFs from reference (unchanged)
ref_orfs = initialize_orf(reference, name="gag", start=790, end=2292, ...)

# 2. Align each ORF to query using mappy
orf_aligner = OrfAligner()
orf_alignment = orf_aligner.align_orf(reference, query, ref_orfs)

# 3. Extract query ORF using alignment result
query_orf_seq = query.seq[orf_alignment.query_start:orf_alignment.query_end]

# 4. Create query ORF object (similar to current)
query_orf = OriginalORF(
    name="gag",
    start=orf_alignment.query_start,
    end=orf_alignment.query_end,
    nucleotides=query_orf_seq,
    ...
)

# 5. Create MappedORF for analysis
mapped_orf = MappedORF(
    reference=ref_orfs,
    query=query_orf,
    orientation="forward" if orf_alignment.strand == 1 else "reverse",
    cigar_hit=orf_alignment.cigar_hit  # NEW: include for coordinate mapping
)

# 6. Use cigar_hit for indel analysis
indel_impact = analyze_indels(mapped_orf.cigar_hit)
```

### Updated MappedORF

```python
@dataclass(frozen=True)
class MappedORF:
    reference: OriginalORF
    query: OriginalORF
    orientation: str
    cigar_hit: CigarHit  # NEW: aligntools CigarHit for coordinate mapping
    
    @cached_property
    def amino_alignment(self) -> Alignment:
        # Keep existing BioPython alignment for amino acid comparison
        return detailed_aligner.align(self.reference.protein, self.query.protein)
    
    @cached_property
    def coordinate_mapping(self):
        # Use aligntools coordinate mapping
        return self.cigar_hit.coordinate_mapping
    
    @cached_property
    def indel_impact(self) -> int:
        # Use CigarHit to analyze indels
        deletions = list(self.cigar_hit.deletions())
        insertions = list(self.cigar_hit.insertions())
        # Calculate impact...
```

---

## Performance Expectations

### Before (MAFFT):
- **Whole sequence alignment**: ~5-10 seconds per 9kb sequence
- **Memory**: ~200-300 MB per alignment
- **Throughput**: ~10-20 sequences/minute

### After (Mappy):
- **Per-ORF alignment**: ~50-100ms for 9 ORFs = ~0.5-1 second total
- **Small region alignment**: ~10ms each for 3 regions = ~30ms
- **Total per sequence**: ~0.5-1.5 seconds (**5-10x faster**)
- **Memory**: ~50-100 MB (**50-75% reduction**)
- **Throughput**: ~60-120 sequences/minute (**6-10x improvement**)

---

## Testing Strategy

### Unit Tests

1. **ORF Aligner Tests** (`test_orf_aligner.py`):
   - Test alignment of known ORFs
   - Test handling of insertions/deletions
   - Test reverse strand detection
   - Test no-match scenarios

2. **Region Aligner Tests** (`test_region_aligner.py`):
   - Test MSD detection
   - Test PSI alignment
   - Test RRE alignment

3. **Integration Tests**:
   - Run full analysis on test sequences
   - Compare defect detection with previous results
   - Validate ORF positions
   - Check coordinate mappings

### Validation

Run against existing test data and ensure:
- Same defects detected
- Same ORF boundaries (within tolerance)
- Same or better accuracy
- Significantly faster performance

---

## Migration Checklist

### Code Changes

- [ ] Create `orf_aligner.py` with OrfAligner class
- [ ] Create `region_aligner.py` with RegionAligner class
- [ ] Update `find_orf.py` to use OrfAligner
- [ ] Update `mapped_orf.py` to include cigar_hit
- [ ] Update `intact.py` MSD/PSI/RRE detection
- [ ] Remove `aligned_sequence.py`
- [ ] Remove `wrappers.py` mafft function
- [ ] Update all imports

### Dependency Changes

- [ ] Add `mappy>=2.24` to `pyproject.toml`
- [ ] Remove MAFFT from `Dockerfile`
- [ ] Remove MAFFT from `.github/workflows/main.yml`

### Documentation

- [ ] Update `docs/workflow.md`
- [ ] Update `docs/installation.md`
- [ ] Update `README.md`
- [ ] Update code docstrings

### Testing

- [ ] Write unit tests for new components
- [ ] Update integration tests
- [ ] Run full test suite
- [ ] Performance benchmarking
- [ ] Validation against known data

---

## Advantages of This Design

1. ✅ **Uses minimap2's strength** - local alignment for everything
2. ✅ **Uses aligntools properly** - leverages existing CigarHit, CoordinateMapping
3. ✅ **Unified approach** - mappy for all alignments (ORFs and small regions)
4. ✅ **Simpler codebase** - single alignment tool, removes MAFFT complexity
5. ✅ **Better performance** - 5-10x faster
6. ✅ **Fewer dependencies** - no BioPython PairwiseAligner needed
7. ✅ **Consistent coordinate mapping** - same mechanism for all features

---

## Risk Mitigation

### Risk: Different alignment results

**Mitigation:**
- Extensive testing against existing data
- Validate that same defects are detected
- Adjust minimap2 presets if needed (asm5 vs asm20)

### Risk: Small region detection less accurate

**Mitigation:**
- Use BioPython PairwiseAligner (already proven in detailed_aligner.py)
- Can tune alignment parameters if needed
- Test extensively on known MSD mutations

### Risk: Edge cases not handled

**Mitigation:**
- Comprehensive unit tests
- Handle no-alignment cases gracefully
- Provide clear error messages

---

## Conclusion

This redesign properly leverages minimap2 as a **local aligner** for all alignment tasks, uses **aligntools** library as intended, and completely removes MAFFT. The result is faster, simpler, more maintainable code with a unified alignment approach.

**Key Innovation:** Align what you need at the right granularity using a single tool (mappy):
- Large regions (ORFs): align the ORF directly
- Small regions (MSD/PSI/RRE): align with context window (±100bp)

---

**Next Steps:** Implement Phase 1 (Core Infrastructure) and validate with unit tests.
