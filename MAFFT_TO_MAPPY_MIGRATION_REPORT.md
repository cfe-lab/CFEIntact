# CFEIntact: MAFFT to Mappy Migration Report

**Date:** January 27, 2026  
**Prepared by:** GitHub Copilot AI Assistant  
**Project:** CFEIntact - HIV-1 Proviral Intactness Checker

---

## Executive Summary

This report provides a comprehensive analysis of the CFEIntact codebase with a focus on migrating from MAFFT (Multiple Alignment using Fast Fourier Transform) to mappy (minimap2 Python bindings) for sequence alignment operations. The analysis covers current MAFFT usage patterns, mappy capabilities, migration challenges, and a detailed implementation plan.

**Key Finding:** While mappy offers significant performance advantages, the migration presents substantial technical challenges, particularly in maintaining compatibility with critical features like Major Splice Donor (MSD) site detection. A hybrid approach is recommended.

---

## 1. Project Overview

### 1.1 CFEIntact Purpose
CFEIntact is a bioinformatics tool designed to analyze HIV-1 consensus sequences for proviral intactness. It identifies various defects including:
- ORF (Open Reading Frame) mutations
- Major Splice Donor site mutations
- Packaging Signal (PSI) deletions
- Rev Response Element (RRE) deletions
- Hypermutation signatures
- Internal inversions and scrambling
- Non-HIV sequences

### 1.2 Current Architecture
The tool employs a two-stage alignment strategy:
1. **BLAST**: Region-based alignment for subtype determination and structural analysis
2. **MAFFT**: Global multiple sequence alignment for detailed feature detection

---

## 2. Current MAFFT Usage Analysis

### 2.1 Usage Locations

MAFFT is currently used in **one primary location**:

**File:** `src/cfeintact/aligned_sequence.py`
```python
@cached_property
def alignment(self) -> MultipleSeqAlignment:
    return wrappers.mafft([self.reference, self.this])
```

**Wrapper Implementation:** `src/cfeintact/wrappers.py`
```python
def mafft(sequences: Iterable[SeqRecord]) -> MultipleSeqAlignment:
    with tempfile.NamedTemporaryFile() as alignment_input, tempfile.NamedTemporaryFile() as alignment_output:
        SeqIO.write(sequences, alignment_input.name, "fasta")
        subprocess.run(["mafft", "--quiet", alignment_input.name],
                       shell=False, stdout=alignment_output, check=True)
        alignment: MultipleSeqAlignment = AlignIO.read(alignment_output.name, "fasta")
        return alignment
```

### 2.2 Methodology

**Alignment Type:** Global pairwise alignment  
**Input:** Two sequences (reference and query)  
**Output:** MultipleSeqAlignment object from BioPython  
**Caching:** Results cached via `@cached_property` decorator

**Call Hierarchy:**
```
AlignedSequence.alignment (property)
  └─> wrappers.mafft()
       └─> subprocess.run(["mafft", ...])
```

### 2.3 Usage Patterns

The `AlignedSequence` class is instantiated in multiple contexts:

1. **Subtype-to-HXB2 alignment** (Line 773 of intact.py):
   ```python
   aligned_subtype = AlignedSequence(this=reference, reference=st.HXB2())
   ```

2. **Query-to-Subtype alignment** (Lines 775-776):
   ```python
   forward_aligned_sequence = AlignedSequence(this=sequence, reference=aligned_subtype.this)
   reverse_aligned_sequence = forward_aligned_sequence.reverse()
   ```

3. **Coordinate mapping** for feature localization (PSI, RRE, MSD)

### 2.4 Downstream Dependencies

The alignment result is consumed by:

1. **Cigar string generation** (aligned_sequence.py:36):
   ```python
   return Cigar.from_msa(reference=seq_reference, query=seq_query)
   ```

2. **Coordinate mapping** for translating HXB2 positions to query positions

3. **Feature detection**:
   - Major Splice Donor site (MSD) - Lines 268-281 of intact.py
   - Packaging Signal (PSI) - Lines 284-338
   - Rev Response Element (RRE) - Lines 341-377
   - Hypermutation analysis - Lines 188-241

4. **ORF analysis** via `check_reading_frame()`

### 2.5 Purpose

MAFFT serves three critical functions:

1. **Global Alignment**: Provides end-to-end alignment with gap handling
2. **Coordinate Translation**: Maps HXB2 reference coordinates to query positions
3. **Feature Preservation**: Maintains structural features for regulatory element detection

### 2.6 Current Advantages

✅ **Accurate global alignment** with sophisticated gap penalties  
✅ **Standardized output** compatible with BioPython ecosystem  
✅ **Well-tested** with established behavior in production  
✅ **Multiple sequence alignment** support (though only 2 sequences used)  
✅ **Established installation** via apt/brew packages

### 2.7 Current Disadvantages

❌ **External dependency** requiring system installation  
❌ **Subprocess overhead** with file I/O for each alignment  
❌ **Performance**: Slower than modern lightweight aligners  
❌ **Complexity**: Full MSA algorithm for pairwise tasks  
❌ **Deployment**: Docker dependency for MAFFT installation

---

## 3. Mappy (Minimap2) Analysis

### 3.1 Overview

**Mappy** is a Python binding for minimap2, a versatile sequence alignment tool designed by Heng Li. It's optimized for:
- Long read alignment (PacBio, Nanopore)
- Assembly-to-assembly alignment
- Short read mapping
- Splice-aware alignment

### 3.2 Technical Capabilities

**Strengths:**
- ⚡ **Performance**: 10-50x faster than traditional aligners for long sequences
- 📦 **Pure Python binding**: No external binaries required
- 🔧 **Flexible presets**: Configurable for various alignment scenarios
- 💾 **Low memory footprint**: Efficient k-mer indexing
- 🎯 **Native PAF/CIGAR output**: Direct alignment representation

**Limitations:**
- ⚠️ **Heuristic-based**: May miss optimal alignment in highly divergent regions
- ⚠️ **Limited MSA support**: Primarily designed for pairwise alignment
- ⚠️ **Different output format**: Requires adaptation layer
- ⚠️ **Gap scoring**: Less configurable than traditional aligners

### 3.3 API Overview

```python
import mappy

# Index reference sequence
aligner = mappy.Aligner(seq=reference_seq, preset="map-ont")  # or "ava-ont", "asm5", etc.

# Align query
for hit in aligner.map(query_seq):
    print(hit.q_st, hit.q_en)  # Query start/end
    print(hit.r_st, hit.r_en)  # Reference start/end
    print(hit.cigar_str)       # CIGAR string
    print(hit.strand)          # Strand orientation
```

### 3.4 Preset Options

Relevant presets for CFEIntact:
- `map-ont`: Long noisy reads (Oxford Nanopore)
- `asm5`: Assembly-to-assembly (≤5% divergence)
- `asm10`: Assembly-to-assembly (≤10% divergence)
- `asm20`: Assembly-to-assembly (≤20% divergence)
- `splice`: Splice-aware alignment
- `sr`: Short reads

**Recommendation**: `asm10` or `asm20` for HIV sequence alignment (considering ~5-15% divergence)

---

## 4. Critical Challenge: Major Splice Donor Detection

### 4.1 Current Implementation

The MSD detection relies on precise alignment-based coordinate mapping:

```python
# From intact.py:268-281
def has_mutated_major_splice_donor_site(alignment, splice_donor_start_pos, splice_donor_end_pos, splice_donor_sequence):
    # Find position in alignment where non-gap characters reach the reference position
    sd_begin = [m.start() for m in re.finditer(r"[^-]", str(alignment[0].seq))][splice_donor_start_pos]
    sd_end = [m.start() for m in re.finditer(r"[^-]", str(alignment[0].seq))][splice_donor_end_pos]
    
    # Extract corresponding sequence from query
    sd = alignment[1].seq[sd_begin:(sd_end + 1)]
    context = alignment[1].seq[(sd_begin - 10):(sd_end + 1 + 10)]
    
    if sd.upper() != splice_donor_sequence.upper():
        return Defect(...)
```

**Key Requirements:**
1. Precise positional mapping from reference to query
2. Access to aligned sequences with gap characters
3. Context extraction around target positions

### 4.2 Challenge with Mappy

Mappy outputs CIGAR strings without explicit gap-aligned sequences:
```
Example: 100M2D50M3I75M
```

This requires:
- ✅ CIGAR parsing to reconstruct alignment
- ✅ Manual position tracking through CIGAR operations
- ⚠️ Additional implementation complexity

### 4.3 Solution Strategy

**Option A: CIGAR-based reconstruction**
```python
def map_position_via_cigar(cigar_str, ref_pos):
    ref_idx, query_idx = 0, 0
    for length, op in parse_cigar(cigar_str):
        if op in 'M=X':  # Match/mismatch
            if ref_idx + length > ref_pos:
                return query_idx + (ref_pos - ref_idx)
            ref_idx += length
            query_idx += length
        elif op == 'D':  # Deletion
            ref_idx += length
        elif op == 'I':  # Insertion
            query_idx += length
    return None
```

**Option B: Hybrid approach** (RECOMMENDED)
- Use mappy for initial alignment and orientation detection
- Fall back to BioPython's PairwiseAligner for precise feature mapping
- Best of both worlds: performance + accuracy

---

## 5. Migration Plan

### 5.1 Recommended Approach: **Hybrid Implementation**

**Phase 1: Parallel Implementation**
1. Add mappy as dependency
2. Implement `MappyAlignedSequence` class
3. Add feature flag for mappy vs MAFFT
4. Run both aligners in test suite

**Phase 2: Feature Parity**
1. Implement CIGAR-based coordinate mapping
2. Validate MSD detection accuracy
3. Benchmark performance differences
4. Validate all 932 lines of intact.py logic

**Phase 3: Gradual Migration**
1. Use mappy for orientation detection
2. Use mappy for ORF finding (less position-critical)
3. Keep MAFFT/BioPython for precise feature detection
4. Monitor production metrics

**Phase 4: Full Migration (Optional)**
1. Once confidence is high, make mappy default
2. Keep MAFFT as fallback option
3. Update documentation

### 5.2 Implementation Components

#### Component 1: Mappy Wrapper
```python
# src/cfeintact/mappy_aligner.py
import mappy

class MappyAligner:
    def __init__(self, preset="asm20"):
        self.preset = preset
    
    def align(self, reference: SeqRecord, query: SeqRecord) -> MappyAlignment:
        aligner = mappy.Aligner(seq=str(reference.seq), preset=self.preset)
        hits = list(aligner.map(str(query.seq)))
        
        if not hits:
            raise ValueError("No alignment found")
        
        # Take best hit
        best_hit = max(hits, key=lambda h: h.mlen)
        return MappyAlignment(best_hit, reference, query)
```

#### Component 2: Coordinate Mapper
```python
class MappyCoordinateMapper:
    def __init__(self, cigar_str, ref_start, query_start):
        self.cigar = self.parse_cigar(cigar_str)
        self.ref_start = ref_start
        self.query_start = query_start
    
    def ref_to_query(self, ref_pos):
        # Implementation using CIGAR operations
        pass
```

#### Component 3: Compatibility Layer
```python
class MappyAlignedSequence(AlignedSequence):
    """Drop-in replacement for MAFFT-based AlignedSequence"""
    
    @cached_property
    def alignment(self) -> MultipleSeqAlignment:
        # Convert mappy hit to MSA format for compatibility
        pass
    
    @cached_property
    def coordinate_mapping(self):
        # Return aligntools-compatible mapping
        pass
```

### 5.3 Testing Strategy

1. **Unit Tests**
   - Test coordinate mapping accuracy
   - Validate CIGAR parsing
   - Test edge cases (no alignment, multiple hits)

2. **Integration Tests**
   - Run entire test suite with mappy
   - Compare outputs: `expected-results-*` directories
   - Validate bit-for-bit compatibility where needed

3. **Performance Tests**
   - Benchmark alignment speed
   - Memory usage comparison
   - Large dataset validation

4. **Regression Tests**
   - Ensure MSD detection accuracy
   - Validate PSI and RRE detection
   - Check all defect types

---

## 6. Risk Assessment

### 6.1 High Risks

| Risk | Impact | Mitigation |
|------|--------|------------|
| **MSD detection inaccuracy** | High - Core feature failure | Implement exhaustive validation; keep MAFFT fallback |
| **Changed alignment quality** | High - Different defect detection | A/B testing; gradual rollout |
| **CIGAR parsing bugs** | Medium - Wrong coordinates | Extensive unit tests; edge case coverage |
| **Performance regression** | Low - Unlikely but possible | Benchmark before merge |

### 6.2 Migration Complexity

- **Code changes required**: ~5 new files, ~3 modified files
- **Test updates**: All integration tests need validation
- **Documentation updates**: Installation instructions, architecture docs
- **Estimated effort**: 40-60 hours (including testing)

---

## 7. Advantages and Disadvantages Summary

### 7.1 Advantages of Migration

✅ **Performance**: 10-50x faster alignment  
✅ **Reduced dependencies**: No external binary (pip installable)  
✅ **Modern codebase**: Active development, better maintained  
✅ **Simplified deployment**: Easier Docker builds  
✅ **Memory efficiency**: Better for large-scale processing  
✅ **Future-proof**: Industry-standard aligner

### 7.2 Disadvantages of Migration

❌ **Reimplementation effort**: Coordinate mapping logic  
❌ **Different alignment characteristics**: May affect edge cases  
❌ **Testing burden**: Need to validate all features  
❌ **Learning curve**: Team needs to understand mappy API  
❌ **Potential accuracy trade-offs**: Heuristic vs. exact alignment  
❌ **Breaking change risk**: Could affect existing workflows

---

## 8. Recommendations

### 8.1 Primary Recommendation: **Hybrid Approach**

Implement mappy as the primary aligner while retaining MAFFT/BioPython capabilities for critical features:

```python
USE_MAPPY = True  # Feature flag

if USE_MAPPY:
    aligned = MappyAlignedSequence(this=sequence, reference=reference)
else:
    aligned = AlignedSequence(this=sequence, reference=reference)
```

**Benefits:**
- ✅ Performance gains where safe
- ✅ Maintains accuracy for critical features
- ✅ Gradual validation possible
- ✅ Fallback option always available

### 8.2 Alternative: **Status Quo**

Keep MAFFT if:
- Performance is already acceptable
- Risk tolerance is low
- Development resources are limited

**Note**: MAFFT is battle-tested and working well. Migration should be driven by clear performance needs, not just novelty.

### 8.3 Timeline Recommendation

- **Week 1-2**: Implement mappy wrapper and coordinate mapper
- **Week 3-4**: Integration and initial testing
- **Week 5-6**: Validation against test suite
- **Week 7-8**: Performance benchmarking and optimization
- **Week 9-10**: Documentation and gradual rollout

---

## 9. Technical Specifications for Implementation

### 9.1 Dependencies to Add

```toml
# pyproject.toml
dependencies = [
    "biopython>=1.83",
    "click>=8.0",
    "scipy>=1.1",
    "numpy==2.4.1",
    "aligntools==1.2.2",
    "mappy>=2.24",  # NEW
]
```

### 9.2 Files to Create

1. `src/cfeintact/mappy_aligner.py` - Core mappy wrapper
2. `src/cfeintact/mappy_aligned_sequence.py` - Compatibility layer
3. `src/cfeintact/cigar_parser.py` - CIGAR string utilities
4. `tests/test_mappy_alignment.py` - Unit tests
5. `tests/test_mappy_coordinate_mapping.py` - Coordinate mapping tests

### 9.3 Files to Modify

1. `src/cfeintact/aligned_sequence.py` - Add factory method
2. `src/cfeintact/intact.py` - Add mappy usage option
3. `src/cfeintact/main.py` - Add CLI flag
4. `Dockerfile` - Remove mafft installation (optional)
5. `.github/workflows/main.yml` - Remove mafft from CI

### 9.4 Configuration Options

Add CLI flags:
```bash
--use-mappy / --use-mafft  # Aligner selection
--mappy-preset PRESET      # asm5, asm10, asm20, etc.
```

---

## 10. Conclusion

The migration from MAFFT to mappy is **technically feasible** but presents **significant implementation challenges**, particularly around precise feature detection like the Major Splice Donor site. 

**Key Takeaway**: A **hybrid approach** offers the best balance:
- Leverage mappy's performance for general alignment tasks
- Retain proven methods for position-critical features
- Minimize risk while enabling future optimization

The CFEIntact codebase is well-structured for this migration, with clear abstraction boundaries (`AlignedSequence` class) that facilitate a gradual, safe transition.

**Final Recommendation**: Proceed with hybrid implementation, validate extensively, and migrate incrementally based on production data and performance metrics.

---

## Appendix A: MAFFT vs Mappy Comparison Table

| Feature | MAFFT | Mappy |
|---------|-------|-------|
| **Speed (10kb seqs)** | ~5-10 seconds | ~0.1-0.5 seconds |
| **Installation** | System package | pip install |
| **Output format** | MSA (FASTA) | CIGAR/PAF |
| **Gap handling** | Sophisticated | Heuristic |
| **Accuracy** | Very high | High |
| **Memory usage** | Moderate | Low |
| **Multiple alignment** | Native | Pairwise only |
| **Maintenance** | Mature | Active |
| **Python integration** | Subprocess | Native binding |

---

## Appendix B: Key Code Locations Reference

- **MAFFT wrapper**: `src/cfeintact/wrappers.py:14-38`
- **AlignedSequence class**: `src/cfeintact/aligned_sequence.py:13-56`
- **MSD detection**: `src/cfeintact/intact.py:268-281`
- **Coordinate mapping usage**: Throughout `src/cfeintact/intact.py`
- **ORF finding**: `src/cfeintact/find_orf.py`
- **Secondary aligner**: `src/cfeintact/detailed_aligner.py` (BioPython PairwiseAligner)

---

**Report End**
