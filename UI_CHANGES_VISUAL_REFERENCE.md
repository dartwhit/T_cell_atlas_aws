# UI Changes - Visual Reference

## Before Changes

```
┌─────────────────────────────────────┐
│ Differential Expression setup       │
├─────────────────────────────────────┤
│                                     │
│ Cell Cluster: [Fibroblasts ▼]      │
│                                     │
│ ☐ Compare by disease                │
│                                     │
│ [Show gene in UMAP space]           │
│                                     │
└─────────────────────────────────────┘
```

**Problems:**
- Label "Compare by disease" is confusing for Khanna (which uses timepoints)
- No indication of what comparison is actually being shown
- Users don't know if they're seeing "SSc vs Healthy" or "After vs Before Treatment"

---

## After Changes - Disease Datasets (Tabib, Gur, Ma, TMKMH)

### When Checkbox is UNCHECKED (showing cluster markers):
```
┌─────────────────────────────────────┐
│ Differential Expression setup       │
├─────────────────────────────────────┤
│                                     │
│ Cell Cluster: [Fibroblasts ▼]      │
│                                     │
│ ☐ Compare by disease                │
│                                     │
│ ┌─────────────────────────────────┐ │
│ │ 🟠 Showing:                     │ │
│ │    Cluster markers (all cells)  │ │
│ └─────────────────────────────────┘ │
│                                     │
│ [Show gene in UMAP space]           │
│                                     │
└─────────────────────────────────────┘
```

### When Checkbox is CHECKED (showing disease comparison):
```
┌─────────────────────────────────────┐
│ Differential Expression setup       │
├─────────────────────────────────────┤
│                                     │
│ Cell Cluster: [Fibroblasts ▼]      │
│                                     │
│ ☑ Compare by disease                │
│                                     │
│ ┌─────────────────────────────────┐ │
│ │ 🔵 Showing:                     │ │
│ │    SSc vs Healthy               │ │
│ └─────────────────────────────────┘ │
│                                     │
│ [Show gene in UMAP space]           │
│                                     │
└─────────────────────────────────────┘
```

---

## After Changes - Khanna Dataset

### When Checkbox is UNCHECKED (showing cluster markers):
```
┌─────────────────────────────────────┐
│ Differential Expression setup       │
├─────────────────────────────────────┤
│                                     │
│ Cell Cluster: [Fibroblasts ▼]      │
│                                     │
│ ☐ Compare by timepoint              │
│                                     │
│ ┌─────────────────────────────────┐ │
│ │ 🟠 Showing:                     │ │
│ │    Cluster markers (all cells)  │ │
│ └─────────────────────────────────┘ │
│                                     │
│ [Show gene in UMAP space]           │
│                                     │
└─────────────────────────────────────┘
```

### When Checkbox is CHECKED (showing timepoint comparison):
```
┌─────────────────────────────────────┐
│ Differential Expression setup       │
├─────────────────────────────────────┤
│                                     │
│ Cell Cluster: [Fibroblasts ▼]      │
│                                     │
│ ☑ Compare by timepoint              │
│                                     │
│ ┌─────────────────────────────────┐ │
│ │ 🔵 Showing:                     │ │
│ │    After vs Before Treatment    │ │
│ └─────────────────────────────────┘ │
│                                     │
│ [Show gene in UMAP space]           │
│                                     │
└─────────────────────────────────────┘
```

---

## Key UI Improvements

1. **Dynamic Checkbox Label**
   - Changes from "Compare by disease" to "Compare by timepoint" for Khanna
   - Accurately describes what the comparison represents

2. **Comparison Information Box**
   - **Orange (🟠)** when showing cluster markers
     - Indicates standard differential expression between clusters
   - **Blue (🔵)** when showing condition comparison
     - Shows exact comparison being performed
     - "SSc vs Healthy" for disease datasets
     - "After vs Before Treatment" for Khanna

3. **Clarity**
   - Users immediately know what data they're viewing
   - No confusion about what "disease" means in different contexts
   - Clear distinction between cluster markers and condition comparisons

4. **Visual Hierarchy**
   - Info box uses color coding for quick recognition
   - Positioned logically below the checkbox that controls it
   - Styled with border and padding for clear separation

---

## Color Scheme

- **Orange (#FF9800 / #F57C00)**: Standard view (cluster markers)
- **Blue (#2196F3 / #1976D2)**: Comparison view (condition-based DE)
- Background colors provide subtle distinction
- Border-left provides strong visual indicator

---

## Responsive to Dataset Changes

When user changes dataset:
1. Checkbox label updates immediately
2. Comparison info box updates to show new comparison label
3. If comparison data not available, info box is hidden
4. Smooth transition without page reload
