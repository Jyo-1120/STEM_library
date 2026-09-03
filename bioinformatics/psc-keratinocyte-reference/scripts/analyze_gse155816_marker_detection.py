#!/usr/bin/env python3

import csv
import gzip
import math
from pathlib import Path

import numpy as np


ROOT = Path.cwd()
DATA = ROOT / "data" / "public_reference" / "GSE155816" / "extracted"
OUT = ROOT / "output" / "psc_kc_reference" / "GSE155816_marker_validation"
OUT.mkdir(parents=True, exist_ok=True)

SAMPLES = {
    "K82_P2": DATA / "GSM4712539_K82_counts_matrix.txt.gz",
    "K86_P2": DATA / "GSM4712540_K86_counts_matrix.txt.gz",
    "K82_P5": DATA / "GSM4824577_K82_t2_counts_matrix.txt.gz",
    "K86_P5": DATA / "GSM4824578_K86_t2_counts_matrix.txt.gz",
}

META = {
    "K82_P2": ("K82", "P2", 3367),
    "K86_P2": ("K86", "P2", 3978),
    "K82_P5": ("K82", "P5", 4292),
    "K86_P5": ("K86", "P5", 3576),
}

MODULES = {
    "Basal": ["TP63", "KRT5", "KRT14", "ITGA6", "ITGB4", "COL17A1"],
    "Holoclone_cycle": ["FOXM1", "ANLN", "AURKB", "CCNA2", "CKAP2L", "HMGB2", "LMNB1"],
    "General_cycle": ["MKI67", "TOP2A", "TYMS", "UBE2C", "CENPF"],
    "Differentiating": ["KRT1", "KRT10", "IVL", "DSG1", "TGM1", "CLDN1", "SPINK5"],
    "Terminal": ["FLG", "LOR", "ABCA12", "SPRR1A", "SPRR1B"],
    "Mesenchymal": ["VIM", "COL1A1", "COL1A2", "DCN", "LUM"],
}
TARGET_GENES = sorted({g for genes in MODULES.values() for g in genes})


def mad(x):
    med = np.median(x)
    return np.median(np.abs(x - med))


def process_sample(sample, path):
    print(f"Streaming {sample}: {path.name}", flush=True)
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cells = header[1:]
        n = len(cells)
        total = np.zeros(n, dtype=np.int64)
        detected = np.zeros(n, dtype=np.int32)
        mitochondrial = np.zeros(n, dtype=np.int64)
        target = {gene: np.zeros(n, dtype=np.int32) for gene in TARGET_GENES}

        for line_no, line in enumerate(handle, start=2):
            gene, values = line.rstrip("\n").split("\t", 1)
            arr = np.fromstring(values, sep="\t", dtype=np.int32)
            if arr.size != n:
                raise ValueError(
                    f"{path.name}:{line_no} has {arr.size} values; expected {n}"
                )
            total += arr
            detected += arr > 0
            if gene.startswith("MT-"):
                mitochondrial += arr
            if gene in target:
                target[gene] = arr.copy()

    mt_pct = np.divide(
        mitochondrial * 100.0,
        total,
        out=np.zeros(n, dtype=float),
        where=total > 0,
    )

    # The paper used sample-specific outlier removal but did not publish numeric
    # cutoffs. This transparent screen is for sensitivity analysis only.
    upper_count = np.median(total) + 4 * 1.4826 * mad(total)
    upper_feature = np.median(detected) + 4 * 1.4826 * mad(detected)
    qc_keep = (
        (detected > 500)
        & (detected < max(7500, upper_feature))
        & (total > 500)
        & (total < upper_count)
        & (mt_pct < 20)
    )

    donor, passage, paper_cells = META[sample]
    qc_row = {
        "sample": sample,
        "donor": donor,
        "passage": passage,
        "matrix_barcodes": n,
        "paper_reported_cells": paper_cells,
        "transparent_qc_cells": int(qc_keep.sum()),
        "median_features": float(np.median(detected)),
        "median_UMIs": float(np.median(total)),
        "median_mt_pct": float(np.median(mt_pct)),
        "upper_count_cutoff": float(upper_count),
        "upper_feature_cutoff": float(max(7500, upper_feature)),
    }

    gene_rows = []
    for gene in TARGET_GENES:
        values = target[gene]
        for label, keep in (("all_matrix_barcodes", np.ones(n, dtype=bool)),
                            ("transparent_qc", qc_keep)):
            denom = int(keep.sum())
            gene_rows.append({
                "sample": sample,
                "donor": donor,
                "passage": passage,
                "cell_set": label,
                "gene": gene,
                "n_cells": denom,
                "detected_cells": int(((values > 0) & keep).sum()),
                "detected_pct": float(((values > 0) & keep).sum() * 100 / denom)
                if denom else math.nan,
                "total_UMIs": int(values[keep].sum()),
                "mean_UMI_per_cell": float(values[keep].mean()) if denom else math.nan,
            })

    module_rows = []
    normalized = {
        gene: np.log1p(
            np.divide(values * 10000.0, total,
                      out=np.zeros(n, dtype=float), where=total > 0)
        )
        for gene, values in target.items()
    }
    module_scores = {}
    for module, genes in MODULES.items():
        module_scores[module] = np.mean([normalized[g] for g in genes], axis=0)
        present = np.sum([target[g] > 0 for g in genes], axis=0)
        for label, keep in (("all_matrix_barcodes", np.ones(n, dtype=bool)),
                            ("transparent_qc", qc_keep)):
            denom = int(keep.sum())
            module_rows.append({
                "sample": sample,
                "donor": donor,
                "passage": passage,
                "cell_set": label,
                "module": module,
                "n_cells": denom,
                "module_detected_pct": float((present[keep] > 0).mean() * 100)
                if denom else math.nan,
                "mean_genes_detected": float(present[keep].mean()) if denom else math.nan,
                "mean_log_normalized_module_score": float(module_scores[module][keep].mean())
                if denom else math.nan,
            })

    cell_rows = []
    for i, cell in enumerate(cells):
        row = {
            "cell": f"{sample}_{cell}",
            "sample": sample,
            "donor": donor,
            "passage": passage,
            "nFeature": int(detected[i]),
            "nUMI": int(total[i]),
            "percent_mt": float(mt_pct[i]),
            "transparent_qc_keep": bool(qc_keep[i]),
        }
        row.update({f"module_{k}": float(v[i]) for k, v in module_scores.items()})
        cell_rows.append(row)
    return qc_row, gene_rows, module_rows, cell_rows


all_qc, all_genes, all_modules, all_cells = [], [], [], []
for sample, path in SAMPLES.items():
    qc, genes, modules, cells = process_sample(sample, path)
    all_qc.append(qc)
    all_genes.extend(genes)
    all_modules.extend(modules)
    all_cells.extend(cells)


def write_csv(path, rows, compressed=False):
    opener = gzip.open if compressed else open
    kwargs = {"mode": "wt", "newline": ""} if compressed else {"mode": "w", "newline": ""}
    with opener(path, **kwargs) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


write_csv(OUT / "sample_qc_sensitivity.csv", all_qc)
write_csv(OUT / "marker_detection_by_sample.csv", all_genes)
write_csv(OUT / "module_scores_by_sample.csv", all_modules)
write_csv(OUT / "cell_qc_and_module_scores.csv.gz", all_cells, compressed=True)
print(f"Finished: {OUT}", flush=True)
