import csv
import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# SETTINGS
# ============================================================

FOLDER = Path(__file__).resolve().parent
OUTPUT = FOLDER / "final_results"
OUTPUT.mkdir(exist_ok=True)

FILES = [
    FOLDER / f"output-Set{i}.txt"
    for i in range(1, 11)
]

EXPECTED_EVENTS = 5_000_000
K = 10

# Names and IDs supplied in the course.
SPECIES = {
    211: "Carbon monoxide",
    -211: "Carbon-13 monoxide",
    321: "Nitric oxide",
    -321: "Ionised NO",
    2212: "Water",
    -2212: "Heavy water",
    3122: "Methane",
    -3122: "Methyl ion",
    3312: "Ethylene",
    -3312: "Ionised ethylene",
    3334: "Ozone",
    -3334: "Superoxide anion",
}

PAIRS = [211, 321, 2212, 3122, 3312, 3334]
IDS = list(SPECIES)
INDEX = {code: i for i, code in enumerate(IDS)}

# Momentum is expressed in the units supplied by the course:
# 10^-23 kg m/s.
EDGES = np.array([
    0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 6, 8, np.inf
])

LABELS = [
    "0–0.5", "0.5–1", "1–1.5", "1.5–2", "2–2.5",
    "2.5–3", "3–4", "4–6", "6–8", "≥8"
]


# ============================================================
# FUNCTIONS
# ============================================================

def momentum(px, py, pz):
    """Return the magnitude of the momentum vector."""
    return math.hypot(px, py, pz)


def asymmetry_statistics(positive, negative):
    """
    Calculate the pooled asymmetry and subsampling uncertainty.

    Each input contains one count per subsample.

    SD measures variation between subsample asymmetries.
    SE estimates uncertainty on the full-sample asymmetry.
    """
    positive = np.asarray(positive, dtype=float)
    negative = np.asarray(negative, dtype=float)

    combined = positive + negative

    if combined.sum() == 0:
        return math.nan, math.nan, math.nan, math.nan

    asymmetry = (
        (positive.sum() - negative.sum()) / combined.sum()
    )

    # Do not treat undefined subsample asymmetries as zero.
    if np.any(combined == 0):
        return asymmetry, math.nan, math.nan, math.nan

    subsample_asymmetries = (
        (positive - negative) / combined
    )

    sd = np.std(subsample_asymmetries, ddof=1)
    se = sd / math.sqrt(len(subsample_asymmetries))

    z = asymmetry / se if se > 0 else math.nan

    return asymmetry, sd, se, z


def save_csv(filename, rows):
    """Save a list of dictionaries as a CSV table."""
    with (OUTPUT / filename).open(
        "w", newline="", encoding="utf-8"
    ) as file:
        writer = csv.DictWriter(
            file, fieldnames=list(rows[0])
        )
        writer.writeheader()
        writer.writerows(rows)


def markdown_table(headers, rows):
    """Create a Markdown table for the README."""
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]

    for row in rows:
        lines.append(
            "| " + " | ".join(str(value) for value in row) + " |"
        )

    return "\n".join(lines)


def passes_three_se(z):
    """Use a strict three-standard-error screening threshold."""
    return math.isfinite(z) and abs(z) > 3


# ============================================================
# READ ALL TEN FILES
# ============================================================

def read_data():
    for path in FILES:
        if not path.is_file():
            raise FileNotFoundError(
                f"Cannot find {path}\n"
                "Keep output-Set1.txt through output-Set10.txt "
                "in the same folder as this script."
            )

    # Dimensions: subsample, species, momentum bin.
    counts = np.zeros((K, len(IDS), len(LABELS)), dtype=np.int64)

    # A second partition uses each input file as one subsample.
    # This checks sensitivity to the choice of subsamples.
    file_counts = np.zeros_like(counts)

    subsample_sizes = np.zeros(K, dtype=np.int64)
    file_sizes = np.zeros(K, dtype=np.int64)

    events = 0
    molecules = 0
    empty_events = 0

    for file_index, path in enumerate(FILES):
        print(f"\nReading {path.name}...", flush=True)

        with path.open("r") as file:
            numbered_lines = enumerate(file, start=1)

            for line_number, line in numbered_lines:
                if not line.strip():
                    continue

                # Each event begins with: event_ID molecule_count.
                header = line.split()

                if len(header) != 2:
                    raise ValueError(
                        f"{path.name}, line {line_number}: "
                        "expected an event header with two values."
                    )

                event_id = int(header[0])
                number_of_molecules = int(header[1])

                if number_of_molecules < 0:
                    raise ValueError(
                        f"{path.name}: negative molecule count."
                    )

                # Round-robin assignment keeps complete events together.
                group = events % K
                read = 0

                while read < number_of_molecules:
                    try:
                        line_number, line = next(numbered_lines)
                    except StopIteration:
                        raise ValueError(
                            f"{path.name}: event {event_id} is incomplete."
                        ) from None

                    if not line.strip():
                        continue

                    values = line.split()

                    if len(values) != 4:
                        raise ValueError(
                            f"{path.name}, line {line_number}: "
                            "expected px py pz molecule_ID."
                        )

                    px, py, pz = map(float, values[:3])
                    molecule_id = int(values[3])

                    if not all(
                        math.isfinite(value) for value in (px, py, pz)
                    ):
                        raise ValueError(
                            f"{path.name}, line {line_number}: "
                            "invalid momentum."
                        )

                    # Select only the twelve requested species.
                    if molecule_id in INDEX:
                        p = momentum(px, py, pz)

                        if not math.isfinite(p):
                            raise ValueError("Invalid momentum magnitude.")

                        bin_index = int(
                            np.searchsorted(EDGES, p, side="right") - 1
                        )

                        species_index = INDEX[molecule_id]

                        counts[group, species_index, bin_index] += 1
                        file_counts[
                            file_index, species_index, bin_index
                        ] += 1

                    read += 1

                # Empty events are included in the event denominator.
                events += 1
                molecules += number_of_molecules
                empty_events += number_of_molecules == 0

                subsample_sizes[group] += 1
                file_sizes[file_index] += 1

                if events % 500_000 == 0:
                    print(
                        f"Processed {events:,} events...",
                        flush=True
                    )

    if events != EXPECTED_EVENTS:
        raise ValueError(
            f"Read {events:,} events; expected {EXPECTED_EVENTS:,}. "
            "Check that all ten complete files are present."
        )

    if not np.all(subsample_sizes == subsample_sizes[0]):
        raise ValueError("The primary subsamples are not equal in size.")

    if not np.all(file_sizes == file_sizes[0]):
        raise ValueError(
            "The file-partition comparison requires equal event counts."
        )

    if not np.array_equal(
        counts.sum(axis=0), file_counts.sum(axis=0)
    ):
        raise ValueError("The two counting methods disagree.")

    return (
        counts, file_counts, subsample_sizes,
        events, molecules, empty_events
    )


# ============================================================
# CALCULATE RESULTS
# ============================================================

def analyse(counts, file_counts, sizes, events):
    species_rows = []
    integrated_rows = []
    bin_rows = []

    totals = counts.sum(axis=2)
    file_totals = file_counts.sum(axis=2)

    for code in IDS:
        i = INDEX[code]

        subsample_means = totals[:, i] / sizes
        sd = np.std(subsample_means, ddof=1)
        se = sd / math.sqrt(K)

        species_rows.append({
            "ID": code,
            "Species": SPECIES[code],
            "Total": int(totals[:, i].sum()),
            "Mean_per_event": totals[:, i].sum() / events,
            "Subsample_SD": sd,
            "Full_sample_SE": se,
        })

    for code in PAIRS:
        ip = INDEX[code]
        im = INDEX[-code]

        a, sd, se, z = asymmetry_statistics(
            totals[:, ip], totals[:, im]
        )

        _, file_sd, file_se, file_z = asymmetry_statistics(
            file_totals[:, ip], file_totals[:, im]
        )

        integrated_rows.append({
            "Positive_ID": code,
            "Negative_ID": -code,
            "Asymmetry": a,
            "Subsample_SD": sd,
            "Full_sample_SE": se,
            "Approx_Z": z,
            "A_over_SD": a / sd if sd > 0 else math.nan,
            "File_partition_SD": file_sd,
            "File_partition_SE": file_se,
            "File_partition_Z": file_z,
        })

        for b, label in enumerate(LABELS):
            positive = counts[:, ip, b]
            negative = counts[:, im, b]

            a, sd, se, z = asymmetry_statistics(positive, negative)

            bin_rows.append({
                "Positive_ID": code,
                "Negative_ID": -code,
                "Momentum_bin": label,
                "Positive_count": int(positive.sum()),
                "Negative_count": int(negative.sum()),
                "Asymmetry": a,
                "Subsample_SD": sd,
                "Full_sample_SE": se,
                "Approx_Z": z,
                "Minimum_pair_count_per_subsample": int(
                    (positive + negative).min()
                ),
            })

    return species_rows, integrated_rows, bin_rows


# ============================================================
# CREATE THE REQUIRED MOMENTUM-ASYMMETRY PLOT
# ============================================================

def plot_asymmetries(integrated_rows, bin_rows):
    fig, axes = plt.subplots(
        3, 2, figsize=(14, 13), layout="constrained"
    )

    x = np.arange(len(LABELS))

    for ax, code, overall in zip(
        axes.flat, PAIRS, integrated_rows
    ):
        rows = [
            row for row in bin_rows
            if row["Positive_ID"] == code
        ]

        # Multiply by 100 to display asymmetry as a percentage.
        a = np.array([r["Asymmetry"] for r in rows]) * 100
        se = np.array([r["Full_sample_SE"] for r in rows]) * 100

        valid = np.isfinite(a) & np.isfinite(se)

        overall_a = overall["Asymmetry"] * 100
        overall_se = overall["Full_sample_SE"] * 100

        ax.axhline(
            0, color="gray", linestyle="--", linewidth=1
        )

        ax.axhspan(
            overall_a - overall_se,
            overall_a + overall_se,
            color="orange",
            alpha=0.18,
            label="Integrated A ± 1 SE",
        )

        ax.axhline(
            overall_a, color="darkorange", linewidth=1
        )

        ax.errorbar(
            x[valid],
            a[valid],
            yerr=se[valid],
            fmt="o",
            color="tab:blue",
            capsize=3,
            markersize=4,
            label="Bin A ± 1 SE",
        )

        ax.set_title(
            f"{SPECIES[code]} / {SPECIES[-code]}\n"
            f"IDs {code} / {-code}"
        )

        ax.set_xticks(x, LABELS, rotation=40, ha="right")
        ax.set_xlabel("Momentum interval (10⁻²³ kg m/s)")
        ax.set_ylabel("Asymmetry (%)")
        ax.grid(axis="y", alpha=0.2)
        ax.legend(fontsize=8)
        ax.margins(y=0.15)

    fig.suptitle(
        "Molecular asymmetry by momentum\n"
        "5 million events · 10 subsamples · "
        "independent vertical scales",
        fontsize=15,
    )

    fig.savefig(OUTPUT / "asymmetry_vs_momentum.png", dpi=200)
    fig.savefig(OUTPUT / "asymmetry_vs_momentum.pdf")
    plt.close(fig)


# ============================================================
# SUPPORTING MOMENTUM SPECTRA
# ============================================================

def plot_spectra(counts, sizes):
    finite_edges = EDGES[:-1]
    widths = np.diff(finite_edges)
    centres = (finite_edges[:-1] + finite_edges[1:]) / 2

    fig, axes = plt.subplots(
        3, 2, figsize=(14, 12), layout="constrained"
    )

    for ax, code in zip(axes.flat, PAIRS):
        for molecule_id, colour, marker in [
            (code, "tab:blue", "o"),
            (-code, "tab:orange", "s"),
        ]:
            # Exclude only the open-ended bin from the density plot:
            # its width is not finite.
            density = (
                counts[:, INDEX[molecule_id], :-1]
                / sizes[:, None]
                / widths
            )

            mean = density.mean(axis=0)
            se = density.std(axis=0, ddof=1) / math.sqrt(K)

            valid = mean > 0

            ax.errorbar(
                centres[valid],
                mean[valid],
                xerr=widths[valid] / 2,
                yerr=se[valid],
                fmt=marker,
                color=colour,
                capsize=2,
                markersize=4,
                label=f"ID {molecule_id}",
            )

        ax.set_title(f"{SPECIES[code]} / {SPECIES[-code]}")
        ax.set_xlabel("Momentum magnitude (10⁻²³ kg m/s)")
        ax.set_ylabel("Mean count / event / momentum unit")
        ax.set_yscale("log")
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8)

    fig.suptitle(
        "Momentum spectra for 0 ≤ p < 8\n"
        "Vertical bars: ±1 SE; horizontal bars: bin widths",
        fontsize=15,
    )

    fig.savefig(OUTPUT / "momentum_spectra.png", dpi=200)
    plt.close(fig)


# ============================================================
# WRITE THE README FROM THE CALCULATED RESULTS
# ============================================================

def create_readme(
    species_rows, integrated_rows, file_checks,
    events, molecules, empty_events
):
    selected = sum(row["Total"] for row in species_rows)

    mean_table = markdown_table(
        ["ID", "Species", "Total", "Mean/event", "SD", "SE"],
        [
            [
                r["ID"],
                r["Species"],
                f'{r["Total"]:,}',
                f'{r["Mean_per_event"]:.7f}',
                f'{r["Subsample_SD"]:.6g}',
                f'{r["Full_sample_SE"]:.6g}',
            ]
            for r in species_rows
        ],
    )

    asymmetry_table = markdown_table(
        [
            "Pair", "A (%)", "SE (percentage points)",
            "SD (percentage points)", "Approx. Z", "A/SD",
            "File-partition Z"
        ],
        [
            [
                f'{r["Positive_ID"]} / {r["Negative_ID"]}',
                f'{100 * r["Asymmetry"]:.5f}',
                f'{100 * r["Full_sample_SE"]:.5f}',
                f'{100 * r["Subsample_SD"]:.5f}',
                f'{r["Approx_Z"]:.2f}',
                f'{r["A_over_SD"]:.2f}',
                f'{r["File_partition_Z"]:.2f}',
            ]
            for r in integrated_rows
        ],
    )

    stable = []
    sensitive = []
    inconclusive = []

    for r in integrated_rows:
        pair = f'{r["Positive_ID"]}/{r["Negative_ID"]}'
        primary_pass = passes_three_se(r["Approx_Z"])
        file_pass = passes_three_se(r["File_partition_Z"])

        if primary_pass and file_pass:
            stable.append(pair)
        elif primary_pass != file_pass:
            sensitive.append(pair)
        else:
            inconclusive.append(pair)

    def pair_list(values):
        return ", ".join(values) if values else "none"

    file_table = markdown_table(
        ["File", "Selected molecules", "Fraction with p ≥ 8 (%)"],
        [
            [
                r["File"],
                f'{r["Selected"]:,}',
                f'{100 * r["Fraction_p_ge_8"]:.4f}',
            ]
            for r in file_checks
        ],
    )

    report = f"""# PRA2003 — Molecular asymmetry analysis

**Author:** Marie Piekarek  
**Final deliverable: part 2**

## 1. Aim

Analyse all five million events, select the twelve course species,
calculate their mean counts per event and integrated asymmetries,
and plot asymmetry against momentum with statistical uncertainties.

The species names follow the course table. Positive and negative IDs
identify paired species; they do not necessarily represent electric charge.

## 2. Dataset

| Quantity | Value |
|---|---:|
| Input files | 10 |
| Events analysed | {events:,} |
| Empty events retained | {empty_events:,} |
| Molecule rows read | {molecules:,} |
| Selected molecule rows | {selected:,} |
| Subsamples | {K} |
| Events per subsample | {events // K:,} |

The input files are `output-Set1.txt` through `output-Set10.txt`.

Each event contains a header with its ID and molecule count, followed
by rows containing `px py pz molecule_ID`.

All events, including empty events and events without a selected
species, remain in the denominator for mean counts per event.

## 3. Statistical method

### Mean counts

For each species:

`Mean count per event = total species count / total events`

Complete events are assigned round-robin to ten equal subsamples.
All molecules belonging to an event stay together.

For ten subsample estimates x_j:

`SD = sqrt(sum((x_j - average(x_j))²) / (K - 1))`

`SE = SD / sqrt(K)`

SD measures the variation between subsample estimates.
SE estimates the statistical uncertainty on the full-sample result.

Both quantities are reported because the course notes refer to
subsample SD. **The plotted error bars and the ± uncertainties
discussed here use full-sample SE.**

### Paired asymmetry

`A = (N_positive - N_negative) / (N_positive + N_negative)`

The central result uses counts pooled over all events.
It is not an average of individual-event ratios.

The asymmetry is also calculated separately within each subsample.
Its uncertainty is estimated as the SD of these ten asymmetries
divided by sqrt(10). This retains correlations between members of
the pair and is an approximate uncertainty for the pooled ratio.

- Positive A: more positive-ID molecules.
- Negative A: more negative-ID molecules.
- A = 0: equal counts.

An asymmetry of 1% means a 50.5% versus 49.5% split of the pair total.

### Sigma screening rule

`Z = A / SE(A)`

The screening rule used here is `abs(Z) > 3`.
Sigma refers to the full-sample SE.

These scores are approximate standardized departures from zero,
not exact Gaussian significances or p-values.
Only ten subsamples are used, and no multiple-testing correction
is applied.

A second calculation treats each input file as one subsample.
This checks sensitivity to the partition; it does not automatically
correct for differences in simulation conditions.

## 4. Mean molecules per event

{mean_table}

## 5. Integrated asymmetries

These results include all momenta.

{asymmetry_table}

The A/SD column is also provided for comparison with the literal
subsample-SD convention in the course notes. It must not be confused
with A/SE.

### Do we see molecular asymmetry?

Pairs passing the approximate three-SE threshold under both
partitions: **{pair_list(stable)}**.

Pairs whose threshold classification changes with the partition:
**{pair_list(sensitive)}**.

Pairs not exceeding the threshold under either partition:
**{pair_list(inconclusive)}**.

Failure to exceed a threshold is not proof that the true asymmetry
is exactly zero. These statements describe statistical evidence
within the supplied simulation, not a demonstrated physical cause.

## 6. Asymmetry as a function of momentum

For every selected molecule:

`p = sqrt(px² + py² + pz²)`

Momentum is expressed in units of 10^-23 kg m/s.

The bin edges are:

`0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 6, 8, infinity`

Each finite bin includes its lower edge and excludes its upper edge.
The last bin contains all p ≥ 8.

Within each bin, the pooled asymmetry and its subsampling uncertainty
are calculated using the same procedure as for the integrated result.

![Asymmetry versus momentum](asymmetry_vs_momentum.png)

**Figure 1.** Blue points show binned asymmetries with ±1 SE error bars.
The orange band shows the integrated asymmetry ±1 SE.
The dashed line indicates zero asymmetry.

The horizontal axis shows ordered momentum intervals, not an equally
spaced numerical momentum scale. The final point represents an
open-ended interval. Each panel has its own vertical scale.

The integrated band and bin estimates share observations and are
correlated. Their overlap is not an independent statistical test.
A nonzero bin asymmetry alone does not demonstrate momentum dependence.

All numerical bin results are saved in `momentum_asymmetry.csv`.
Undefined subsample ratios are reported as NaN instead of zero.

## 7. Supporting momentum spectra

![Momentum spectra](momentum_spectra.png)

**Figure 2.** Mean counts per event per unit momentum.
Vertical error bars show ±1 SE; horizontal bars show bin widths.
The vertical axis is logarithmic.

The open-ended ≥8 bin is excluded only from this density plot because
it has no finite width. It remains included in the asymmetry plot and
integrated results.

Overlapping spectra can still have a measurable asymmetry:
the relative differences may be small compared with the total counts.

## 8. Input-file momentum comparison

{file_table}

The supplied data show a substantially smaller fraction of selected
molecules with p ≥ 8 in file 10 than in files 1–9.
The cause is not established from the supplied information.

All ten files are retained. Their momentum distributions should not
be assumed identical without checking the simulation settings and
units. Round-robin subsamples mix files, so the primary error bars
describe sampling variation of this particular combined dataset;
they do not account for an unknown difference in generation conditions.

Consequently, an interpretation of momentum trends as a physical
effect remains provisional. The file-partition comparison is a
sensitivity diagnostic, not a complete correction.

## 9. Limitations

- Subsampling uncertainties rely on suitable independence assumptions.
- Only ten subsamples are available, so uncertainty estimates fluctuate.
- The pooled asymmetry is a nonlinear ratio; its SE is approximate.
- Rare species have less precise estimates.
- Six pairs and multiple bins are examined without a multiple-testing
  correction.
- No exact Gaussian discovery significance is claimed.
- Systematic effects are not estimated.
- A formal test of a momentum trend is not performed here.
- The input-file momentum difference requires clarification.

## 10. Conclusion

The combined dataset provides evidence of integrated molecular
asymmetry for **{pair_list(stable)}**, with the approximate three-SE
threshold exceeded under both partitions.

The results for **{pair_list(sensitive)}** depend on the subsampling
partition. The results for **{pair_list(inconclusive)}** do not give
strong evidence under this screening rule.

The required momentum-resolved asymmetries and statistical error bars
are shown for all six pairs. Their physical interpretation must account
for the observed difference between input-file momentum distributions.

## 11. Results files

- `mean_counts.csv`: totals, means, SDs and SEs.
- `integrated_asymmetry.csv`: integrated results and partition comparison.
- `momentum_asymmetry.csv`: results and uncertainties for every bin.
- `file_momentum_check.csv`: file-level momentum comparison.
- `asymmetry_vs_momentum.png`: required figure.
- `asymmetry_vs_momentum.pdf`: vector version of the required figure.
- `momentum_spectra.png`: supporting spectra.

Keep the images beside this README so they display on GitHub.
"""

    (OUTPUT / "README.md").write_text(report, encoding="utf-8")


# ============================================================
# RUN EVERYTHING
# ============================================================

def main():
    (
        counts, file_counts, sizes,
        events, molecules, empty_events
    ) = read_data()

    species_rows, integrated_rows, bin_rows = analyse(
        counts, file_counts, sizes, events
    )

    # Inspect the high-momentum fraction separately for each file.
    file_checks = []

    for i, path in enumerate(FILES):
        selected = int(file_counts[i].sum())
        high_momentum = int(file_counts[i, :, -1].sum())

        file_checks.append({
            "File": path.name,
            "Selected": selected,
            "Count_p_ge_8": high_momentum,
            "Fraction_p_ge_8": (
                high_momentum / selected if selected else math.nan
            ),
        })

    save_csv("mean_counts.csv", species_rows)
    save_csv("integrated_asymmetry.csv", integrated_rows)
    save_csv("momentum_asymmetry.csv", bin_rows)
    save_csv("file_momentum_check.csv", file_checks)

    plot_asymmetries(integrated_rows, bin_rows)
    plot_spectra(counts, sizes)

    create_readme(
        species_rows,
        integrated_rows,
        file_checks,
        events,
        molecules,
        empty_events,
    )

    print("\nANALYSIS COMPLETE")
    print(f"Events analysed: {events:,}")
    print(f"Molecule rows read: {molecules:,}")

    print("\nINTEGRATED ASYMMETRIES")
    print("Pair           A          SE          Approx. Z")

    for row in integrated_rows:
        pair = f'{row["Positive_ID"]}/{row["Negative_ID"]}'
        print(
            f'{pair:<12} '
            f'{row["Asymmetry"]: .7f}  '
            f'{row["Full_sample_SE"]:.7f}  '
            f'{row["Approx_Z"]:.2f}'
        )

    print(f"\nREADME and figures saved in:\n{OUTPUT}")


if __name__ == "__main__":
    main()
    