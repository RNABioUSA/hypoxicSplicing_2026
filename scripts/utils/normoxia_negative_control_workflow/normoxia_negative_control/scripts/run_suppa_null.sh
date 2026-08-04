#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: $0 CONFIG_ENV" >&2
    exit 2
}

[[ $# -eq 1 ]] || usage
config_file=$1
[[ -s "$config_file" ]] || { echo "Config not found: $config_file" >&2; exit 1; }

# shellcheck disable=SC1090
source "$config_file"

: "${WORKFLOW_ROOT:?Set WORKFLOW_ROOT in the config}"
: "${OUTPUT_ROOT:?Set OUTPUT_ROOT in the config}"
: "${GTF:?Set GTF in the config}"
: "${SUPPA_PY:=suppa.py}"

analysis_root="${OUTPUT_ROOT}/analysis"
membership="${analysis_root}/permutation_membership.tsv"
master_tpm="${analysis_root}/suppa/all_samples_iso_tpm.txt"
suppa_root="${analysis_root}/suppa"
per_event_dir="${suppa_root}/event_annotation/per_event"
ioe_master="${suppa_root}/event_annotation/events_all.ioe"

[[ -s "$membership" ]] || {
    echo "Permutation membership not found; run run_null_statistics.R first: $membership" >&2
    exit 1
}
[[ -s "$master_tpm" ]] || {
    echo "Master isoform TPM file not found: $master_tpm" >&2
    exit 1
}
[[ -s "$GTF" ]] || { echo "GTF not found: $GTF" >&2; exit 1; }

if [[ "$SUPPA_PY" == */* ]]; then
    [[ -x "$SUPPA_PY" ]] || { echo "SUPPA executable not found: $SUPPA_PY" >&2; exit 1; }
else
    command -v "$SUPPA_PY" >/dev/null 2>&1 || {
        echo "SUPPA command not found: $SUPPA_PY" >&2
        exit 1
    }
fi

mkdir -p "$per_event_dir"

if [[ -s "$ioe_master" ]]; then
    echo "[SUPPA2] event annotation present: $ioe_master"
else
    echo "[SUPPA2] generating event annotation from: $GTF"
    "$SUPPA_PY" generateEvents \
        -i "$GTF" \
        -o "${per_event_dir}/events" \
        -f ioe \
        -e SE SS MX RI FL

    mapfile -t ioe_files < <(find "$per_event_dir" -maxdepth 1 -type f -name '*.ioe' -print | sort)
    ((${#ioe_files[@]} > 0)) || {
        echo "SUPPA2 did not generate any .ioe files." >&2
        exit 1
    }
    awk 'FNR == 1 && NR != 1 { next } { print }' "${ioe_files[@]}" >"$ioe_master"
    [[ -s "$ioe_master" ]] || { echo "Failed to create: $ioe_master" >&2; exit 1; }
fi

subset_suppa_matrix() {
    local input_file=$1
    local sample_csv=$2
    local output_file=$3

    awk -F '\t' -v OFS='\t' -v wanted="$sample_csv" '
        BEGIN {
            n_wanted = split(wanted, requested, ",")
        }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                header_pos[$i] = i
            }
            header = ""
            for (k = 1; k <= n_wanted; k++) {
                if (!(requested[k] in header_pos)) {
                    print "Sample not found in SUPPA matrix header: " requested[k] > "/dev/stderr"
                    exit 3
                }
                data_pos[k] = header_pos[requested[k]] + 1
                header = (k == 1 ? requested[k] : header OFS requested[k])
            }
            print header
            next
        }
        {
            out = $1
            for (k = 1; k <= n_wanted; k++) {
                out = out OFS $(data_pos[k])
            }
            print out
        }
    ' "$input_file" >"$output_file"
}

mapfile -t perm_ids < <(awk -F '\t' 'NR > 1 {print $1}' "$membership" | sort -u)
((${#perm_ids[@]} > 0)) || { echo "No permutations found in $membership" >&2; exit 1; }

for perm_id in "${perm_ids[@]}"; do
    perm_dir="${suppa_root}/${perm_id}"
    mkdir -p "$perm_dir"
    group_a_tpm="${perm_dir}/group_A_iso_tpm.txt"
    group_b_tpm="${perm_dir}/group_B_iso_tpm.txt"
    group_a_psi="${perm_dir}/group_A.psi"
    group_b_psi="${perm_dir}/group_B.psi"
    out_prefix="${perm_dir}/diffSplice"
    dpsi_file="${out_prefix}.dpsi"

    if [[ -s "$dpsi_file" ]]; then
        echo "[SUPPA2] present: $perm_id"
        continue
    fi

    group_a_csv=$(
        awk -F '\t' -v p="$perm_id" '
            NR > 1 && $1 == p && $3 == "A" {
                out = out (out == "" ? "" : ",") $2
            }
            END { print out }
        ' "$membership"
    )
    group_b_csv=$(
        awk -F '\t' -v p="$perm_id" '
            NR > 1 && $1 == p && $3 == "B" {
                out = out (out == "" ? "" : ",") $2
            }
            END { print out }
        ' "$membership"
    )
    [[ -n "$group_a_csv" && -n "$group_b_csv" ]] || {
        echo "Could not resolve both pseudo-groups for $perm_id" >&2
        exit 1
    }

    echo "[SUPPA2] $perm_id"
    subset_suppa_matrix "$master_tpm" "$group_a_csv" "$group_a_tpm"
    subset_suppa_matrix "$master_tpm" "$group_b_csv" "$group_b_tpm"

    "$SUPPA_PY" psiPerEvent \
        -i "$ioe_master" \
        -e "$group_a_tpm" \
        -o "${perm_dir}/group_A"
    "$SUPPA_PY" psiPerEvent \
        -i "$ioe_master" \
        -e "$group_b_tpm" \
        -o "${perm_dir}/group_B"

    [[ -s "$group_a_psi" && -s "$group_b_psi" ]] || {
        echo "SUPPA2 PSI calculation failed for $perm_id" >&2
        exit 1
    }

    "$SUPPA_PY" diffSplice \
        -m empirical \
        -gc \
        -s \
        -i "$ioe_master" \
        -p "$group_a_psi" "$group_b_psi" \
        -e "$group_a_tpm" "$group_b_tpm" \
        -o "$out_prefix"

    [[ -s "$dpsi_file" ]] || {
        echo "SUPPA2 diffSplice output missing for $perm_id" >&2
        exit 1
    }
done

echo "[SUPPA2] completed: $suppa_root"
