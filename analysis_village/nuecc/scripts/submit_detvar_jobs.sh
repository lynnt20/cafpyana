#!/usr/bin/env bash
# Submit grid jobs for all detector variation samples.

NGRID=25
CONFIG="analysis_village/nuecc/configs/nuecc_mc_detvar.py"
SAMPLE_DIR="data/sample_lists/sbnd/prod_2025B/v10_06_00_10"

for sample_list in "${SAMPLE_DIR}"/mc_MCP2025B_1e20_10_*_SystVar_*.txt; do
    # Extract the variation name from the filename, e.g. SystVar_2xSCE -> 2xSCE
    fname=$(basename "${sample_list}" .txt)
    variant=$(echo "${fname}" | sed 's/.*SystVar_\(.*\)_caf_flat_caf_sbnd/\1/')
    output_tag="detvar_$(echo "${variant}" | tr '[:upper:]' '[:lower:]')"

    echo "Submitting: variant=${variant}  output=${output_tag}"
    python run_df_maker.py \
        -ngrid "${NGRID}" \
        -c "${CONFIG}" \
        -l "${sample_list}" \
        -o "${output_tag}"
done
