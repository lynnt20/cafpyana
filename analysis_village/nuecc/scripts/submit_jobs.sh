#!/usr/bin/env bash
# Submit grid jobs for a user-defined set of (sample_list, config, output_tag, ngrid) quads.
# Each job is 4 consecutive entries in JOBS: sample_list, config, output_tag, ngrid.

JOBS=(
    # CV
    # "data/sample_lists/sbnd/prod_2025B/v10_06_00_10/mc_MCP2025B_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_CV_caf_flat_caf_sbnd.txt"
    # "analysis_village/nuecc/configs/nuecc_mc_detvar.py"
    # "detvar_cv"
    # "25"
    # # 2xSCE
    # "data/sample_lists/sbnd/prod_2025B/v10_06_00_10/mc_MCP2025B_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_2xSCE_caf_flat_caf_sbnd.txt"
    # "analysis_village/nuecc/configs/nuecc_mc_detvar.py"
    # "detvar_2xsce"
    # "25"
    # # 0xSCE
    # "data/sample_lists/sbnd/prod_2025B/v10_06_00_10/mc_MCP2025B_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_0xSCE_caf_flat_caf_sbnd.txt"
    # "analysis_village/nuecc/configs/nuecc_mc_detvar.py"
    # "detvar_0xsce"
    # "25"
    # # PMTGainFluct
    # "data/sample_lists/sbnd/prod_2025B/v10_06_00_10/mc_MCP2025B_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_PMTGainFluct_caf_flat_caf_sbnd.txt"
    # "analysis_village/nuecc/configs/nuecc_mc_detvar.py"
    # "detvar_pmtgainfluct"
    # "25"
    # # PMTHighNoise
    # "data/sample_lists/sbnd/prod_2025B/v10_06_00_10/mc_MCP2025B_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_PMTHighNoise_caf_flat_caf_sbnd.txt"
    # "analysis_village/nuecc/configs/nuecc_mc_detvar.py"
    # "detvar_pmthighnoise"
    # "25"
    # # PMTLowEff
    # "data/sample_lists/sbnd/prod_2025B/v10_06_00_10/mc_MCP2025B_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_PMTLowEff_caf_flat_caf_sbnd.txt"
    # "analysis_village/nuecc/configs/nuecc_mc_detvar.py"
    # "detvar_pmtloweff"
    # "25"
    # # WireMod_XThetaXW
    # "data/sample_lists/sbnd/prod_2025B/v10_06_00_10/mc_MCP2025B_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_WireMod_XThetaXW_caf_flat_caf_sbnd.txt"
    # "analysis_village/nuecc/configs/nuecc_mc_detvar.py"
    # "detvar_wiremodxw"
    # "25"
    # # WireMod_YZ
    # "data/sample_lists/sbnd/prod_2025B/v10_06_00_10/mc_MCP2025B_1e20_10_prodgenie_corsika_proton_rockbox_sbnd_SystVar_WireMod_YZ_caf_flat_caf_sbnd.txt"
    # "analysis_village/nuecc/configs/nuecc_mc_detvar.py"
    # "detvar_wiremodyz"
    # "25"
    # mc low energy
    # "analysis_village/nuecc/file_lists/mc_lowenergy.paths"
    # "analysis_village/nuecc/configs/nuecc_mc.py"
    # "mc_lowenergy"
    # "100"
    # # offbeamlight
    # "analysis_village/nuecc/file_lists/dt_offbeamlight.paths"
    # "analysis_village/nuecc/configs/nuecc.py"
    # "offbeamlight"
    # "25"
    # # intime MC
    # "analysis_village/nuecc/file_lists/mc_intime.paths"
    # "analysis_village/nuecc/configs/nuecc_mc.py"
    # "mc_intime"
    # "200"
    

    # MC no systs 
    "analysis_village/nuecc/file_lists/mc_nominal_ar23.paths"
    "analysis_village/nuecc/configs/nuecc_mc_detvar.py"
    "mc_nosyst"
    "25"
    # 
analysis_village/nuecc/configs/nuecc_mc_detvar.py -l analysis_village/nuecc/file_lists/mc_wiremodyz_add.paths -o detvar_wiremodyz_add
)

if (( ${#JOBS[@]} % 4 != 0 )); then
    echo "Error: JOBS has ${#JOBS[@]} entries, which is not a multiple of 4. Check for a missing or extra line." >&2
    exit 1
fi

for ((i=0; i<${#JOBS[@]}; i+=4)); do
    sample_list="${JOBS[i]}"
    config="${JOBS[i+1]}"
    output_tag="${JOBS[i+2]}"
    ngrid="${JOBS[i+3]}"

    if [[ -z "${sample_list}" || -z "${config}" || -z "${output_tag}" || -z "${ngrid}" ]]; then
        echo "Error: empty field in job at index ${i} (sample_list='${sample_list}', config='${config}', output_tag='${output_tag}', ngrid='${ngrid}')" >&2
        exit 1
    fi

    echo "Submitting: config=${config}  list=${sample_list}  output=${output_tag}  ngrid=${ngrid}"
    python run_df_maker.py \
        -ngrid "${ngrid}" \
        -c "${config}" \
        -l "${sample_list}" \
        -o "${output_tag}"
done
