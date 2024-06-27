#!/usr/bin/env python
import os
import re
import subprocess
import multiprocessing
import pandas as pd
import pysam


# Global settings/variables
RTGPATH     = "./starfish_analysis/rtg-tools-3.12.1"
STARFISH_PY = "./starfish_analysis/starfish/starfish.py"
REFPATH     = "./starfish_analysis/GRCh38.sdf"
GLOBAL_THREADS = 4
# Setting up starfish environment
STARFISH_ENV = os.environ.copy()
STARFISH_ENV["PATH"] = f"{RTGPATH}:{STARFISH_ENV['PATH']}"
STARFISH_ENV["PATH"] = f"{RTGPATH}/jre/bin:{STARFISH_ENV['PATH']}"
    

def run_intersections_parallel(arglist, intersection_func, maxprocesses=1, subsetn=0):
    # Subset the supplied paired_paths to subsetn elements if positive integer supplied
    arglist = arglist[:subsetn] if subsetn > 0 else arglist
    # Create threadpool and map the filter/intersect function to the paired path list (ppl)
    with multiprocessing.Pool(maxprocesses) as threadpool:
        intersection_stats = threadpool.starmap(intersection_func, arglist)
    # Return collected stats on filtering and intersections (order is non-deterministic)
    return intersection_stats


def get_preset_slivar_filter(presetname, bedpath=None):
    presets = {
        "chiprange":             ["--info", "INFO.AF > 0.02 && INFO.AF < 0.2 && INFO.VD > 5"],
        "chiprange_bed":         ["--info", "INFO.AF > 0.02 && INFO.AF < 0.2 && INFO.VD > 5", "--region", bedpath],
        "chiprange_pass":        ["--info", "INFO.AF > 0.02 && INFO.AF < 0.2 && INFO.VD > 5 && variant.FILTER == 'PASS'"],
        "chiprange_pass_bed":    ["--info", "INFO.AF > 0.02 && INFO.AF < 0.2 && INFO.VD > 5 && variant.FILTER == 'PASS'", "--region", bedpath],
        "chiprange_noupper":     ["--info", "INFO.AF > 0.02 && INFO.VD > 5"],
        "chiprange_noupper_bed": ["--info", "INFO.AF > 0.02 && INFO.VD > 5", "--region", bedpath],
        "germline":              ["--info", "INFO.AF > 0.30 && INFO.VD > 5"],
        "overzero":              ["--info", "INFO.AF > 0.00"],
        "overzero_bed":          ["--info", "INFO.AF > 0.00", "--region", bedpath],
    }
    return presets[presetname] if presetname in presets else []


def pair_bld_sal_vcfs(pipeline_outdir, directory_whitelist, metadatadf):
    # Getting relative paths to vcf directories of interest
    vcfdirs = [os.path.join(pipeline_outdir, d, "vcfs") for d in os.listdir(pipeline_outdir) if d in directory_whitelist]

    # Creating a dictionary mapping sample names (of non-controls) to relative vcf paths
    sample_vcf_dict = {}
    for vd in vcfdirs:
        vcfs = [v for v in os.listdir(vd) if v.endswith(".vardict.norm.vcf.gz")]
        samples = [v.split('_')[0] if v.startswith("BT") else v.split('.')[0] for v in vcfs]
        for vcf, sample in zip(vcfs, samples):
            sample_vcf_dict[sample] = os.path.join(vd, vcf)

    # Getting metadata from dispatch list
    dispatchids = set(metadatadf["Spec ID"])

    # Checking ids from dispatch list not in the vcf lookup (should only be controls)
    assert all([True if k.startswith("Control") else False for k in sample_vcf_dict if k not in dispatchids])

    # Now pair ABC/saliva specids using UPN
    pairdf = metadatadf[metadatadf["Spec ID"].isin(sample_vcf_dict) & metadatadf["Collection Type"].isin({"ABC_BLD", "ABC_SALIVA"})] \
        .sort_values(["UPN", "Collection Type"])
    pairdf = pairdf.assign(vcfpath=pairdf["Spec ID"].map(sample_vcf_dict))
    paired_paths = pairdf[["UPN", "Collection Type", "vcfpath"]].pivot(index="UPN", columns="Collection Type", values="vcfpath")
    return paired_paths


def calculate_intersection_variant_stats(intid, a_vcf, b_vcf, ab_vcf, ba_vcf, a_label, b_label):
    # For each variant extract/calculate:
    # intid, source, intersecting, varid (string of chr:pos:ref:alt), DP, VD, VAF
    var_header = ("comparison_id", "source", "intersecting", "varid", "TYPE", "DP", "VD", "VAF", "PMEAN", "filter")
    var_records = []
    vcfpaths = (a_vcf, b_vcf, ab_vcf, ba_vcf)
    sources = (a_label, b_label, a_label, b_label)
    isintersecting = (False, False, True, True)
    for vcfpath, source, intersecting in zip(vcfpaths, sources, isintersecting):
        try:
            varfile = pysam.VariantFile(vcfpath).fetch()
        except ValueError:
            # ValueError indicates empty vcf, skip
            continue
        for var in varfile:
            varid = f"{var.chrom}:{var.pos}:{var.ref}:{var.alts[0]}"
            type = var.info["TYPE"]
            dp = var.info["DP"]
            vd = var.info["VD"]
            af = var.info["AF"][0]
            pmean = var.info["PMEAN"]
            filterstring = '|'.join(var.filter)
            rec = (intid, source, intersecting, varid, type, dp, vd, af, pmean, filterstring)
            var_records.append(rec)
    return pd.DataFrame.from_records(var_records, columns=var_header)


def filter_intersect_starfish_acrossmedia(inputs, slivar_filter):
    upn, vcfs = inputs
    # TODO: Define intermediate output location in main function / config
    outputpath = os.path.join("starfish_outputs", upn.replace(".", "-"))
    # Regular expression to extract filter counts from slivar stderr
    varcounts_re = re.compile("evaluated (?P<totalvars>\d+) total variants and wrote (?P<passfilter>\d+) variants that passed your slivar expressions")
    if not os.path.exists(outputpath):
        os.mkdir(outputpath)
    bldvcf_unsorted = os.path.join(outputpath, f"{os.path.basename(vcfs.ABC_BLD).rstrip('.vcf.gz')}.bld.flt.uns.vcf.gz")
    salvcf_unsorted = os.path.join(outputpath, f"{os.path.basename(vcfs.ABC_SALIVA).rstrip('.vcf.gz')}.sal.flt.uns.vcf.gz")
    bldvcf = os.path.join(outputpath, f"{os.path.basename(vcfs.ABC_BLD).rstrip('.vcf.gz')}.bld.flt.vcf.gz")
    salvcf = os.path.join(outputpath, f"{os.path.basename(vcfs.ABC_SALIVA).rstrip('.vcf.gz')}.sal.flt.vcf.gz")

    # Blood vcf filtering
    slivar_cmd_bld = ["slivar", "expr", "-v", vcfs.ABC_BLD, "-o", bldvcf_unsorted] + slivar_filter
    slivar_bld_stderr = str(subprocess.run(slivar_cmd_bld, capture_output=True).stderr)
    bld_sort_out = subprocess.run(["bcftools", "sort", bldvcf_unsorted, "-o", bldvcf], capture_output=True)  # Sorting blood vcf (required)
    subprocess.run(["bcftools", "index", "--tbi", bldvcf])  # Indexing blood vcf (required)

    # Blood filtering stats
    bld_counts = varcounts_re.search(slivar_bld_stderr)
    bld_total = int(bld_counts.group("totalvars"))
    bld_pass  = int(bld_counts.group("passfilter"))

    # Saliva vcf filtering
    slivar_cmd_sal = ["slivar", "expr", "-v", vcfs.ABC_SALIVA, "-o", salvcf_unsorted] + slivar_filter
    slivar_sal_stderr = str(subprocess.run(slivar_cmd_sal, capture_output=True).stderr)
    sal_sort_out = subprocess.run(["bcftools", "sort", salvcf_unsorted, "-o", salvcf], capture_output=True)  # Sorting saliva vcf (required)
    subprocess.run(["bcftools", "index", "--tbi", salvcf])  # Indexing saliva vcf (required)

    # Saliva filtering stats
    sal_counts = varcounts_re.search(slivar_sal_stderr)
    sal_total = int(sal_counts.group("totalvars"))
    sal_pass  = int(sal_counts.group("passfilter"))

    # Starfish intersection stats
    starfish_cmd = [STARFISH_PY, "-t", REFPATH, "-V", bldvcf, salvcf, "-O", outputpath,
                    "--threads", str(GLOBAL_THREADS), "--squash-ploidy", "--all-records"]
    subprocess.run(starfish_cmd, env=STARFISH_ENV)
    a_vcf   = os.path.join(outputpath, "A.vcf.gz")
    b_vcf   = os.path.join(outputpath, "B.vcf.gz")
    ab_vcf  = os.path.join(outputpath, "A_and_B.vcf.gz")
    ba_vcf  = os.path.join(outputpath, "B_and_A.vcf.gz")
    varstats = calculate_intersection_variant_stats(upn, a_vcf, b_vcf, ab_vcf, ba_vcf, "BLD", "SAL")
    # varstats.to_csv(os.path.join(outputpath, f"{upn.replace('.', '-')}_varstats.csv"), index=False)

    # Getting raw counts from intersection vcfs
    both_count = len(list(pysam.VariantFile(os.path.join(outputpath, "AB.vcf.gz")).fetch()))
    bld_count  = len(list(pysam.VariantFile(os.path.join(outputpath, "A.vcf.gz")).fetch()))
    sal_count  = len(list(pysam.VariantFile(os.path.join(outputpath, "B.vcf.gz")).fetch()))

    print(f"{upn} complete.")
    # Returning collected stats on intersection (and variant stats dataframe?)
    # - WARN: Keep an eye on memory usage as these dataframes collect if being run in parallel
    # UPN, BLD_VCF, SAL_VCF, totalvars_BLD, filteredvars_BLD, totalvars_SAL, filteredvars_SAL, BLD_only, SAL_only, BLD_SAL
    return (varstats, (upn, os.path.basename(vcfs.ABC_BLD), os.path.basename(vcfs.ABC_SALIVA),
            bld_total, bld_pass, sal_total, sal_pass, bld_count, sal_count, both_count))


def main():
    # Input values and paths
    bedpath    = "./ten_gene_bed/CHIP_ten_genes_padded_20bp.bed"
    filtername = "chiprange_bed"
    panelname  = "CHIP-only"

    toplevel_vcfdir  = "./pipeline_outputs"
    output_directory = "./bldsal_refactor_testout"
    output_basename  = "ARCH_bld_sal_chip_only_tengenes_270624"
    maxintersections = 0  # 0 For unlimited
    stats_csvname    = os.path.join(output_directory, f"{output_basename}_{panelname}_flt-{filtername}.csv")
    varstats_csvname = os.path.join(output_directory, f"{output_basename}_{panelname}_flt-{filtername}_allvarstats.csv")

    paneldirs = {
        "CHIP-only": {'ARCH_run04_plate1', 'ARCH_run04_plate2'},
    }
    directory_whitelist = paneldirs[panelname]

    # Pair vcfs for comparison using dispatch list metadata
    metadata_path = "CHIP_phase 2_plate_layouts_CHIP samples.xlsx"
    metadatadf    = pd.read_excel(metadata_path, sheet_name="Dispatchlist")
    paired_paths  = pair_bld_sal_vcfs(toplevel_vcfdir, directory_whitelist, metadatadf)

    # Prepare args for intersection function with required presets
    slivar_filter = get_preset_slivar_filter(filtername, bedpath)
    arglist = [(row, slivar_filter) for row in paired_paths.iterrows()]

    dfs_and_stats = run_intersections_parallel(arglist, filter_intersect_starfish_acrossmedia, maxprocesses=4, subsetn=maxintersections)

    # Aggregate variant stats
    allvarstats = pd.concat([t[0] for t in dfs_and_stats if not t[0].empty])
    allvarstats.to_csv(varstats_csvname, index=False)

    # Aggregate intersection stats
    stats_header = ("UPN", "BLD_VCF", "SAL_VCF", "totalvars_BLD", "filteredvars_BLD",
                    "totalvars_SAL", "filteredvars_SAL", "BLD_only", "SAL_only", "BLD_SAL")
    intersection_stats = [t[1] for t in dfs_and_stats]
    statdf = pd.DataFrame.from_records(intersection_stats, columns=stats_header)

    # Have csv with raw counts of variant stats, want to normalise the intersections to percentages/fractions and write out
    int_totals = statdf.BLD_only + statdf.SAL_only + statdf.BLD_SAL
    final_statdf = statdf.assign(
        pct_BLD_only=statdf.BLD_only / int_totals,
        pct_SAL_only=statdf.SAL_only / int_totals,
        pct_BLD_SAL=statdf.BLD_SAL / int_totals,
    )
    final_statdf.to_csv(stats_csvname, index=False)


if __name__ == "__main__":
    main()