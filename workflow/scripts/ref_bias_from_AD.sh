#!/usr/bin/env bash
# ref_bias_from_ad.sh
#
# Computes reference-allele read balance at heterozygous sites, per sample
# and per population, using the AD (allelic depth) FORMAT field of a bialVCF.
# A value near 0.5 = no detectable mapping bias. Consistently >0.5 = reads
# carrying the reference allele are being favored at het sites.
#
# Requires: bcftools, awk, sort, coreutils
#
# Usage:
#   1. Edit the CONFIG section below (or export the variables before running).
#   2. Create pop_map.txt: two tab-separated columns, no header:
#        sample1<TAB>population1
#        sample2<TAB>population1
#        sample3<TAB>population2
#   3. ./ref_bias_from_ad.sh
#
# Outputs:
#   per_sample/<sample>.ref_frac.txt   - one line per het site: ref_frac, z-score
#   per_sample_summary.tsv             - mean/median ref balance per sample
#   population_summary.tsv             - mean ref balance averaged across samples per population

set -euo pipefail

# ---------------- CONFIG ----------------
VCF="${VCF:-/mnt/sda/saad/WW_Trimmed_Reads_WWPaper/Lsinapis-vgMapping/results/genotyping_rescaled/merged.all.DTOLREF.sitefilt.bQ20.mq30.snps5.noIndel.Q30.dp6-30.AB.jointCall.rescaled.biallelic.fmiss0.1.bcf}" #should be biallelic and contain AD field
POP_MAP="${POP_MAP:-/mnt/sda/saad/WW_Trimmed_Reads_WWPaper/Lsinapis-vgMapping/Pops_GTS.txt}"
MIN_DEPTH="${MIN_DEPTH:-6}"       # minimum ref+alt depth to trust a site
OUTDIR="${OUTDIR:-.}"
# -----------------------------------------

mkdir -p "${OUTDIR}/per_sample"
cd "${OUTDIR}"

echo "[1/4] Checking AD field is present..."
if ! bcftools view -h "${VCF}" | grep -q '^##FORMAT=<ID=AD'; then
  echo "ERROR: AD FORMAT field not declared in VCF header. Regenerate calls with a caller/setting that outputs AD." >&2
  exit 1
fi

#echo "[2/4] Restricting to biallelic SNPs..."
#bcftools view -m2 -M2 -v snps "${VCF}" -Oz -o biallelic.vcf.gz
#bcftools index -f biallelic.vcf.gz

echo "[3/4] Extracting per-sample het-site ref fraction + z-score..."
for s in $(bcftools query -l "${VCF}"); do
  bcftools query -s "${s}" -f '[%GT\t%AD]\n' "${VCF}" \
    | awk -v OFS='\t' -v mind="${MIN_DEPTH}" '
        $1=="0/1" || $1=="0|1" || $1=="1/0" || $1=="1|0" {
          n = split($2, ad, ",")
          if (n < 2 || ad[1] == "." || ad[2] == ".") next
          ref = ad[1]; alt = ad[2]; depth = ref + alt
          if (depth >= mind) {
            frac = ref / depth
            z = (ref - depth * 0.5) / sqrt(depth * 0.25)
            print frac, z
          }
        }' > "per_sample/${s}.ref_frac.txt"
done

echo "[4/4] Summarizing per sample and per population..."
{
  echo -e "sample\tn_het_sites\tmean_ref_frac\tmedian_ref_frac\tpct_sites_skewed_absZ_gt_1.96"
  for f in per_sample/*.ref_frac.txt; do
    s=$(basename "${f}" .ref_frac.txt)
    n=$(wc -l < "${f}")
    if [ "${n}" -eq 0 ]; then
      echo -e "${s}\t0\tNA\tNA\tNA"
      continue
    fi
    mean=$(awk '{sum+=$1} END{printf "%.4f", sum/NR}' "${f}")
    median=$(sort -k1,1n "${f}" | awk '{a[NR]=$1} END{
      if (NR % 2 == 1) printf "%.4f", a[(NR+1)/2]
      else printf "%.4f", (a[NR/2] + a[NR/2+1]) / 2
    }')
    pct_skew=$(awk '{if ($2 < -1.96 || $2 > 1.96) c++} END{printf "%.2f", (c/NR)*100}' "${f}")
    echo -e "${s}\t${n}\t${mean}\t${median}\t${pct_skew}"
  done
} > per_sample_summary.tsv

# Population-level: mean of per-sample means, weighted equally by sample
awk -F'\t' '
  NR==FNR { pop[$1] = $2; next }          # pop_map.txt: sample, population
  FNR==1 { next }                          # skip per_sample_summary.tsv header
  $3 != "NA" {
    p = pop[$1]
    if (p == "") { print "WARNING: sample " $1 " not found in pop_map.txt" > "/dev/stderr"; next }
    sum[p] += $3
    cnt[p]++
  }
  END {
    print "population\tn_samples\tmean_ref_frac_across_samples"
    for (p in sum) printf "%s\t%d\t%.4f\n", p, cnt[p], sum[p] / cnt[p]
  }
' "${POP_MAP}" per_sample_summary.tsv > population_summary.tsv

echo ""
echo "Done. Results in ${OUTDIR}/:"
echo "  per_sample_summary.tsv   - per-sample ref allele balance at het sites"
echo "  population_summary.tsv   - averaged per population"
echo ""
column -t population_summary.tsv
