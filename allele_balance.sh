# script to calculate the proportion of reference allele and alternative allele at heterozygous biallelic positions in VCF for every sample indivisually
# consider only sites with FORMAT/DP >= 6
#!/usr/bin/env bash
set -euo pipefail

VCF="${1:?Usage: bash allele_balance.sh input.vcf.gz [output.tsv]}"
OUT="${2:-allele_balance.tsv}"

command -v bcftools >/dev/null 2>&1 || {
    echo "ERROR: bcftools is not available." >&2
    exit 1
}

[[ -f "$VCF" ]] || {
    echo "ERROR: Input file not found: $VCF" >&2
    exit 1
}

# Check required FORMAT fields.
HEADER=$(bcftools view -h "$VCF")
for TAG in GT DP AD; do
    if [[ "$HEADER" != *"##FORMAT=<ID=${TAG},"* ]]; then
        echo "ERROR: FORMAT/$TAG is missing from the VCF header." >&2
        exit 1
    fi
done

# Supply sample names first so samples without qualifying sites
# are also included in the output.
{
    bcftools query -l "$VCF"

    bcftools view -m2 -M2 -Ou "$VCF" |
        bcftools query -f '[%SAMPLE\t%GT\t%DP\t%AD\n]'
} |
awk '
BEGIN {
    FS = OFS = "\t"
}

# Sample list, in original VCF order.
NF == 1 {
    samples[++ns] = $1
    next
}

{
    sample = $1
    gt = $2
    dp = $3

    # Accept phased and unphased diploid heterozygotes.
    if (gt !~ /^(0[\/|]1|1[\/|]0)$/)
        next

    if (dp !~ /^[0-9]+$/ || dp + 0 <= 6)
        next

    # Skip missing, invalid, or zero-total allele depths.
    if (split($4, ad, ",") != 2)
        next

    if (ad[1] !~ /^[0-9]+$/ || ad[2] !~ /^[0-9]+$/)
        next

    total = ad[1] + ad[2]
    if (total <= 0)
        next

    ref_sum[sample] += ad[1] / total
    alt_sum[sample] += ad[2] / total
    n[sample]++
}

END {
    print "sample", "mean_frequency_reference", "mean_frequency_alternative"

    for (i = 1; i <= ns; i++) {
        sample = samples[i]

        if (n[sample] > 0)
            printf "%s\t%.8f\t%.8f\n", sample, \
                ref_sum[sample] / n[sample], \
                alt_sum[sample] / n[sample]
        else
            print sample, "NA", "NA"
    }
}
' > "$OUT"

echo "Saved: $OUT" >&2
