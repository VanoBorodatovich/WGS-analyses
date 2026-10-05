#!/usr/bin/env bash
set -euo pipefail

# Usage: bash individual_heterozygosity.sh input.vcf.gz [DP_EXCLUSIVE=5] \
#          [MIN_MINOR_READS=2] [MIN_MINOR_FRACTION=0] > heterozygosity.tsv
# Requires bcftools and awk. Streams the VCF; no index or temporary files needed.
# Input: one record per base, including invariant bases; not gVCF blocks.
# Counts invariant single bases and biallelic SNPs only, using diploid GT.
# Failed heterozygotes are excluded from BOTH numerator and denominator.
# Homozygous calls are not reclassified using AD.

if (( $# < 1 || $# > 4 )); then
    echo "Usage: $0 input.vcf[.gz] [DP_EXCLUSIVE=5] [MIN_MINOR_READS=2] [MIN_MINOR_FRACTION=0]" >&2
    exit 1
fi
vcf=$1
dp_min=${2:-6}
minor_reads=${3:-2}
minor_fraction=${4:-0}
out="${vcf%.*}_het_dpmin${dp_min}_admin${minor_reads}.tsv"
[[ -r "$vcf" ]] || { echo "Cannot read: $vcf" >&2; exit 1; }
command -v bcftools >/dev/null || { echo 'bcftools is required.' >&2; exit 1; }
[[ "$dp_min" =~ ^[0-9]+$ && "$minor_reads" =~ ^[0-9]+$ ]] || {
    echo 'Depth and minor-read thresholds must be nonnegative integers.' >&2; exit 1;
}
awk -v x="$minor_fraction" 'BEGIN {
    # BEGIN runs once without reading input. The expression checks numeric syntax
    # and the allowed range; ~ means matches a regular expression, && means AND.
    # ! reverses the result: exit 0 means valid, exit 1 makes Bash report an error.
    exit !(x ~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)$/ && x >= 0 && x <= 0.5)
    # The closing brace below ends this one-time validation block.
}' || { echo 'Minor fraction must be between 0 and 0.5.' >&2; exit 1; }

header=$(bcftools view -h "$vcf")
for tag in GT DP AD; do
    [[ "$header" == *"##FORMAT=<ID=${tag},"* ]] || {
        echo "Required FORMAT/$tag is absent from the VCF header." >&2; exit 1;
    }
done
samples=$(bcftools query -l "$vcf")
[[ -n "$samples" ]] || { echo 'No samples found.' >&2; exit 1; }

bcftools query -f '%REF\t%ALT[\t%GT\t%DP\t%AD]\n' "$vcf" |
awk -v names="$samples" -v dp_min="$dp_min" \
    -v min_reads="$minor_reads" -v min_frac="$minor_fraction" '
# BEGIN runs once before any site records are read; -v above passes Bash values to awk.
BEGIN {
    # FS splits input fields at tabs; OFS separates output columns with tabs.
    FS = OFS = "\t"
    # Split the newline-separated sample names into sample[1], sample[2], etc.; n is their count.
    n = split(names, sample, "\n")
    # End the one-time setup block.
}
# This unnamed block runs once for every input line, representing one VCF site.
{
    # $1 is REF and $2 is ALT. Require a single A/C/G/T REF and either ALT=. or one A/C/G/T ALT.
        if ($1 !~ /^[ACGTacgt]$/ || ($2 != "." && $2 !~ /^[ACGTacgt]$/)) next
    # Process samples 1 through n; i++ increments the sample index after each iteration.
    for (i = 1; i <= n; i++) {
        # REF and ALT occupy fields 1 and 2; each sample then occupies three fields: GT, DP, AD.
        # Thus its GT column j is 3 for sample 1, 6 for sample 2, 9 for sample 3, etc.
        j = 3 + (i-1)*3
        # Read the fields numbered j, j+1 and j+2 into variables; semicolons separate assignments.
        gt = $j; dp = $(j+1); ad = $(j+2)
        # Reject missing/noninteger depth or DP <= threshold; default 5 therefore requires DP >= 6.
        # [0-9]+ means one or more digits; +0 forces numeric interpretation; continue skips this sample/site.
        if (dp !~ /^[0-9]+$/ || dp+0 < dp_min) continue
        # Split GT at either / or | into alleles[1] and alleles[2]; require exactly two parts (diploid).
        if (split(gt, alleles, /[\/|]/) != 2) continue
        # Each allele must be exactly 0 or 1; reject missing alleles and allele indices 2 or higher.
        if (alleles[1] !~ /^[01]$/ || alleles[2] !~ /^[01]$/) continue
        # When ALT is absent (.), accept only 0/0 or 0|0; a nonreference allele would be inconsistent.
        if ($2 == "." && (alleles[1] != 0 || alleles[2] != 0)) continue
        # Count this genotype after site, genotype and depth checks, but BEFORE allele-support checks.
        # Array counters are separate for each sample; uninitialized counters start numerically at zero.
        depth_pass[i]++
        # Different allele indices mean a heterozygote; only these calls enter the AD-filtering block.
        if (alleles[1] != alleles[2]) {
            # Split AD at commas; require exactly two counts: counts[1]=REF and counts[2]=ALT.
            if (split(ad, counts, ",") != 2 ||
                # Continue that condition: reject either count if it is not a nonnegative integer.
                counts[1] !~ /^[0-9]+$/ || counts[2] !~ /^[0-9]+$/) {
                # Count the invalid/missing AD, then skip this sample/site without counting it as callable.
                bad_ad[i]++; continue
            # End the invalid-AD conditional block.
            }
            # Sum reference and alternative read counts; this AD sum can differ from FORMAT/DP.
            total = counts[1] + counts[2]
            # A zero AD sum cannot define an allele fraction: count invalid AD and skip the sample/site.
            if (total <= 0) { bad_ad[i]++; continue }
            # Ternary syntax condition ? value_if_true : value_if_false selects the smaller read count.
            # This is the less-supported allele within this individual, not the population minor allele.
            minor = (counts[1]+0 < counts[2]+0 ? counts[1]+0 : counts[2]+0)
            # Reject if the smaller count OR its fraction of the AD sum is below its threshold.
            # Equality passes; min_reads=2 requires both alleles to have at least two supporting reads.
            if (minor < min_reads || minor/total < min_frac) {
                # Record failed support and skip; the rejected heterozygote enters neither H numerator nor denominator.
                failed[i]++; continue
            # End the allele-support rejection block.
            }
            # All heterozygote checks passed: add one to the heterozygous-site numerator for this sample.
            het[i]++
        # End the heterozygote-only block; homozygotes bypass all AD checks.
        }
        # Count every accepted homozygote or accepted heterozygote in the callable-site denominator.
        callable[i]++
    # End the loop over samples at this site.
    }
 # End processing of this site; awk proceeds to the next input line.
}
# END runs once after all input sites have been processed.
END {
    # Print the first four column headings; the trailing backslash continues the same statement below.
    # The next code line supplies the final three headings; OFS inserts tabs between all headings.
    print "Sample", "Heterozygosity", "Heterozygous_sites", "Callable_sites", \
          "Depth_passing_sites", "Het_failed_support", "Het_missing_or_invalid_AD"
    # Output one summary row for each sample, retaining the original VCF sample order.
    for (i = 1; i <= n; i++) {
        # If callable sites exist, calculate het/callable; otherwise output NA (division by zero is avoided).
        # sprintf formats the ratio as text with up to 10 significant digits (%.10g).
        h = (callable[i] > 0 ? sprintf("%.10g", het[i]/callable[i]) : "NA")
        # Print the name, ratio and counters; +0 makes never-incremented counters print as 0 rather than blank.
        print sample[i], h, het[i]+0, callable[i]+0, depth_pass[i]+0, failed[i]+0, bad_ad[i]+0
    # End the output loop.
    }
# Close END; the single quote below also closes the awk program passed by Bash.
}' > "$out"
