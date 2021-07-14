require(VariantAnnotation)
require(GenomicRanges)

get_cis_geno = function(tab, chrom, left, right, genome_build = "GRCh38") {
  gr = GRanges(chrom, IRanges(left, right))
  sp = ScanVcfParam(which = gr)
  vcf = readVcf(tab, genome_build, param = sp)
  gt = geno(vcf)$GT
  if (nrow(gt) == 0) return(NULL)
  allele1 = substr(gt, 1, 1)
  class(allele1) = "numeric"
  allele2 = substr(gt, 3, 3)
  class(allele2) = "numeric"
  allele1 + allele2
}


