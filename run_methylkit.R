suppressPackageStartupMessages(library(methylKit))
suppressPackageStartupMessages(library(genomation))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(ggrepel))


message("Running methylkit")

# set up input and output directories
genome = snakemake@config[["genome"]]
genome_dir = file.path(snakemake@config[["genome_root_dir"]], genome)
project_dir = snakemake@config[["output_dir"]]
message(paste0("project dir is ",project_dir))
processed_dir = file.path(project_dir, "processed")

message(paste0("Looking for processed sequences in ",processed_dir))


rdata_dir = file.path(project_dir, "rdata")
exit_code <- dir.create(rdata_dir, showWarnings=FALSE, recursive=TRUE)
if (exit_code) {
    message(sprintf("Created r data directory %s", rdata_dir))
}

output_dir = file.path(project_dir, "differential_output")
exit_code <- dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
if (exit_code) {
    message(sprintf("Created output directory %s", output_dir))
}
setwd(output_dir)

# sample data

message("loading sample info")
sample_data = snakemake@config[["samples"]]
sample_df = ldply(sample_data, data.frame)
colnames(sample_df) <- c("condition", "sample_string")
s = str_split(sample_df$sample_string," ")
sample_df = data.frame(condition = rep(sample_df$condition, sapply(s, length)), sample = unlist(s)) %>%
    mutate(name = paste(condition, row_number(condition), sep="_"))


conditions = sample_df$condition
treatment_vector = as.numeric(as.factor(conditions))
sample_df$treatment = treatment_vector
all_samples = sample_df$sample
sample_names = as.character(sample_df$name)

message("Running on samples")
message(paste(all_samples, sep=", "))
message("in conditions")
message(paste(conditions, sep=", "))
message("with names")
message(paste(sample_names, sep=", "))


write.table(sample_df, file.path(rdata_dir, "sample_info.txt"), sep="\t", row.names=FALSE, quote=FALSE)

# load annotation

gtf.file = file.path(genome_dir, snakemake@config[["annotation_filename"]])
message(sprintf("Loading gene annotation from %s", gtf.file))
transcript.db = makeTxDbFromGFF(gtf.file, format="gtf", dataSource="Ensembl")

genenames.file = snakemake@input[[1]]
message(sprintf("Loading gene names from %s", genenames.file))
geneid_to_genename = read.table(genenames.file, sep="\t", fill=TRUE)
colnames(geneid_to_genename) <- c("gene_id","gene_name")

exons = sort(exons(transcript.db))
cds = sort(cds(transcript.db))
genes = sort(genes(transcript.db))
transcript = sort(transcripts(transcript.db))
intergenic = gaps(genes)
introns = sort(unlist(intronsByTranscript(transcript.db)))
promoters = sort(promoters(transcript.db, upstream=1500, downstream=1000))
tss = sort(promoters(transcript.db, upstream=0, downstream=1))

annotation.all.features=GRangesList(genes=genes,
                                    intergenic=intergenic,
                                    transcript=transcript,
                                    cds=cds,
                                    exons=exons,
                                    introns=introns,
                                    promoters=promoters,
                                    TSSes=tss)
rm(exons, cds, genes, transcript, promoters, tss, introns, intergenic)

# process bismark output into methylRawList object, calculate some basic stats and do some basic filtering.

message("Loading methylation calls")
processed_file_suffix = "_trimmed_bismark_bt2.CpG_report.txt.gz"
processed_file_list = as.list(paste(processed_dir,"/",all_samples,processed_file_suffix,sep=""))

methRaw = methRead(location = processed_file_list,
                sample.id   = as.list(sample_names),
                treatment   = treatment_vector,
                assembly    = genome,
                context     = "CpG",
                mincov      = 10,
                pipeline    = "bismarkCytosineReport",
                header      = FALSE,
                skip        = 0,
                sep         = "\t",
                dbtype      = "tabix",
                dbdir       = rdata_dir
                )

message("Computing basic methylation stats")
outfile = file.path(output_dir, "MethylationStats.pdf")
pdf(outfile)
par(pty="s", mfrow=c(1,1))
for (i in seq_along(all_samples)) {
        getMethylationStats(methRaw[[i]],plot=TRUE,both.strands=FALSE)
        getCoverageStats(methRaw[[i]],plot=TRUE,both.strands=FALSE)
    }
dev.off()

message("Filtering samples based on read coverage")
methFiltered=filterByCoverage(methRaw, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.9, save.db=FALSE)

message("Normalizing samples based on read coverage")
methNormalized = normalizeCoverage(methFiltered)

message("Merging samples to improve coverage")
meth=methylKit::unite(methNormalized, destrand=FALSE)
rm(methFiltered, methNormalized)

message("Between-sample correlation analysis")

outfile = file.path(output_dir, "sample_correlation.pdf")
pdf(outfile, height=10, width=10)
    getCorrelation(meth,plot=TRUE)
dev.off()

outfile = file.path(output_dir, "sample_clusters.pdf")
pdf(outfile, height=10, width=10)
    clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
dev.off()

outfile = file.path(output_dir, "sample_PCA.pdf")
# pdf(outfile, height=10, width=10)
#     PCASamples(meth)
# dev.off()

# do my own better PCA plot
mat = getData(meth)
# remove rows containing NA values, they might be introduced at unite step
mat = mat[ rowSums(is.na(mat))==0, ] 
meth.mat = mat[, meth@numCs.index]/ (mat[,meth@numCs.index] + mat[,meth@numTs.index] )                                      
names(meth.mat) = meth@sample.ids
pca_obj = prcomp(meth.mat)
pca_df = data.frame(pca_obj$rotation) %>%
    rownames_to_column(var = "name") %>%
    left_join(sample_df)
p = ggplot(pca_df, aes(x=PC1,y=PC2, colour=condition)) +
    geom_point() +
    geom_text_repel(aes(label=name)) +
    theme_bw()
ggsave(outfile, p)
rm(mat)

################################################################

message("Performing tiling analysis")
tiles = tileMethylCounts(methRaw, win.size=200, step.size=200)
tiles.united = methylKit::unite(tiles)

######

message("Running pairwise differential methylation analysis")

unique_conditions = unique(conditions)
pairs_df = ldply(combn(unique_conditions, 2, simplify=FALSE))
datalist = list()
tiledatalist = list()
for (row in 1:nrow(pairs_df)) {
    condition_1 = pairs_df$V1[row]
    condition_2 = pairs_df$V2[row]
    message(condition_1)
    message(condition_2)
    df_subset = sample_df %>% filter(condition == condition_1 | condition == condition_2)
    meth_subset = reorganize(meth, df_subset$name, df_subset$treatment)
    # tiles_subset = reorganize(tiles.united, df_subset$name, df_subset$treatment)

    message("Calculating differentially methylated bases")
    myDiff = calculateDiffMeth(meth_subset, overdispersion="MN", test="Chisq", mc.cores=8, chunk.size=1e5)
    diffAnn = annotateWithGeneParts(as(myDiff,"GRanges"), annotation.all.features)

    myDiff25p = getMethylDiff(myDiff,difference=25,qvalue=0.01)
    diffAnn25p = annotateWithGeneParts(as(myDiff25p,"GRanges"), annotation.all.features)

    pdf(sprintf("annotationPercentages_%s_vs_%s_.p25q0.01.pdf", condition_1, condition_2))
        genomation::plotTargetAnnotation(diffAnn25p,precedence=TRUE, main=sprintf("%s vs %s annotation (25% changed, q<0.01)", condition_1, condition_2))
    dev.off()

    message("Calculating differentially methylated tiles")
    tileDiff = calculateDiffMeth(tiles_subset, overdispersion="MN", test="Chisq")
    tileDiffAnn = annotateWithGeneParts(as(tileDiff,"GRanges"), annotation.all.features, intersect.chr = TRUE)

    message("Exporting differential analysis data")
    
    # raw data
    diffAnn.data = getAssociationWithTSS(diffAnn) %>%
        rownames_to_column(var = "rowid") %>%
        mutate(transcript_id = gsub("\\..*","",rowid))

    myDiff.withrows = getData(myDiff) %>%
        rowid_to_column(var = "rowid")

    transcript_ids = diffAnn.data$transcript_id

    transcriptid_to_geneid = AnnotationDbi::select(transcript.db, keys = unique(transcript_ids), columns="GENEID", keytype="TXNAME") %>%
        distinct()
    all_gene_data = geneid_to_genename %>%
        inner_join(transcriptid_to_geneid, by=c("gene_id" = "GENEID"))

    all_meth_data = diffAnn.data %>%
        left_join(all_gene_data, by=c("transcript_id" = "TXNAME")) %>%
        left_join(myDiff.withrows, by=c("target.row" = "rowid")) %>%
        filter(qvalue < 0.01 & abs(meth.diff) > 25) %>%
        dplyr::select(chr, start, end, gene_id, transcript_id, gene_name, dist_to_feature = dist.to.feature, pvalue, qvalue, meth.diff)

    all_meth_output_file = sprintf(sprintf("differential_methylation_%s_vs_%s_p25q0.01.txt", condition_1, condition_2))
    message(sprintf("Writing output to %s", all_meth_output_file))
    write.table(all_meth_data, all_meth_output_file, sep="\t", row.names=FALSE, quote=FALSE)

    all_meth_data$comparison = sprintf("%s_vs_%s", condition_1, condition_2)
    datalist[[row]] = all_meth_data

    # tiling data
    message("Exporting tiled data")
    tileDiffAnn.data = getAssociationWithTSS(tileDiffAnn) %>%
        rownames_to_column(var = "rowid") %>%
        mutate(transcript_id = gsub("\\..*","",rowid))

    tileDiff.withrows = getData(tileDiff) %>%
        rowid_to_column(var = "rowid")

    transcript_ids = tileDiffAnn.data$transcript_id

    transcriptid_to_geneid = AnnotationDbi::select(transcript.db, keys = unique(transcript_ids), columns="GENEID", keytype="TXNAME") %>%
        distinct()
    all_gene_data = geneid_to_genename %>%
        inner_join(transcriptid_to_geneid, by=c("gene_id" = "GENEID"))

    all_tiling_data = tileDiffAnn.data %>%
        left_join(all_gene_data, by=c("transcript_id" = "TXNAME")) %>%
        left_join(tileDiff.withrows, by=c("target.row" = "rowid")) %>%
        filter(qvalue < 0.01) %>%
        dplyr::select(chr, start, end, gene_id, transcript_id, gene_name, dist_to_feature = dist.to.feature, pvalue, qvalue, meth.diff)

    tiled_output_file = sprintf(sprintf("differential_methylation_tiled_%s_vs_%s_p25q0.01.txt", condition_1, condition_2))
    message(sprintf("Writing tiled output to %s", tiled_output_file))
    write.table(all_tiling_data, tiled_output_file, sep="\t", row.names=FALSE, quote=FALSE)

    all_tiling_data$comparison = sprintf("%s_vs_%s", condition_1, condition_2)
    tiledatalist[[row]] = all_tiling_data

}

all_comparison_data <- dplyr::bind_rows(datalist)
all_tiling_comparison_data <- dplyr::bind_rows(tiledatalist)


compare_comparisons = all_comparison_data %>%
    dplyr::select(chr,start,end,gene_id,transcript_id,gene_name,dist_to_feature,qvalue,meth.diff,comparison) %>%
    pivot_wider(names_from=c(comparison),values_from=c(qvalue,meth.diff))

write.table(compare_comparisons, "all_comparisons_p25q0.01.txt", sep="\t", row.names=FALSE, quote=FALSE)

compare_tiling_comparisons = all_tiling_comparison_data %>%
    dplyr::select(chr,start,end,gene_id,transcript_id,gene_name,dist_to_feature,qvalue,meth.diff,comparison) %>%
    pivot_wider(names_from=c(comparison),values_from=c(qvalue,meth.diff))

write.table(compare_tiling_comparisons, "all_tiling_comparisons_q0.01.txt", sep="\t", row.names=FALSE, quote=FALSE)


message("Saving workspace")
save.image("R_workspace.Rdata")




