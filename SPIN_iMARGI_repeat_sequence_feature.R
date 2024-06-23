# test whetehr SPIN DNA sequence and SPIN specific RNA sequence have similar features 

require(GenomicRanges)
require(tidyverse)
require(tibble)
require(data.table)
require(magrittr)
require(InteractionSet)
require(regioneR)

`%ni%` = Negate(`%in%`)

# ---------- load in SPIN data ----------

H1_spin <- fread("/mnt/extraids/SDSC_NFS/wenxingzhao/database/4DN/SPIN/H1_new.SPIN.JAWG.25kb.9_state.bed") %>% 
  dplyr::filter(!grepl("NAN",V4)) %>%
  makeGRangesFromDataFrame(.,seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = T)

HFF_spin <- fread("/mnt/extraids/SDSC_NFS/wenxingzhao/database/4DN/SPIN/HFF_25kb_hg38.4DN_Sheng.bed") %>% 
  filter(!grepl("NAN",V4)) %>% 
  makeGRangesFromDataFrame(.,seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = T)



STATES <- c("Speckle",  "Interior_Act1", "Interior_Act2", "Interior_Act3",
            "Interior_Repr1","Interior_Repr2","Near_Lm1","Near_Lm2","Lamina")

color_ <- c("#8b254a","#c14e4c","#ec7b57","#f2b579","#dbd291","#a8d29f", "#5fbba2","#7d9a98","#54508b")


## ---- Read in bedpe and transform into .gi ----------

bedpe_2_gi <- function(bedpe_file){
  
  bedpe <- fread(bedpe_file)
  # print(colnames(bedpe))
  
  margi_ctrl_RNA <- bedpe %>% dplyr::select(`#chrom1`,start1,end1,strand1) %>%
    dplyr::mutate(start1 = start1 + 1) %>%
    dplyr::mutate(strand1 = ifelse(strand1=="+","-","+")) %>%
    GenomicRanges::makeGRangesFromDataFrame(.,
                                            keep.extra.columns = F,
                                            seqnames.field = "#chrom1",
                                            start.field = "start1",
                                            end.field = "end1",
                                            strand.field = "strand1" )
  margi_ctrl_DNA <- bedpe %>% dplyr::select(chrom2,start2,end2,strand2) %>%
    dplyr::mutate(start2 = start2 + 1) %>%
    GenomicRanges::makeGRangesFromDataFrame(.,
                                            keep.extra.columns = F,
                                            seqnames.field = "chrom2",
                                            start.field = "start2",
                                            end.field = "end2"
    )
  
  
  margi_gi <- GInteractions(h1_margi_ctrl_RNA, h1_margi_ctrl_DNA)
  mcols(margi_gi)$anchor1.name <- NULL
  mcols(margi_gi)$anchor2.name <- NULL
  
  # filter 1k pairs
  pairdist_ = pairdist(margi_gi)
  margi_gi_1k <- margi_gi[unique(c(which(is.na(pairdist_)),which(pairdist_>1000)))]
  margi_gi_1k_nochrM = margi_gi_1k[which(as.logical(seqnames(anchors(margi_gi_1k)$first)!="chrM") & as.logical(seqnames(anchors(margi_gi_1k)$second)!="chrM"))]
  
  return(intra_margi_gi_1k_nochrM)

}

H1_gi <- bedpe_2_gi("/std_results/merged_iMARGI/iMARGI_H1_control/iMARGI_H1_control.mapq30.1k.bedpe.gz")
HFF_gi <- bedpe_2_gi("/std_results/merged_iMARGI/iMARGI_H1_control/iMARGI_HFF_control.mapq30.1k.bedpe.gz")


# ---------- load in iMARGI data ----------

get10Mgi <- function(gi){
  
  # return RNA-DNA pairs, where the genomic distances between RNA and DNA ends are longer than 10Mb
  pairdist <- pairdist(gi)
  ctrl_gi_longdist <- gi[c(which(is.na(pairdist)), which(pairdist>10000000))]
  return(ctrl_gi_longdist)
  
}

getintergi <- function(gi){
  
  # return RNA-DNA pairs, where the RNA and DNA ends mapped to different chromosomes
  pairdist <- pairdist(gi)
  ctrl_gi_longdist <- gi[which(is.na(pairdist))]
  return(ctrl_gi_longdist)
  
}

H1_10M_gi <- H1_gi %>% get10Mgi
HFF_10M_gi <- HFF_gi %>% get10Mgi

H1_inter_gi <- H1_gi %>% getintergi
HFF_inter_gi <- HFF_gi %>% getintergi

# ---------- load in repeats ----------

hg38_human_rpt <- fread("/mnt/extraids/SDSC_NFS/wenxingzhao/database/Repeat_masker_human/UCSC_tbl_hg38_repeatmasker_all_fields.bed", header = T) %>% 
  filter(repFamily %in% c("Alu", "srpRNA", "snRNA", "SVA", 
                          "L1", "L2",  "ERVL","ERV1",
                          "LTR", "Low complexity")) %>% 
  makeGRangesFromDataFrame(., seqnames.field = "genoName", start.field = "genoStart", end.field = "genoEnd", 
                           strand.field = "strand", keep.extra.columns = T)



## ---- Observed RNA repeat abundance in each SPIN ----------

rpt_expect_over_obs <- function(gi_longdist, spin_gr){
  
  
  dna_ovlp <- findOverlaps(anchors(gi_longdist)$second, spin_gr, ignore.strand=T, minoverlap = 20) %>% 
    as.data.frame() %>% set_colnames(c("readID", "SPIN_ID")) 
  
  rna_ovlp <- findOverlaps(anchors(gi_longdist)$first, hg38_human_rpt, ignore.strand=F, minoverlap = 20) %>% 
    as.data.frame() %>% set_colnames(c("readID", "rpt_ID")) 
  
  spin_rna_ct <- dna_ovlp %>% left_join(rna_ovlp) %>% filter(!is.na(rpt_ID)) %>% 
    mutate(SPIN_name = spin_gr$V4[SPIN_ID], rpt_name = hg38_human_rpt$repFamily[rpt_ID]) %>% 
    group_by(SPIN_name, rpt_name) %>% summarise(ct = n())
  
  spin_rna_ct_wider <- spin_rna_ct %>% pivot_wider(., names_from = "rpt_name", values_from = "ct")
  
  
  
  ## ---- Expect Calculation1: total repeats RNA received in each SPIN state ----------
  
  RNA_is_rpt_gi <- gi_longdist[queryHits(findOverlaps(anchors(gi_longdist)$first, hg38_human_rpt, ignore.strand=F, minoverlap = 20)) %>% unique]
  
  total_repeats_each_SPIN <- findOverlaps(anchors(RNA_is_rpt_gi)$second, spin_gr, ignore.strand=T, minoverlap = 20) %>% as.data.frame() %>% 
    set_colnames(c("readID","SPIN_ID")) %>%                
    mutate(SPIN = spin_gr$V4[SPIN_ID]) %>% 
    group_by(SPIN) %>% 
    summarise(total_reads = n()) %>% 
    mutate(total_DNA_SPIN_norm = total_reads/sum(total_reads))
  
  
  ## ---- Expect Calculation2: overall abundance of each repeat category in iMARGI RNA ends ----------
  DNA_is_SPIN_gi <- gi_longdist[queryHits(findOverlaps(anchors(gi_longdist)$second, spin_gr, ignore.strand=T, minoverlap = 20)) %>% unique]
  
  rna_ovlp <- findOverlaps(anchors(DNA_is_SPIN_gi)$first, hg38_human_rpt, ignore.strand=F, minoverlap = 20) %>% 
    as.data.frame() %>% set_colnames(c("readID", "rpt_ID")) %>%  
    mutate(rpt_type = hg38_human_rpt$repFamily[rpt_ID])
  
  rna_rpt_abundance <- rna_ovlp$rpt_type %>% table() %>% as.data.frame() %>% set_colnames(c("rpt_name", "abundance"))
  
  
  ## ---- Expect Calculation3: expected rpt abundance at each SPIN state ----------
  
  expected_rpt_at_SPIN <- rna_rpt_abundance$abundance %o% total_repeats_each_SPIN$total_DNA_SPIN_norm
  rownames(expected_rpt_at_SPIN) <- rna_rpt_abundance$rpt_name
  colnames(expected_rpt_at_SPIN) <- total_repeats_each_SPIN$SPIN
  expected_rpt_at_SPIN_df <- expected_rpt_at_SPIN %>% t() %>% as.data.frame()
  
  
  ## ---- observed vs expected ----------
  
  obs_mtx <- spin_rna_ct_wider %>% column_to_rownames("SPIN_name")
  final <- (obs_mtx/expected_rpt_at_SPIN_df) %>% rownames_to_column("SPIN_name") %>% 
    arrange(factor(SPIN_name, levels=STATES))
}


H1_rpt_final_inter <- rpt_expect_over_obs(gi_longdist = H1_inter_gi, spin_gr = H1_spin)
HFF_rpt_final_inter <- rpt_expect_over_obs(gi_longdist = HFF_inter_gi, spin_gr = HFF_spin)

H1_rpt_final_10M <- rpt_expect_over_obs(gi_longdist = H1_10M_gi, spin_gr = H1_spin)
HFF_rpt_final_10M <- rpt_expect_over_obs(gi_longdist = HFF_10M_gi, spin_gr = HFF_spin)

write.csv(x = H1_rpt_final_inter, 
          file = "/dataOS/wenxingzhao/project/Rproj/4DN_marker/II_SPIN_states_specific_RNA/Joint_analysis_SPIN-RNA/data/H1_repeat_SPIN_mtx_use_interchromosomal_pairs_only.csv",
          quote = F, row.names = T, col.names = T)

write.csv(x = HFF_rpt_final_inter, 
          file = "/dataOS/wenxingzhao/project/Rproj/4DN_marker/II_SPIN_states_specific_RNA/Joint_analysis_SPIN-RNA/data/HFF_repeat_SPIN_mtx_use_interchromosomal_pairs_only.csv",
          quote = F, row.names = T, col.names = T)

write.csv(x = H1_rpt_final_10M, 
          file = "/dataOS/wenxingzhao/project/Rproj/4DN_marker/II_SPIN_states_specific_RNA/Joint_analysis_SPIN-RNA/data/H1_repeat_SPIN_mtx_use_RD_dist_gt_10M_pairs_only.csv",
          quote = F, row.names = T, col.names = T)

write.csv(x = HFF_rpt_final_10M, 
          file = "/dataOS/wenxingzhao/project/Rproj/4DN_marker/II_SPIN_states_specific_RNA/Joint_analysis_SPIN-RNA/data/HFF_repeat_SPIN_mtx_use_RD_dist_gt_10M_pairs_only.csv",
          quote = F, row.names = T, col.names = T)


## ---- plottt ----------
require(ComplexHeatmap)

make_rna_ept_ht <- function(Obs_exp_mtx, TITLE){
  
  final_plt <- Obs_exp_mtx %>% column_to_rownames("SPIN_name") %>% 
    dplyr::select_at(c("Alu", "srpRNA","SVA","snRNA", "L1", "L2",  "ERVL","ERV1","LTR")) %>% 
    log2(.)
  
  col_fun = circlize::colorRamp2(c(-0.3, -0.1,  0, 0.1,  0.3), c("#3e76af", "#cedee9", "#ffffff","#db9ca3", "#bc332f"))
  
  ht <- Heatmap(final_plt, name = "log2(Obs/Exp)",
                heatmap_legend_param = list(
                  legend_direction = "horizontal"
                ),
                row_names_side = "left", 
                column_names_side = "bottom", cluster_rows = F, cluster_columns = F,
                col = col_fun,
                rect_gp = gpar(col = "black", lwd = 2),
                column_names_rot = 45,
                row_names_gp = gpar(fontsize=14),
                column_names_gp = gpar(fontsize=14),
                column_title = TITLE,
                column_title_gp = gpar(fontsize=14)
  )
  
  return(ht)
}

ht_H1_RNA_rpt_10M <- make_rna_ept_ht(H1_rpt_final_10M, TITLE = "RNA DNA end dist. >10M")
ht_H1_RNA_rpt_inter <- make_rna_ept_ht(H1_rpt_final_inter, TITLE = "Interchromosomal pairs")
ht_HFF_RNA_rpt_10M <- make_rna_ept_ht(HFF_rpt_final_10M, TITLE = "RNA DNA end dist. >10M")
ht_HFF_RNA_rpt_inter <- make_rna_ept_ht(HFF_rpt_final_inter, TITLE = "Interchromosomal pairs")


## ---- generate the correponding DNA SPIN repeat abundance features ----------

dna_ept_feature_matrix <- function(spin_gr){
  
  dna_SPIN_ovlp <- findOverlaps(hg38_human_rpt, spin_gr, ignore.strand=T, minoverlap = 20) %>% 
    as.data.frame() %>% set_colnames(c("rpt_ID", "SPIN_ID")) 
  
  spin_rpt_ct <- dna_SPIN_ovlp %>% 
    mutate(SPIN_name = spin_gr$V4[SPIN_ID], rpt_name = hg38_human_rpt$repFamily[rpt_ID]) %>% 
    group_by(SPIN_name, rpt_name) %>% summarise(ct = n())
  
  spin_rpt_ct_wider <- spin_rpt_ct %>% pivot_wider(., names_from = "rpt_name", values_from = "ct")
  
  SPIN_len_Mb <- split(spin_gr,spin_gr$V4) %>% lapply(.,function(x){sum(width(x))}) %>% as.numeric() %/% 1000000
  
  spin_rpt_ct_wider_norm <- spin_rpt_ct_wider %>% column_to_rownames("SPIN_name") %>% sweep(., 1, STATS = SPIN_len_Mb,FUN = "/")
  
  spin_rpt_ct_wider_norm_scale <- spin_rpt_ct_wider_norm %>% scale 
  
  spin_rpt_final_plt <- spin_rpt_ct_wider_norm_scale[STATES,] %>% as.data.frame() %>% 
    dplyr::select_at(c("Alu", "srpRNA","SVA","snRNA", "L1", "L2",  "ERVL","ERV1","LTR")) 
  
}

H1_DNA_spt_mtx <- dna_ept_feature_matrix(H1_spin)
HFF_DNA_spt_mtx <- dna_ept_feature_matrix(HFF_spin)

write.csv(x = H1_DNA_spt_mtx, 
          file = "/dataOS/wenxingzhao/project/Rproj/4DN_marker/II_SPIN_states_specific_RNA/Joint_analysis_SPIN-RNA/data/H1_repeat_SPIN_mtx.csv",
          quote = F, row.names = T, col.names = T)

write.csv(x = HFF_DNA_spt_mtx, 
          file = "/dataOS/wenxingzhao/project/Rproj/4DN_marker/II_SPIN_states_specific_RNA/Joint_analysis_SPIN-RNA/data/HFF_repeat_SPIN_mtx.csv",
          quote = F, row.names = T, col.names = T)



col_fun2 = circlize::colorRamp2(c(-2, -0.6,  0, 0.6,  2), c("#3e76af", "#cedee9", "#ffffff","#db9ca3", "#bc332f"))

plot_dna_rpt_spin <- function(dna_rpt_mtx){
  h <- Heatmap(dna_rpt_mtx, 
               heatmap_legend_param = list(
                 legend_direction = "horizontal"
               ),
               name = "zscore",
               row_names_side = "left", 
               column_names_side = "bottom", cluster_rows = F, cluster_columns = F,
               col = col_fun2,
               rect_gp = gpar(col = "black", lwd = 2),
               column_names_rot = 45,
               row_names_gp = gpar(fontsize=14),
               column_names_gp = gpar(fontsize=14),
               column_title = " \nDNA RE abundance",
               column_title_gp = gpar(fontsize=16)
  )
  return(h)
}

H1_DNA_ht <- plot_dna_rpt_spin(H1_DNA_spt_mtx)
HFF_DNA_ht <- plot_dna_rpt_spin(HFF_DNA_spt_mtx)

draw(H1_DNA_ht, heatmap_legend_side = "bottom")
draw(HFF_DNA_ht, heatmap_legend_side = "bottom")


## ---- combine all three plots for both cell line  ----------

H1_list <- ht_H1_RNA_rpt_10M + ht_H1_RNA_rpt_inter
HFF_list <- ht_HFF_RNA_rpt_10M + ht_HFF_RNA_rpt_inter

draw(H1_list, 
     column_title = "                  Repeating elements containing caRNA enrichment",
     heatmap_legend_side = "bottom",
     column_title_gp = gpar(fontsize = 16),
     gap = unit(1, "cm"))

draw(HFF_list, 
     column_title = "                  Repeating elements containing caRNA enrichment",
     heatmap_legend_side = "bottom",
     column_title_gp = gpar(fontsize = 16),
     gap = unit(1, "cm"))


