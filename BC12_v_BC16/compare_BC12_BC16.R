# compare_BC12_BC16.R
# --- Dylan H. Ross
# --- 2021/09/13
# --- --- Comparison of DEGs from BC12 and BC16


# venn_diagram
# --- compares DEGs from BC12 and BC16 analyses and sorts out unique/shared genes
# --- considers only up/down regulated genes (determined by the direction parameter)
# --- LFCs are reversed because in this context UP/DOWN is meant to be treatment v. control
# --- --- and the LFCs were computed in the opposite direction (control v. treatment)
# --- returns the list of genes -> list(bc12, shared, bc16)
venn_diagram <- function(bc12_lfc, bc16_lfc, direction) {
    bc12_ <- c()
    shared_ <- c()
    bc16_ <- c()
    for (i in 1 : length(bc12_lfc$Geneid)) {
        gene <- bc12_lfc$Geneid[i]
        lfc <- bc12_lfc$log2FoldChange[i]
        # test whether the gene is up/down regulated as specified
        if (((lfc < 0) & (direction == "UP")) | ((lfc > 0) & (direction == "DOWN"))) {
            # test whether the gene belongs to the opposite set and has the same direction
            if (is.element(gene, bc16_lfc$Geneid)) {
                lfc2 <- bc16_lfc$log2FoldChange[which(bc16_lfc$Geneid == gene)]
                if (((lfc2 < 0) & (direction == "UP")) | ((lfc2 > 0) & (direction == "DOWN"))) {
                    shared_ <- c(shared_, gene)
                }
            } else {
                bc12_ <- c(bc12_, gene)
            }
        }
    }
    for (i in 1 : length(bc16_lfc$Geneid)) {
        gene <- bc16_lfc$Geneid[i]
        lfc <- bc16_lfc$log2FoldChange[i]
        # test whether the gene is up/down regulated as specified
        if (((lfc < 0) & (direction == "UP")) | ((lfc > 0) & (direction == "DOWN"))) {
            # test whether the gene belongs to the opposite set
            if (!(is.element(gene, bc12_lfc$Geneid))) {
                bc16_ <- c(bc16_, gene)
            }
        }
    }
    return(list(bc12 = bc12_, shared = shared_, bc16 = bc16_))
}


# main
# --- main execution sequence
main <- function() {

    # load the LFC data
    bc12_lfc <- data.frame(read.csv("BC12_v_CTRL_fliv_deseq2_results_shrinkage_padj<0.05.csv"))
    bc16_lfc <- data.frame(read.csv("BC16_v_CTRL_fliv_deseq2_results_shrinkage_padj<0.05.csv"))

    # sort shared/unique genes
    vd_up <- venn_diagram(bc12_lfc, bc16_lfc, direction = "UP")
    vd_dn <- venn_diagram(bc12_lfc, bc16_lfc, direction = "DOWN")

    # print results of venn diagram
    message(paste("BC12 genes", length(bc12_lfc$Geneid)))
    message(paste("BC16 genes", length(bc16_lfc$Geneid)))
    message("UP")
    message(paste("BC12", length(vd_up$bc12)))
    message(paste("shared", length(vd_up$shared)))
    message(paste("BC16", length(vd_up$bc16)))
    message("DOWN")
    message(paste("BC12", length(vd_dn$bc12)))
    message(paste("shared", length(vd_dn$shared)))
    message(paste("BC16", length(vd_dn$bc16)))
    message("UP shared")
    print(sort(vd_up$shared))
    message("DOWN shared")
    print(sort(vd_dn$shared))
}

main()
