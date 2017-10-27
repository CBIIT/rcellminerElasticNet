#--------------------------------------------------------------------------------------------------
# Select drugs used for elastic net analyses.
#--------------------------------------------------------------------------------------------------
library(rcellminer)

# Must include: Topotecan, Etoposide, Cisplatin

drugAnnot <- getFeatureAnnot(rcellminerData::drugData)[["drug"]]

fdaDrugAnnot <- drugAnnot[which(drugAnnot$FDA_STATUS == "FDA approved"), ]
ctDrugAnnot  <- drugAnnot[which(drugAnnot$FDA_STATUS == "Clinical trial"), ]

write.table(x = fdaDrugAnnot[, "NSC", drop=FALSE], file = "inst/extdata/cm_fda_approved.txt",
            row.names = FALSE, quote=FALSE)

write.table(x = fdaDrugAnnot[1:60, "NSC", drop=FALSE], file = "inst/extdata/cm_fda_approved_discovery.txt",
            row.names = FALSE, quote=FALSE)

write.table(x = fdaDrugAnnot[61:120, "NSC", drop=FALSE], file = "inst/extdata/cm_fda_approved_vr_home.txt",
            row.names = FALSE, quote=FALSE)

write.table(x = fdaDrugAnnot[121:140, "NSC", drop=FALSE], file = "inst/extdata/cm_fda_approved_al_mac.txt",
            row.names = FALSE, quote=FALSE)

write.table(x = fdaDrugAnnot[141:158, "NSC", drop=FALSE], file = "inst/extdata/cm_fda_approved_al_linux.txt",
            row.names = FALSE, quote=FALSE)

write.table(x = ctDrugAnnot[, "NSC", drop=FALSE], file = "inst/extdata/cm_clinical_trial.txt",
            row.names = FALSE, quote=FALSE)


cmPriority <- c("609699", "757804", "141540", "757036", "326231", "134727",
                "71423", "27640", "760087", "721517", "698037", "226080",
                "759857", "747973", "83142")
cmPriorityDrugAnnot <- drugAnnot[cmPriority, ]

write.table(x = cmPriorityDrugAnnot[, "NSC", drop=FALSE], 
            file = "inst/extdata/cm_priority.txt",
            row.names = FALSE, quote=FALSE)
#--------------------------------------------------------------------------------------------------

