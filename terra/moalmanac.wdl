version 1.0

workflow MolecularOncologyAlmanac {
    input {
        String patientId
        String tumorType = "Unknown"
        String stage = "Unknown"
        String config = "/moalmanac/config.ini"
        String dbs = "/moalmanac/annotation-databases.ini"
        String preclinicalDbs = "/moalmanac/preclinical-databases.ini"
        File? snvHandle
        File? indelHandle
        File? segHandle
        File? calledCNHandle
        File? fusionHandle
        File? burdenHandle
        File? germlineHandle
        File? validationHandle
        File? mutationalSignatures
        String purity = "Unknown"
        String ploidy = "Unknown"
        String microsatelliteStatus = "unk"
        Boolean wgd = false
        Boolean disableMatchmaking = false

        Int RAM = 8
        Int SSD = 50
        Int preemptible = 3

        String docker_tag = "0.9.0_v.2025-02-07"
    }

    meta {
        author: "Brendan Reardon"
        email: "breardon@broadinstitute.org"
        laboratory: "Van Allen Lab"
        institution: "Dana-Farber Cancer Institute, Broad Institute of MIT & Harvard"
        github: "https://github.com/vanallenlab/moalmanac"
        license: "GPL-2.0"
    }

    call almanacTask {
        input:
            patientId = patientId,
            tumorType = tumorType,
            stage = stage,
            config = config,
            dbs = dbs,
            preclinicalDbs = preclinicalDbs,
            snvHandle = snvHandle,
            indelHandle = indelHandle,
            segHandle = segHandle,
            calledCNHandle = calledCNHandle,
            fusionHandle = fusionHandle,
            burdenHandle = burdenHandle,
            germlineHandle = germlineHandle,
            validationHandle = validationHandle,
            mutationalSignatures = mutationalSignatures,
            purity = purity,
            ploidy = ploidy,
            microsatelliteStatus = microsatelliteStatus,
            wgd = wgd,
            disableMatchmaking = disableMatchmaking,
            RAM = RAM,
            SSD = SSD,
            preemptible = preemptible,
            docker_tag = docker_tag
    }

    output {
        File actionable = almanacTask.actionable
        File somaticScored = almanacTask.somaticScored
        File somaticFiltered = almanacTask.somaticFiltered
        File germlineACMG = almanacTask.germlineACMG
        File germlineCancer = almanacTask.germlineCancer
        File germlineHereditary = almanacTask.germlineHereditary
        File inputMetadata = almanacTask.inputMetadata
        File integrated = almanacTask.integrated
        File log = almanacTask.log
        File matchmaker = almanacTask.matchmaker
        File execution = almanacTask.execution
        File msiVariants = almanacTask.msiVariants
        File mutationalBurden = almanacTask.mutationalBurden
        File preclinicalEfficacy = almanacTask.preclinicalEfficacy
        File validationOverlap = almanacTask.validationOverlap
        File report = almanacTask.report
        File tarGz = almanacTask.tarGz
    }
}

task almanacTask {
    input {
        String patientId
        String tumorType
        String stage
        String config
        String dbs
        String preclinicalDbs
        File? snvHandle
        File? indelHandle
        File? segHandle
        File? calledCNHandle
        File? fusionHandle
        File? burdenHandle
        File? germlineHandle
        File? validationHandle
        File? mutationalSignatures
        String purity
        String ploidy
        String microsatelliteStatus
        Boolean wgd
        Boolean disableMatchmaking

        Int RAM
        Int SSD
        Int preemptible

        String docker_tag
    }

    command <<<
        python /moalmanac/moalmanac.py \
        --patient_id ~{patientId} \
        --tumor_type ~{tumorType} \
        --stage ~{stage} \
        --config ~{config} \
        --dbs ~{dbs} \
        --preclinical-dbs ~{preclinicalDbs} \
        ~{"--snv_handle " + snvHandle} \
        ~{"--indel_handle " + indelHandle} \
        ~{"--cnv_handle " + segHandle} \
        ~{"--called_cn_handle " + calledCNHandle} \
        ~{"--fusion_handle " + fusionHandle} \
        ~{"--bases_covered_handle " + burdenHandle} \
        ~{"--germline_handle " + germlineHandle} \
        ~{"--validation_handle " + validationHandle} \
        ~{"--mutational_signatures " + mutationalSignatures} \
        ~{"--purity " + purity} \
        ~{"--ploidy " + ploidy} \
        ~{"--ms_status " + microsatelliteStatus} \
        ~{if wgd then "--wgd" else ""} \
        ~{if disableMatchmaking then "--disable_matchmaking" else ""}

        touch ~{patientId}.validation_overlap.png
        touch ~{patientId}.matchmaker.txt

        mkdir docs
        cp -r /docs/* docs/

        tar -zcf ~{patientId}.almanac.tar.gz ~{patientId}* docs
    >>>

    output {
        File actionable = "~{patientId}.actionable.txt"
        File somaticScored = "~{patientId}.somatic.scored.txt"
        File somaticFiltered = "~{patientId}.somatic.filtered.txt"
        File germlineACMG = "~{patientId}.germline.acmg.txt"
        File germlineCancer = "~{patientId}.germline.cancer_related.txt"
        File germlineHereditary = "~{patientId}.germline.hereditary_cancers.txt"
        File inputMetadata = "~{patientId}.input-metadata.txt"
        File integrated = "~{patientId}.integrated.summary.txt"
        File log = "~{patientId}.log"
        File matchmaker = "~{patientId}.matchmaker.txt"
        File execution = "~{patientId}.moalmanac-execution.json"
        File msiVariants = "~{patientId}.msi_variants.txt"
        File mutationalBurden = "~{patientId}.mutational_burden.txt"
        File preclinicalEfficacy = "~{patientId}.preclinical.efficacy.txt"
        File validationOverlap = "~{patientId}.validation_overlap.png"
        File report = "~{patientId}.report.html"
        File tarGz = "~{patientId}.almanac.tar.gz"
    }

    runtime {
        disks: "local-disk ~{SSD} SSD"
        docker: "vanallenlab/moalmanac:~{docker_tag}"
        memory: "~{RAM} GB"
        preemptible: preemptible
    }
}
