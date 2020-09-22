workflow MolecularOncologyAlmanac {
    String patientId
    String tumorType
    String? stage
    File? snvHandle
    File? indelHandle
    File? segHandle
    File? fusionHandle
    File? burdenHandle
    File? germlineHandle
    File? validationHandle
    String? purity
    String? ploidy
    String? microsatellite_status
    String? whole_genome_doubling
    String? disable_matchmaking

    Int? RAM = 6
    Int? SSD = 25
    Int? preemptible = 3

    String? docker_tag = "0.4.21"

    meta {
        author: "Brendan Reardon"
        email: "breardon@broadinstitute.org"
        laboratory: "Van Allen Lab"
        institution: "Dana-Farber Cancer Institute, Broad Institute of MIT & Harvard"
        github: "https://github.com/vanallenlab/moalmanac"
        license: "Apache-2.0 License"
    }

    call almanacTask {
        input: patientId=patientId,
            tumorType=tumorType,
            stage=stage,
            snvHandle=snvHandle,
            indelHandle=indelHandle,
            segHandle=segHandle,
            fusionHandle=fusionHandle,
            burdenHandle=burdenHandle,
            germlineHandle=germlineHandle,
            validationHandle=validationHandle,
            purity=purity,
            ploidy=ploidy,
            microsatellite_status=microsatellite_status,
            whole_genome_doubling=whole_genome_doubling,
            disable_matchmaking=disable_matchmaking,
            RAM=RAM,
            SSD=SSD,
            preemptible=preemptible,
            docker_tag=docker_tag
    }
}

task almanacTask {
    String patientId
    String tumorType
    String? stage
    File? snvHandle
    File? indelHandle
    File? segHandle
    File? fusionHandle
    File? burdenHandle
    File? germlineHandle
    File? validationHandle
    String? purity
    String? ploidy
    String? microsatellite_status
    String? whole_genome_doubling
    String? disable_matchmaking

    Int? RAM
    Int? SSD
    Int? preemptible

    String? docker_tag

    command {
        if [ "${whole_genome_doubling}" == "True" ]; then
            wgd_arg="--wgd"; else
            wgd_arg="";
        fi

        if [ "${disable_matchmaking}" == "True" ]; then
            matchmaking_arg="--disable_matchmaking"; else
            matchmaking_arg="";
        fi

        python /moalmanac/moalmanac.py --patient_id ${patientId} --tumor_type ${tumorType} ${"--stage " + stage} \
        ${"--snv_handle " + snvHandle} ${"--indel_handle " + indelHandle} ${"--cnv_handle " + segHandle} \
        ${"--fusion_handle " + fusionHandle} ${"--bases_covered_handle " + burdenHandle} \
        ${"--germline_handle " + germlineHandle} ${"--validation_handle " + validationHandle} \
        ${"--purity " + purity} ${"--ploidy " + ploidy} ${"--ms_status " + microsatellite_status} \
        $wgd_arg $matchmaking_arg

        mv /moalmanac/build/index.html ${patientId}.report.html

        touch ${patientId}.sigs.context.txt ${patientId}.sigs.cosmic.txt
        touch ${patientId}.sigs.tricontext.counts.png ${patientId}.sigs.tricontext.normalized.png
        touch ${patientId}.validation_overlap.png
        touch ${patientId}.matchmaker.txt

        tar -zcf ${patientId}.almanac.tar.gz ${patientId}*
    }

    output  {
        File actionable = "${patientId}.actionable.txt"
        File somaticScored = "${patientId}.somatic.scored.txt"
        File somaticFiltered = "${patientId}.somatic.filtered.txt"
        File germlineACMG = "${patientId}.germline.acmg.txt"
        File germlineCancer = "${patientId}.germline.cancer_related.txt"
        File germlineHereditary = "${patientId}.germline.hereditary_cancers.txt"
        File integrated = "${patientId}.integrated.summary.txt"
        File matchmaker = "${patientId}.matchmaker.txt"
        File msiVariants = "${patientId}.msi_variants.txt"
        File mutationalBurden = "${patientId}.mutational_burden.txt"
        File preclinicalEfficacy = "${patientId}.preclinical.efficacy.txt"
        File deconstructSigsContext = "${patientId}.sigs.context.txt"
        File deconstructSigsCosmic = "${patientId}.sigs.cosmic.txt"
        File triContextCounts = "${patientId}.sigs.tricontext.counts.png"
        File triContextNormalized = "${patientId}.sigs.tricontext.normalized.png"
        File validationOverlap = "${patientId}.validation_overlap.png"
        File report = "${patientId}.report.html"
        File tarGz = "${patientId}.almanac.tar.gz"
    }

    runtime {
        disks: "local-disk " + SSD + " SSD"
        docker: "vanallenlab/moalmanac:" + docker_tag
        memory: RAM + " GB"
        preemptible: preemptible
    }
}
