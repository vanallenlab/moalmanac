import logger
from config import COLNAMES
import json


class Writer:
    section = 'outputs'
    score_bin = COLNAMES[section]['score_bin']
    almanac_bin = COLNAMES[section]['almanac_bin']
    cancerhotspots_bin = COLNAMES[section]['cancerhotspots_bin']
    cancerhotspots3D_bin = COLNAMES[section]['cancerhotspots3d_bin']
    cgc_bin = COLNAMES[section]['cgc_bin']
    gsea_pathways_bin = COLNAMES[section]['gsea_pathways_bin']
    gsea_cm_bin = COLNAMES[section]['gsea_cm_bin']
    cosmic_bin = COLNAMES[section]['cosmic_bin']
    clinvar_bin = COLNAMES[section]['clinvar_bin']
    acmg_bin = COLNAMES[section]['acmg_bin']
    hereditary_bin = COLNAMES[section]['hereditary_bin']
    msi_bin = COLNAMES[section]['msi_bin']

    sensitive_bin = COLNAMES[section]['sensitive_bin']
    resistance_bin = COLNAMES[section]['resistance_bin']
    prognostic_bin = COLNAMES[section]['prognostic_bin']
    sensitive_implication = COLNAMES[section]['sensitive_implication']
    resistance_implication = COLNAMES[section]['resistance_implication']
    prognostic_implication = COLNAMES[section]['prognostic_implication']
    sensitive_map = COLNAMES[section]['sensitive_implication_map']
    resistance_map = COLNAMES[section]['resistance_implication_map']
    prognostic_map = COLNAMES[section]['prognostic_implication_map']
    sensitive_therapy = COLNAMES[section]['sensitive_therapy']
    resistance_therapy = COLNAMES[section]['resistance_therapy']
    sensitive_therapy_strategy = COLNAMES[section]['sensitive_therapy_strategy']
    resistance_therapy_strategy = COLNAMES[section]['resistance_therapy_strategy']
    sensitive_therapy_type = COLNAMES[section]['sensitive_therapy_type']
    resistance_therapy_type = COLNAMES[section]['resistance_therapy_type']
    favorable_prognosis = COLNAMES[section]['favorable_prognosis']
    sensitive_oncotree_code = COLNAMES[section]['sensitive_oncotree_code']
    resistance_oncotree_code = COLNAMES[section]['resistance_oncotree_code']
    prognostic_oncotree_code = COLNAMES[section]['prognostic_oncotree_code']
    sensitive_description = COLNAMES[section]['sensitive_description']
    resistance_description = COLNAMES[section]['resistance_description']
    prognostic_description = COLNAMES[section]['prognostic_description']
    sensitive_url = COLNAMES[section]['sensitive_url']
    resistance_url = COLNAMES[section]['resistance_url']
    prognostic_url = COLNAMES[section]['prognostic_url']
    sensitive_citation = COLNAMES[section]['sensitive_citation']
    resistance_citation = COLNAMES[section]['resistance_citation']
    prognostic_citation = COLNAMES[section]['prognostic_citation']

    feature_type = COLNAMES[section]['feature_type']
    feature = COLNAMES[section]['feature']
    alt_type = COLNAMES[section]['alt_type']
    alt = COLNAMES[section]['alt']
    chr = COLNAMES[section]['chr']
    start = COLNAMES[section]['start']
    end = COLNAMES[section]['end']
    ref = COLNAMES[section]['ref']
    allele1 = COLNAMES[section]['allele1']
    allele2 = COLNAMES[section]['allele2']
    tumor_f = COLNAMES[section]['tumor_f']
    coverage = COLNAMES[section]['coverage']
    clinvar = COLNAMES[section]['clinvar']
    number_germline_mutations = COLNAMES[section]['number_germline_mutations']

    exac_common = COLNAMES[section]['exac_common']
    exac_af = COLNAMES[section]['exac_af']
    exac_ac = COLNAMES[section]['exac_ac']
    exac_an = COLNAMES[section]['exac_an']
    exac_afr_ac = COLNAMES[section]['exac_afr_ac']
    exac_amr_ac = COLNAMES[section]['exac_amr_ac']
    exac_eas_ac = COLNAMES[section]['exac_eas_ac']
    exac_fin_ac = COLNAMES[section]['exac_fin_ac']
    exac_nfe_ac = COLNAMES[section]['exac_nfe_ac']
    exac_sas_ac = COLNAMES[section]['exac_sas_ac']
    exac_oth_ac = COLNAMES[section]['exac_oth_ac']
    exac_afr_an = COLNAMES[section]['exac_afr_an']
    exac_amr_an = COLNAMES[section]['exac_amr_an']
    exac_eas_an = COLNAMES[section]['exac_eas_an']
    exac_fin_an = COLNAMES[section]['exac_fin_an']
    exac_nfe_an = COLNAMES[section]['exac_nfe_an']
    exac_sas_an = COLNAMES[section]['exac_sas_an']
    exac_oth_an = COLNAMES[section]['exac_oth_an']

    spanningfrags = COLNAMES[section]['spanningfrags']
    left_gene = COLNAMES[section]['left_gene']
    left_chr = COLNAMES[section]['left_chr']
    left_start = COLNAMES[section]['left_start']
    right_gene = COLNAMES[section]['right_gene']
    right_chr = COLNAMES[section]['right_chr']
    right_start = COLNAMES[section]['right_start']

    validation_tumor_f = COLNAMES[section]['validation_tumor_f']
    validation_coverage = COLNAMES[section]['validation_coverage']
    validation_detection_power = COLNAMES[section]['validation_detection_power']

    feature_display = COLNAMES[section]['feature_display']
    preclinical_efficacy = COLNAMES[section]['efficacy_obs']
    connections = COLNAMES[section]['connections']
    rationale = COLNAMES[section]['rationale']
    patient_id = COLNAMES[section]['patient_id']
    tumor = COLNAMES[section]['tumor']
    normal = COLNAMES[section]['normal']

    @staticmethod
    def create_output_name(folder, patient_id, output_suffix):
        return f'{folder}/{patient_id}.{output_suffix}'

    @staticmethod
    def export_series(series, output_name):
        series.to_csv(output_name, sep='\t', index=True)

    @staticmethod
    def export_dataframe(df, output_name):
        df.to_csv(output_name, sep='\t', index=False)

    @staticmethod
    def export_dataframe_indexed(df, output_name, index_label):
        df.to_csv(output_name, sep='\t', index_label=index_label)

    @staticmethod
    def log(label, filename, dataframe, add_line_break=False):
        logger.Messages.dataframe_size(
            label=f"Writing {label} to {filename}",
            dataframe=dataframe,
            add_line_break=add_line_break
        )

    @staticmethod
    def save_figure(figure, output_name):
        figure.savefig(output_name, bbox_inches='tight')

    @staticmethod
    def sort_columns(df, columns, ascending_boolean):
        return df.sort_values(columns, ascending=ascending_boolean)

    @staticmethod
    def return_nonzero_bin_idx(series):
        return series[series.astype(float) != 0].index


class Actionable:
    sort_columns = [Writer.almanac_bin, Writer.sensitive_map,  Writer.resistance_map, Writer.prognostic_map]
    output_columns = [Writer.score_bin,
                      Writer.sensitive_implication, Writer.resistance_implication, Writer.prognostic_implication,
                      Writer.feature_type, Writer.feature, Writer.alt_type, Writer.alt,
                      Writer.tumor_f, Writer.coverage, Writer.exac_af, Writer.exac_common, Writer.clinvar,
                      Writer.sensitive_bin,
                      Writer.sensitive_therapy, Writer.sensitive_therapy_strategy, Writer.sensitive_therapy_type,
                      Writer.sensitive_oncotree_code,
                      Writer.sensitive_description, Writer.sensitive_citation, Writer.sensitive_url,
                      Writer.resistance_bin,
                      Writer.resistance_therapy, Writer.resistance_therapy_strategy, Writer.resistance_therapy_type,
                      Writer.resistance_oncotree_code,
                      Writer.resistance_description, Writer.resistance_citation, Writer.resistance_url,
                      Writer.prognostic_bin, Writer.favorable_prognosis,
                      Writer.prognostic_oncotree_code,
                      Writer.prognostic_description, Writer.prognostic_citation, Writer.prognostic_url,
                      Writer.number_germline_mutations,
                      Writer.validation_coverage, Writer.validation_tumor_f, Writer.validation_detection_power,
                      Writer.feature_display, Writer.preclinical_efficacy,
                      Writer.patient_id, Writer.tumor, Writer.normal]

    output_suffix = 'actionable.txt'

    @classmethod
    def write(cls, df, patient_id, folder):
        df[Writer.patient_id] = patient_id
        df_sorted = Writer.sort_columns(df=df, columns=cls.sort_columns, ascending_boolean=False)
        output_dataframe = df_sorted.loc[:, cls.output_columns]
        output_name = Writer.create_output_name(folder, patient_id, cls.output_suffix)
        Writer.log(label="Clinically relevant observations", filename=output_name, dataframe=output_dataframe)
        Writer.export_dataframe(df=output_dataframe, output_name=output_name)
        return df_sorted


class GermlineACMG:
    sort_columns = [Writer.clinvar_bin, Writer.feature, Writer.feature_type]
    output_columns = [Writer.feature_type, Writer.feature, Writer.alt_type, Writer.alt,
                      Writer.chr, Writer.start, Writer.end, Writer.ref, Writer.allele1, Writer.allele2,
                      Writer.tumor_f, Writer.coverage,
                      Writer.clinvar, Writer.exac_common, Writer.exac_af, Writer.exac_ac, Writer.exac_an,
                      Writer.exac_afr_ac, Writer.exac_amr_ac, Writer.exac_eas_ac, Writer.exac_fin_ac,
                      Writer.exac_nfe_ac, Writer.exac_sas_ac, Writer.exac_oth_ac,
                      Writer.exac_afr_an, Writer.exac_amr_an, Writer.exac_eas_an, Writer.exac_fin_an,
                      Writer.exac_nfe_an, Writer.exac_sas_an, Writer.exac_oth_an,
                      Writer.patient_id, Writer.tumor, Writer.normal]

    output_suffix = 'germline.acmg.txt'

    bin = Writer.acmg_bin

    @classmethod
    def write(cls, df, patient_id, folder):
        df[Writer.patient_id] = patient_id
        df_sorted = Writer.sort_columns(df=df, columns=cls.sort_columns, ascending_boolean=False)
        idx = Writer.return_nonzero_bin_idx(df.loc[:, cls.bin])
        output_dataframe = df_sorted.loc[idx, cls.output_columns]
        output_name = Writer.create_output_name(folder, patient_id, cls.output_suffix)
        Writer.log(label="Germline variants in ACMG genes", filename=output_name, dataframe=output_dataframe)
        Writer.export_dataframe(df=output_dataframe, output_name=output_name)


class GermlineCancer:
    sort_columns = [Writer.almanac_bin, Writer.cancerhotspots_bin, Writer.cancerhotspots3D_bin,
                    Writer.cgc_bin, Writer.gsea_pathways_bin, Writer.gsea_cm_bin, Writer.cosmic_bin,
                    Writer.exac_common, Writer.exac_af]
    sort_ascending = [False, False, False, False, False, False, False, True, True]
    output_columns = [Writer.score_bin, Writer.sensitive_bin, Writer.resistance_bin, Writer.prognostic_bin,
                      Writer.feature_type, Writer.feature, Writer.alt_type, Writer.alt,
                      Writer.chr, Writer.start, Writer.end, Writer.ref, Writer.allele1, Writer.allele2,
                      Writer.tumor_f, Writer.coverage,
                      Writer.clinvar, Writer.exac_common, Writer.exac_af, Writer.exac_ac, Writer.exac_an,
                      Writer.exac_afr_ac, Writer.exac_amr_ac, Writer.exac_eas_ac, Writer.exac_fin_ac,
                      Writer.exac_nfe_ac, Writer.exac_sas_ac, Writer.exac_oth_ac,
                      Writer.exac_afr_an, Writer.exac_amr_an, Writer.exac_eas_an, Writer.exac_fin_an,
                      Writer.exac_nfe_an, Writer.exac_sas_an, Writer.exac_oth_an,
                      Writer.patient_id, Writer.tumor, Writer.normal]

    output_suffix = 'germline.cancer_related.txt'

    almanac_bin = Writer.almanac_bin
    hotspots_bin = Writer.cancerhotspots_bin
    cgc_bin = Writer.cgc_bin

    @classmethod
    def get_cancer_idx(cls, df):
        idx_almanac = Writer.return_nonzero_bin_idx(df.loc[:, cls.almanac_bin])
        idx_hotspot = Writer.return_nonzero_bin_idx(df.loc[:, cls.hotspots_bin])
        idx_cgc = Writer.return_nonzero_bin_idx(df.loc[:, cls.cgc_bin])
        return idx_almanac.union(idx_hotspot).union(idx_cgc)

    @classmethod
    def write(cls, df, patient_id, folder):
        df[Writer.patient_id] = patient_id
        df_sorted = Writer.sort_columns(df=df, columns=cls.sort_columns, ascending_boolean=cls.sort_ascending)
        idx = cls.get_cancer_idx(df)
        output_dataframe = df_sorted.loc[idx, cls.output_columns]
        output_name = Writer.create_output_name(folder, patient_id, cls.output_suffix)
        Writer.log(label="Germline variants in cancer related genes", filename=output_name, dataframe=output_dataframe)
        Writer.export_dataframe(df=output_dataframe, output_name=output_name)


class GermlineHereditary:
    sort_columns = [Writer.clinvar_bin, Writer.feature, Writer.feature_type]
    output_columns = [Writer.feature_type, Writer.feature, Writer.alt_type, Writer.alt,
                      Writer.chr, Writer.start, Writer.end, Writer.ref, Writer.allele1, Writer.allele2,
                      Writer.tumor_f, Writer.coverage,
                      Writer.clinvar, Writer.exac_common, Writer.exac_af, Writer.exac_ac, Writer.exac_an,
                      Writer.exac_afr_ac, Writer.exac_amr_ac, Writer.exac_eas_ac, Writer.exac_fin_ac,
                      Writer.exac_nfe_ac, Writer.exac_sas_ac, Writer.exac_oth_ac,
                      Writer.exac_afr_an, Writer.exac_amr_an, Writer.exac_eas_an, Writer.exac_fin_an,
                      Writer.exac_nfe_an, Writer.exac_sas_an, Writer.exac_oth_an,
                      Writer.patient_id, Writer.tumor, Writer.normal]

    output_suffix = 'germline.hereditary_cancers.txt'

    bin = Writer.hereditary_bin

    @classmethod
    def write(cls, df, patient_id, folder):
        df[Writer.patient_id] = patient_id
        df_sorted = Writer.sort_columns(df=df, columns=cls.sort_columns, ascending_boolean=False)
        idx = Writer.return_nonzero_bin_idx(df.loc[:, cls.bin])
        output_dataframe = df_sorted.loc[idx, cls.output_columns]
        output_name = Writer.create_output_name(folder, patient_id, cls.output_suffix)
        Writer.log(
            label="Germline variants in genes related to hereditary cancers",
            filename=output_name,
            dataframe=output_dataframe
        )
        Writer.export_dataframe(df=output_dataframe, output_name=output_name)


class Illustrations:
    @classmethod
    def write(cls, figure, folder, profile_id, output_suffix):
        output_name = Writer.create_output_name(folder, profile_id, output_suffix)
        Writer.save_figure(figure, output_name)


class Integrated:
    section = 'integrative'
    somatic = COLNAMES[section]['somatic']
    copynumber = COLNAMES[section]['copynumber']
    fusion = COLNAMES[section]['fusion']
    germline = COLNAMES[section]['germline']

    output_columns = [somatic, copynumber, fusion, germline]

    output_suffix = 'integrated.summary.txt'

    @classmethod
    def write(cls, df, patient_id, folder):
        df_sorted = df.sort_index()
        output_dataframe = df_sorted.loc[:, cls.output_columns]
        output_name = Writer.create_output_name(folder, patient_id, cls.output_suffix)
        Writer.log(label="Integrated summary", filename=output_name, dataframe=output_dataframe)
        Writer.export_dataframe(df=output_dataframe, output_name=output_name)


class MSI:
    sort_columns = [Writer.almanac_bin, Writer.cancerhotspots_bin, Writer.cancerhotspots3D_bin,
                    Writer.cgc_bin, Writer.gsea_pathways_bin, Writer.gsea_cm_bin, Writer.cosmic_bin,
                    Writer.exac_common, Writer.exac_af]
    sort_ascending = [False, False, False, False, False, False, False, True, True]
    output_columns = [Writer.score_bin, Writer.sensitive_bin, Writer.resistance_bin, Writer.prognostic_bin,
                      Writer.feature_type, Writer.feature, Writer.alt_type, Writer.alt,
                      Writer.chr, Writer.start, Writer.end, Writer.ref, Writer.allele1, Writer.allele2,
                      Writer.tumor_f, Writer.coverage, Writer.exac_af, Writer.exac_common, Writer.clinvar,
                      Writer.number_germline_mutations,
                      Writer.spanningfrags,
                      Writer.left_gene, Writer.left_chr, Writer.left_start,
                      Writer.right_gene, Writer.right_chr, Writer.right_start,
                      Writer.rationale, Writer.patient_id, Writer.tumor, Writer.normal,
                      Writer.almanac_bin, Writer.cancerhotspots_bin, Writer.cancerhotspots3D_bin,
                      Writer.cgc_bin, Writer.gsea_pathways_bin, Writer.gsea_cm_bin, Writer.cosmic_bin]

    output_suffix = 'msi_variants.txt'

    bin = Writer.msi_bin

    @classmethod
    def write(cls, df, patient_id, folder):
        df[Writer.patient_id] = patient_id
        df_sorted = Writer.sort_columns(df, cls.sort_columns, False)
        output_dataframe = df_sorted.loc[:, cls.output_columns]
        output_name = Writer.create_output_name(folder, patient_id, cls.output_suffix)
        Writer.log(
            label="Variants in genes related to microsatellite instability",
            filename=output_name,
            dataframe=output_dataframe
        )
        Writer.export_dataframe(df=output_dataframe, output_name=output_name)


class MutationalBurden:
    section = 'burden'
    patient = COLNAMES[section]['patient']
    tumor_type = COLNAMES[section]['tumor_type']
    ontology = COLNAMES[section]['ontology']
    code = COLNAMES[section]['code']
    bases_covered = COLNAMES[section]['bases_covered']
    n_nonsyn_mutations = COLNAMES[section]['n_nonsyn_mutations']
    mutational_burden = COLNAMES[section]['mutational_burden']
    percentile_tcga = COLNAMES[section]['percentile_tcga']
    percentile_tcga_tissue = COLNAMES[section]['percentile_tcga_tissue']
    high_burden_boolean = COLNAMES[section]['high_burden_boolean']

    output_columns = [patient, tumor_type, ontology, code,
                      bases_covered, n_nonsyn_mutations, mutational_burden,
                      percentile_tcga, percentile_tcga_tissue, high_burden_boolean]

    output_suffix = 'mutational_burden.txt'

    @classmethod
    def write(cls, df, patient_id, folder):
        df[Writer.patient_id] = patient_id
        output_dataframe = df.loc[:, cls.output_columns]
        output_name = Writer.create_output_name(folder, patient_id, cls.output_suffix)
        Writer.log(label="Tumor mutational burden", filename=output_name, dataframe=output_dataframe)
        Writer.export_dataframe(df=output_dataframe, output_name=output_name)


class PreclinicalEfficacy:
    output_suffix = 'preclinical.efficacy.txt'

    @classmethod
    def write(cls, df, patient_id, folder):
        output_name = Writer.create_output_name(folder, patient_id, cls.output_suffix)
        Writer.log(
            label="Preclinical efficacy of clinically relevant relationships",
            filename=output_name,
            dataframe=df
        )
        Writer.export_dataframe(df=df, output_name=output_name)


class PreclinicalMatchmaking:
    output_suffix = 'matchmaker.txt'

    @classmethod
    def write(cls, df, patient_id, folder):
        output_name = Writer.create_output_name(folder, patient_id, cls.output_suffix)
        Writer.log(
            label="Genomic similarity to cancer cell lines",
            filename=output_name,
            dataframe=df,
            add_line_break=True
        )
        Writer.export_dataframe(df=df, output_name=output_name)


class SomaticFiltered:
    sort_columns = [Writer.feature, Writer.feature_type]
    output_columns = [Writer.feature_type, Writer.feature, Writer.alt_type, Writer.alt,
                      Writer.chr, Writer.start, Writer.end, Writer.ref, Writer.allele1, Writer.allele2,
                      Writer.tumor_f, Writer.coverage,
                      Writer.spanningfrags,
                      Writer.left_gene, Writer.left_chr, Writer.left_start,
                      Writer.right_gene, Writer.right_chr, Writer.right_start,
                      Writer.rationale, Writer.patient_id, Writer.tumor, Writer.normal]

    output_suffix = 'somatic.filtered.txt'

    @classmethod
    def write(cls, df, patient_id, folder):
        df[Writer.patient_id] = patient_id
        df_sorted = Writer.sort_columns(df=df, columns=cls.sort_columns, ascending_boolean=False)
        output_dataframe = df_sorted.loc[:, cls.output_columns]
        output_name = Writer.create_output_name(folder, patient_id, cls.output_suffix)
        Writer.log(label="Somatic variants that were filtered", filename=output_name, dataframe=output_dataframe)
        Writer.export_dataframe(df=output_dataframe, output_name=output_name)


class SomaticScored:
    sort_columns = [Writer.almanac_bin, Writer.cancerhotspots_bin, Writer.cancerhotspots3D_bin,
                    Writer.cgc_bin, Writer.gsea_pathways_bin, Writer.gsea_cm_bin, Writer.cosmic_bin,
                    Writer.validation_detection_power, Writer.validation_coverage, Writer.number_germline_mutations,
                    Writer.exac_common, Writer.exac_af]
    sort_ascending = [False, False, False, False, False, False, False, False, False, False, True, True]
    output_columns = [Writer.score_bin, Writer.sensitive_bin, Writer.resistance_bin, Writer.prognostic_bin,
                      Writer.feature_type, Writer.feature, Writer.alt_type, Writer.alt,
                      Writer.chr, Writer.start, Writer.end, Writer.ref, Writer.allele1, Writer.allele2,
                      Writer.tumor_f, Writer.coverage, Writer.exac_af, Writer.exac_common, Writer.clinvar,
                      Writer.number_germline_mutations,
                      Writer.spanningfrags,
                      Writer.left_gene, Writer.left_chr, Writer.left_start,
                      Writer.right_gene, Writer.right_chr, Writer.right_start,
                      Writer.validation_coverage, Writer.validation_tumor_f, Writer.validation_detection_power,
                      Writer.rationale, Writer.patient_id, Writer.tumor, Writer.normal,
                      Writer.almanac_bin, Writer.cancerhotspots_bin, Writer.cancerhotspots3D_bin,
                      Writer.cgc_bin, Writer.gsea_pathways_bin, Writer.gsea_cm_bin, Writer.cosmic_bin]

    output_suffix = 'somatic.scored.txt'

    @classmethod
    def write(cls, df, patient_id, folder):
        df[Writer.patient_id] = patient_id
        df_sorted = Writer.sort_columns(df=df, columns=cls.sort_columns, ascending_boolean=cls.sort_ascending)
        output_dataframe = df_sorted.loc[:, cls.output_columns]
        output_name = Writer.create_output_name(folder, patient_id, cls.output_suffix)
        Writer.log(label="Somatic variants that were scored", filename=output_name, dataframe=output_dataframe)
        Writer.export_dataframe(df=output_dataframe, output_name=output_name)


class Strategies:
    output_suffix = 'therapeutic_strategies.txt'

    @classmethod
    def write(cls, df, patient_id, folder):
        output_name = Writer.create_output_name(folder, patient_id, cls.output_suffix)
        Writer.log(label="Therapeutic strategies that were highlighted", filename=output_name, dataframe=df)
        Writer.export_dataframe_indexed(df=df, output_name=output_name, index_label="Assertion / Strategy")


class Json:
    @staticmethod
    def write(handle, dictionary):
        with open(handle, 'w') as json_handle:
            json.dump(dictionary, json_handle, sort_keys=True, indent=4)
