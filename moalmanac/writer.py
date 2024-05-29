from config import COLNAMES
import json


class Writer:
    section = "outputs"

    def __init__(self, strings):
        self.strings = strings

    """
    Defining properties for each datasource bin
    """
    @property
    def score_bin(self):
        return self.strings[self.section]['score_bin']

    @property
    def almanac_bin(self):
        return self.strings[self.section]['almanac_bin']

    @property
    def cancer_hotspots_bin(self):
        return self.strings[self.section]['cancer_hotspots_bin']

    @property
    def cgc_bin(self):
        return self.strings[self.section]['cgc_bin']

    @property
    def gsea_pathways_bin(self):
        return self.strings[self.section]['gsea_pathways_bin']

    @property
    def gsea_cm_bin(self):
        return self.strings[self.section]['gsea_cm_bin']

    @property
    def cosmic_bin(self):
        return self.strings[self.section]['cosmic_bin']

    @property
    def clinvar_bin(self):
        return self.strings[self.section]['clinvar_bin']

    @property
    def acmg_bin(self):
        return self.strings[self.section]['acmg_bin']

    @property
    def hereditary_bin(self):
        return self.strings[self.section]['hereditary_bin']

    @property
    def msi_bin(self):
        return self.strings[self.section]['msi_bin']

    """
    Defining properties for moalmanac specific annotations
    """

    @property
    def sensitive_bin(self):
        return self.strings[self.section]['sensitive_bin']

    @property
    def resistance_bin(self):
        return self.strings[self.section]['resistance_bin']

    @property
    def prognostic_bin(self):
        return self.strings[self.section]['prognostic_bin']

    @property
    def sensitive_implication(self):
        return self.strings[self.section]['sensitive_implication']

    @property
    def resistance_implication(self):
        return self.strings[self.section]['resistance_implication']

    @property
    def prognostic_implication(self):
        return self.strings[self.section]['prognostic_implication']

    @property
    def sensitive_map(self):
        return self.strings[self.section]['sensitive_map']

    @property
    def resistance_map(self):
        return self.strings[self.section]['resistance_map']

    @property
    def prognostic_map(self):
        return self.strings[self.section]['prognostic_map']

    @property
    def sensitive_therapy(self):
        return self.strings[self.section]['sensitive_therapy']

    @property
    def resistance_therapy(self):
        return self.strings[self.section]['resistance_therapy']

    @property
    def sensitive_therapy_strategy(self):
        return self.strings[self.section]['sensitive_therapy_strategy']

    @property
    def resistance_therapy_strategy(self):
        return self.strings[self.section]['resistance_therapy_strategy']

    @property
    def sensitive_therapy_type(self):
        return self.strings[self.section]['sensitive_therapy_type']

    @property
    def resistance_therapy_type(self):
        return self.strings[self.section]['resistance_therapy_type']

    @property
    def favorable_prognosis(self):
        return self.strings[self.section]['favorable_prognosis']

    @property
    def sensitive_oncotree_code(self):
        return self.strings[self.section]['sensitive_oncotree_code']

    @property
    def resistance_oncotree_code(self):
        return self.strings[self.section]['resistance_oncotree_code']

    @property
    def prognostic_oncotree_code(self):
        return self.strings[self.section]['prognostic_oncotree_code']

    @property
    def sensitive_description(self):
        return self.strings[self.section]['sensitive_description']

    @property
    def resistance_description(self):
        return self.strings[self.section]['resistance_description']

    @property
    def prognostic_description(self):
        return self.strings[self.section]['prognostic_description']

    @property
    def sensitive_url(self):
        return self.strings[self.section]['sensitive_url']

    @property
    def resistance_url(self):
        return self.strings[self.section]['resistance_url']

    @property
    def prognostic_url(self):
        return self.strings[self.section]['prognostic_url']

    @property
    def sensitive_citation(self):
        return self.strings[self.section]['sensitive_citation']

    @property
    def resistance_citation(self):
        return self.strings[self.section]['resistance_citation']

    @property
    def prognostic_citation(self):
        return self.strings[self.section]['prognostic_citation']

    """
    Defining properties for describing biomarkers
    """

    @property
    def feature_type(self):
        return self.strings[self.section]['feature_type']

    @property
    def feature(self):
        return self.strings[self.section]['feature']

    @property
    def alt_type(self):
        return self.strings[self.section]['alt_type']

    @property
    def alt(self):
        return self.strings[self.section]['alt']

    @property
    def chr(self):
        return self.strings[self.section]['chr']

    @property
    def start(self):
        return self.strings[self.section]['start']

    @property
    def end(self):
        return self.strings[self.section]['end']

    @property
    def ref(self):
        return self.strings[self.section]['ref']

    @property
    def allele1(self):
        return self.strings[self.section]['allele1']

    @property
    def allele2(self):
        return self.strings[self.section]['allele2']

    @property
    def tumor_f(self):
        return self.strings[self.section]['tumor_f']

    @property
    def coverage(self):
        return self.strings[self.section]['coverage']

    @property
    def clinvar(self):
        return self.strings[self.section]['clinvar']

    @property
    def number_germline_mutations(self):
        return self.strings[self.section]['number_germline_mutations']

    """
    Defining properties for ExAC annotations
    """

    @property
    def exac_common(self):
        return self.strings[self.section]['exac_common']

    @property
    def exac_af(self):
        return self.strings[self.section]['exac_af']

    @property
    def exac_ac(self):
        return self.strings[self.section]['exac_ac']

    @property
    def exac_an(self):
        return self.strings[self.section]['exac_an']

    @property
    def exac_afr_ac(self):
        return self.strings[self.section]['exac_afr_ac']

    @property
    def exac_amr_ac(self):
        return self.strings[self.section]['exac_eas_ac']

    @property
    def exac_eas_ac(self):
        return self.strings[self.section]['exac_eas_ac']

    @property
    def exac_fin_ac(self):
        return self.strings[self.section]['exac_nfe_ac']

    @property
    def exac_nfe_ac(self):
        return self.strings[self.section]['exac_nfe_ac']

    @property
    def exac_sas_ac(self):
        return self.strings[self.section]['exac_sas_ac']

    @property
    def exac_oth_ac(self):
        return self.strings[self.section]['exac_oth_ac']

    @property
    def exac_afr_an(self):
        return self.strings[self.section]['exac_afr_an']

    @property
    def exac_amr_an(self):
        return self.strings[self.section]['exac_amr_an']

    @property
    def exac_eas_an(self):
        return self.strings[self.section]['exac_eas_an']

    @property
    def exac_fin_an(self):
        return self.strings[self.section]['exac_fin_an']

    @property
    def exac_nfe_an(self):
        return self.strings[self.section]['exac_nfe_an']

    @property
    def exac_sas_an(self):
        return self.strings[self.section]['exac_sas_an']

    @property
    def exac_oth_an(self):
        return self.strings[self.section]['exac_oth_an']

    """
    Defining properties to describe fusions
    """

    @property
    def spanningfrags(self):
        return self.strings[self.section]['spanningfrags']

    @property
    def left_gene(self):
        return self.strings[self.section]['left_gene']

    @property
    def left_chr(self):
        return self.strings[self.section]['left_chr']

    @property
    def left_start(self):
        return self.strings[self.section]['left_start']

    @property
    def right_gene(self):
        return self.strings[self.section]['right_gene']

    @property
    def right_chr(self):
        return self.strings[self.section]['right_chr']

    @property
    def right_start(self):
        return self.strings[self.section]['right_start']

    """
    Define properties for annotations to annotate somatic variants with those observed in validation sequencing 
    """

    @property
    def validation_tumor_f(self):
        return self.strings[self.section]['validation_tumor_f']

    @property
    def validation_coverage(self):
        return self.strings[self.section]['validation_coverage']

    @property
    def validation_detection_power(self):
        return self.strings[self.section]['validation_detection_power']

    """
    Defining remaining properties
    """
    @property
    def feature_display(self):
        return self.strings[self.section]['feature_display']

    @property
    def preclinical_efficacy(self):
        return self.strings[self.section]['preclinical_efficacy']

    @property
    def connections(self):
        return self.strings[self.section]['connections']

    @property
    def rationale(self):
        return self.strings[self.section]['rationale']

    @property
    def patient_id(self):
        return self.strings[self.section]['patient_id']

    @property
    def tumor(self):
        return self.strings[self.section]['tumor']

    @property
    def normal(self):
        return self.strings[self.section]['normal']

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
    def write(cls, df, patient_id, strings, folder):
        df[Writer.patient_id] = patient_id
        df_sorted = Writer.sort_columns(df, cls.sort_columns, False)
        output_name = Writer.create_output_name(folder, patient_id, cls.output_suffix)
        Writer.export_dataframe(df_sorted.loc[:, cls.output_columns], output_name)
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
        df_sorted = Writer.sort_columns(df, cls.sort_columns, False)
        idx = Writer.return_nonzero_bin_idx(df.loc[:, cls.bin])
        output_name = Writer.create_output_name(folder, patient_id, cls.output_suffix)
        Writer.export_dataframe(df_sorted.loc[idx, cls.output_columns], output_name)


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
        df_sorted = Writer.sort_columns(df, cls.sort_columns, cls.sort_ascending)
        idx = cls.get_cancer_idx(df)
        output_name = Writer.create_output_name(folder, patient_id, cls.output_suffix)
        Writer.export_dataframe(df_sorted.loc[idx, cls.output_columns], output_name)


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
        df_sorted = Writer.sort_columns(df, cls.sort_columns, False)
        idx = Writer.return_nonzero_bin_idx(df.loc[:, cls.bin])
        output_name = Writer.create_output_name(folder, patient_id, cls.output_suffix)
        Writer.export_dataframe(df_sorted.loc[idx, cls.output_columns], output_name)


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
        output_name = Writer.create_output_name(folder, patient_id, cls.output_suffix)
        Writer.export_dataframe_indexed(df_sorted.loc[:, cls.output_columns], output_name, Writer.feature)


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
        output_name = Writer.create_output_name(folder, patient_id, cls.output_suffix)
        Writer.export_dataframe(df_sorted.loc[:, cls.output_columns], output_name)


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
        output_name = Writer.create_output_name(folder, patient_id, cls.output_suffix)
        Writer.export_dataframe(df.loc[:, cls.output_columns], output_name)


class PreclinicalEfficacy:
    output_suffix = 'preclinical.efficacy.txt'

    @classmethod
    def write(cls, df, patient_id, folder):
        output_name = Writer.create_output_name(folder, patient_id, cls.output_suffix)
        Writer.export_dataframe(df, output_name)


class PreclinicalMatchmaking:
    output_suffix = 'matchmaker.txt'

    @classmethod
    def write(cls, df, patient_id, folder):
        output_name = Writer.create_output_name(folder, patient_id, cls.output_suffix)
        Writer.export_dataframe(df, output_name)


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
        df_sorted = Writer.sort_columns(df, cls.sort_columns, False)
        output_name = Writer.create_output_name(folder, patient_id, cls.output_suffix)
        Writer.export_dataframe(df_sorted.loc[:, cls.output_columns], output_name)


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
        df_sorted = Writer.sort_columns(df, cls.sort_columns, cls.sort_ascending)
        output_name = Writer.create_output_name(folder, patient_id, cls.output_suffix)
        Writer.export_dataframe(df_sorted.loc[:, cls.output_columns], output_name)


class Strategies:
    output_suffix = 'therapeutic_strategies.txt'

    @classmethod
    def write(cls, df, patient_id, folder):
        output_name = Writer.create_output_name(folder, patient_id, cls.output_suffix)
        Writer.export_dataframe_indexed(df, output_name, 'Assertion / Strategy')


class Json:
    @staticmethod
    def write(handle, dictionary):
        with open(handle, 'w') as json_handle:
            json.dump(dictionary, json_handle, sort_keys=True, indent=4)
