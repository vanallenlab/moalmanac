from datasources import Oncotree


class OntologyMapper(object):
    ontology = Oncotree.ontology
    code = Oncotree.code

    @staticmethod
    def get_index_match(df, col, tumor_type):
        return df.index[df[col].str.lower() == tumor_type][0]

    @classmethod
    def get_mapped_ontology(cls, df, tumor_type):
        ontology_match = tumor_type in df[cls.ontology].str.lower().tolist()
        code_match = tumor_type in df[cls.code].str.lower().tolist()

        if ontology_match:
            idx = cls.get_index_match(df, cls.ontology, tumor_type)
        elif code_match:
            idx = cls.get_index_match(df, cls.code, tumor_type)
        else:
            return {cls.ontology: tumor_type, cls.code: 'NA'}

        return {
            cls.ontology: df.loc[idx, cls.ontology],
            cls.code: df.loc[idx, cls.code]
        }

    @classmethod
    def map(cls, dbs, tumor_type):
        oncotree = Oncotree.import_ds(dbs)
        return cls.get_mapped_ontology(oncotree, str(tumor_type).lower())
