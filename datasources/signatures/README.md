# Databases: Mutational signature aetiologies

The Molecular Oncology Almanac can optionally utilize a list of Mutational Signature `id`s and their associated aetiologies to annotate [provided mutational signature contributions](../../docs/description-of-inputs.md#mutational-signatures).

We recommend generating a file from [the Sanger Institute's website](https://cancer.sanger.ac.uk/signatures/) that follows the format of the [example file](#example).

By default, the Molecular Oncology Almanac reads an [empty version of this datasource](./aetiologies.empty.tsv).

## Example

| id | aetiology |
| --- | --- |
| SBS1 | Example aetiology A |
| SBS2 | Example aetiology B |
| SBS3 | Unknown |

## References

1. [Alexandrov, L.B., Kim, J., Haradhvala, N.J. *et al*. The repertoire of mutational signatures in human cancer. *Nature* **578**, 94–101 (2020).](https://www.nature.com/articles/s41586-020-1943-3)
2. [Díaz-Gay, M. *et al*. Assigning mutational signatures to individual samples and individual somatic mutations with SigProfilerAssignment. *Bioinformatics* **39**, btad756 (2023).](https://academic.oup.com/bioinformatics/article/39/12/btad756/7473371)
