# Nextclade sort minimizer creator

Script for generating nextclade minimizer indices from INSDC reference sequences with custom levels of specificity. Nextclade sort performs fast local alignment of sequences to a reference based on k-mers of the reference that are stored in a minimizer index.

Read more about nextclade sort, the parameters you can use when running it and how to run it with a custom minimizer [here](https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextclade-cli/reference.html#nextclade-sort).

The user can configure the number of k-mers to keep in the minimizer index - a higher level will increase match likelihood and is required to classify shorter segments, or segments with higher divergence from the reference. However, keeping more k-mers leads to increased run time.

In Pathoplexus we use nextclade sort to classify rsv sub-types and in Genspectrum it is used to classify influenza segments and HA and NA subtypes.

## Creating the minimizer

First activate your environment using
```
micromamba create -f environment.yaml
micromamba activate minimizer-creator
```
Then create a config.yaml file (like the example) with the list of all INSDC references you would like to include in your minimizer index. For example:

```yaml
references:
  seg4_H1: U08903.1
  seg4_H2: CY005413.1
  seg3__H5N1: NC_007359.1
threshold: 28 # default used by nextclade to generate minimizer indices
```

Note that the key e.g. `seg4_H1` will be the `dataset` name that is returned in the nextclade sort results.tsv.

```tsv
index	seqName	dataset	score	numHits
7	MW874350.1	seg3__H5N1	0.4829120176662018	521
7	MW874350.1	seg3__H1N1	0.34937439846005774	363
7	MW874350.1	seg3__H3N2	0.22552301255230126	252
```

You can also set the specificity of your minimizer index as a value `threshold` from 1 to 31. This specifies the fraction of kmers that are kept. (To be exact we keep all k-mers by checking if their hash is below the threshold `2^{threshold}`, where the maximum hash value is `2^32`).

Run the creation script:

```bash 
snakemake
```

## Using your minimizer in Loculus

### Ingest

The minimizer index can be used in ingest to identify segments, e.g.:

```yaml 
segment_identification:
  method: minimizer
  minimizer_index: https://anna-parker.github.io/InfluenzaAReferenceDB/results/influenza_segments.minimizer.json
  minimizer_parser:
  - segment
  - subtype
  - annotation
```
Note that you can also define the structure of the reference keys in the `minimizer_parser` to add metadata fields based on the nextclade sort results. For example the `minimizer_parser` above will split the nextclade sort key `seg6_N2_H9N2` into `segment: seg6`, `subtype: N2` and `annotation: H9N2` which can all be added to the output metadata.

### Preprocessing

In preprocessing it can be used to validate that the uploaded sequence has the highest sort score of a dataset in your `accepted_dataset_matches` list.

```yaml
require_nextclade_sort_match: true
minimizer_url: "https://anna-parker.github.io/InfluenzaAReferenceDB/results/rsv_segments.minimizer.json"
accepted_dataset_matches: 
- rsv-a
```

