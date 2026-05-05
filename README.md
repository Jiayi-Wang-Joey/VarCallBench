## Benchmarking Variant Calling Methods Using LongBench Datasets

This benchmark evaluates variant calling methods using the LongBench datasets ([LongBench](https://github.com/mritchielab/LongBench.io)), 
including lrRNA-seq data from the H211 and H526 cell lines generated with PacBio Kinnex (Mas-Seq), ONT cDNA R10, and ONT dRNA004 technologies.

The workflow is implemented using [Omnibenchmark](https://omnibenchmark.org/).

### Reproduction
To reproduce the results, you first need to create the apptainer image (`sif`) from the `envs/*.def` files. 
Then you run:

```bash
ob run benchmark.yaml
```

### Citation
```bash
@ARTICLE{Wang2026-lv,
  title       = "Systematic benchmarking of small variant calling pipelines for
                 long-read {RNA} sequencing data",
  author      = "Wang, Jiayi and Robinson, Mark D",
  journal     = "bioRxiv",
  institution = "bioRxiv",
  month       =  may,
  year        =  2026
}

```

### Data availability

All LongBench raw sequencing files (FASTQ) are publicly available. Please refer the paper for the links. 
The output VCFs can be found on [Zenodo](https://zenodo.org/records/19857089).

### Contact
If you have any question regarding the code, please contact jiayi.wang2@uzh.ch.
