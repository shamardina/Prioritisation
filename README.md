Disclaimer
==========

This set of scripts is not intended for general reuse as is. It is shared for reference, adaption and discussion. The input files are the result of the previous steps in the WGS-processing pipeline.

The prioritisation script is written by Olga Shamardina based on previous work and experience of Karyn Megy, Sri Deevi and Rutendo Mapeta.


Prioritisation of variants
==========================

The main script is [prioritisation.sh](prioritisation.sh), it has one required argument: the configuration file (see [example.conf](example.conf)).

```sh
./prioritisation.sh some_project.conf
```


[parse_gene_list.py](parse_gene_list.py) can be used to export from XLSX to tab-separated CSV, saved in the same location:

```sh
./parse_gene_list.py /path/to/genelist.xlsx
```


All other scripts are called from [prioritisation.sh](prioritisation.sh).

