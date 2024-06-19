# introns-selection

- Code for subselecting introns from reference genome and finding nested introns in splice-junction (SJ) data.

- For running the python code, please install the requirement by ```pip install -r requirements.txt```
  You also need to have [bedtools installed](https://bedtools.readthedocs.io/en/latest/content/installation.html) and
  added to PATH.
- Note that .bed files are 0-indexed with right bound of the interval being exclusive, while the .gtf files and SJ files
  from STAR (see page 14 of the [manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf))
  are 1-indexed with right bound of the interval being inclusive. Bedtools (and PyBedTools) are able to correctly parse
  different formats, but we need to keep the
  differences in mind when processing the files by different tools (custom aws / python / etc. scripts).

## Creating .bed file with introns from reference genome

- The code is adapted from Ricfrid's [intership project](https://github.com/dirfcir/Pol_II_RvdM), which also refers to
  the original [Poll II speed project](https://github.com/beyergroup/ElongationRate/tree/main).
- The scripts [get_intron_from_gtf_file.sh](./get_introns_from_gtf_file.sh)
  and [get_introns_from_gtf_file.py](./get_introns_from_gtf_file.py) should yield equivalent results.
- Please note that reference .gtf file may differ in naming convention - for example, the UTRs may be recorded as 'three_prime_utr' and 'five_prime_utr' instead of 'UTR'.
If that's the case, the code need to be modified accordingly.

## Finding nested introns

- Finding nested introns is implemented in the notebook [find_nested_introns.ipynb](./find_nested_introns.ipynb)
- Example of SJ data can be downloaded from
  this [link](https://drive.google.com/file/d/1MnhSRSMya5N33H6FdDk2OI97UUNcao1m/view?usp=sharing).
