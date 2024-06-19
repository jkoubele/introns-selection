# introns-selection
Code for subselecting introns from reference genome and finding nested introns in splice-junction (SJ) data.

For running the python code, please install the requirement by ```pip install -r requirements.txt```
You also need to have [bedtools installed](https://bedtools.readthedocs.io/en/latest/content/installation.html) and added to PATH.
## Creating .bed file with introns from reference genome
 - The code is adapted from Ricfrid's [intership project](https://github.com/dirfcir/Pol_II_RvdM), which also refers to the original [Poll II speed project](https://github.com/beyergroup/ElongationRate/tree/main).
 - The scripts get_intron_from_gtf_file.sh and get_intron_from_gtf_file.py should yield equivalent results.
 - Please note that reference .gtf file may differ in naming convention - for example, the UTR may be recorded as '' 



