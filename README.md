[![Powered by RDKit](https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC)](https://www.rdkit.org/)
    
![Pandas](https://img.shields.io/badge/pandas-%23150458.svg?style=for-the-badge&logo=pandas&logoColor=white)

### Description:
This script can accept SMILES-string as input, either from a .smi-file or as a string. The script will then "charge" the input SMILES according to the chosen pH-value and range, producing either an output file or return a string with the charged structure. This script was originally written by prof. Ruth Brenk (University of Bergen), and was adapted to work with RDKit and .csv-files. 

We offer no guarantees that this script will work with all or any valid SMILES-structures.
We acknowledge the existence of a lot of spaghetti code within this script, and its dire need of refactoring.


### charger_RDKit_moldesc.py usage

if input file is named input.smi, output will by default saved in a file called output.smi:

`
charger_RDKit_moldesc.py -i ./csvs/input.smi -r 1.0 -p 4 -o
`

if input is string, and ouput should be printed:

`
charger_rdkit_moldesc.py -str CCN -R 1.0 - p 4 -o
`

Expected output: CC[NH3+], CCN 
