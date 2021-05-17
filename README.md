# devaplot
## Dependencies
matplotlib  
numpy  
pandas  
## Usage

```
usage: devaplot.py [-h] [-M FLOAT] [-m FLOAT] [-e INT] [-f STR] [-F] [-l] [-D INT] [-s FLOAT,FLOAT]
                   [-g INT,INT[,INT,INT[,INT,INT,...]]] [-x X_TICK] [-t STR] [-T STR] [-d INT] vcf_file

Plot genome depth with nucleotide polymorphisms

positional arguments:
  vcf_file              VCF file with AD format

optional arguments:
  -h, --help            show this help message and exit
  -M FLOAT, --major FLOAT
                        Threshold to include variant at position (default: 0.1)
  -m FLOAT, --minor FLOAT
                        Threshold to include base variant (default: 0.1)
  -e INT, --extend INT  Extend variant bar for INT position to both sides (default: 4)
  -f STR, --figure STR  Output figure (default: None)
  -F, --force           Force overwite output (default: False)
  -l, --log             Plot depth in log scale (default: False)
  -D INT, --dpi INT     Image DPI (default: 300)
  -s FLOAT,FLOAT, --size FLOAT,FLOAT
                        Size of output (default: 16,9)
  -g INT,INT[,INT,INT[,INT,INT,...]], --gap INT,INT[,INT,INT[,INT,INT,...]]
                        Gap position and size (default: None)
  -x X_TICK, --x-tick X_TICK
                        Tick interval (default: 500)
  -t STR, --table-relative STR
                        Save relative variant table to STR. Use '' for STDOUT (default: None)
  -T STR, --table-absolute STR
                        Save aboslute variant table to STR. Use '' for STDOUT (default: None)
  -d INT, --depth INT   Depth of position to report variant (default: 20)
```
## License
[GNU General Public License 3.0](https://www.gnu.org/licenses/gpl-3.0.en.html)
