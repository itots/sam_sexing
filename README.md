# sam_sexing.py

The script estimates the sex of individuals based on the ratio of the counts
of reads mapped on Y and X chromosomes. The script uses the alignment
summary statistics reported by 
[samtools idxstats](http://www.htslib.org/doc/samtools-idxstats.html). 

The script was tested with the outputs produced by samtools v1.9, in the
 environment of Python 3.7.
## Dependencies
Python 3
```
os
argparse
pandas
seaborn
matplotlib
```
## Inputs
You should specify the path to the directory containing idxstats output
 files. <br>
The idxstats output is
> TAB-delimited with each line consisting of reference 
>sequence name, sequence length, # mapped read-segments and # 
>unmapped read-segments

like this:
```
NW_021160369.1	81779	38	0
NW_021160370.1	89493	0	0
NW_021160371.1	90411	2	0
NW_021160372.1	134439	86	0
```
The files are assumed to have the extension of ".txt". If the files have
 different extensions, please specify that by using `--ext` option. If files
  have no extension, specify that by `--ext ""`.

## Usage
```
chmod +x sam_sexing.py
```
```
usage: sam_sexing.py [-h] -i IDXSTATS_DIR -y YCHROM -x XCHROM [-t THRESHOLD]
                     [-o OUT_DIR] [-p] [--min_y MIN_Y] [--min_x MIN_X]
                     [--two_thresholds] [--upper_threshold UPPER_THRESHOLD]
                     [--lower_threshold LOWER_THRESHOLD] [--ext EXT]
                     [--nbins N_BINS] [--dpi DPI] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -i IDXSTATS_DIR, --input IDXSTATS_DIR
                        Path to the directory containing idxstats files
  -y YCHROM, --ychrom YCHROM
                        Name of Y chromosome in reference sequence
  -x XCHROM, --xchrom XCHROM
                        Name of X chromosome in reference sequence
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold for sexing
  -o OUT_DIR, --output OUT_DIR
                        Name of output directory. It is current directory by default.
  -p, --plot            Output scatter plot and histogram, disabled by default
  --min_y MIN_Y         Minimum number of Y chromosome reads required for sexing
  --min_x MIN_X         Minimum number of X chromosome reads required for sexing
  --two_thresholds      Two thresholds will be used for sexing. If you use this option, you should specify upper and lower thresholds.
  --upper_threshold UPPER_THRESHOLD
                        Upper threshold for sexing
  --lower_threshold LOWER_THRESHOLD
                        Lower threshold for sexing
  --ext EXT             Extension
  --nbins N_BINS        Number of bins for histogram
  --dpi DPI             DPI for plots
  -v, --version         show program's version number and exit

```

## Examples
```
./sam_sexing.py -i ./input_dir -y NC_027914.1 -x NC_041774.1 -p
```
This provides tab-delimited file with each line consisting of sample (file) 
name, estimated sex, # reads on Y chromosome, # reads on X chromosome, and
 Y/X ratio, in current directory. The scatter plot of the number of Y and X
  chromosome reads and the histogram of Y/X ratio are also produced if `-p` option is specified.
### 


## Author
ITO, Tsuyoshi

## License
This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

