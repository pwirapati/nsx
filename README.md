# nsx: nanostring expression data analysis

Install the package:

```
library(devtools)
install_github("pwirapati/nsx")
```

Load the package

```
library(nsx)
```

Prepare the input data by unzipping the RCC files into a directory. The files can also be individually compressed (such as `*.RCC.gz` typically provided in NCBI GEO).

As an example, we use a publicly available dataset from NCBI GEO.

```
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE143382&format=file","GSE143382_RAW.tar")
untar("GSE143382_RAW.tar",exdir="RCC")    # extract the RCC files into directory RCC
list.files("RCC")
```

Read in the raw data:

```
raw_data <- readRCCset("RCC")
```

Plot data distribution:

```
rawQCplot(raw_data)
```
(Subsets can be shown by supplying a set of indices. For example, adding `o=1:12` will show only
the first 12. Adding `o=12:1` will show the same set in reversed order.)


Normalize the data:

```
norm_data <- bgsc_norm(raw_data)
```

The normalized data matrix can be found in `norm_data$z`.

Show normalization curve of an individual sample:

```
bgsc_plot( norm_data, 1)
```
(The number refers to the sample position in the normalized data matrix `norm_data$z`).
