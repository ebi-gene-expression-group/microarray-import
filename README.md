# microarray-import

This submodule stores scripts associated with Micorarray processing for Atlas data production.
Following microarray technologies are supported in data production.
- Affymetrix - single channel arrays
- Agilent - 1 and 2 color arrays
- Illumina - single channel arrays

# Utils

## Check CEL files in directory

Aims to check one by one the CEL files in a directory, showing any individual failure:

```
./check_cel_files_in_dir.R <dir-with-.CEL-files>
```

it will exit with non-zero status if any of the CEL files fail to be loaded with `affyio` Bioconductor package. Individual error messages are printed per file, if they come occur.
